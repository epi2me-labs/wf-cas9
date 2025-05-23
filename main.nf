#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/ingress'


process getVersions {
   label "cas9"
   publishDir "${params.out_dir}", mode: 'copy', pattern: "versions.txt"
    cpus 1
    memory "2 GB"
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    bedtools --version | sed 's/^bedtools /betools,/' >> versions.txt
    """
}


process getParams {
   label "cas9"
   cache false
   publishDir "${params.out_dir}", mode: 'copy', pattern: "params.json"
    cpus 2
    memory "2 GB"
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

process make_tiles {
    label 'cas9'
    cpus 2
    memory "2 GB"
    input:
        path "chrom.sizes"
        path "targets.bed"
    output:
        path 'tiles.bed', emit: tiles
        path 'tiles_int_targets.bed', emit: tiles_inter_targets
    script:
    """
    bedtools makewindows -g chrom.sizes -w 100 -i 'srcwinnum' | gzip > tiles.bed
    bedtools intersect -a tiles.bed -b targets.bed -wb > tiles_int_targets.bed
    """
}

process build_index{
    /*
    Build minimap index from reference genome
    */
    label "cas9"
    cpus params.threads
    memory "15 GB"
    input:
        path "reference"
    output:
        path "genome_index.mmi", emit: index
        path "chrom.sizes", emit: chrom_sizes
    script:
    """
        minimap2 -t $task.cpus -I 16G -x map-ont -d genome_index.mmi reference
        samtools faidx reference
        cut -f 1,2 reference.fai >> chrom.sizes
    """
}

process align_reads {
    label "cas9"
    cpus Math.min(params.threads, 20)
    memory "15 GB"
    input:
        path "genome_index.mmi"
        path reference_fasta
        tuple val(meta), path("reads.fastq")
    output:
        path "${meta.alias}_aln_stats.csv", emit: aln_stats
        tuple val(meta), path("${meta.alias}_fastq_pass.bed"), emit: bed
        tuple val(meta), path("${meta.alias}.bam"), emit: bam
    script:
    def mm2_threads = Math.max(task.cpus - 4, 1)
    """
    minimap2 -t ${mm2_threads} -K 20M -ax map-ont "${reference_fasta}" "reads.fastq" | \
        samtools sort -m 400M -O bam -@ 2 - | tee "${meta.alias}.bam" | \
        bedtools bamtobed -i stdin | sort -k 1,1 -k2,2n > "${meta.alias}_fastq_pass.bed"
    # Get a csv with columns: [read_id, alignment_accuracy]
    samtools index "${meta.alias}.bam"
    bamstats "${meta.alias}.bam" | \
        # Add sample id column
        sed "s/\$/\t${meta.alias}/" | \
        # Fix header
        sed '1s/${meta.alias}/sample_id/' > "${meta.alias}_aln_stats.csv"
    """
}

process target_coverage {
    /* Call the python processing script and get back CSVs that will be used in the report
    emits
        target_coverage: tiled csv for creating plots
    # NOTE
    use \W\+\W as strand may move columns in future versions

    emits tsv with these columns
        chr start end target cov_f cov_r sample_id

     */
    label "cas9"
    cpus 3
    memory "2 GB"
    input:
        path 'targets.tsv'
        path 'tiles.tsv'
        path 'tile_target_intersection.tsv'
        tuple val(meta),
              path('align.bed')
    output:
        tuple val(meta),
              path("*_target_cov.tsv"),
              emit: target_coverage


    script:
    """
    # Get alignment coverage at tiles per strand, and add sample_id column
    echo "chr\tstart\tend\ttarget\tcoverage\tstrand\tsample_id" > "${meta.alias}_target_cov.tsv"

    if grep -q "\\W+" align.bed
      then
        cat align.bed | grep "\\W+" | bedtools coverage -a tile_target_intersection.tsv -b - -wa | \
            cut -f 1,2,3,8,9 | sed "s/\$/\t+\t${meta.alias}/" >> "${meta.alias}_target_cov.tsv"
      else
        echo "_\t0\t1\ttest_id\t0\t+" > p.bed
        cat p.bed| grep "\\W+" | bedtools coverage -a tile_target_intersection.tsv -b - | \
            cut -f 1,2,3,8,9 | sed "s/\$/\t+\t${meta.alias}/" >> "${meta.alias}_target_cov.tsv"
    fi

    if grep -q "\\W-" align.bed
      then
        cat align.bed | grep "\\W-" | bedtools coverage -a tile_target_intersection.tsv -b - | \
            cut -f 1,2,3,8,9 | sed "s/\$/\t-\t${meta.alias}/" >> "${meta.alias}_target_cov.tsv"
      else
        echo "_\t0\t1\ttest_id\t0\t-\n" > n.bed;
        cat n.bed | grep "\\W-" | bedtools coverage -a tile_target_intersection.tsv -b - | \
            cut -f 1,2,3,8,9 | sed "s/\$/\t-\t${meta.alias}/" >> "${meta.alias}_target_cov.tsv"
    fi
    """

}

process target_summary {
    /*
    Make a target summary bed file with a row per target. Columns:
        chr,
        start,
        end,
        target,
        number of reads,
        num bases covered,
        target length,
        fracTargAln,
        medianCov,
        num positive
        strand reads,
        num negative,
        strand reads
    */
    label "cas9"
    cpus 2
    memory "2 GB"
    input:
        path "targets.bed"
        path "tiles.bed"
        path "tiles_inter_targets.bed"
        path "chrom.sizes"
        tuple val(meta),
              path("align.bed")
    output:
        path('*_target_summary.bed'), emit: table
    script:
    """
    # Map targets to aln.
    cat "align.bed" | bedtools intersect -a - -b "targets.bed" -wb  > aln_targets.bed

    # chr, start, stop, target, overlaps, covered_bases, len(target), frac_covered
    # This forms first few columns of output table. Sort by target
    bedtools coverage -a "targets.bed" -b "align.bed" | sort -k 4 > target_summary_temp.bed

    # Get alignment coverage at tiles per strand
    cat "align.bed" | bedtools coverage -a "tiles_inter_targets.bed" -b -  > target_cov.bed

    # Get median coverage (col 9) by target (col 8)
    bedtools groupby -i target_cov.bed -g 8 -c 9 -o median | sort -k 1 | cut -f 2  > median_coverage.bed

    # Strand bias
    # First add strand column to the targets BED
    cat targets.bed | awk -v OFS="\t" '{print \$0, "0", "+"}' > target_strand.bed
    # Then get coverage per strand
    cat aln_targets.bed | bedtools coverage -s -b - -a target_strand.bed |  sort -k 4 | cut -f 7  > pos.bed  || true
    cat aln_targets.bed | bedtools coverage -S -b - -a target_strand.bed |  sort -k 4 | cut -f 7  > neg.bed || true

    paste target_summary_temp.bed \
        median_coverage.bed \
        pos.bed \
        neg.bed > tmp1

    # Add sample_id column
    sed "s/\$/\t${meta.alias}/" tmp1 > tmp2

    # Add run_id column
    sed "s/\$/\t${meta.run_ids.join(',')}/" tmp2 > ${meta.alias}_target_summary.bed

    rm median_coverage.bed pos.bed neg.bed tmp1 tmp2
    """
}

process coverage_summary {
    label "cas9"
    cpus 2
    memory "2 GB"
    input:
        path 'targets.bed'
        tuple val(meta),
              path('align.bed')
    output:
        path("${meta.alias}_coverage_summary.csv"), emit: summary
        path("${meta.alias}_read_to_target.bed"), emit: read_to_target
        tuple val(meta),
              path("*_on_target.bed"),
              emit: on_target_bed
    script:
    """
    # Get all non-intersecting reads from aln/targets using '-v'
     bedtools intersect -a "align.bed" -b 'targets.bed' -wa -wb -v \
        | cut -f 1-4 \
        | awk -F '\\t' -v OFS='\\t' '{print \$0,"off_target"}' > off.bed

    bedtools intersect -a 'align.bed' -b 'targets.bed' -wa -wb | cut -f 1-4,10  > "${meta.alias}_on_target.bed"

    numread_on=\$(cat "${meta.alias}_on_target.bed" | wc -l | tr -d ' ')
    numread_off=\$(cat off.bed | wc -l | tr -d ' ')

    cat "${meta.alias}_on_target.bed" off.bed > tmp

    bases_on=\$(cat "${meta.alias}_on_target.bed"   | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}')
    bases_off=\$(cat off.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}')
    rm off.bed

    echo "\${numread_on}\t\${numread_off}\n\${bases_on}\t\${bases_off}" > ${meta.alias}_coverage_summary.csv

    # Add sample id columns
    sed "s/\$/\t${meta.alias}/" tmp > "${meta.alias}_read_to_target.bed"
    """
}


process background {
    label "cas9"
    cpus 2
    memory "4 GB"
    input:
        path 'targets.tsv'
        path 'tiles.tsv'
        path 'chrom_sizes.tsv'
        tuple val(meta),
              path('align.bed')
    output:
        path('off_target_hotspots.bed'), emit: hotspots
        path('coverage_at_tiles.tsv'), emit: tiles_coverage
    script:
    """
    # For each tile that does not intersect with a slopped target (intersect -v) get the coverage.
    bedtools slop -i targets.tsv -g chrom_sizes.tsv -b 1000 \
        | tee  targets_padded.bed \
        | bedtools intersect -v -a align.bed -b - -wa \
        | bedtools coverage -a tiles.tsv -b - \
        | awk '\$5 > 0 {print \$5}' \
        | sed "s/\$/\t${meta.alias}\toff_target/" > off_target_cov_at_tiles.tsv

    # write header
    echo "cov\tsample_id\ttarget_status" > coverage_at_tiles.tsv

    # For each tile does intersect with a slopped target get the coverage.
    bedtools slop -i targets.tsv -g chrom_sizes.tsv -b 1000 \
        | bedtools intersect -a align.bed -b - -wa \
        | bedtools coverage -a tiles.tsv -b - \
        | awk '\$5 > 0 {print \$5}' \
        | sed "s/\$/\t${meta.alias}\ton_target/"  >> coverage_at_tiles.tsv

    # Coverage_at tiles is a TSV with cols: cov, sample_id, target_status (on_target/off_target)
    cat off_target_cov_at_tiles.tsv >> coverage_at_tiles.tsv

    # Get all contiguous regions of background alignments
    cat targets_padded.bed | bedtools intersect -a align.bed -b - -v  | \
        bedtools merge -i - | bedtools coverage -a - -b align.bed | \
        cut -f 1-4 > off_target_hotspots.bed

    # Cols: chrom. start, stop, coverage, sample_id
    sed "s/\$/\t${meta.alias}/" off_target_hotspots.bed > tmp2
    mv tmp2 off_target_hotspots.bed
    """
}


process get_on_target_reads {
    label "cas9"
    cpus 2
    memory "2 GB"
    input:
        tuple val(meta),
              path("input.fastq"),
              path("on_target.bed")
    output:
         tuple val(meta),
               path("${meta.alias}_ontarget.fastq"),
               emit: ontarget_fastq
    script:
    """
    cat "on_target.bed" | cut -f 4 > seqids
    cat "input.fastq"| seqkit grep -f seqids -o "${meta.alias}_ontarget.fastq"
    """
}


process get_on_target_bams {
    label "cas9"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta),
              path("on_target.bed"),
              path("input.bam")
    output:
        tuple val(meta),
              path("${meta.alias}_on_target.bam"),
              emit: on_target_bam

    script:
    """
    samtools view "input.bam" -L "on_target.bed" \
        -O bam > ${meta.alias}_on_target.bam
    """
}


process build_tables {
    label "cas9"
    cpus 2
    input:
        path 'read_to_target.tsv'
        path 'aln_summary.tsv'
        path 'target_summary.tsv'
    output:
        path 'target_summary.csv', emit: target_summary
        path 'sample_summary.csv', emit: sample_summary
        path 'read_target_summary.tsv', emit: read_target_summary
    script:
    """
    workflow-glue build_tables \
        --target_summary target_summary.tsv \
        --read_to_target read_to_target.tsv \
        --aln_summary aln_summary.tsv
    """
}

process makeReport {
   label "wf_common"
   publishDir "${params.out_dir}", mode: 'copy', pattern: "wf-cas9-*.html"
   cpus 2
   memory "4 GB"
    input:
        path "versions/*"
        path "params.json"
        path 'per-read-stats.tsv'
        path 'tile_coverage.tsv'
        path target_coverage
        path 'target_summary_table.tsv'
        path 'coverage_summary.tsv'
        path off_target_hotspots
        val wf_version

    output:
        path "wf-cas9-*.html", emit: report
    script:
        report_name = "wf-cas9-report.html"
        def opttcov = target_coverage.name.startsWith('OPTIONAL_FILE') ? '' : "--target_coverage ${target_coverage}"
        def optbghot = off_target_hotspots.name.startsWith('OPTIONAL_FILE') ? '' : "--off_target_hotspots ${off_target_hotspots}"

    """
    workflow-glue report $report_name \
        --read_stats per-read-stats.tsv \
        --versions versions \
        --params params.json \
        --tile_coverage tile_coverage.tsv \
        --coverage_summary coverage_summary.tsv \
        --target_summary target_summary_table.tsv \
        --wf_version $wf_version \
        $opttcov \
        $optbghot
    """
}


process combine_stats {
    label "cas9"
    cpus 2
    memory "2 GB"
    input:
        tuple val(meta),
              path('stats.tsv.gz')
        output:
            path('stats.tsv')
    """
    gunzip -c stats.tsv.gz |
        # Add sample_id column
        sed "s/\$/\t${meta.alias}/" > stats.tsv
    """
}

process pack_files_into_sample_dirs {
    label "cas9"
    cpus 2
    memory "2 GB"
    input:
        tuple val(meta),
              path(sample_files)
    output:
        path "${meta.alias}", emit: results_dir
    """
    mkdir "${meta.alias}"
    for file in $sample_files; do
        mv \$file "${meta.alias}"
    done;
    """
}

// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process publish {
    // publish inputs to output directory
    label "cas9"
    cpus 2
    memory "4 GB"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}


// workflow module
workflow pipeline {
    take:
        input_reads
        ref_genome
        targets
    main:

        build_index(ref_genome)

        per_read_stats = input_reads
            .map { meta, reads, stats_dir -> [meta, stats_dir.resolve('per-read-stats.tsv.gz')] }

        // remove fastcat stats from reads channel
        reads = input_reads.map { meta, reads, stats_dir -> [meta, reads] }

        //summariseReads(reads)
        software_versions = getVersions()
        workflow_params = getParams()

        align_reads(
            build_index.out.index,
            ref_genome,
            reads)

        make_tiles(
            build_index.out.chrom_sizes,
            targets)

        coverage_summary(
            targets,
            align_reads.out.bed)

        get_on_target_reads(
            reads
            .join(coverage_summary.out.on_target_bed))

        get_on_target_bams(
            coverage_summary.out.on_target_bed
            .join(align_reads.out.bam))

        target_summary(targets,
            make_tiles.out.tiles,
            make_tiles.out.tiles_inter_targets,
            build_index.out.chrom_sizes,
            align_reads.out.bed)


        target_coverage(targets,
            make_tiles.out.tiles,
            make_tiles.out.tiles_inter_targets,
            align_reads.out.bed)

        background(targets,
            make_tiles.out.tiles,
            build_index.out.chrom_sizes,
            align_reads.out.bed)

        tar_cov_tsv = target_coverage.out.target_coverage
            .map {meta, target_cov -> target_cov}
            .collectFile(name: 'target_coverage', keepHeader: true)

        tile_cov = background.out.tiles_coverage.collectFile(name: 'tile_cov', keepHeader: true)
        bg_hotspots = background.out.hotspots.collectFile(name: 'hotspots')

        build_tables(
            coverage_summary.out.read_to_target.collectFile(name: 'read_to_target'),
            align_reads.out.aln_stats.collectFile(name: 'aln_stats', keepHeader: true),
            target_summary.out.table.collectFile(name: 'target_summary')
        )

        report = makeReport(
                    software_versions,
                    workflow_params,
                    per_read_stats  | combine_stats | collectFile(keepHeader: true),
                    tile_cov,
                    tar_cov_tsv,
                    build_tables.out.target_summary,
                    build_tables.out.read_target_summary,
                    bg_hotspots,
                    workflow.manifest.version
        )

        pack_files_into_sample_dirs(
            coverage_summary.out.on_target_bed.concat(
            get_on_target_reads.out.ontarget_fastq,
            get_on_target_bams.out.on_target_bam,
            target_coverage.out.target_coverage,
            per_read_stats).groupTuple())

        results = makeReport.out.report
             .concat(build_tables.out.sample_summary,
             build_tables.out.target_summary,
             pack_files_into_sample_dirs.out.results_dir)

    emit:
        results
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    Pinguscript.ping_start(nextflow, workflow, params)

    ref_genome = file(params.reference_genome, type: "file")
    if (!ref_genome.exists()) {
        error "--reference_genome: File doesn't exist, check path."
    }
    targets = file(params.targets, type: "file")
    if (!targets.exists()) {
        error "--targets: File doesn't exist, check path."
    }

    def line
    targets.withReader { line = it.readLine() }
    if (line.split("\t").size() != 4){
        error 'Target file should have 4 cols: chr start end target_name'
    }

    samples = fastq_ingress([
        "input":params['fastq'],
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "analyse_unclassified":params.analyse_unclassified,
        "stats": true,
        "per_read_stats": true,
        "fastcat_extra_args": ""])

    pipeline(samples, ref_genome, targets)
    publish(pipeline.out.results)
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)

}