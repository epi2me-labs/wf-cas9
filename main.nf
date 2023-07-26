#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress'


process getVersions {
   label "cas9"
    cpus 1
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
    cpus 1
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
    cpus 1
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
    input:
        path "reference"
    output:
        path "genome_index.mmi", emit: index
        path "chrom.sizes", emit: chrom_sizes
    script:
    """
        minimap2 -t $task.cpus -x map-ont -d genome_index.mmi reference
        samtools faidx reference
        cut -f 1,2 reference.fai >> chrom.sizes
    """
}

process align_reads {
    label "cas9"
    cpus params.threads
    memory params.minimap2_max_memory
    input:
        path "genome_index.mmi"
        path "reference"
        tuple val(meta), path("reads.fastq")
    output:
        path "${meta.alias}_aln_stats.csv", emit: aln_stats
        tuple val(meta), path("${meta.alias}_fastq_pass.bed"), emit: bed
        tuple val(meta), path("${meta.alias}.bam"), emit: bam
    script:
    """
    minimap2 -t $task.cpus -ax map-ont "genome_index.mmi" "reads.fastq" | \
        samtools sort -O bam -@ $task.cpus - | tee "${meta.alias}.bam" | \
        bedtools bamtobed -i stdin | sort -k 1,1 -k2,2n --parallel $task.cpus > "${meta.alias}_fastq_pass.bed"
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
    cpus 1
    input:
        path 'targets.tsv'
        path 'tiles.tsv'
        path 'tile_target_intersection.tsv'
        tuple val(meta),
              path('align.bed')
    output:
        path('target_cov.bed'), emit: target_coverage


    script:
    """
    # Get alignment coverage at tiles per strand

    echo "chr\tstart\tend\ttarget\tcoverage\tstrand\tsample_id" > target_cov.bed

    if grep -q "\\W+" align.bed
      then
        cat align.bed | grep "\\W+" | bedtools coverage -a tile_target_intersection.tsv -b - -wa | \
            cut -f 1,2,3,8,9 | sed "s/\$/\t+\t${meta.alias}/" >> target_cov.bed
      else
        echo "_\t0\t1\ttest_id\t0\t+" > p.bed
        cat p.bed| grep "\\W+" | bedtools coverage -a tile_target_intersection.tsv -b - | \
            cut -f 1,2,3,8,9 | sed "s/\$/\t+\t${meta.alias}/" >> target_cov.bed
    fi

    # Add the strand sample_id columns
    if grep -q "\\W-" align.bed
      then
        cat align.bed | grep "\\W-" | bedtools coverage -a tile_target_intersection.tsv -b - | \
            cut -f 1,2,3,8,9 | sed "s/\$/\t-\t${meta.alias}/" >> target_cov.bed
      else
        echo "_\t0\t1\ttest_id\t0\t-\n" > n.bed;
        cat n.bed | grep "\\W-" | bedtools coverage -a tile_target_intersection.tsv -b - | \
            cut -f 1,2,3,8,9 | sed "s/\$/\t-\t${meta.alias}/" >> target_cov.bed
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
    cpus 1
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
    cat "align.bed" | bedtools intersect -a - -b "targets.bed" -wb > aln_targets.bed

    # chr, start, stop, target, overlaps, covered_bases, len(target), frac_covered
    # This forms first few columns of output table
    bedtools coverage -a "targets.bed" -b "align.bed" > target_summary_temp.bed

    # Get alignment coverage at tiles per strand
    cat "align.bed" | bedtools coverage -a "tiles_inter_targets.bed" -b -  > target_cov.bed

    # Get median coverage (col 9) by target (col 8)
    bedtools groupby -i target_cov.bed -g 8 -c 9 -o median | cut -f 2  > median_coverage.bed

    # Strand bias
    cat aln_targets.bed | grep "\\W+\\W" | bedtools coverage -b - -a "targets.bed" | cut -f 5  > pos.bed  || true
    cat aln_targets.bed | grep "\\W-\\W" | bedtools coverage -b - -a "targets.bed" | cut -f 5  > neg.bed || true

    paste target_summary_temp.bed \
        median_coverage.bed \
        pos.bed \
        neg.bed > ${meta.alias}_target_summary.bed

    # Add sample_id column
    sed -i "s/\$/\t${meta.alias}/" ${meta.alias}_target_summary.bed

    rm median_coverage.bed pos.bed neg.bed
    """
}

process coverage_summary {
    label "cas9"
    cpus 1
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

    cat "${meta.alias}_on_target.bed" off.bed > "${meta.alias}_read_to_target.bed"

    bases_on=\$(cat "${meta.alias}_on_target.bed"   | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}')
    bases_off=\$(cat off.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}')
    rm off.bed

    echo "\${numread_on}\t\${numread_off}\n\${bases_on}\t\${bases_off}" > ${meta.alias}_coverage_summary.csv

    # Add sample id columns
    sed "s/\$/\t${meta.alias}/" "${meta.alias}_coverage_summary.csv"

    sed -i "s/\$/\t${meta.alias}/" "${meta.alias}_read_to_target.bed"
    """
}


process background {
    label "cas9"
    cpus 1
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
    cpus 1
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
    cpus 1
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
   label "cas9"
   cpus 1
    input:
        path "versions/*"
        path "params.json"
        path 'seq_summaries'
        path 'tile_coverage.tsv'
        path target_coverage
        path 'target_summary_table.tsv'
        path 'coverage_summary.tsv'
        path off_target_hotspots

    output:
        path "wf-cas9-*.html", emit: report
    script:
        report_name = "wf-cas9-report.html"
        def opttcov = target_coverage.name.startsWith('OPTIONAL_FILE') ? '' : "--target_coverage ${target_coverage}"
        def optbghot = off_target_hotspots.name.startsWith('OPTIONAL_FILE') ? '' : "--off_target_hotspots ${off_target_hotspots}"

    """
    workflow-glue report $report_name \
        --read_stats seq_summaries \
        --versions versions \
        --params params.json \
        --tile_coverage tile_coverage.tsv \
        --coverage_summary coverage_summary.tsv \
        --target_summary target_summary_table.tsv \
        $opttcov \
        $optbghot
    """
}

process pack_files_into_sample_dirs {
    label "cas9"
    cpus 1
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
process output {
    // publish inputs to output directory
    label "cas9"
    cpus 1
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

        // put fastcat stats into results channels
        stats = input_reads
        .map { meta, reads, stats_dir -> [meta, stats_dir.resolve('per-read-stats.tsv')] }

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

        tar_cov_tsv = target_coverage.out.target_coverage.collectFile(name: 'target_coverage', keepHeader: true)
        tile_cov = background.out.tiles_coverage.collectFile(name: 'tile_cov', keepHeader: true)
        bg_hotspots = background.out.hotspots.collectFile(name: 'hotspots')

        read_stats = stats
                        .map {it -> it[1]}
                        .collectFile(keepHeader:true, name: 'stats')

        build_tables(
            coverage_summary.out.read_to_target.collectFile(name: 'read_to_target'),
            align_reads.out.aln_stats.collectFile(name: 'aln_stats', keepHeader: true),
            target_summary.out.table.collectFile(name: 'target_summary')
        )
        report = makeReport(
                    software_versions,
                    workflow_params,
                    read_stats,
                    tile_cov,
                    tar_cov_tsv,
                    build_tables.out.target_summary,
                    build_tables.out.read_target_summary,
                    bg_hotspots,
        )
      
        pack_files_into_sample_dirs(
            coverage_summary.out.on_target_bed.concat(
            get_on_target_reads.out.ontarget_fastq,
            get_on_target_bams.out.on_target_bam,
            stats).groupTuple())

        results = makeReport.out.report
             .concat(build_tables.out.sample_summary,
             build_tables.out.target_summary,
             pack_files_into_sample_dirs.out.results_dir)

    emit:
        results
        telemetry = workflow_params
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }

    ref_genome = file(params.reference_genome, type: "file")
    if (!ref_genome.exists()) {
        println("--ref_genome: File doesn't exist, check path.")
        exit 1
    }
    targets = file(params.targets, type: "file")
    if (!targets.exists()) {
        println("--targets: File doesn't exist, check path.")
        exit 1
    }
    def line
    targets.withReader { line = it.readLine() }
    if (line.split("\t").size() != 4){
        println('Target file should have 4 cols: chr start end target_name')
        exit 1
    }

    // reads = fastq_ingress([
    //     "input":params.fastq,
    //     "sample":params.sample,
    //     "sample_sheet":params.sample_sheet])
    //     .map {it -> [it[1], it[0]]}
    
    reads = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "fastcat_stats": true
    ])

    pipeline(reads, ref_genome, targets)
    output(pipeline.out.results)
}

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }

    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }

}