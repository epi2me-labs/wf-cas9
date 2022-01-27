#!/usr/bin/env nextflow

// Developer notes
// 
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended pratices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion.
//   iii) a second concreate, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress'
include { start_ping; end_ping } from './lib/ping'


process summariseReads {
    // concatenate fastq and fastq.gz in a dir

   label "cas9"
    cpus 1
    input:
        tuple path(directory), val(sample_id), val(type)
    output:
        tuple val(sample_id), path("${sample_id}.stats"), emit: stats
        tuple val(sample_id), path("${sample_id}.fastq"), emit: reads

    shell:
    """
    fastcat -s ${sample_id} -r ${sample_id}.stats -x ${directory} > ${sample_id}.fastq
    """
}


process getVersions {
   label "cas9"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
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

process build_index{
    /*
    Build minimap index from reference genome
    */
    label "cas9"
    cpus params.threads

    input:
        file reference
    output:
        path "genome_index.mmi", emit: index
    script:
    """
        minimap2 -t $params.threads -x map-ont -d genome_index.mmi $reference
    """
}

process align_reads {
    /*
    TODO: The number of off-target is quite a lot higher than in the tutorial
    Tutorial uses mini_align rather than minimap2 directly.

    mini_align \
    -r "$reference_genome" -i "$input_file" \
    -p "$output_folder/alignments" \
    -t 4 -m
    */
    label "cas9"
    input:
        path index
        path reference
        tuple val(sample_id), file(fastq_reads)
    output:
        tuple val(sample_id), path("${sample_id}.sam"), emit: sam
        tuple val(sample_id), path("${sample_id}_fastq_pass.bed"), emit: bed
    script:
    """
    minimap2 -t $params.threads -m 4 -ax map-ont $index $fastq_reads > ${sample_id}.sam
    bedtools bamtobed -i ${sample_id}.sam | bedtools sort > ${sample_id}_fastq_pass.bed
    """
}

process overlaps {
    /* Call the python processing script and get back CSVs that will be used in the report
    emits
        target_coverage: tiled csv for creating plots
        coverage summary: csv for table with stats for all targets aggregated
        target_summary: csv with summary info for each plot


     */
    label "cas9"
    input:
        path targets_bed
        path genome
        tuple val(sample_id),
              path(alignment_bed)
    output:
        tuple val(sample_id), path('*coverage_summary.csv'), emit: coverage_summary
        tuple val(sample_id), path('*target_coverage.csv'), emit: target_coverage
        tuple val(sample_id), path('*target_summary.csv'), emit: target_summary
    script:
    """
    target_overlaps.py \
    $targets_bed \
    $alignment_bed \
    $genome \
    $sample_id
    """
}

process makeReport {
   label "cas9"
    input:
        path "versions/*"
        path "params.json"
        tuple val(sample_ids),
              path(seq_summaries),
              path(coverage_summary),
              path(target_coverage)
    output:
        path "wf-cas9-*.html", emit: report
    script:
        report_name = "wf-cas9-" + params.report_name + '.html'
    """
    report.py $report_name \
        --summaries $seq_summaries \
        --versions versions \
        --params params.json \
        --coverage_summary $coverage_summary \
        --target_coverage $target_coverage \
        --sample_ids $sample_ids
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory

    publishDir "${params.out_dir}/output/${sample_id}", mode: 'copy', pattern: "*"
    input:
        tuple val(sample_id), path(fname)
    output:
        path fname
    """
    echo "Writing output files"
    echo $fname
    """
}

process output_report {
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*report.html"

    input:
        path fname
    output:
        path fname
    """
    echo "Copying report"
    """
}


// workflow module
workflow pipeline {
    take:
        reads
        ref_genome
        targets
    main:
        println(reads)
        build_index(ref_genome)
        summariseReads(reads)
//         sample_ids = summariseConcatReads.out.summary.flatMap({it -> it[0]})
        software_versions = getVersions()
        workflow_params = getParams()

        align_reads(
            build_index.out.index,
            ref_genome,
            summariseReads.out.reads
        )
        overlaps(targets,
                ref_genome,
                align_reads.out.bed
        )

        report = makeReport(software_versions.collect(),
                        workflow_params,
                        summariseReads.out.stats
                        .join(overlaps.out.coverage_summary)
                        .join(overlaps.out.target_coverage)
                        )
    emit:
        results = summariseReads.out.stats
        report
        // TODO: use something more useful as telemetry
        telemetry = workflow_params
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    start_ping()

    ref_genome = file(params.ref_genome, type: "file")
    if (!ref_genome.exists()) {
        println("--ref_genome: File doesn't exist, check path.")
        exit 1
    }
    targets = file(params.targets, type: "file")
    if (!targets.exists()) {
        println("--targets: File doesn't exist, check path.")
        exit 1
    }

    samples = fastq_ingress(
        params.fastq, params.out_dir, params.sample, params.sample_sheet, params.sanitize_fastq)


    pipeline(samples, ref_genome, targets)
    output(pipeline.out.results)
    output_report(pipeline.out.report)
    end_ping(pipeline.out.telemetry)
}
