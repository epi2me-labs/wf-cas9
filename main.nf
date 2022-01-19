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
        tuple val(sample_id), path("${sample_id}.sam"), emit: aligned_reads
    script:
    """
    minimap2 -t $params.threads -m 4 -ax map-ont $index $fastq_reads > ${sample_id}.sam
    """
}




process makeReport {
   label "cas9"
    input:
        path "versions/*"
        path "params.json"
        tuple val(sample_ids), path(seq_summaries)
    output:
        path "wf-cas9-*.html", emit: report
    script:
        report_name = "wf-cas9-" + params.report_name + '.html'
    """
    report.py $report_name \
        --versions versions \
        --params params.json
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
   label "cas9"
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

//         report = makeReport(software_versions.collect(),
//                         workflow_params
//                         summariseReads.out.stats
//                         )
    emit:
        results = summariseReads.out.stats
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
//     output(pipeline.out.results)
    end_ping(pipeline.out.telemetry)
}
