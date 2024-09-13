# Cas9 enrichment workflow

Summarise the results of Cas9 enrichment sequencing.



## Introduction

<!---This section of documentation typically contains a list of things the workflow can perform also any other intro.--->
The ONT Cas9 sequencing kit allows the enrichment of genomic
regions of interest by amplifying target regions from adapters ligated to Cas9 cleavage sites.
The purpose of this workflow is to assess the effectiveness of such Cas9 enrichment, 
but it can be applied to other enrichment approaches. The workflow outputs
help assess the effectiveness of the enrichement strategy and can be used to diagnose issues such as poorly performing probes.

This workflow can be used for the following:

+ To obtain simple Cas9 enrichment sequencing summaries.
+ Statistics of the coverage of each target.
+ Plots of coverage across each target.
+ Identification of off-target hot-spots



## Compute requirements

Recommended requirements:

+ CPUs = 16
+ Memory = 64 GB

Minimum requirements:

+ CPUs = 6
+ Memory = 16 GB

Approximate run time: Approximately 30 minutes to process a single sample of 100K reads wit the minimum requirements.

ARM processor support: True




## Install and run


These are instructions to install and run the workflow on command line.
You can also access the workflow via the
[EPI2ME Desktop application](https://labs.epi2me.io/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage
compute and software resources,
therefore Nextflow will need to be
installed before attempting to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop)
or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)
to provide isolation of the required software.
Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.
This is controlled by the
[`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles)
parameter as exemplified below.

It is not required to clone or download the git repository
in order to run the workflow.
More information on running EPI2ME workflows can
be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow.
This will pull the repository in to the assets folder of
Nextflow and provide a list of all parameters
available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-cas9 --help
```
To update a workflow to the latest version on the command line use
the following command:
```
nextflow pull epi2me-labs/wf-cas9
```

A demo dataset is provided for testing of the workflow.
It can be downloaded and unpacked using the following commands:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-cas9/wf-cas9-demo.tar.gz
tar -xzvf wf-cas9-demo.tar.gz
```
The workflow can then be run with the downloaded demo data using:
```
nextflow run epi2me-labs/wf-cas9 \
	--fastq 'wf-cas9-demo/fastq/sample_1' \
	--full_report \
	--reference_genome 'wf-cas9-demo/grch38/grch38_chr19_22.fa.gz' \
	--targets 'wf-cas9-demo/targets.bed' \
	-profile standard
```

For further information about running a workflow on
the command line see https://labs.epi2me.io/wfquickstart/




## Related protocols

<!---Hyperlinks to any related protocols that are directly related to this workflow, check the community for any such protocols.--->

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).



## Input example

<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts FASTQ files as input.

The FASTQ input parameters for this workflow accept one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second cases (i and ii), a sample name can be supplied with `--sample`. In the last case (iii), the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`.

```
(i)                     (ii)                 (iii)
input_reads.fastq   ─── input_directory  ─── input_directory
                        ├── reads0.fastq     ├── barcode01
                        └── reads1.fastq     │   ├── reads0.fastq
                                             │   └── reads1.fastq
                                             ├── barcode02
                                             │   ├── reads0.fastq
                                             │   ├── reads1.fastq
                                             │   └── reads2.fastq
                                             └── barcode03
                                              └── reads0.fastq
```



## Input parameters

### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| reference_genome | string | FASTA reference file. | Full path to a FASTA reference genome file conating the target regions of interest. |  |
| targets | string | A tab-delimited BED file of target regions. | Each row should contain the following fields: chromosome/contig name, start and end of target region, and a name to identify the target. An example row would look like this: `chr19 13204400 13211100 SCA6` |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |
| full_report | boolean | Select this option to write a full report that contains plots giving a graphical representation of coverage at each target region. | In cases where there are many target to visualise, the report loading time can be slow and so it's is recommended to set `full_report` to false in such cases. | False |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Number of CPU threads to use per workflow task. | The total CPU resource used by the workflow is constrained by the executor configuration. | 8 |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | ./wf-cas9-report.html | Report for all samples | aggregated |
| Sample summary table | ./sample_summary.csv | Summary statistics for each sample. | aggregated |
| Target summary table | ./target_summary.csv | Summary statistics for each sample-target combination. | aggregated |
| Per file read stats | ./{{ alias }}/{{ alias }}_per-read-stats.tsv.gz | A TSV with per-file read statistics, including all samples. | per-sample |
| On-target BAM | .//{{ alias }}/{{ alias }}_on_target.bam | BAM file containing alignments that map to one of the given targets. | per-sample |
| On-target BED | .//{{ alias }}/{{ alias }}_on_target.bed | BED file containing summarising alignments that map to one of the given targets. | per-sample |
| On-target FASTQ | .//{{ alias }}/{{ alias }}_on_target.fastq | FASTQ file containing reads that map to one of the given targets. | per-sample |




## Pipeline overview

<!---High level numbered list of main steps of the workflow and hyperlink to any tools used. If multiple workflows/different modes perhaps have subheadings and numbered steps. Use nested numbering or bullets where required.--->
### 1. Concatenate input files and generate per read stats.

The [fastcat/bamstats](https://github.com/epi2me-labs/fastcat) tool is used to concatenate multifile samples to be processed by the workflow.
### 2. Align read to the reference genome.
The reads are then aligned to the reference genome supplied by the user using [minimap2](https://github.com/lh3/minimap2).

### 3. Generate target coverage data.
This stage of the workflow generates coverage data for each of targets that are used to make
the plots and tables in the report.

First, the reference genome is split into consecutive windows of 100bp. The coverage statistics are 
calcaulted across these windows.

The user must supply a tab-delinted BED file (`targets`) detailing the genomic locations of the targets of interest.
Here is an example file describing two targets.

```
chr19	13204400	13211100	SCA6
chr22	45791500	45799400	SCA10
```
The columns are:
+ chromosome
+ start position
+ end position
+ target name

*Note*: the file does not contain column names. 

With the target locations defined, [bedtools](https://bedtools.readthedocs.io/en/latest/) is used to generate
information regarding the alignment coverage at each of the targets and also of 
background reads (see section 4) per strand. 

### 4. Background identification
In a sequencing enrichment experiment, it can useful to know if there are any off-target genomic regions (regions that are not defined in the `targets` file) that are being preferentially encountered. This information can be used to inform primer design.

We define off-target regions here as any region not within 1kb of a defined target.
Hot-spots are further defined as contiguous regions of off-target alignemnts containing at least 10 reads. 





## Troubleshooting

<!---Any additional tips.--->
+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).



## FAQ's

<!---Frequently asked questions, pose any known limitations as FAQ's.--->

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-template/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).



## Related blog posts


See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.



