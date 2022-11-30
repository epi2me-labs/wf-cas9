# wf-cas9

This repository contains a [nextflow](https://www.nextflow.io/) workflow
for the multiplexed analysis of Cas9-enrichment sequencing. 
## Introduction
This workflow generated a report summarising the results of Cas9 enrichment sequencing.

Users provide a reference genome, FASTQ ONT reads, and a BED file containing enrichment regions.
The reads are first mapped the reference genome using [minimap2](https://github.com/lh3/minimap2), and
various plots and tables are generated summarising the enrichment results.
## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker or singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-cas9 --help
```
to see the options for the workflow.

### Workflow inputs:
* Folder of fastq reads.
* Genome reference file.
* Target BED file with 4 columns:
  * chromosome
  * start
  * end
  * target_name


To test on a small dataset with two targets and two chromosomes:
```shell
cd wf-cas9
nextflow run . --fastq test_data/fastq/ --ref_genome \
test_data/grch38/grch38_chr19_22.fa.gz --targets test_data/targets.bed
```

To evaluate on a larger dataset, use the evaluation script:
```
evaluation run_evaluation.sh <out_dir> [optional_nexflow config]
```
**Workflow outputs**

The primary outputs of the workflow include:

* A folder per sample containing:
  * BAM file filtered to contain reads overlapping with targets (*_on_target.bam).
  * BED file with alignment information for on-target reads (*on_target.bed).
  * A simple text file providing a summary of sequencing reads (*.stats).
* sample_summary.csv - read and alignment summary for each sample.
* target_summary.csv - read and alignment summary for reads overlapping each target.
* A combined HTML report detailing the primary findings of the workflow across all samples.
By default, the report contains sequencing quality plots and two tables that summarise targeted sequencing results:
* On/off-target reads per sample.
* Summaries of each sample/target pair.

Using `--full_report`, the report will also contain the following elements that may be useful for
diagnosing issues with the experiment. These are turned off by default as they can lead to slow loading of the
HTML report:
* Plots of stranded coverage at each target.
* Histograms of on and off-target coverage for each sample.
* Off-target hotspot region tables.
## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html)