# wf-cas9

wf-cas9 is a [nextflow](https://www.nextflow.io/) workflow
for the multiplexed analysis of Oxford Nanopore Cas9 enrichment sequencing. 




## Introduction
The ONT Cas9 sequencing kit allows the enrichment of genomic
regions of interest by amplifying target regions from adapters ligated to Cas9 cleavage sites.
The purpose of this workflow is to assess the effectiveness of such Cas9 enrichment, 
but it can be applied to other enrichment approaches. The workflow outputs
help assess the effectiveness 
of the enrichement strategy and can be used to diagnose issues such as poorly performing probes.

Inputs to the workflow are: a reference genome file, FASTQ reads from enrichment sequencing,
and a BED file detailing the regions of interests (targets).
The main outputs are a report containing summary statistics and plots which give an overview of 
the enrichment, and a BAM file with target-overlapping reads.

The main steps of the workflow are alignemnt of reads to the genome using 
[minimap2](https://github.com/lh3/minimap2) and the analaysis
of read-target overlap with [bedtools](https://github.com/arq5x/bedtools2).








## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[singularity](https://docs.sylabs.io/guides/latest/user-guide/) to provide isolation of
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

The main inputs are:
* Folder of FASTQ reads.
* Genome reference file.
* Target BED file with 4 columns:
  * chromosome
  * start
  * end
  * target_name


To test on a small dataset with two targets and two chromosomes:
first download and unpack the demo data
```shell
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-cas9/wf-cas9-demo.tar.gz \
  && tar -xvf wf-cas9-demo.tar.gz

```shell
nextflow run epi2me-labs/wf-cas9 \
  --fastq wf-cas9-demo/fastq/ \
  --reference_genome wf-cas9-demo/grch38/grch38_chr19_22.fa.gz \
  --targets wf-cas9-demo/targets.bed
```

**Workflow outputs**

The primary outputs of the workflow include:

* A folder per sample containing:
  * BAM file filtered to contain reads overlapping with targets (*_on_target.bam).
  * BED file with alignment information for on-target reads (*on_target.bed).
  * BED file containing windowed coverage for each target (*target_cov.bed).
  * A simple text file providing a summary of sequencing reads (*.stats).
* sample_summary.csv - read and alignment summary for each sample.
* target_summary.csv - read and alignment summary for reads overlapping each target.
* A combined HTML report detailing the primary findings of the workflow across all samples including:
  * Sequencing quality plots. 
  * Tables summarising targeted sequencing results.
  * Plots of stranded coverage at each target.
  * Histograms of on and off-target coverage for each sample.
  * Off-target hotspot region tables.




## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html)