# Workflow template

This repository contains a [nextflow](https://www.nextflow.io/) workflow
for the multiplexed analysis of Cas9-targeted sequencing. 

## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker of conda is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-template --help
```
to see the options for the workflow.

### Workflow inputs:
* folder of fastq reads.
* genome reference file.
* target bed file with 4 columns:
  * chromosome
  * start
  * end
  * target_name
 

To test on a small dataset with two targets and two chromosomes:
```shell
cd wf-cas9
nextflow run . --fastq test_data/fastq/ --ref_genome \
test_data/grch38/grch38_chr19_22.fa.gz --targets test_data/targets.bed \
-profile conda -resume
```

To evaluate on a larger dataset, use the evaluation script:
``` 
evaluation run_evaluation.sh <out_dir> [optional_nexflow config]
```
**Workflow outputs**

The primary outputs of the workflow include:

* a per-sample on-target reads fastq file. 
* a per-sample simple text file providing a summary of sequencing reads.
* a combined HTML report document detailing the primary findings of the workflow across all samples.

By default, the report contains sequencing quality plots and two tables that summarize targeted sequencing results:
* on/off-target reads per sample.
* summaries of each sample/target pair. 

Using the `--debug_mode`, the report will also contain the following elements that may be useful for 
diagnosing issues with the experiment. These are turned off by default as they can lead to slow loading of the 
html report:
* plots of stranded coverage at each target.
* histograms of on and off-target coverage for each sample.
* off-target hotspot region tables.






## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)
