## Introduction
The ONT Cas9 sequencing kit allows the enrichment of gemomic
regions of interest by amplifying targets regions using adapters bound to Cas9 cleaved regions.
This purpose of this workflow is to assess the effectiveness of such Cas9 enrichment, 
but can be applied to other approaches. 

This workflow provides a report that is useful to help understand the effectiveness and issues
of the enrichement straetegy. 

Inputs to the workflow are: a reference genome file, fastq reads from enrichment sequencing ,
and a BED file detailing the regions of interests (targets).
The main output is a report with summary statistics and plots which illustrate 
the effectiveness of enrichment. Other outputs inclide a BAM alignment file covering the
target regions and a BED file with similar imformation.

The workflow innvolves the genome alignment of input reads using 
[minimap2](https://github.com/lh3/minimap2) and the analaysis 
of read-target overlap with [bedtools](https://github.com/arq5x/bedtools2).


