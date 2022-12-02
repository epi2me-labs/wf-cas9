## Introduction
The ONT Cas9 sequencing kit allows the enrichment of gemomic
regions of interest by amplifying target regions from adapters ligated to Cas9 cleavage sites.
The purpose of this workflow is to assess the effectiveness of such Cas9 enrichment, 
but it can be applied to other enrichment approaches. The workflow outputs
help assess the effectiveness 
of the enrichement strategy and can be used to diagnose issues such as poorly performing probles.

Inputs to the workflow are: a reference genome file, fastq reads from enrichment sequencing,
and a BED file detailing the regions of interests (targets).
The main outputs are a report containing summary statistics and plots which give an overview of 
the enrichment, and a BAM file with target-overlapping reads.

The main steps of the workflow are alignemnt of reds to the genome using 
[minimap2](https://github.com/lh3/minimap2) and the analaysis
of read-target overlap with [bedtools](https://github.com/arq5x/bedtools2).




