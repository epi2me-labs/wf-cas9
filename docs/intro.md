## Introduction
Cas9 enrichemnt sequencing 
This workflow is suitable for analysisng the effectiveness of Cas9 enrighment strategies
but can be applied to other appeaches. The ONT Cas9 sequencing kit allows the enrichment of
regions of interest using by cleaving and apapter-binding with Cas9 flanking regions of.

This workflow provides a report that is useful to help understand the effectiveness and issues
of the enrichement straetegy. 

The user needs to supply a reference genome file, fastq reads from their enrichment sequencing 
experiment and a BED file detailing the regions of interests (targets) that they were aiming to enrich.
The main output is a report with summary statistics and plots which illustrate 
the effectiveness of enrichment. Other outputs inclide a BAM alignemnt file covering the
target regions and a BED file with simiar imformation.

The workflow innvolves the genome alignment of input reads using 
[minimap2](https://github.com/lh3/minimap2) and the analaysis 
of read-target overlap with [bedtools](https://github.com/arq5x/bedtools2).


