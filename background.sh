#!/usr/bin/env zsh

aln="/Users/Neil.Horner/work/testing/cas9/test_data/fastq_pass.bed"
targets="/Users/Neil.Horner/work/workflow_outputs/cas9/targets.bed"
OUTDIR="/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/outtest/background"
stats="/Users/Neil.Horner/work/workflow_outputs/cas9/seqstats.csv"
chr_sizes="/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/chrom.sizes"


# Make tiles bed - TODO: compress as it can get big
#bedtools makewindows -g $chr_sizes -w 100 -i 'srcwinnum' > $OUTDIR/windows.bed

bedtools slop -i $targets -g $chr_sizes -b 1000 | \
# remove reads that overlap slopped targets
bedtools intersect -v -a $aln -b - -wa | \
bedtools coverage -a $OUTDIR/windows.bed -b - > $OUTDIR/tiles_background_cov.bed




