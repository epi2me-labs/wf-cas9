#!/usr/bin/env bash

aln="/Users/Neil.Horner/work/testing/cas9/test_data/fastq_pass.bed"
targets="/Users/Neil.Horner/work/workflow_outputs/cas9/targets.bed"
windows="/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/windows"

OUTDIR="/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/outtest"
chr_sizes="/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/chrom_small.txt"

# Question1
# Depth of coverag

# Make tiles bed
bedtools makewindows -g $chr_sizes -w 100 -i 'srcwinnum' > $OUTDIR/windows.bed

# Bed file for mapping tile to target
bedtools intersect -a $OUTDIR/windows.bed -b $targets -wb > $OUTDIR/tiles_int_targets.bed

# Get alignment coverage at tiles per strand
#header="chr_sizes start end target coverage #_bases_covered tile_size fracTileCovered\n"
#printf $header > $OUTDIR/positive_target_cov.bed
cat $aln | grep "\+\$" | bedtools coverage -a $OUTDIR/tiles_int_targets.bed -b - -wa -wb  \
> $OUTDIR/positive_target_cov.bed
#$OUTDIR/positive_target_cov.bed

cat $aln | grep "\-$" | bedtools coverage -a $OUTDIR/tiles_int_targets.bed -b - -wa | \
cut -f 1,2,3,8,9,10,11,12 > $OUTDIR/negative_target_cov.bed


#cat $aln | grep "\-$" | bedtools coverage -a $OUTDIR/tiles_int_targets.bed -b - > $OUTDIR/negative_target_cov.bed


