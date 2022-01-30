#!/usr/bin/env bash

aln="/Users/Neil.Horner/work/testing/cas9/test_data/fastq_pass.bed"
targets="/Users/Neil.Horner/work/workflow_outputs/cas9/targets.bed"
OUTDIR="/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/outtest"
stats="/Users/Neil.Horner/work/workflow_outputs/cas9/seqstats.csv"
chr_sizes="/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/chrom.sizes"

# chr, start, stop (target), target, overlaps, covered_bases, len(target), frac_covered
bedtools coverage -a $targets -b $aln | bedtools sort > $OUTDIR/target_summary_temp.bed

# Need to add following columns
# - kbases - kbases of coverage - DONE
# - median coverage   DONE
# mean_read_len DONE
# mean_accuracy - use aln_tagets.bed and seq stats from fastcat in python
# strand bias - I think we can get this from target_summary process?

# Make tiles bed - TODO: compress as it can get big
#bedtools makewindows -g $chr_sizes -w 100 -i 'srcwinnum' > $OUTDIR/windows.bed

# Bed file for mapping tile to target
bedtools intersect -a $OUTDIR/windows.bed -b $targets -wb > $OUTDIR/tiles_int_targets.bed

# Get alignment coverage at tiles per strand
cat $aln | bedtools coverage -a $OUTDIR/tiles_int_targets.bed -b - > $OUTDIR/target_cov.bed
bedtools groupby -i $OUTDIR/target_cov.bed -g 1 -c 9 -o median | cut -f 2  > $OUTDIR/median_coverage.bed

# Map targets to aln
alntargets=$OUTDIR/aln_tagets.bed
cat $aln | bedtools intersect -a - -b $targets -wb > $alntargets

#cat $alntargets | bedtools coverage -a - -b $targets -wb > $OUTDIR/test.bed
#
#cat $alntargets | bedtools coverage -a - -b $targets -wb |bedtools > test.bed

# Strand bias
cat $alntargets | grep '\W+\W' | bedtools coverage -a - -b $targets -wb | \
bedtools sort | bedtools groupby -g 10 -c 1 -o count    > $OUTDIR/pos.bed

cat $alntargets | grep '\W-\W' | bedtools coverage -a - -b $targets -wb \
| bedtools groupby -g 10 -c 1 -o count | cut -f 2 > $OUTDIR/neg.bed

# Mean read len
cat $alntargets | bedtools coverage -a - -b $targets -wb | bedtools groupby -g 10 -c 13 -o mean >  $OUTDIR/mean_read_len.bed

# Kbases of coverage
cat $alntargets | bedtools coverage -a - -b $targets -wb | bedtools groupby -g 10, -c 12 -o sum | cut -f 2 > $OUTDIR/kbases.bed

paste $OUTDIR/target_summary_temp.bed \
      $OUTDIR/mean_read_len.bed \
      $OUTDIR/kbases.bed \
      $OUTDIR/median_coverage.bed \
      $OUTDIR/pos.bed \
      $OUTDIR/neg.bed > $OUTDIR/target_summary.bed


# Strand bias


