#!/usr/bin/env bash

aln="/Users/Neil.Horner/work/testing/cas9/test_data/fastq_pass.bed"
targets="/Users/Neil.Horner/work/workflow_outputs/cas9/targets.bed"
OUTDIR="/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/outtest"
stats="/Users/Neil.Horner/work/workflow_outputs/cas9/seqstats.csv"
chr_sizes="/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/chrom_small.txt"

# chr, start, stop (target), target, overlaps, covered_bases, len(target), frac_covered
bedtools coverage -a $targets -b $aln > $OUTDIR/target_summary.bed

# Need to add following columns
# - kbases - kbases of coverage - DONE
# - median coverage   DONE
# mean_read_len DONE
# mean_accuracy - Need to output a mapping of read_to_target so can be done in PYTHON DONE
# strand bias - I think we can get this from target_summary process?

# Make tiles bed
#bedtools makewindows -g $chr_sizes -w 100 -i 'srcwinnum' > $OUTDIR/windows.bed

# Bed file for mapping tile to target
bedtools intersect -a $OUTDIR/windows.bed -b $targets -wb > $OUTDIR/tiles_int_targets.bed

# Get alignment coverage at tiles per strand
cat $aln | bedtools coverage -a $OUTDIR/tiles_int_targets.bed -b - > $OUTDIR/target_cov.bed
bedtools groupby -i $OUTDIR/target_cov.bed -g 1 -c 9 -o median > $OUTDIR/median_coverage.bed


# Mean read len
cat $aln | bedtools coverage -a - -b $targets -wb | bedtools groupby -g 1 -c 9 -o mean > $OUTDIR/mean_read_len.bed

# Kbases of coverage
cat $aln | bedtools coverage -a - -b $targets -wb | bedtools groupby -g 1, -c 9 -o sum > $OUTDIR/kbases.bed


# Mapping of read to target to add t summary table in python
cat $aln | bedtools intersect -a $aln -b $targets -wa -wb | \
cut -f 1,4,10  > $OUTDIR/read_target_map.bed


# Strand bias


