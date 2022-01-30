#!/usr/bin/env bash

aln="/Users/Neil.Horner/work/testing/cas9/test_data/fastq_pass.bed"
targets="/Users/Neil.Horner/work/workflow_outputs/cas9/targets.bed"
OUTDIR="/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/outtest"
stats="/Users/Neil.Horner/work/workflow_outputs/cas9/seqstats.csv"

bedtools coverage -a targets -b aln > $OUTDIR/target_summary.bed