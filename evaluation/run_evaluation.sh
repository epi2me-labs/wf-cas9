#!/usr/bin/env bash

site="https://ont-exd-int-s3-euwst1-epi2me-labs.s3-eu-west-1.amazonaws.com"


if [[ "$#" -lt 1 ]]; then
    echo "usage: run_evaluation_.sh <outdir> [nextflow.config]"
    exit 1
fi

if [[ "$#" -eq 1 ]]; then
    CONFIG="-c ../nextflow.config"
fi

if [[ "$#" -eq 2 ]]; then
  CONFIG="-c $2";
fi


PROJDIR=$1
#INPUT="$PROJDIR/sample_cas9/fastq_pass"
INPUT="$PROJDIR/multisample_cas9"

if [[ ! -d "$PROJDIR/sample_cas9" ]];
then
# The reads
  wget -P $PROJDIR "$site/cas9_tutorial/sample_cas9.tar.gz"
  tar -xzvf "$PROJDIR/sample_cas9.tar.gz"
#  # Cas9 targets
  wget -P $PROJDIR $site/cas9_tutorial/targets.bed
#  # Reference genome
  wget -P $PROJDIR $site/grch38.tar.gz
  tar -xzvf $PROJDIR/grch38.tar.gz -C $PROJDIR
fi

nextflow run ../ --fastq $INPUT --ref_genome $PROJDIR/grch38/grch38.fasta.gz \
 --targets $PROJDIR/targets.bed --out_dir $PROJDIR/output \
 -w $PROJDIR/workspace -profile standard -resume $CONFIG