Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | ./wf-cas9-report.html | Report for all samples | aggregated |
| Sample summary table | ./sample_summary.csv | Summary statistics for each sample. | aggregated |
| Target summary table | ./target_summary.csv | Summary statistics for each sample-target combination. | aggregated |
| Per file read stats | ./{{ alias }}/{{ alias }}_per-read-stats.tsv.gz | A TSV with per-file read statistics, including all samples. | per-sample |
| On-target BAM | .//{{ alias }}/{{ alias }}_on_target.bam | BAM file containing alignments that map to one of the given targets. | per-sample |
| On-target BED | .//{{ alias }}/{{ alias }}_on_target.bed | BED file containing summarising alignments that map to one of the given targets. | per-sample |
| On-target FASTQ | .//{{ alias }}/{{ alias }}_on_target.fastq | FASTQ file containing reads that map to one of the given targets. | per-sample |
