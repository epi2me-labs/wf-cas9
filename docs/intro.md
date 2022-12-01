## Introduction
This workflow generated a report summarising the results of Cas9 enrichment sequencing.

Users provide a reference genome, FASTQ ONT reads, and a BED file containing enrichment regions.
The reads are first mapped to the reference genome using [minimap2](https://github.com/lh3/minimap2), and
various plots and tables are generated summarising the enrichment results.
