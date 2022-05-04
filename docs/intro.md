## Introduction
This workflow generated a report summarising the results of Cas9 enrichment sequencing.

Users provide a reference genome, fastq ONT reads, and a bed file containing enrichment regions. 
The reads are first mapped the reference genome using [minimap2](https://github.com/lh3/minimap2), and 
various plots and tables are generated summarizing the enrichment results.
