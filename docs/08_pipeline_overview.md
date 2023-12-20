<!---High level numbered list of main steps of the workflow and hyperlink to any tools used. If multiple workflows/different modes perhaps have subheadings and numbered steps. Use nested numbering or bullets where required.--->
### 1. Concatenate input files and generate per read stats.

The [fastcat/bamstats](https://github.com/epi2me-labs/fastcat) tool is used to concatenate multifile samples to be processed by the workflow.
### 2. Align read to the reference genome.
The reads are then aligned to the reference genome supplied by the user using [minimap2](https://github.com/lh3/minimap2).

### 3. Generate target coverage data.
This stage of the workflow generates coverage data for each of targets that are used to make
the plots and tables in the report.

First, the reference genome is split into consecutive windows of 100bp. The coverage statistics are 
calcaulted across these windows.

The user must supply a tab-delinted BED file (`targets`) detailing the genomic locations of the targets of interest.
Here is an example file describing two targets.

```
chr19	13204400	13211100	SCA6
chr22	45791500	45799400	SCA10
```
The columns are:
+ chromosome
+ start position
+ end position
+ target name

*Note*: the file does not contain column names. 

With the target locations defined, [bedtools](https://bedtools.readthedocs.io/en/latest/) is used to generate
information regarding the alignment coverage at each of the targets and also of 
background reads (see section 4) per strand. 

### 4. Background identification
In a sequencing enrichment experiment, it can useful to know if there are any off-target genomic regions (regions that are not defined in the `targets` file) that are being preferentially encountered. This information can be used to inform primer design.

We define off-target regions here as any region not within 1kb of a defined target.
Hot-spots are further defined as contiguous regions of off-target alignemnts containing at least 10 reads. 

