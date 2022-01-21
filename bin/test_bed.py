from pathlib import Path
import sys
from pybedtools import BedTool
import pandas as pd
import pysam
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# TODO
# - Should the alignemt file be limited to primary alignments. I guess it should

targets_file="/Users/Neil.Horner/work/testing/cas9/test_data/targets.bed"
aln_file="/Users/Neil.Horner/work/testing/cas9/test_data/fastq_pass_m4.sam"
# aln_file="/Users/Neil.Horner/work/testing/cas9/test_data/fastq_pass_m4_10percent.sam"  # 0/1 subsampled file for testing
ref_genome = "/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/grch38.fasta.gz"
ref = pysam.FastaFile(ref_genome)


# create bam: samtools view -b fastq_pass.sam > fastq_pass.bam
# bam to bed: bedtools bamtobed -i fastq_pass.bam | bedtools sort > fastq_pass.bed
sorted_bed = "/Users/Neil.Horner/work/testing/cas9/test_data/fastq_pass.bed"

targets = BedTool(targets_file)
df_targets = targets.to_dataframe()
aln = BedTool(sorted_bed)
aln_df = aln.to_dataframe()

# Note: may have to change this if we are looking for background genome-wide
chr_used = np.unique(df_targets.chrom)
sizes = {k: v for (k, v) in zip(ref.references, ref.lengths) if k in chr_used  }


# How much coverage to be counted as an overlap (default is 1bp)
on_target_depth = targets.coverage(aln).to_dataframe().sort_index()

on_target_depth.rename(columns={
    'score': 'num_overlaps',
    'strand': 'target_bases_with_aln',
    'thickStart': 'target_len',
    'thickEnd': 'fraction_target_with_aln'
    }, inplace=True)

# Will probably not use this table. It's all on the third tutorial table
on_target_depth.to_csv('/Users/Neil.Horner/work/testing/cas9/output/on_target_depth.csv')


on_off = aln.coverage(targets).to_dataframe().sort_index()

on_off.rename(columns={
    'blockCount': 'frac_overlap',
    'thickEnd': 'target_overlap',
    'itemRgb': 'read_len'},
    inplace=True)

on = on_off[on_off['frac_overlap'] > 0]
off = on_off[on_off['frac_overlap'] == 0]

df_on_off = pd.DataFrame(
    [[len(on), on.read_len.sum() / 1000, on.read_len.mean()],
    [len(off), off.read_len.sum() / 1000, off.read_len.mean()],
    [len(on_off), on_off.read_len.sum() / 1000, on_off.read_len.mean()]],
    index=['Reads', 'KBs', 'Mean_read_length'],
    columns=[['On-target', 'Non-target', 'All']]).T


df_on_off.to_csv('/Users/Neil.Horner/work/testing/cas9/output/coverage_summary.csv')


# Plots

# Intersect loses strand information, so do intersection on each strand
rev = aln.to_dataframe()
rev = rev[rev.strand == '-']
rev_bed = BedTool.from_dataframe(rev)

# for_ = aln.to_dataframe()
# for_ = for_[for_.strand == '+']
# for_bed = BedTool.from_dataframe(for_)
#
# bt = BedTool()
# rev_overlaps = targets.intersect(rev_bed).to_dataframe()
# rev_overlaps['strand'] = '-'
# for_overalps = aln.intersect(targets).to_dataframe()
# for_overalps['strand'] = '+'

# overlaps = pd.concat([rev_overlaps, for_overalps])

# Tiling operation

dfs_fwd = []
dfs_rev  = []
tile_size = 100

# make some tiles
tile_dfs = []
for chrom, size in sizes.items():
    starts = list(range(0, size, tile_size))
    df = pd.DataFrame.from_dict({'start': starts})
    df['end'] = df.start + tile_size - 1
    df['chrom'] = chrom
    df = df[['chrom', 'start', 'end']]
    tile_dfs.append(df)

tiles_bed = BedTool().from_dataframe(pd.concat(tile_dfs))

df_all_target_tiles = targets.intersect(tiles_bed).to_dataframe().groupby('name')

# I'm assuming I can do this in pybedtools. But just use pandas for now
fwd_aln = BedTool().from_dataframe(aln_df[aln_df.strand == '+'])
rev_aln = BedTool().from_dataframe(aln_df[aln_df.strand == '-'])


def get_target_overlaps(df_target_tiles_):
    t = BedTool().from_dataframe(df_target_tiles_)
    fwd_reads_int_tiles = t.intersect(fwd_aln)
    rev_reads_int_tiles = t.intersect(rev_aln)
    # These can be variable length. Zero hits are no included
    df_fwd_tile_cov = \
        fwd_reads_int_tiles.to_dataframe().groupby('start').count()[['chrom']]
    df_fwd_tile_cov.rename(columns={'chrom': 'overlaps'}, inplace=True)
    df_rev_tile_cov = \
        rev_reads_int_tiles.to_dataframe().groupby('start').count()[['chrom']]
    df_rev_tile_cov.rename(columns={'chrom': 'overlaps'}, inplace=True)

    ## Merge back to the tiles so we don't lose uncovered tiles
    f = df_target_tiles_.merge(df_fwd_tile_cov[['overlaps']], left_on='start',
                           right_index=True, how='left')
    r = df_target_tiles_.merge(df_rev_tile_cov[['overlaps']], left_on='start',
                               right_index=True, how='left')
    result = f.merge(r[['start', 'overlaps']], left_on='start', right_on='start',
            suffixes=['_f', '_r']).fillna(0)
    return result


 # in production version we may want to do each target in seperate process
 # for now do all here
results = []


for target, df_target_tiles in df_all_target_tiles:
    df_target_tiles.sort_values(by='start', inplace=True)
    result = get_target_overlaps(df_target_tiles)
    results.append(result)


result_df = pd.concat(results).reset_index(drop=True)
result_df.rename(columns={'name': 'target'}, inplace=True)
result_df.to_csv("/Users/Neil.Horner/work/testing/cas9/output/target_coverage.csv")








