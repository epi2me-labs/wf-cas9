from pathlib import Path
import sys
from pybedtools import BedTool
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


targets_file="/Users/Neil.Horner/work/workflow_outputs/cas9/targets.bed"
aln_file="/Users/Neil.Horner/work/testing/cas9/report_test_data/fastq_pass_m4.sam"

# create bam: samtools view -b fastq_pass.sam > fastq_pass.bam
# bam to bed: bedtools bamtobed -i fastq_pass.bam | bedtools sort > fastq_pass.bed
sorted_bed = "/Users/Neil.Horner/work/testing/cas9/report_test_data/fastq_pass.bed"

targets = BedTool(targets_file)
df_targets = targets.to_dataframe()
aln = BedTool(sorted_bed)


# How much coverage to be counted as an overlap (default is 1bp)
on_target_depth = targets.coverage(aln).to_dataframe().sort_index()

on_target_depth.rename(columns={
    'score': 'num_overlaps',
    'strand': 'target_bases_with_aln',
    'thickStart': 'target_len',
    'thickEnd': 'fraction_target_with_aln'
    }, inplace=True)

# Will probably not use this table. It's all on the third tutorial table
on_target_depth.to_csv('output/on_target_depth.csv')


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


df_on_off.to_csv('output/on_off_targets.csv')


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
aln_df = aln.to_dataframe()
dfs = []
tile_size = 100
for (chrom, strand), df in aln_df.groupby(['chrom', 'strand']):
    starts = list(range(df.start.min(), df.end.max(), tile_size))
    df = pd.DataFrame.from_dict({'start': starts})
    df['end'] = df.start + tile_size - 1
    df['chrom'] = chrom
    df['strand'] = strand
    df = df[['chrom', 'start', 'end', 'strand']]
    dfs.append(df)
    break
tiles = BedTool().from_dataframe(pd.concat(dfs))

target_tiles = targets.intersect(tiles).to_dataframe().groupby('name')

for target, df in target_tiles:
    t = BedTool().from_dataframe(df)
    reads_int_tiles = t.intersect(aln)
    df_reads_target = reads_int_tiles.to_dataframe().groupby('start').count()[['chrom']]
    df_reads_target.rename(columns={'chrom': 'overlaps'}, inplace=True)
    sns.lineplot(df_reads_target.index, df_reads_target.overlaps)
    plt.show()
    print
print('p')

    # reads_int_tiles = target_tiles.intersect(aln)
    # df_reads_target = reads_int_tiles.to_dataframe()
    # df_reads_target.rename(columns={'score': 'strand'}, inplace=True)
    #
    # for strand, dfstr in df_reads_target.groupby(['strand']):
    #     dfc = dfstr.groupby('start').count()[['chrom']]
    #     sns.lineplot(dfc.index, dfc.chrom)
    #     plt.show()
    #     print


# hits = []
# for target in targets:
#     b = BedTool(str(target), from_string=True)
#     target_hits = aln.intersect(b).to_dataframe()
#     target_hits['tname'] = target.name
#     hits.append(target_hits)







