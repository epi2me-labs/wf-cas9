#!/usr/bin/env python

import argparse
from functools import reduce
from pathlib import Path
import math
import sys
from pybedtools import BedTool
import pandas as pd
import numpy as np
import pysam
import seaborn as sns
import matplotlib.pyplot as pl


def get_target_overlaps(df_target_tiles_, fwd_aln, rev_aln):
    """Get dataframe of read coverage across a target.

    :returns
        pd.Dataframe for a target with columns
            chrom
            start: start of tle
            end: end of tile
            name: target name
            overlaps_f: num positive strand overlaps
            overlaps_r: num negative strand overlaps
    """
    t = BedTool().from_dataframe(df_target_tiles_)
    fwd_reads_int_tiles = t.intersect(fwd_aln)
    rev_reads_int_tiles = t.intersect(rev_aln)
    # These can be variable length. Zero hits are no included
    df_fwd_tile_cov = \
        fwd_reads_int_tiles.to_dataframe().groupby('start').count()[
            ['chrom']]
    df_fwd_tile_cov.rename(columns={'chrom': 'overlaps'}, inplace=True)
    df_rev_tile_cov = \
        rev_reads_int_tiles.to_dataframe().groupby('start').count()[
            ['chrom']]
    df_rev_tile_cov.rename(columns={'chrom': 'overlaps'}, inplace=True)

    ## Merge back to the tiles so we don't lose uncovered tiles
    f = df_target_tiles_.merge(df_fwd_tile_cov[['overlaps']],
                               left_on='start',
                               right_index=True, how='left')
    r = df_target_tiles_.merge(df_rev_tile_cov[['overlaps']],
                               left_on='start',
                               right_index=True, how='left')
    result = f.merge(r[['start', 'overlaps']], left_on='start',
                     right_on='start',
                     suffixes=['_f', '_r']).fillna(0)
    return result


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("targets", help="bed file of single target region")
    parser.add_argument("alignment_bed", help="alignment")
    parser.add_argument("ref_genome", help='Reference genome fasta')
    parser.add_argument("sample_id", help="alignment")
    parser.add_argument("seq_stats", help="Sequence summary stats")

    args = parser.parse_args()

    ref = pysam.FastaFile(args.ref_genome)
    stats = pd.read_csv(args.seq_stats, index_col=0, sep='\t')
    # create bam: samtools view -b fastq_pass.sam > fastq_pass.bam
    # bam to bed: bedtools bamtobed -i fastq_pass.bam | bedtools sort > fastq_pass.bed
    # sorted_bed = "/Users/Neil.Horner/work/testing/cas9/test_data/fastq_pass.bed"

    # Get
    targets = BedTool(args.targets)
    df_targets = targets.to_dataframe()
    aln = BedTool(args.alignment_bed)
    aln_cov = aln.coverage(targets)
    aln_cov_df = aln_cov.to_dataframe().sort_index()

    # Map target names onto aln
    target_read_map = targets.intersect(aln, wo=True).to_dataframe()
    target_read_map = target_read_map[['name', 'thickEnd']]
    target_read_map.rename(columns={'name':'target', 'thickEnd': 'read_id'}, inplace=True)
    # aln_off =

    # Note: may have to change this if we are looking for background genome-wide
    chr_used = np.unique(df_targets.chrom)
    sizes = {k: v for (k, v) in zip(ref.references, ref.lengths) if
             k in chr_used}

    # How much coverage to be counted as an overlap (default is 1bp)

    aln_cov_df.rename(columns={
        'name': 'seq_id',
        'blockCount': 'frac_overlap',
        'start': 'read_start',
        'end': 'read_end',
        'itemRgb': 'read_len',
        'thickStart': 'has_overlap',
        'thickEnd': 'bases_aligning'
    }, inplace=True)

    # Will probably not use this table. It's all on the third tutorial table
    # on_target_depth.to_csv(
    #     '{}_on_target_depth.csv'.format(args.sample_id))

    # Try to map target name to aln
    # aln_test = aln.intersect(targets, wb=True).to_dataframe().sort_index()
    ###

    # header = ['chrom', 'start', 'stop', 'seq_id', '2', 'strand', '3', 't_start', 't_end',
    #           'target', 'overlap_bases', '6', '7', 'read_len', 'frac_overlap']
    aln_cov_df = aln_cov_df.merge(stats[['mean_quality']], left_on='seq_id', right_on='read_id')

    # on_off.drop(columns=[x for x in header if x.isnumeric()], inplace=True)
    # on_off.rename(columns={
    #     'blockCount': 'frac_overlap',
    #     'thickEnd': 'target_overlap',
    #     'itemRgb': 'read_len'},
    #     inplace=True)

    on = aln_cov_df[aln_cov_df['frac_overlap'] > 0]
    on = on.merge(target_read_map, left_on='seq_id', right_on='read_id', how='left')
    # target_read_overlaps = targets.coverage(BedTool(on))
    off = aln_cov_df[aln_cov_df['frac_overlap'] == 0]

    # Move this up?
    # Map target names back the coverage df
    # on['target_name'] = on.apply(lambda row: )
    # median_cov = pd.DataFrame(on.groupby(['tile_start', 'target']).count())
    # for k in median_cov:
    #     print(k)

    df_on_off = pd.DataFrame(
        [[len(on), on.read_len.sum() / 1000, on.read_len.mean()],
         [len(off), off.read_len.sum() / 1000, off.read_len.mean()],
         [len(aln_cov_df), aln_cov_df.read_len.sum() / 1000, aln_cov_df.read_len.mean()]],
        columns=['Reads', 'KBs', 'Mean_read_length'],
        index=[['On-target', 'Non-target', 'All']])

    df_on_off.to_csv(
        '{}_coverage_summary.csv'.format(args.sample_id))

    # Plots
    # Intersect loses strand information, so do intersection on each strand
    tile_size = 100

    # make some tiles
    tile_dfs = []
    for chrom, size in sizes.items():
        if chrom != 'chr1':
            continue
        starts = list(range(0, size, tile_size))
        df = pd.DataFrame.from_dict({'start': starts})
        df['end'] = df.start + tile_size - 1
        df['chrom'] = chrom
        df = df[['chrom', 'start', 'end']]
        tile_dfs.append(df)

    tiles_bed = BedTool().from_dataframe(pd.concat(tile_dfs))

    df_all_target_tiles = targets.intersect(tiles_bed).to_dataframe().groupby(
        'name')

    # I'm assuming I can do this in pybedtools. But just use pandas for now
    fwd_aln = BedTool().from_dataframe(aln_cov_df[aln_cov_df.strand == '+'])
    rev_aln = BedTool().from_dataframe(aln_cov_df[aln_cov_df.strand == '-'])


    # in production version we may want to do each target in seperate process
    # for now do all here
    results = []

    for target, df_target_tiles in df_all_target_tiles:
        df_target_tiles.sort_values(by='start', inplace=True)
        result = get_target_overlaps(df_target_tiles, fwd_aln, rev_aln)
        results.append(result)

    result_df = pd.concat(results).reset_index(drop=True)
    result_df.rename(columns={'name': 'target'}, inplace=True)
    result_df.to_csv(
        "{}_target_coverage.csv".format(args.sample_id))
    on_target_depth = targets.coverage(aln).to_dataframe().sort_index()
    on_target_depth = on_target_depth.rename(columns={
        'name': 'target',
        'thickStart': 'tSize',
        'strand':  'basesCovererd',
        'thickEnd': 'fracTAln'}).set_index('target', drop=True)
    on_target_depth.drop(columns=['score'], inplace=True)
    # target summaries table
    f = on[on['strand'] == '+'].groupby(
        ['target']).count()['chrom']
    r = on[on['strand'] == '-'].groupby(
        ['target']).count()['chrom']
    bias = (f - r) / (f + r)
    bias.columns = ['strand_bias']
    mean_read_len = on.groupby(['target']).mean()[['read_len']]
    result_df['all_overlaps'] = result_df.overlaps_f + result_df.overlaps_r
    median_coverage = result_df.groupby(['target']).median()[['all_overlaps']]
    median_coverage.columns = ['median_coverage']
    # on_target_depth['kbases'] = on.bes + result_df.overlaps_r
    kbases = on.groupby(['target']).sum()[['bases_aligning']] / 1000
    kbases.columns = ['kbases']
    mean_quality = on.groupby(['target']).mean()[['mean_quality']]

    frames = [on_target_depth, kbases,  median_coverage, mean_quality,
              mean_read_len, bias]
    on_target_depth = reduce(lambda left, right: pd.merge(
        left, right, on='target', how='left'), frames)
    on_target_depth = on_target_depth.rename(
        columns={'read_len': 'mean_read_len', 'chrom_x': 'chrom'})\
        .sort_values(by=['chrom', 'start'])\
        .drop(columns=['chrom_y'])
    on_target_depth.reset_index(inplace=True, drop=False)

    on_target_depth.to_csv('target_summary.csv')

if __name__ == '__main__':
    target_file = "/Users/Neil.Horner/work/workflow_outputs/cas9/targets.bed"
    aln_bed = "/Users/Neil.Horner/work/testing/cas9/test_data/fastq_pass.bed"
    genome_file = "/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/grch38.fasta.gz"
    stats_file = "/Users/Neil.Horner/work/workflow_outputs/cas9/seqstats.csv"
    sys.argv.extend([target_file, aln_bed, genome_file, 'test', stats_file])
    main()