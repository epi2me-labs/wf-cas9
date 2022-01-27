#!/usr/bin/env python

import argparse
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
    stats = pd.read_csv(args.seq_stats)
    # create bam: samtools view -b fastq_pass.sam > fastq_pass.bam
    # bam to bed: bedtools bamtobed -i fastq_pass.bam | bedtools sort > fastq_pass.bed
    # sorted_bed = "/Users/Neil.Horner/work/testing/cas9/test_data/fastq_pass.bed"

    # Get
    targets = BedTool(args.targets)
    df_targets = targets.to_dataframe()
    aln = BedTool(args.alignment_bed)
    aln = aln.intersect(targets, wo=True)
    aln_df = aln.to_dataframe()

    # Note: may have to change this if we are looking for background genome-wide
    chr_used = np.unique(df_targets.chrom)
    sizes = {k: v for (k, v) in zip(ref.references, ref.lengths) if
             k in chr_used}

    # How much coverage to be counted as an overlap (default is 1bp)
    on_target_depth = targets.coverage(aln).to_dataframe().sort_index()
    on_target_depth.rename(columns={
        'score': 'num_overlaps',
        'strand': 'target_bases_with_aln',
        'thickStart': 'target_len',
        'thickEnd': 'fraction_target_with_aln'
    }, inplace=True)

    # Will probably not use this table. It's all on the third tutorial table
    on_target_depth.to_csv(
        '{}_on_target_depth.csv'.format(args.sample_id))

    # Try to map target name to aln
    # aln_test = aln.intersect(targets, wb=True).to_dataframe().sort_index()
    ###

    header = ['chrom', 'start', 'stop', '1', '2', 'strand', '3', 't_start', 't_end',
              'target', 'overlap_bases', '6', '7', 'read_len', 'frac_overlap']
    on_off = aln.coverage(targets).to_dataframe(names=header).sort_index()

    on_off.drop(columns=[x for x in header if x.isnumeric()], inplace=True)
    on_off.rename(columns={
        'blockCount': 'frac_overlap',
        'thickEnd': 'target_overlap',
        'itemRgb': 'read_len'},
        inplace=True)

    on = on_off[on_off['frac_overlap'] > 0]
    off = on_off[on_off['frac_overlap'] == 0]

    # Move this up?
    # Map target names back the coverage df
    # on['target_name'] = on.apply(lambda row: )
    # median_cov = pd.DataFrame(on.groupby(['tile_start', 'target']).count())
    # for k in median_cov:
    #     print(k)

    df_on_off = pd.DataFrame(
        [[len(on), on.read_len.sum() / 1000, on.read_len.mean()],
         [len(off), off.read_len.sum() / 1000, off.read_len.mean()],
         [len(on_off), on_off.read_len.sum() / 1000, on_off.read_len.mean()]],
        index=['Reads', 'KBs', 'Mean_read_length'],
        columns=[['On-target', 'Non-target', 'All']]).T

    df_on_off.to_csv(
        '{}_coverage_summary.csv'.format(args.sample_id))
    f = on_off[on_off['strand'] == '-'].groupby(
                        ['target']).count()['chrom']
    r = on_off[on_off['strand'] == '+'].groupby(
                        ['target']).count()['chrom']
    bias = (f - r) / (f + r)
    mean_read_len = on_off.groupby(['target']).mean()['read_len']


    # Plots
    # Intersect loses strand information, so do intersection on each strand
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

    df_all_target_tiles = targets.intersect(tiles_bed).to_dataframe().groupby(
        'name')

    # I'm assuming I can do this in pybedtools. But just use pandas for now
    fwd_aln = BedTool().from_dataframe(aln_df[aln_df.strand == '+'])
    rev_aln = BedTool().from_dataframe(aln_df[aln_df.strand == '-'])


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

    # target summaries table

    df['unstranded_overlaps'] = df.overlaps_f + df.overlaps_r
    median_coverage = df.groupby(['target']).median()['unstranded_overlaps']
    kbases = df.groupby(['target']).sum()[['total']]
    # mean score
    #strand_bias = bias
    # mean read length -> mean_read_len

    print(p)



if __name__ == '__main__':
    target_file = "/Users/Neil.Horner/work/workflow_outputs/cas9/targets.bed"
    aln_bed = "/Users/Neil.Horner/work/testing/cas9/test_data/fastq_pass.bed"
    genome_file = "/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/grch38.fasta.gz"
    stats_file = "/Users/Neil.Horner/work/workflow_outputs/cas9/seqstats.csv"
    sys.argv.extend([target_file, aln_bed, genome_file, stats_file, 'test'])
    main()