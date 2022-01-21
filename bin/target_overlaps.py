#!/usr/bin/env python

import argparse
from pathlib import Path
import sys
from pybedtools import BedTool
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as pl


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("target", help="bed file of single target region",
                        dest='target', required=True)
    parser.add_argument("alignment", help="alignment", dest="aln",
                        required=True)
    parser.add_argument("sample_id", help="alignment", dest="sample_id",
                        required=True)
    args = parser.parse_args()

    targets = BedTool(args.target)
    df_targets = targets.to_dataframe()

    aln = BedTool(args.aln)
    # How much coverage to be counted as an overlap (default is 1bp)
    on_target_depth = targets.coverage(aln).to_dataframe().sort_index()

    on_target_depth.rename(columns={
        'score': 'num_overlaps',
        'strand': 'target_bases_with_aln',
        'thickStart': 'target_len',
        'thickEnd': 'fraction_target_with_aln'
    }, inplace=True)

    # Will probably not use this table. It's all on the third tutorial table
    on_target_depth.to_csv('{}_on_target_depth.csv'.format(args.sample_id))


    ###### Table 2 - I think
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

    df_on_off.to_csv('{}_on_off_targets.csv'.format(args.csv))


    ########## data for plots
    # Tiling operation
    aln_df = aln.to_dataframe()
    dfs = []
    tile_size = 100
    i = 0 # for testing
    for (chrom, strand), df in aln_df.groupby(['chrom', 'strand']):
        starts = list(range(df.start.min(), df.end.max(), tile_size))
        df = pd.DataFrame.from_dict({'start': starts})
        df['end'] = df.start + tile_size - 1
        df['chrom'] = chrom
        df['strand'] = strand
        df = df[['chrom', 'start', 'end', 'strand']]
        dfs.append(df)
        if i == 1:
            break
    tiles = BedTool().from_dataframe(pd.concat(dfs))

    target_tiles = targets.intersect(tiles).to_dataframe().groupby('name')

    for target, df_f in target_tiles:
        tf = BedTool().from_dataframe(df_f)
        reads_int_tiles = tf.intersect(aln)
        df_reads_target = \
        reads_int_tiles.to_dataframe().groupby('start').count()[['chrom']]
        df_reads_target.rename(columns={'chrom': 'overlaps'}, inplace=True)
        df_reads_target.to_csv()
        # sns.lineplot(df_reads_target.index, df_reads_target.overlaps)
        # plt.show()
        # print
    print('p')



if __name__ == '__main__':
    target_file = "/Users/Neil.Horner/work/workflow_outputs/cas9/targets.bed"
    aln_file = "/Users/Neil.Horner/work/testing/cas9/report_test_data/fastq_pass.sam"
    sys.argv.extend([target_file, aln_file])
    main()