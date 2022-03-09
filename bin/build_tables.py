#!/usr/bin/env python
"""build summary tables for report and output CSVs."""

import argparse

from natsort import natsort_keygen
import pandas as pd


def main(target_summary, on_off_bed, aln_sum):
    """Entry point."""
    header = [
        'chr', 'start', 'end', 'target', 'nreads', 'nbases',
        'tsize', 'coverage_frac', 'median_cov', 'p', 'n', 'sample_id']

    frames = []

    df_ono_ff = pd.read_csv(
        on_off_bed, sep='\t',
        names=['chr', 'start', 'end', 'read_id', 'target', 'sample_id'],
        index_col=False)

    stats_df = pd.read_csv(aln_sum, sep='\t', index_col=False)
    stats_df.rename(columns={stats_df.columns[-1]: 'sample_id'}, inplace=True)
    df_on_off = df_ono_ff.merge(
        stats_df[['name', 'read_length', 'acc']],
        left_on='read_id', right_on='name')

    main_df = pd.read_csv(
        target_summary, sep='\t', names=header, index_col=False)

    for id_, df in main_df.groupby('sample_id'):
        df = df.drop(columns=['sample_id'])
        if len(df) == 0:
            continue
        df_on_off = df_on_off.astype({
            'start': int,
            'end': int,
            'read_length': int,
            'acc': float
        })
        read_len = df_on_off.groupby(['target']).mean()[['read_length']]
        read_len.columns = ['mean_read_length']
        if len(read_len) > 0:
            df = df.merge(read_len, left_on='target', right_index=True)
        else:
            df['mean_read_length'] = 0

        kbases = df_on_off.groupby(['target']).sum()[['read_length']] / 1000
        kbases.columns = ['kbases']
        if len(kbases) > 0:
            df = df.merge(kbases, left_on='target', right_index=True)
        else:
            df['kbases'] = 0

        acc = df_on_off.groupby(['target']).mean()[['acc']]
        acc.columns = ['mean_acc']
        df = df.merge(acc, left_on='target', right_index=True)

        df['strand_bias'] = (df.p - df.n) / (df.p + df.n)
        df.drop(columns=['p', 'n'], inplace=True)
        df.insert(0, 'sample', id_)
        frames.append(df)

    if len(frames) > 0:
        df_all = pd.concat(frames)
        df_all = df_all.astype({
            'start': int,
            'end': int,
            'nreads': int,
            'nbases': int,
            'tsize': int
        })

        df_all = df_all.round({
            'strand_bias': 2,
            'coverage_frac': 2,
            'kbases': 2,
            'mean_read_length': 1,
            'mean_acc': 2})

        df_all = df_all[[
            'sample', 'chr', 'start', 'end', 'target', 'tsize',
            'kbases', 'coverage_frac', 'median_cov', 'nreads',
            'mean_read_length', 'mean_acc', 'strand_bias']]
        df_all.sort_values(
            by=["sample", "chr", "start"], key=natsort_keygen(), inplace=True)
    else:
        df_all = pd.DataFrame()
    df_all.to_csv('target_summary.csv', index=False)

    # Get the sample summary
    stats_df = stats_df.astype({
        'read_length': int,
        'acc': float})

    # reads_per_sample = df_all.groupby('sample').sum()['nreads']
    gb = df_all.groupby('sample')
    dfs = []
    for sid, df in gb:
        df['s_kbases'] = df['kbases'] * (
            df['nreads'] / df['nreads'].sum())
        df['s_mean_read_length'] = df['mean_read_length'] * (
                    df['nreads'] / df['nreads'].sum())
        df['s_mean_acc'] = df['mean_acc'] * (
                    df['nreads'] / df['nreads'].sum())
        df['s_strand_bias'] = df['strand_bias'] * (
                    df['nreads'] / df['nreads'].sum())
        sample_df = df[['s_kbases', 's_mean_read_length', 's_mean_acc',
                      's_strand_bias']]
        sample_df = sample_df.sum()
        sample_df['sample_id'] = sid
        dfs.append(sample_df)

    sample_summary = pd.concat(dfs, axis=1).T
    sample_summary.set_index('sample_id', drop=True, inplace=True)

    sample_summary.to_csv('sample_summary.csv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--target_summary", help="Target summary bed.")
    parser.add_argument(
        "--aln_summary", help="Alignment summary from pomoxis/stats_from_bam.")
    parser.add_argument(
        "--on_off", help="bed file of xx .")
    args = parser.parse_args()
    main(args.target_summary, args.on_off, args.aln_summary)
