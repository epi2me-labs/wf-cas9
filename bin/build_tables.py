#!/usr/bin/env python
"""build summary tables for report and output CSVs."""

import argparse
import pandas as pd
from natsort import natsort_keygen


def main(target_summary, on_off_bed, aln_sum):
    """Entry point."""
    header = ['chr', 'start', 'end', 'target', 'nreads', 'nbases',
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

    reads_per_sample = stats_df.groupby('sample_id')['read_length'].count()
    s_kbases = df_on_off.groupby('sample_id')['read_length'].sum() / 1000
    s_mean_read_len = stats_df.groupby('sample_id')['read_length'].mean()
    s_acc = stats_df.groupby('sample_id')['acc'].mean()

    # Merge the target name onto the stats to get reads that map to a target
    stats_df = stats_df.merge(df_on_off[pd.notna(df_on_off.target)][['name', 'target']], left_on='name', right_on='name', how='right')
    stats_df['p'] = stats_df.direction.map(lambda x: 1 if x == '+' else 0)
    stats_df['n'] = stats_df.direction.map(lambda x: 1 if x == '-' else 0)
    s_bias = stats_df.groupby(['sample_id']).sum()[['p', 'n']]
    s_bias['bias'] = (s_bias.p - s_bias.n) / (s_bias.p + s_bias.n)

    # Build the summary cs
    sample_summary = pd.DataFrame({
        'reads': reads_per_sample,
        'kbases': s_kbases,
        'mean_read_length': s_mean_read_len,
        'mean_accuracy': s_acc,
        'strand_bias': s_bias['bias']
    })
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
    # from pathlib import Path
    # dir_ = Path('/Volumes/Groups/custflow/active/nhorner/wf-cas9-stuff/work/e8/35adbc9e77237d9d9d63a317cb8980')
    # on_off = dir_ / 'on_off'
    # summary = dir_ / 'aln_stats'
    # table = dir_ / 'target_summary'
    # main(table, on_off, summary)