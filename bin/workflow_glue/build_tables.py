#!/usr/bin/env python
"""build summary tables for report and output CSVs."""

from natsort import natsort_keygen
import pandas as pd
from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("build_tables")
    parser.add_argument(
        "--target_summary", help="Target summary bed.")
    parser.add_argument(
        "--aln_summary", help="Alignment summary from bamstats.")
    parser.add_argument(
        "--read_to_target", help="bed file including read_id, target, sample_id")
    return parser


def read_target_summary(df):
    """Build table summarising on target/off-target read status."""
    df.loc[df.target != 'off_target', 'target'] = 'on_target'

    agg = df.groupby(['sample_id', 'target']).agg(
        mean_len=('read_length', 'mean'),
        num_reads=('read_id', 'count'),
        kbases_mapped=('read_length', 'sum'))

    agg.kbases_mapped /= 1000
    agg = agg.astype('int')
    agg.reset_index(inplace=True)
    result = agg.pivot(
        index='sample_id',
        columns=['target'],
        values=['mean_len', 'num_reads', 'kbases_mapped'])

    # Create emptpy columns if there are no on-target reads
    if ('num_reads', 'on_target') not in result:
        result[[('mean_len', 'on_target')]] = 0
        result[[('num_reads', 'on_target')]] = 0
        result[[('kbases_mapped', 'on_target')]] = 0
    return result


def main(args):
    """Entry point."""
    header = [
        'chr', 'start', 'end', 'target', 'nreads', 'nbases',
        'tsize', 'coverage_frac', 'median_cov', 'p', 'n', 'sample_id', 'run_id']

    frames = []

    df_read_to_target = pd.read_csv(
        args.read_to_target, sep='\t',
        names=['chr', 'start', 'end', 'read_id', 'target', 'sample_id'],
        index_col=False)

    read_stats_df = pd.read_csv(args.aln_summary, sep='\t', index_col=False)

    df_read_to_target = df_read_to_target.merge(
        read_stats_df[['sample_id', 'name', 'read_length']],
        left_on=['sample_id', 'read_id'], right_on=['sample_id', 'name'])

    df_read_to_target['align_len'] = (
            df_read_to_target['end'] - df_read_to_target['start']
    )

    df_target_summary = pd.read_csv(
        args.target_summary, sep='\t', names=header, index_col=False)

    for sample_id, df in df_target_summary.groupby('sample_id'):
        df = df.drop(columns=['sample_id'])
        if len(df) == 0:
            continue

        read_len = (
            df_read_to_target.loc[df_read_to_target.sample_id == sample_id]
            [['target', 'read_length']]
            .groupby(['target'])
            .agg(mean_read_length=('read_length', 'mean'))
        )

        if len(read_len) > 0:
            df = df.merge(read_len, left_on='target', right_index=True)
        else:
            df['mean_read_length'] = 0

        # Kbases is the approximate number of bases mapping to a target.
        # Deletions and insertions within the reads will mean the actual value may
        # vary slightly
        kbases = (
            df_read_to_target[
                df_read_to_target.sample_id == sample_id][['target', 'align_len']]
            .groupby(['target']).sum() / 1000
        )
        kbases.columns = ['kbases']
        if len(kbases) > 0:
            df = df.merge(kbases, left_on='target', right_index=True)
        else:
            df['kbases'] = 0

        df['strand_bias'] = (df.p - df.n) / (df.p + df.n)
        df.drop(columns=['p', 'n'], inplace=True)
        df.insert(0, 'sample', sample_id)
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
            'mean_read_length': 1})

        df_all = df_all[[
            'sample', 'run_id', 'chr', 'start', 'end', 'target', 'tsize',
            'kbases', 'coverage_frac', 'median_cov', 'nreads',
            'mean_read_length', 'strand_bias']]
        df_all.sort_values(
            by=["sample", "chr", "start"], key=natsort_keygen(), inplace=True)
    else:
        df_all = pd.DataFrame()
    df_all.to_csv('target_summary.csv', index=False)

    # Convert some of the target summary data to sample summary data.
    gb = df_all.groupby('sample')
    dfs = []
    for sid, df in gb:
        sdf = pd.DataFrame()
        sdf['kbases'] = df['kbases']
        sdf['mean_read_length'] =\
            df['mean_read_length'] * (df['nreads'] / df['nreads'].sum())
        sdf['strand_bias'] =\
            df['strand_bias'] * (df['nreads'] / df['nreads'].sum())
        sample_df = sdf.sum()
        sample_df['nreads'] = df['nreads'].sum()
        sample_df = sample_df.round(2)
        sample_df['sample_id'] = sid
        # the `run_id` column should contain the same value for all rows for this
        # sample
        sample_df["run_id"], = df["run_id"].unique()
        dfs.append(sample_df)

    if dfs:
        sample_summary = pd.concat(dfs, axis=1).T
        sample_summary.set_index('sample_id', drop=True, inplace=True)
        # move the `run_id` column to the beginning of the dataframe
        sample_summary.insert(0, "run_id", sample_summary.pop("run_id"))
    else:
        sample_summary = pd.DataFrame()

    sample_summary.to_csv('sample_summary.csv')

    read_target_summary_table = read_target_summary(df_read_to_target)
    read_target_summary_table.to_csv('read_target_summary.tsv', sep='\t')
