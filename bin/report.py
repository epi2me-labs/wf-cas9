#!/usr/bin/env python
"""Create workflow report."""

import argparse
from pathlib import Path
from typing import List

from aplanat import hist, lines
from aplanat.components import fastcat
from aplanat.components import simple as scomponents
from aplanat.report import WFReport
from bokeh.layouts import gridplot
from bokeh.models import Legend, Panel, Tabs
from natsort import natsort_keygen, natsorted
import pandas as pd


def plot_target_coverage(report: WFReport, target_coverages: Path):
    """Make coverage plots of each target.

    Detailing positive and negative strand coverage.

    :param target_coverages
    """
    section = report.add_section()
    section.markdown('''
    ### Target coverage

    Each of the following plots show the amount of coverage, per strand
    in discretized bins of 100 bp.
    ''')

    header = ["chr", "start", "end", "target", "coverage_f", 'coverage_r',
              "sample_id"]
    tabs = []
    main_df = pd.read_csv(target_coverages, names=header, sep='\t')
    for id_, df in main_df.groupby('sample_id'):
        dfg = df.groupby('target')

        ncols = dfg.ngroups if dfg.ngroups < 5 else 4
        plots = []

        for i, (target, df) in enumerate(dfg):
            chrom = df.loc[df.index[0], 'chr']
            ymax = max(df.coverage_f.max(), df.coverage_r.max())
            ylim = [0, ymax * 1.05]  # a bit of space at top of plot

            # if no coverage across target, force some height.
            if not any([df.coverage_f.any(), df.coverage_r.any()]):
                ylim = [0, 10]

            p = lines.line(
                [df.start.values, df.start.values],  # x-values
                [df.coverage_f, df.coverage_r],      # y-values
                title="{}".format(target),
                x_axis_label='{}'.format(chrom),
                y_axis_label='',
                colors=['#1A85FF', '#D41159'],
                ylim=ylim,
                height=200, width=300
            )
            p.xaxis.formatter.use_scientific = False
            p.xaxis.major_label_orientation = 3.14 / 6

            plots.append([chrom, df.start.values[0], p])

        sorted_plots = [p[2] for p in natsorted(plots, key=lambda x: x[0])]

        legend_plot = sorted_plots[ncols - 1]
        legend_plot.width = legend_plot.width + 80
        legend = Legend(
            items=[("+", legend_plot.renderers[0:1]),
                   ("-", legend_plot.renderers[1:])])
        legend_plot.add_layout(legend, 'right')

        grid = gridplot(sorted_plots, ncols=ncols)
        tabs.append(Panel(child=grid, title=id_))

        # Extract target coverage
        cov = pd.DataFrame(df.coverage_f + df.coverage_r)
        cov.columns = ['coverage']
    cover_panel = Tabs(tabs=tabs)
    section.plot(cover_panel)


def make_coverage_summary_table(report: WFReport,
                                sample_ids: List[str],
                                table_file: Path,
                                seq_stats: List[Path],
                                on_offs: Path):
    """
    Summary table detailing all on and off target reads.

    On target here means:
    >=1bp overlap with target and off target the rest. Do we need to change the
    definition here to exclude proximal hits from the off-targets as is done
    later.

    Merge the data from all samples into a single table.

    :param seq_stats:  the summary from fastcat
    """
    section = report.add_section()
    section.markdown('''
        ### Summary of on-target and off-target reads.
        On target reads are defined here as any read that contains at least 1pb
        overlap with a target region and off target reads have 0 overlapping
        bases.
        ''')
    sample_frames = []

    id_stats = {k: v for k, v in zip(sample_ids, seq_stats)}
    df_onoff = pd.read_csv(
        on_offs,
        sep='\t',
        names=['chr', 'start', 'end', 'read_id', 'target', 'sample_id'],
        index_col=False)

    main_df = pd.read_csv(
        table_file, sep='\t', names=['on target', 'off target', 'sample_id'])

    for id_, df in main_df.groupby('sample_id'):
        df = df.drop(columns=['sample_id'])
        df['all'] = df['on target'] + df['off target']

        df = df.T
        df.columns = ['num_reads', 'kbases mapped']
        df['kbases mapped'] = df['kbases mapped'] / 1000

        df_stats = pd.read_csv(id_stats[id_], sep='\t')

        df_m = df_onoff[df_onoff['sample_id'] == id_]
        df_m = df_m.merge(
            df_stats[['read_id', 'read_length']],
            left_on='read_id', right_on='read_id')
        df_m['target'] = df_m['target'].fillna('OFF')

        mean_read_len = [df_m[df_m.target != 'OFF'].read_length.mean(),
                         df_m[df_m.target == 'OFF'].read_length.mean(),
                         df_m.read_length.mean()]

        df['mean read length'] = mean_read_len
        df.fillna(0, inplace=True)
        df = df.astype('int')

        # Melt the dataframe into a single row
        dfu = df.unstack().to_frame().sort_index(level=1).T
        dfu.insert(0, 'sample', id_)
        sample_frames.append(dfu)

    df_all_samples = pd.concat(sample_frames)

    df_all_samples.sort_values(
        by=["sample"],
        key=natsort_keygen(),
        inplace=True)
    # Sort the multiindex columns
    df_all_samples = df_all_samples.T.sort_index(ascending=False).T

    section.table(df_all_samples, searchable=True, paging=True, index=False)


def make_target_summary_table(report: WFReport, sample_ids: List,
                              table_file: Path,
                              seq_stats, on_off_bed):
    """Create a table of target summary statistics.

    TODO: missing mean accuracy column
    """
    section = report.add_section()
    section.markdown('''
        ### Target regions summary

        This table provides a summary of all the target region detailing:

        * chr, start, end: the location of the target region.
        * \\#reads: number of reads mapped to target region.
        * \\#basesCov: number of bases in target with at least 1x coverage.
        * targetLen: length of target region.
        * fracTargAln: proportion of the target with at least 1x coverage.
        * medianCov: median coverage of 100 bp bins.
        * meanReadlen: mean read length of reads mapping to target.
        * strandBias: proportional difference of reads aligning to each strand.
            A value or +1 or -1 indicates complete bias to the foward or
            reverse strand respectively.
        * kbases: kbases of total reads mapped to target.
        ''')

    # Note meanAlnLen: needs to be switched to meanreadLen in next version
    header = ['chr', 'start', 'end', 'target', '#reads', '#basesCov',
              'targetLen', 'fracTargAln', 'medianCov', 'p', 'n', 'sample_id']

    frames = []
    id_stats = {k: v for k, v in zip(sample_ids, seq_stats)}

    df_onoff = pd.read_csv(
        on_off_bed,
        sep='\t',
        names=['chr', 'start', 'end', 'read_id', 'target', 'sample_id'],
        index_col=False)

    main_df = pd.read_csv(
        table_file, sep='\t', names=header, index_col=False)

    for id_, df in main_df.groupby('sample_id'):
        df = df.drop(columns=['sample_id'])
        if len(df) == 0:
            continue

        df_stats = pd.read_csv(id_stats[id_], sep='\t')
        df_on_off = df_onoff.merge(df_stats[['read_id', 'read_length']],
                                   left_on='read_id', right_on='read_id')

        read_len = df_on_off.groupby(['target']).mean()[['read_length']]
        read_len.columns = ['meanReadLen']
        if len(read_len) > 0:
            df = df.merge(read_len, left_on='target', right_index=True)
        else:
            df['meanReadLen'] = 0

        kbases = df_on_off.groupby(['target']).sum()[['read_length']] / 1000
        kbases.columns = ['kbases']
        if len(kbases) > 0:
            df = df.merge(kbases, left_on='target', right_index=True)
        else:
            df['kbases'] = 0

        df['strandBias'] = (df.p - df.n) / (df.p + df.n)
        df.drop(columns=['p', 'n'], inplace=True)
        df.insert(0, 'sample', id_)
        frames.append(df)

    if len(frames) > 0:
        df_all = pd.concat(frames)
        df_all = df_all.astype({
            'start': int,
            'end': int,
            '#reads': int,
            '#basesCov': int,
            'targetLen': int
        })

        df_all = df_all.round({'strandBias': 2,
                               'fracTargAln': 2,
                               'kbases': 2,
                               'meanReadLen': 1})
        df_all.sort_values(
            by=["sample", "chr", "start"],
            key=natsort_keygen(),
            inplace=True)
    else:
        df_all = pd.DataFrame()

    section.table(df_all, searchable=True, paging=True)


def plot_tiled_coverage_hist(report: WFReport, background: List[Path],
                             target_coverage: List[Path]):
    """Coverage histograms.

    Show on-target and off-target (proximal removed) coverage
    over bins.
    """
    section = report.add_section()
    section.markdown('''
            ### Coverage distribution

            This histogram(s) show the coverage distribution of on-target and
            or off-target (background) reads binned by 100bp
            genome tiles. Off-target regions are defined as any region not
            within 1kb of a target.

            The background histogram should naturally be be skewed heavily
            to the left, this noise being expected when many regions in the
            genome have a single read mapping.

            If the targeted sequencing approach has performed well,
            the target histogram should be skewed towards the right
            as there has been a depletion of non-target reads.

            ''')

    header_target = ["chr", "start", "end", "target", "coverage_f",
                     'coverage_r', 'sample_id']
    header_background = ['chr', 'start', 'end', 'tile_name', '#reads',
                         '#bases_cov', 'tileLen', 'fracTileAln', 'sample_id']

    plots = []

    df_target = pd.read_csv(target_coverage, sep='\t', names=header_target)
    df_background = pd.read_csv(background, sep='\t', names=header_background)
    for id_, dfb, in df_background.groupby('sample_id'):
        dft = df_target[df_target.sample_id == id_]
        # Extract target coverage
        tc = pd.DataFrame(dft.coverage_f + dft.coverage_r)
        tc.columns = ['coverage']

        len_bg = len(dfb['#reads'].values)
        len_target = len(tc['coverage'])
        weights = [[1 / len_bg] * len_bg,
                   [1 / len_target] * len_target]

        plot = hist.histogram([dfb['#reads'].values,
                               tc['coverage']],
                              colors=['#1A85FF',
                                      '#D41159'],
                              normalize=True,
                              weights=weights,
                              names=['Background',
                                     'On-target'],
                              x_axis_label='Coverage',
                              y_axis_label='Proportion of reads (normalized'
                              'by class size)',
                              title=id_
                              )
        plots.append(plot)
    grid = gridplot(plots, ncols=3, width=360, height=300)
    section.plot(grid)


def make_offtarget_hotspot_table(report: WFReport, background: Path,
                                 nreads_cutoff=10):
    """Make a table of off-target hotspot regions.

    :param background: fill in
    :param nreads_cutoff: threshold for inclusion in background hotspot table

    Using aplanat.report.FilterableTable to crate a column-searchable table.
    """
    section = report.add_section()
    section.markdown('''
            ### Off-target hotspots

            Off target regions are again defined here as all regions of the
            genome not within 1kb of a target region.

            An off-target hotspot is a off-target region with contiguous
            overlapping reads. These hotspots may indicate incorrectly-
            performing primers. Only regions with {} reads or more are included
            '''.format(nreads_cutoff))
    main_df = pd.read_csv(
        background, sep='\t',
        names=['chr', 'start', 'end', 'numReads', 'sample_id'])

    for id_, df in main_df.groupby('sample_id'):
        df = df[df.numReads >= nreads_cutoff]
        df['hotspotLength'] = df.end - df.start
        df = df[['chr', 'numReads', 'start', 'end', 'hotspotLength']]
        df.sort_values('numReads', ascending=False, inplace=True)
        # Just a few rows for init view until we can use tables in tabs
        tab_params = {'pageLength': 15}
        section.markdown(f'Sample: {id_}')
        section.filterable_table(df, index=False, table_params=tab_params)


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("report", help="Report output file")
    parser.add_argument("--summaries", nargs='+', help="Read summary file.")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--sample_ids", required=True, nargs='+',
        help="List of sample ids")
    parser.add_argument(
        "--coverage_summary", required=True, type=Path,
        help="Contigency table coverage summary csv")
    parser.add_argument(
        "--target_coverage", required=False, default=None,
        type=Path, help="Tiled coverage for each target")
    parser.add_argument(
        "--target_summary", required=True, type=Path,
        help="Summary stats for each target. CSV.")
    parser.add_argument(
        "--background", required=False, default=None, type=Path,
        help="Tiled background coverage")
    parser.add_argument(
        "--off_target_hotspots", required=False, default=None,
        type=Path, help="Tiled background coverage")
    parser.add_argument(
        "--on_off", required=True, type=Path,
        help="Bed file. 5th column containing target or empty for off-target")
    args = parser.parse_args()

    report = WFReport(
        "Workflow for analysis of cas9-targeted sequencing", "wf-cas9",
        revision=args.revision, commit=args.commit)

    intro_section = report.add_section()
    intro_section.markdown('''
    The workflow aids with the quantification of the non-target depletion and
    provides information on mapping characteristics that highlight the cas9
    targeted sequencing protocol performance. The figures plotted include
    depth-of-coverage over the target regions and strand bias over these
    regions. The location and peaks of coverage and local biases in
    strandedness may be used to assess the performance of guide-RNA sequences
    and may highlight guide RNAs that are not performing. A review of likely
    off-target regions over-represented within the sequence collection may
    inform of strategies to refine guide-RNA design.
    ''')

    # Add reads summary section
    section = report.add_section()
    section.markdown("### Read stats")
    for id_, summ in sorted(zip(args.sample_ids, args.summaries)):
        report.add_section(
            section=fastcat.full_report(
                [summ],
                header="{}".format(id_)
            ))

    make_coverage_summary_table(report, args.sample_ids, args.coverage_summary,
                                args.summaries, args.on_off)

    make_target_summary_table(report, args.sample_ids, args.target_summary,
                              args.summaries, args.on_off)

    if args.target_coverage:
        plot_target_coverage(report, args.target_coverage)

    if args.background and args.target_coverage:
        plot_tiled_coverage_hist(report, args.background, args.target_coverage)

    if args.off_target_hotspots:
        make_offtarget_hotspot_table(report, args.off_target_hotspots)

    report.add_section(
        section=scomponents.version_table(args.versions))
    report.add_section(
        section=scomponents.params_table(args.params))

    # write report
    report.write(args.report)


if __name__ == "__main__":
    main()
