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


def plot_target_coverage(report: WFReport, sample_ids,
                         target_coverages: List[Path]):
    """Make coverage plots of each target.

    Detailing positive and negative strand coverage.

    :param target_coverages
    """
    section = report.add_section()
    section.markdown('''
    <br><br>
    ### Target coverage

    Each of the following plots show the amount of coverage, per strand
    in discretized bins of 100 bp.
    ''')

    header = ["chr", "start", "end", 'name_f', "target", "coverage_f",
              'name_r', 'coverage_r']
    tabs = []
    all_cov = []  # merge f + r coverage for use in other functions
    for (id_, t_cov) in zip(sample_ids, target_coverages):
        df = pd.read_csv(t_cov, names=header, sep='\t')
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
        all_cov.append(cov)
    cover_panel = Tabs(tabs=tabs)
    section.plot(cover_panel)
    return all_cov


def make_coverage_summary_table(report: WFReport,
                                sample_ids: List,
                                table_files: List[Path],
                                seq_stats: List[Path],
                                on_offs: List[Path]):
    """
    Summary table detailing all on and off target reads.

    On target here means:
    >=1bp overlap with target and off target the rest. Do we need to change the
    definition here to exclude proximal hits from the off-targets as is done
    later

    :param seq_stats:  the summary from fastcat
    """
    section = report.add_section()
    section.markdown('''
        ### Summary of on-target and off-target reads.
        On target reads are defined here as any read that contains at least 1pb
        overlap with a target region and off target reads have 0 overlapping
        bases.
        ''')

    for id_, table_file, stats, on_off in \
            zip(sample_ids, table_files, seq_stats, on_offs):

        try:  # I'n not sure we need this
            df = pd.read_csv(
                table_file, sep='\t', names=[
                    'on target', 'off target'])
        except pd.errors.EmptyDataError:
            continue

        df['all'] = df['on target'] + df['off target']

        df = df.T
        df.columns = ['num_reads', 'kbases of sequence mapped']
        df['kbases of sequence mapped'] = \
            df['kbases of sequence mapped'] / 1000

        df_stats = pd.read_csv(stats, sep='\t')

        df_onoff = pd.read_csv(
            on_off,
            sep='\t',
            names=[
                'chr',
                'start',
                'end',
                'read_id',
                'target'])

        df_m = df_onoff.merge(df_stats[['read_id', 'read_length']],
                              left_on='read_id', right_on='read_id')

        mean_read_len = [df_m[df_m.target != 'OFF'].read_length.mean(),
                         df_m[df_m.target == 'OFF'].read_length.mean(),
                         df_m.read_length.mean()]

        df['mean read length'] = mean_read_len
        df.fillna(0, inplace=True)
        df = df.astype('int')

        section.markdown(f"Sample id: {id_}")
        section.table(df, searchable=False, paging=False, index=True)


def make_target_summary_table(report: WFReport, sample_ids: List,
                              table_files: List[Path]):
    """Create a table of target summary statistics.

    TODO: missing mean accuracy column
    """
    section = report.add_section()
    section.markdown('''
        <br>
        ### Targeted region summary

        This table provides a summary of all the target region detailing:

        * chr, start, end: the location of the target region
        * \\#reads: number of reads mapped to target region
        * \\#basesCov: number of bases in target with at least 1x coverage
        * targetLen: length of target region
        * fracTargAln: proportion of the target with at least 1x coverage
        * meanAlnlen: mean length alignment of alignment
        * strandBias: proportional difference of reads aligning to each strand.
            A value or +1 or -1 indicates complete bias to the foward or
            reverse strand respectively.
        * kbases: kbases of total reads mapped to target
        ''')

    # Note meanAlnLen: needs to be switched to meanreadLen in next version
    header = ['chr', 'start', 'end', 'target', '#reads', '#basesCov',
              'targetLen', 'fracTargAln', 'meanAlnlen', 'kbases',
              'medianCov', 'p', 'n']

    for (id_, table_file) in zip(sample_ids, table_files):
        df = pd.read_csv(table_file, sep='\t', names=header)
        df.kbases = df.kbases / 1000
        # This bodges a problem with main.nf:target_summary
        df.dropna(inplace=True)

        df['strandBias'] = (df.p - df.n) / (df.p + df.n)
        df.drop(columns=['p', 'n'], inplace=True)
        df.sort_values(
            by=["chr", "start"],
            key=natsort_keygen(),
            inplace=True
        )
        df = df.astype({
            'start': int,
            'end': int,
            '#reads': int,
            '#basesCov': int,
            'targetLen': int,
            'meanAlnlen': int,
            'kbases': int
        })
        df = df.round({'strandBias': 2})

        section.markdown(f"Sample id: {id_}")
        section.table(df, searchable=False, paging=False)


def plot_tiled_coverage_hist(report: WFReport, sample_ids: List,
                             background: List[Path], target_coverage:
                             List[pd.DataFrame]):
    """Coverage histograms.

    Show on-target and off-target (proximal removed) coverage
    over bins.
    """
    section = report.add_section()
    section.markdown('''
            <br>
            ### Coverage distribution

            This histogram(s) show the coverage distribution of on-target and
            or off-target (background) reads binned by 100bp
            genome tiles. Off-target regions are defined as any region not
            within 1kb of a target.

            The background histogram should naturally be be skewed heavily
            to the left, this noise being expected when many regions in the
            genome have a single read mapping.

            The target histogram should be skewed towards the right
            if targeted sequencing approach has enriched for reads at target
            regions.

            ''')
    header = ['chr', 'start', 'end', 'tile_name', '#reads', '#bases_cov',
              'tileLen', 'fracTileAln']

    plots = []
    for id_, bg, tc in zip(sample_ids, background, target_coverage):
        df = pd.read_csv(bg, sep='\t', names=header)

        len_bg = len(df['#reads'].values)
        len_target = len(tc['coverage'])
        weights = [[1 / len_bg] * len_bg,
                   [1 / len_target] * len_target]

        plot = hist.histogram([df['#reads'].values,
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
                              title=id_)
        plots.append(plot)
    grid = gridplot(plots, ncols=3, width=360, height=300)
    section.plot(grid)


def make_offtarget_hotspot_table(report: WFReport, sample_ids: List[str],
                                 background: List[Path],
                                 nreads_cutoff=10):
    """Make a table of off-target hotspot regions.

    :param background: fill in
    :param nreads_cutoff: threshold for inclusion in background hotspot table

    Using aplanat.report.FilterableTable to crate a column-searchable table.
    """
    section = report.add_section()
    section.markdown('''
            <br>
            ### Off-target hotspots

            Off target regions are again defined here as all regions of the
            genome not within 1kb of a target region.

            An off-target hotspot is a off-target region with contiguous
            overlapping reads. These hotspots may indicate incorrectly-
            performing primers. Only regions with {} reads or more are included
            '''.format(nreads_cutoff))
    for (id_, bg) in zip(sample_ids, background):
        df = pd.read_csv(
            bg, sep='\t', names=['chr', 'start', 'end', 'numReads']
        )
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
        "--coverage_summary", required=True, nargs='+', type=Path,
        help="Contigency table coverage summary csv")
    parser.add_argument(
        "--target_coverage", required=True, nargs='+', type=Path,
        help="Tiled coverage for each target")
    parser.add_argument(
        "--target_summary", required=True, nargs='+', type=Path,
        help="Summary stats for each target. CSV.")
    parser.add_argument(
        "--background", required=True, nargs='+', type=Path,
        help="Tiled background coverage")
    parser.add_argument(
        "--off_target_hotspots", required=True, nargs='+', type=Path,
        help="Tiled background coverage")
    parser.add_argument(
        "--on_off", required=True, nargs='+', type=Path,
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
    for id_, summ in zip(args.sample_ids, args.summaries):
        report.add_section(
            section=fastcat.full_report(
                [summ],
                header='#### Read stats: {}'.format(id_)
            ))

    make_coverage_summary_table(report, args.sample_ids, args.coverage_summary,
                                args.summaries, args.on_off)

    make_target_summary_table(report, args.sample_ids, args.target_summary)

    target_coverage = plot_target_coverage(report, args.sample_ids,
                                           args.target_coverage)

    plot_tiled_coverage_hist(report, args.sample_ids, args.background,
                             target_coverage)

    make_offtarget_hotspot_table(report, args.sample_ids,
                                 args.off_target_hotspots)

    report.add_section(
        section=scomponents.version_table(args.versions))
    report.add_section(
        section=scomponents.params_table(args.params))

    # write report
    report.write(args.report)


if __name__ == "__main__":
    main()
