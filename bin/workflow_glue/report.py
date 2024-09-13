#!/usr/bin/env python
"""Create workflow report."""
from math import pi
from pathlib import Path

from dominate.tags import h6, p
from dominate.util import raw
import ezcharts as ezc
from ezcharts.components.ezchart import EZChart
from ezcharts.components.fastcat import SeqSummary
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable, Grid, Tabs
from ezcharts.plots import util
from ezcharts.util import get_named_logger
import pandas as pd

from .util import wf_parser  # noqa: ABS101

# Setup simple globals
Colors = util.Colors


def argparser():
    """Create argument parser."""
    parser = wf_parser("Report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--read_stats",
        help="fastcat read stats file, with multiple samples concatenated")
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
        "--coverage_summary", required=True, type=Path,
        help=("table showing n reads, mean len, "
              "and kbases mapped for on/off-taget reads"))
    parser.add_argument(
        "--target_coverage", required=False, default=None,
        type=Path, help="Tiled coverage for each target for generating target plots")
    parser.add_argument(
        "--target_summary", required=True, type=Path,
        help="Summary stats for each target. TSV.")
    parser.add_argument(
        "--tile_coverage", required=False, default=None, type=Path,
        help="Tiled coverage TSV. columns: cov, sample_id, ontarget/offtarget ")
    parser.add_argument(
        "--off_target_hotspots", required=False, default=None,
        type=Path, help="Tiled background coverage")
    parser.add_argument(
        "--wf_version", default='unknown',
        help="version of the executed workflow")
    return parser


def target_table_and_plots(report, target_coverages, target_summaries):
    """Make coverage plots of each target.

    Detailing positive and negative strand coverage.
    """
    with report.add_section('Targets', 'Targets'):

        h6("Target coverage plots")

        p('''
        Each of the following plots show the amount of coverage per target, for each
        strand in discretized bins of 100 bp.
        ''')

        tabs = Tabs()

        main_df = pd.read_csv(target_coverages, sep='\t')

        for sample_id, df in main_df.groupby('sample_id'):
            if sample_id == 'sample_id':
                continue
            with tabs.add_tab(sample_id):
                with Grid(columns=4):
                    for target, df_target in df.groupby('target'):
                        chrom = df.loc[df.index[0], 'chr']
                        # Show chrom as xAxis label
                        df_target = df_target.rename(columns={'start': chrom})

                        plot = ezc.lineplot(
                            df_target,
                            title=target,
                            x=chrom,
                            y='coverage',
                            hue='strand',
                            marker=False)

                        plot._fig.x_range.start = df_target[chrom].min()
                        plot._fig.xaxis.major_label_orientation = 35 * (pi / 180)
                        plot._fig.xaxis.major_label_standoff = 5
                        plot._fig.add_layout(plot._fig.legend[0], 'right')

                        EZChart(plot, theme='epi2melabs', height='250px')

        h6("Target summaries")
        raw('''
            This table provides summaries for each sample/target combination.
            ''')

        df_all = pd.read_csv(target_summaries, index_col=False)
        df_all = df_all.astype({
            'kbases': int,
            'mean_read_length': int
        })

        tabs = Tabs()
        for sample_id, df in df_all.groupby('sample'):
            with tabs.add_tab(sample_id):
                DataTable.from_pandas(df, use_index=False)

        raw('''
                 <u>Column descriptions</u>:
                 <ul>
                     <li><b>chr</b>, start, end: target location.</li>
                 <li><b>target</b>: the target name. </li>
                 <li><b>nreads</b>: number of reads aligning. </li>
                 <li><b>coverage_frac</b>: fraction of bases within target with non-zero
                   coverage. </li>
                 <li><b>tsize</b>: length of target (in bases). </li>
                 <li><b>median_cov</b>: average read depth across target. </li>
                 <li><b>mean_read_length</b>:  average read length of reads aligning.
                    </li>
                 <li><b>strand_bias</b>: proportional difference of reads aligning to
                 each strand. </li>
                     A value or +1 or -1 indicates complete bias to the forward or
                     reverse strand respectively. </li>
                 <li><b>kbases</b>: number of bases in reads overlapping target. </li>
                 <br>
                 ''')


def coverage_summary_table(report, table_file):
    """
    Summary table detailing all on and off target reads.

    On target here means:
    >=1bp overlap with target and off target the rest. Do we need to change the
    definition here to exclude proximal hits from the off-targets as is done
    later.

    Merge the data from all samples into a single table.
    """
    with report.add_section('Summary of on-target and off-target reads', 'Summary'):
        p(
            """On target reads are defined here as any read that contains at least 1pb
            overlap with a target region. Off target reads contain no target-overlapping
            bases."""
        )
        df_all = pd.read_csv(
            table_file,
            sep='\t',
            header=[0, 1, 2],
            index_col=0)

        # Tidy up the columns index.
        df_all = df_all.droplevel(2, axis=1)
        df_all.index.name = 'sample'

        tabs = Tabs()
        for sample_id, df in df_all.groupby('sample'):
            with tabs.add_tab(sample_id):
                DataTable.from_pandas(df=df)


def tiled_coverage_hist(report, tiled_coverage_tsv):
    """Coverage histograms.

    Show on-target / off-target coverage at 100bp tiled genome bins.
    """
    with report.add_section('Coverage distribution', 'Coverage distribution'):
        p('''These plots show the on-target / off-target coverage distribution
            of genomic regions binned by 100bp.
            Off-target regions are defined as any region not
            within 1kb of a target.

            The background histogram should naturally be be skewed heavily
            to the left, this noise being expected when many regions in the
            genome have a single read mapping.

            If the targeted sequencing approach has performed well,
            the on_target histogram should be skewed towards the right
            indicating a depletion of non-target reads.
            ''')

        df_all = pd.read_csv(tiled_coverage_tsv, sep='\t')

        n_columns = 4
        with Grid(columns=n_columns):
            for i, (sample_id, df) in enumerate(df_all.groupby('sample_id')):
                plot = ezc.histplot(
                    df, title=sample_id, bins=5, hue='target_status', stat='proportion')
                plot._fig.xaxis.axis_label = 'coverage'
                plot._fig.yaxis.axis_label = 'Proportion of reads'
                plot._fig.xaxis.major_label_orientation = 30 * (pi / 180)
                # legend not working
                EZChart(plot, theme='epi2melabs', height='250px')


def make_offtarget_hotspot_table(report, background, nreads_cutoff=10):
    """Make a table of off-target hotspot regions.

    :param background: TSV
        columns:


    :param nreads_cutoff: threshold for inclusion in background hotspot table
    """
    with report.add_section('Off-target hotspots', 'Background'):
        p('''
            Off target regions are again defined here as all regions of the
            genome not within 1kb of a target region.

            An off-target hotspot is a off-target region with contiguous
            overlapping reads. These hotspots may indicate incorrectly-
            performing primers. Only regions with {} reads or more are included
            '''.format(nreads_cutoff))

        main_df = pd.read_csv(
            background, sep='\t',
            names=['chr', 'start', 'end', 'numReads', 'sample_id'])

        tabs = Tabs()
        for sample_id, df in main_df.groupby('sample_id'):
            with tabs.add_tab(sample_id):
                df = df[df.numReads >= nreads_cutoff]
                df.loc[:, 'hotspotLength'] = df.end - df.start
                df = df[['chr', 'numReads', 'start', 'end', 'hotspotLength']]
                df.sort_values('numReads', ascending=False, inplace=True)
                DataTable.from_pandas(df, use_index=False)


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    logger.info('Building report')

    report = LabsReport(
        'Workflow Cas9', 'wf-cas9', args.params, args.versions,
        args.wf_version, head_resources=[*LAB_head_resources])

    with report.add_section('Introduction', 'Intro'):
        p('''
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

    with report.add_section('Read summaries', 'Reads'):
        SeqSummary(args.read_stats)

    coverage_summary_table(report, args.coverage_summary)

    target_table_and_plots(report, args.target_coverage, args.target_summary)

    tiled_coverage_hist(report, args.tile_coverage)

    if args.off_target_hotspots:
        make_offtarget_hotspot_table(report, args.off_target_hotspots)

    report.write(args.report)
