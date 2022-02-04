#!/usr/bin/env python
"""Create workflow report."""

from pathlib import Path
import argparse

from aplanat import hist, lines
from aplanat.components import fastcat
from aplanat.components import simple as scomponents
from aplanat.report import WFReport
from bokeh.layouts import gridplot
from bokeh.models import Legend
from natsort import natsorted, natsort_keygen
import pandas as pd


def _plot_target_coverage(report: WFReport, target_coverage: Path):
    section = report.add_section()
    section.markdown('''
    ### Target coverage 
    ''')

    header = ["chr", "start", "end", 'name_f', "target", "coverage_f",
              'name_r', 'coverage_r']
    df = pd.read_csv(target_coverage, names=header, sep='\t')
    dfg = df.groupby('target')

    ncols = 4
    plots = []
    for i, (target, df) in enumerate(dfg):
        chrom = df.loc[df.index[0], 'chr']
        ymax = max(df.coverage_f.max(), df.coverage_r.max())
        ylim = [0, ymax * 1.05]  # a bit of space at top of plot

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
        items=[("+", legend_plot.renderers[0:1]), ("-", legend_plot.renderers[1:])])
    legend_plot.add_layout(legend, 'right')

    grid = gridplot(sorted_plots, ncols=ncols)

    section.plot(grid)

    # Extract target coverage
    cov = pd.DataFrame(df.coverage_f + df.coverage_r)
    cov.columns = ['coverage']
    return cov


def make_coverage_summary_table(report: WFReport, table_file: Path,
                                 seq_stats: Path, on_off: Path):
    """
    Summary table all on and off target reads. On target here means
    >=1bp overlap with target and off target the rest. Do we need to change the
    definition here to exclude proximal hits from the off-targets as is done
    later
    """
    section = report.add_section()
    section.markdown('''
        ### Summary of on-target and off-target reads. 
        On target
        ''')
    df = pd.read_csv(table_file, sep='\t', names=['on target', 'off target'])
    df['all'] = df['on target'] + df['off target']

    df = df.T
    df.columns = ['num_reads', 'kbases']
    df.kbases = df.kbases / 1000

    df_stats = pd.read_csv(seq_stats, sep='\t')

    df_onoff = pd.read_csv(on_off, sep='\t',
                           names=['chr', 'start', 'end', 'read_id', 'target'])

    df_onoff['target'].fillna('OFF', inplace=True)
    df_m = df_onoff.merge(df_stats[['read_id', 'read_length']],
                   left_on='read_id', right_on='read_id')

    mean_read_len = [df_m[df_m.target != 'OFF'].read_length.mean(),
           df_m[df_m.target == 'OFF'].read_length.mean(),
                     df_m.read_length.mean()]

    df['mean read length'] = mean_read_len
    df = df.astype('int')

    section.table(df, searchable=False, paging=False, index=True)


def _make_target_summary_table(report: WFReport, table_file: Path):
    section = report.add_section()
    section.markdown('''
        ### Summary of each target
        ''')
    header = ['chr', 'start', 'end', 'target', '#reads', '#bases_cov',
              'targetLen', 'fracTargAln', 'meanReadLen', 'kbases',
              'medianCov', 'p', 'n']

    df = pd.read_csv(table_file, sep='\t', names=header)
    df['strand bias'] = (df.p - df.n) / (df.p + df.n)
    df.drop(columns=['p', 'n'], inplace=True)
    df.sort_values(
        by=["chr", "start"],
        key=natsort_keygen(),
        inplace=True
    )
    section.table(df, searchable=False, paging=False)


def _plot_background(report: WFReport, background: Path,
                    target_coverage: pd.DataFrame):
    section = report.add_section()
    section.markdown('''
            ### Coverage distribution
            ''')
    header = ['chr', 'start', 'end', 'tile_name', '#reads', '#bases_cov',
              'tileLen', 'fracTileAln']

    df = pd.read_csv(background, sep='\t', names=header)
    target_weight = len(df) / len(target_coverage)
    weights = [[1] * len(df),
               [target_weight] * len(target_coverage)]

    plot = hist.histogram([df['#reads'].values, target_coverage['coverage']],
                          colors=['#1A85FF', '#D41159'], normalize=True,
                          weights=weights, names=['Background', 'target'])
    section.plot(plot)


def _make_offtarget_hotspot_table(report: WFReport, bg: Path):

    section = report.add_section()
    section.markdown('''
            ### Off-target hotspots
            ''')
    df = pd.read_csv(bg, sep='\t', names=['chr', 'start', 'end', 'num_reads'],
                     )
    df.sort_values('num_reads', ascending=False, inplace=True)
    # t = 'columnDefs: [{ "width": "5%", "targets": [2, 3] }]'
    section.filterable_table(df, index=False, table_params=None)


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
        "--target_coverage", required=True, type=Path,
        help="Tiled coverage for each target")
    parser.add_argument(
        "--target_summary", required=True, type=Path,
        help="Summary stats for each target. CSV.")
    parser.add_argument(
        "--background", required=True, type=Path,
        help="Tiled background coverage")
    parser.add_argument(
        "--off_target_hotspots", required=True, type=Path,
        help="Tiled background coverage")
    parser.add_argument(
        "--on_off", required=True, type=Path,
        help="Bed with 5th column of target name of 'off")
    args = parser.parse_args()

    report = WFReport(
        "Workflow Template Sequencing report", "wf-cas9",
        revision=args.revision, commit=args.commit)

    # Add reads summary section
    for id_, summ in zip(args.sample_ids, args.summaries):
        report.add_section(
            section=fastcat.full_report(
                [summ],
                header='#### Read stats: {}'.format(id_)
            ))

    _make_coverage_summary_table(report, args.coverage_summary, args.summaries[0],
                                 args.on_off)
    _make_target_summary_table(report, args.target_summary)
    target_coverage = _plot_target_coverage(report, args.target_coverage)
    _plot_background(report, args.background, target_coverage)
    _make_offtarget_hotspot_table(report, args.off_target_hotspots)

    report.add_section(
        section=scomponents.version_table(args.versions))
    report.add_section(
        section=scomponents.params_table(args.params))

    # write report
    report.write(args.report)


if __name__ == "__main__":
    import sys
    # sys.argv.extend([
    #     '/Users/Neil.Horner/work/testing/cas9/output/report.html',
    #     '--summaries', '/Users/Neil.Horner/work/testing/cas9/test_data/fastq_pass.stats',
    #     '--versions', '/Users/Neil.Horner/work/testing/cas9/test_data/versions.txt',
    #     '--params', '/Users/Neil.Horner/work/testing/cas9/test_data/params.json',
    #
    #     '--coverage_summary', '/Users/Neil.Horner/work/testing/cas9/output/coverage_summary.csv',
    #     '--target_coverage', '/Users/Neil.Horner/work/testing/cas9/output/target_coverage.csv'
    # ])

    main()
