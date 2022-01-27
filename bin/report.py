#!/usr/bin/env python
"""Create workflow report."""

from pathlib import Path
import argparse

from aplanat import bars, lines
from aplanat.components import fastcat
from aplanat.components import simple as scomponents
from aplanat.report import WFReport
from bokeh.layouts import gridplot
from bokeh.models import Legend
import pandas as pd


def plot_target_coverage(report: WFReport, target_coverage: Path):
    section = report.add_section()
    section.markdown('''
    ### Target coverage 
    ''')

    dfg = pd.read_csv(target_coverage, index_col=0).groupby('target')

    ncols = 5
    plots = []
    for i, (target, df) in enumerate(dfg):
        chrom = df.loc[df.index[0], 'chrom']
        p = lines.line(
            [df.start.values, df.start.values],  # x-values
            [df.overlaps_f, df.overlaps_r],      # y-values
            title="{}".format(target),
            x_axis_label='{}'.format(chrom),
            y_axis_label='',
            colors=['#1A85FF', '#D41159']
            )
        p.xaxis.formatter.use_scientific = False
        p.xaxis.major_label_orientation = 3.14 / 6
        if i == ncols - 1:
            legend = Legend(
                items=[("+", p.renderers[0:1]), ("-", p.renderers[1:])])
            p.add_layout(legend, 'right')

        plots.append(p)

    grid = gridplot(plots, ncols=5, width=250, height=200)

    section.plot(grid)


def make_coverage_summary_table(report: WFReport, table_file: Path):
    section = report.add_section()
    section.markdown('''
        ### TargetSummary 
        ''')
    df = pd.read_csv(table_file)
    df.iloc[:, 1:] = df.iloc[:, 1:].astype(int)
    df.rename(columns={df.columns[0]: ""}, inplace=True)
    section.table(df, searchable=False, paging=False)


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
    report.add_section(
        section=fastcat.full_report(args.summaries))
    report.add_section(
        section=scomponents.version_table(args.versions))
    report.add_section(
        section=scomponents.params_table(args.params))

    make_coverage_summary_table(report, args.coverage_summary)
    plot_target_coverage(report, args.target_coverage)

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
