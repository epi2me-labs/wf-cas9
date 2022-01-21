#!/usr/bin/env python
"""Create workflow report."""

from pathlib import Path
import argparse

from aplanat import bars, lines
from aplanat.components import fastcat
from aplanat.components import simple as scomponents
from aplanat.report import WFReport
from bokeh.layouts import gridplot
import pandas as pd



def target_coverage_plots(report: WFReport, target_coverage: Path):
    section = report.add_section()
    section.markdown('''
    ### Target coverage 
    ''')

    dfg = pd.read_csv(target_coverage, index_col=0).groupby('target')#

    plots = []
    for target, df in dfg:
        chrom = df.loc[df.index[0], 'chrom']
        p = lines.line(
            [df.start.values, df.start.values],  # x-values
            [df.overlaps_f, df.overlaps_r],      # y-values
            title="{}".format(target),
            x_axis_label='{}'.format(chrom),
            y_axis_label='Coverage',
            colors=['blue', 'red'])
        p.xaxis.formatter.use_scientific = False
        p.xaxis.major_label_orientation = 3.14 / 6
        plots.append(p)


    grid = gridplot(plots, ncols=5, plot_width=250, plot_height=200)

    section.plot(grid)


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
        "--coverage_summary", required=True, type=Path,
        help="Contigency table coverage summary csv")
    parser.add_argument(
        "--target_coverage", required=True, type=Path,
        help="Tiled coverage for each target")

    args = parser.parse_args()

    report = WFReport(
        "Workflow Template Sequencing report", "wf-template",
        revision=args.revision, commit=args.commit)

    report.add_section(
        section=fastcat.full_report(args.summaries))
    report.add_section(
        section=scomponents.version_table(args.versions))
    report.add_section(
        section=scomponents.params_table(args.params))

    target_coverage_plots(report, args.target_coverage)

    # write report
    report.write(args.report)


if __name__ == "__main__":
    import sys
    sys.argv.extend([
        '/Users/Neil.Horner/work/testing/cas9/output/report.html',
        '--summaries', '/Users/Neil.Horner/work/testing/cas9/test_data/fastq_pass.stats',
        '--versions', '/Users/Neil.Horner/work/testing/cas9/test_data/versions.txt',
        '--params', '/Users/Neil.Horner/work/testing/cas9/test_data/params.json',
        
        '--coverage_summary', '/Users/Neil.Horner/work/testing/cas9/output/coverage_summary.csv',
        '--target_coverage', '/Users/Neil.Horner/work/testing/cas9/output/target_coverage.csv'
    ])

    main()
