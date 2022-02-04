import pandas as pd
import sys
sys.path.insert(0, '/Users/Neil.Horner/git/aplanat')
from aplanat.report import WFReport
from pathlib import Path
import numpy as np


def _make_offtarget_hotspot_table(report: WFReport, bg: Path):

    section = report.add_section()
    section.markdown('''
            ### Off-target hotspots
            ''')
    df = pd.read_csv(bg, sep='\t', names=['chr', 'start', 'end', 'num_reads'],
                    )
    df = df[['chr', 'start']]
    df['t'] = 1

    # df[0] = 'd'

    # df.drop(columns=['num_reads'], inplace=True)
    # df['+'] = '_'
    # df = df.sort_values('num_reads', ascending=False).reset_index(drop=True)
    df.to_csv('/Users/Neil.Horner/Desktop/t/report_cas9.csv')
    section.filterable_table(df, index=False)


def _make_offtarget_hotspot_table_iso(report: WFReport, bg: Path):
    section = report.add_section()
    section.markdown('''
            ### Off-target hotspots
            ''')
    df = pd.read_csv(bg
                     )
    df.to_csv('/Users/Neil.Horner/Desktop/t/report_cas9.csv')
    section.filterable_table(df, index=False)


if __name__ == '__main__':
    report = WFReport('test', 'test')
    _make_offtarget_hotspot_table(report, '/Users/Neil.Horner/work/workflow_outputs/cas9/workspace/2a/a01b5f36f1974952f95ecaeaeb1e80/fastq_pass_off_target_hotspots.bed')
    # _make_offtarget_hotspot_table_iso(report, '/Users/Neil.Horner/Desktop/t/report_iso.csv')
    report.write('/Users/Neil.Horner/Desktop/t/report_iso.html')

    # report = WFReport('test', 'test')
    # # _make_offtarget_hotspot_table(report, '/Users/Neil.Horner/work/workflow_outputs/cas9/workspace/2a/a01b5f36f1974952f95ecaeaeb1e80/fastq_pass_off_target_hotspots.bed')
    # _make_offtarget_hotspot_table(report,
    #                                   '/Users/Neil.Horner/Desktop/t/report_cas9.csv')
    # report.write('/Users/Neil.Horner/Desktop/t/report_cas9.html')

