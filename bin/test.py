import pandas as pd

import sys
sys.path.insert(0, '/Users/Neil.Horner/git/workflows/wf-cas9/bin/report.py')

from pathlib import Path
from aplanat import report as r
import report


dir_ = Path('/Users/Neil.Horner/work/workflow_outputs/cas9/workspace/cb/66402ed6d345d237048983e97d7856')

rp = r.WFReport('test', 'test')


report._make_coverage_summary_table(rp, dir_ / 'fastq_pass_on_off_summ.csv',
                                    dir_ / 'fastq_pass.stats',
                                    dir_ / 'fastq_pass_on_off.bed')
