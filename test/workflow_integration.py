"""Integration testing of the whole workflow using synthetic data."""

from pathlib import Path

import pandas as pd
from pytest import fixture


@fixture
def wf_out_dir(request):
    """Set workflow directory."""
    return request.config.getoption('--wf_out_dir')

def test_sample_summary(wf_out_dir):
    """Test the sample summary ouput."""
    out_dir = Path(wf_out_dir)

    sample_summary = out_dir / 'sample_summary.csv'
    assert sample_summary.is_file()

    result_df = pd.read_csv(sample_summary, sep=',', index_col=None)
    print(result_df.to_string())
    expected_df = pd.DataFrame([
        ['sample_1', 'run1', 3630, 1728.57, 0.31, 2100],
        ['sample_2', 'run1', 600, 1200, 1.0, 500],
        ], columns=[
        'sample_id', 'run_id', 'kbases', 'mean_read_length', 'strand_bias', 'nreads'])

    (pd.testing.assert_frame_equal
     (expected_df, result_df, check_dtype=False, check_like=True))



def test_target_summary(wf_out_dir):
    """Test the target summary output."""
    out_dir = Path(wf_out_dir)

    target_summary = out_dir / 'target_summary.csv'
    assert target_summary.is_file()

    result_df = pd.read_csv(target_summary, sep=',', index_col=None)

    # The following checks whether the results of the target_summary.tsv are as expected
    # The expected results are known as the reads they are generated from are
    # synthetic reads made to have the given properties.

    expected_df = pd.DataFrame([
        ['sample_1', 'run1', 'chr1', 3000, 6400,
            'SCA3', 3400, 1900.0, 0.94, 750, 1000, 1900.0,  0.23],
        ['sample_1', 'run1', 'chr2', 1000, 3000,
            'SCA10', 2000, 780.0,  1.0,  450, 600,  1300.0, -0.12],
        ['sample_1', 'run1', 'chr3', 3000, 6400,
            'SCA17', 3400, 950.0, 0.94, 375, 500, 1900.0, 1.0],
        ['sample_2', 'run1', 'chr3', 3000, 6400,
            'SCA17', 3400, 600, 0.53, 125, 500, 1200.0, 1.0],
    ], columns=[
        'sample', 'run_id', 'chr', 'start', 'end',
        'target', 'tsize', 'kbases', 'coverage_frac', 'median_cov',
        'nreads', 'mean_read_length', 'strand_bias'])
    (pd.testing.assert_frame_equal
     (expected_df, result_df, check_dtype=False, check_like=True))
