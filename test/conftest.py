#!/usr/bin/env python
"""Pytests argument definitions."""


def pytest_addoption(parser):
    """Define command line arguments for pytest."""
    parser.addoption(
        "--wf_out_dir",
        action='store',
        default='/host/wf-cas9'
    )
    parser.addoption(
        "--sample_id",
        action="store",
        default="sample1"
    )
