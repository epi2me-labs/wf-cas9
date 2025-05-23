# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.2.0]

This workflow is deprecated. It is not longer proactively 
maintained; and reactive maintainance is sporadic and low priority for our
developers. Although the repository is still available, it is now unsupported
and we do not recommend its use. Please contact support@nanoporetech.com for
help with your application.

### Changed
- Update project description.
- Reconcile workflow with wf-template v5.3.4.
- Minor decrease to some memory directives to avoid “Process requirement exceeds available memory” errors when running in WSL.

## [v1.1.3]
### Fixed
- Bug that caused incorrect `median_coverage` and `mean_read_length` in `target_summary.csv`
- `kbases` in `target_summary.csv` now reports total kbases mapped not mean per target.
### Changed
- Reconcile with template 5.3.0 and 5.3.1

## [v1.1.2]
### Changed
- Updated Ezcharts to v0.11.2.

## [v1.1.1]
### Changed
- The name of the column with run IDs from `run_ids` to `run_id`.

## [v1.1.0]
### Added
- A column with sequencing run IDs to the output summary CSV files.

### Changed
- Minor formatting changes of github issue template.

## [v1.0.1]
### Fixed
- Overestimation of kbases mapped values.
- Target coverage plots showing only one sample.

## [v1.0.0]
### Added
- Memory and CPU requirements for each process.

### Changed
- Documentation updated to new format.
- Bumped minimum required Nextflow version to '>=23.04.2'
- Publish target coverage output to output directory.

## [v0.1.12]
### Changed
- Enum choices are enumerated in the --help output.
- Enum choices are enumerated as part of the error message when a user has selected an invalid choice.

### Fixed
- Inputs file names and directories can now contain spaces.
- Replaced --threads option in fastqingress with hardcoded values to remove warning about undefined `param.threads`.
- Handling for file and directory names that contain spaces.

## [v0.1.11]
### Changed
- Bumped minimum required Nextflow version to 22.10.8

### Fixed
- Incompatible default memory options

## [v0.1.10]
### Added
- Configuration for running demo data in AWS

## [v0.1.9]
### Changed
- Generation report with ezcharts

### Fixed
- Memory issue in background calculation with large number of samples

## [v0.1.8]
### Fixed
- sample_sheet format in schema to expect a file

## [v0.1.7]
### Changed
- Updated description in manifest

## [v0.1.6]
### Updated
- Docs
- schema and config for using with new Labs Launcher

## [v0.1.5]
### Changed
- Output target summary csv

## [v0.1.4]
### Changed
- Fastqingress metadata map

## [v0.1.3]
### Changed
- Better help text on cli
- Replace pomoxis/stas_from_bam with fastcat/bamstats

## [v0.1.2]
### Added
- sample_summary.csv with per-sample summary statistics

### Changed
- Added sequence summary plots into tabs
- Minor report formatting
- New docs format

## [v0.1.1]
### Added
- No display of plots by default. Use --full_report flag to display

### Changed
- Merge samples from into single target summary table

### Fixed
- Incorrect values in target summary table
- Errors caused by cases where target had no reads mapping

## [v0.1.0]
### Added
- First version of wf-cas9 workflow

