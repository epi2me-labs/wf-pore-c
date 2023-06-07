# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Changed
- Enum choices are enumerated in the `--help` output
- Enum choices are enumerated as part of the error message when a user has selected an invalid choice
### Fixed
- Replaced `--threads` option in fastqingress with hardcoded values to remove warning about undefined `param.threads`

## [v0.0.7]
### Fixed
- Testing for the cooler tool.
  
## [v0.0.6]
### Added
- Configuration for running demo data in AWS

## [v0.0.5]
### Fixed
- Broken heat map in the pairtools report.
- Meta table repeated tabs.
- Nextflow config example cmd.

### Added
- Cutter parameter help text link to Restriction Enzyme options.

## [v0.0.4]
### Added
- Changed LICENSE to Oxford Nanopore Technologies PLC. Public License Version 1.0.
- Test for Chromunity writer

### Fixed
- Use latest pore-c-py package with fix for the modified bases digest step.

## [v0.0.3]
### Fixed
- Reduce time by using bamindex instead of splitting bam.

### Changed
- Replace input check with fastq ingress.
- Parameters to input fastq or ubam.
- Output a basic report.

## [v0.0.2]
### Fixed
- Create pairs report handling of missing references in pairs file.

### Changed
- Update Pore-c-py package used to v2.0.1
- Improved performance
- Use one pipe for digest, align and annotate processes.

## [v0.0.1]
- First release of Wf-Pore-C

