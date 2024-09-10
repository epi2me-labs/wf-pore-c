# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### [unreleased]
### Fixed
- Set format of `--bam` and `--fastq` in schema to `path`, to enable directories to be selected as input in the EPI2ME application.

### [v1.2.2]
### Fixed
- Capitalised modified base tags additionally removed from monomers if no modified bases for a monomer.

## [v1.2.1]
### Fixed
- bamindex fetch error when running more than one sample and `--chunk_size` is greater than 0.
  
## [v1.2.0]
### Fixed
- `--bed` parameter will now output BED file using the paired_end BAM file.
- Reduce memory usage of BED file creation and sorting.
- Increased memory allocation for `prepare_hic` and `merge_mcool` processes.
- If sample sheet provided and cutter column not present the workflow will instead use `--cutter` parameter.
### Changed
- Bump pore-c-py to v2.1.4 to prevent issues with modified base tags and strip minimap2 tags from inputs.
### Added
- Reduce peak memory usage of minimap2 by adding `--cap-kalloc 100m --cap-sw-mem 50m` to the minimap2 command in the `digest_align_annotate` process.

## [v1.1.0]
### Added
- `--bed` parameter which if set to true will output a BED file that is compatible with downstream tools including the scaffolder [Yahs](https://github.com/c-zhou/yahs).
- `--pairtools_chunksize` parameter which exposes pairtools dedup chunksize parameter, in case peak memory usage of hi_c process needs to be reduced.
- `digest_align_annotate` process uses dedicated pore_c_py container.
- `--max_monomers` parameter, which is set to 250 by default, will filter out any reads that have more than this number of monomers. These reads will not be included in the analysis.
- Output a `filtered_out/{alias}.bam` with any reads that are filtered out due to the max_monomers parameter.

## [v1.0.0]
### Changed
- New documentation.

## [v0.2.0]
### Fixed
- Pairtools merge step single quote the input directory so it will not error with Argument list too long.
- Chromunity parquet files now contain the correct column names.
### Changed
- `--ubam` parameter has been renamed `--bam`
- All other ubam related parameters have been renamed with bam for consistency
- The `--bam_map_threads`, `--digest_annotate_threads` and `bam_bam2fq_threads` threading parameters are now automatically extracted from the `--threads` specifying the maximum number of threads to use for a process.
### Removed
- Default local executor CPU and RAM limits.

## [v0.1.1]
### Changed
- If `--hi_c` parameter set to true the pairs file will be created. 

## [v0.1.0]
### Changed
- GitHub issue templates
- Nextflow minimum version 23.04.2.
- `--sample_id` parameter has been changed to `--sample` for consistency.
- `--summary_json` optional parameter with default set to true, to include an annotation summary json in outputs.
- Remove `--params_sheet` parameter and add all per sample parameters to sample_sheet.

### Added
- `--hi_c` optional parameter with default set to false, to include a `.hic` output file which is compatible with [Juice box](https://www.aidenlab.org/juicebox/).

## [v0.0.8]
* Improve schema parameter explanations and output file descriptions in the README.
* Add a default `--chunk_size` parameter value of 25000.
* Update fastcat which removes need to index ubam.
* Enum choices are enumerated in the `--help` output.
* Enum choices are enumerated as part of the error message when a user has selected an invalid choice.
* Bumped minimum required Nextflow version to 22.10.8.

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
* First release of Wf-Pore-C

