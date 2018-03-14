# Change Log
All notable changes to the bio tools will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/), 
and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]

## [0.2.6] - 2018-03-14
### Added
- Added to bio_programs `bed_to_bigBed`, `run_signalp` and `run_tmhmm`.

### Removed
- Removed deprecated function calls from bio_programs `parse_meme_names`, 
    `create_mygene_report`, `get_all_utrs`, `get_read_length_distribution`, 
    `count_reads` and `count_aligned_reads`.
- Removed deprecated function calls from `bed_utils`.
- Removed deprecated function calls from 
    `fastx_utils`, `bam_utils`.
- Removed deprecated function call from `gtf_utils`.

### Fixed
- Minor changes to `bio_utils.plotting.plot_read_length_distribution`, 
    see [ISSUE #87](https://github.com/dieterich-lab/rp-bp/issues/87).

### Added
- Added *partial support* for GFF3 specifications, in particular this requires removing the
    STOP codons from the CDSs. Reference annotations file with gff extension are now 
    treated as GGF3.
- Handling of additional options passed to star via arguments to `star_utils`. 

## [0.2.5] - 2017-12-08
### Updated
- Removed deprecated function call from `mygene_utils`

## [0.2.4] - 2017-10-26
### Updated
- Version specifications for prereqs

## [0.2.3] - 2017-07-27
### Updated
- Updated download-srr-files to use `misc.shell_utils` and http rather than ftp

## [0.2.2] - 2017-06-14
### Fixed
- References to old gtf functions in `misc.bio`

## [0.2.1] - 2017-06-14
### Fixed
- All references to `misc.bio` and `misc.bio_utils`. Please see
  [Issue #1](https://github.com/bmmalone/pybio-utils/issues/1) in the new repo
  for more details.

## [0.2.0] - 2017-05-31
This is a new version which moves the project from Bitbucket to GitHub.
Additionally, the other utilities (`misc.***`) have been completely
removed. They will be added to a new
[`pymisc-utils`](https://github.com/bmmalone/pymisc-utils) repository.

## [0.1.6] - 2017-05-10
### Fixed
- Missing import in counting alignments for bam files

## [0.1.5] - 2017-05-09
### Updated
- `get-read-length-distribution` script to handle bam, fasta and fastq files.
  See [the bio docs](docs/bio.md#get-read-length-distributions) for more
  details.

## [0.1.4] - 2017-05-09
### Removed
- bed, bam and gtf helpers from `bio.py`. These had already been deprecated for
  quite some time.

## [0.1.3] - 2017-03-30
### Added
- Script to remove duplicate entries from bed files. See
  [the bio docs](docs/bio.md#merge-bed12-files-and-remove-duplicate-entries)
  for more details.

## [0.1.2] and previous versions

The initial versions have not been documented in the change log.


