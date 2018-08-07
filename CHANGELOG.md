# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### `Added`
  - New ViPR project
  - Logo
  - Draft of full documentation and completed Snakemake based testing
  - [#2](https://github.com/nf-core/vipr/pull/17) - Push container to [nfcore/vipr](https://hub.docker.com/r/nfcore/vipr/)
  - [#15](https://github.com/nf-core/vipr/issues/15), [#22](https://github.com/nf-core/vipr/pull/22) - This CHANGELOG
  - [#16](https://github.com/nf-core/vipr/issues/16), [#17](https://github.com/nf-core/vipr/pull/17) - Singularity file
  - [#22](https://github.com/nf-core/vipr/pull/22) - Enable Travis-CI

### `Changed`
  - Resource adjustment based on completed FMDV runs
  - Moved job settings into nextflow.config
  - nf-core submission
  - Enhancement to config, run environment and documentation
  - Adjustments to resource limits after Zika runs and minor cleanup

### `Removed`
  - [#18](https://github.com/nf-core/vipr/pull/18) Local executor
  - [#22](https://github.com/nf-core/vipr/pull/22) Institute logos

### `Fixed`
  - Working and complete version
  - BAI handling
  - Capturing special error to be ignored in gap_fill_assembly
  - pip install in Dockerfile and extended docs
  - [#11](https://github.com/nf-core/vipr/issues/11), [#14](https://github.com/nf-core/vipr/pull/14) -Font problem in logo
  - [#12](https://github.com/nf-core/vipr/issues/12), [#14](https://github.com/nf-core/vipr/pull/14) - Remove groovy to close
  - [#21](https://github.com/nf-core/vipr/pull/21) - Fix containers
  - [#20](https://github.com/nf-core/vipr/pull/20), [#21](https://github.com/nf-core/vipr/pull/21) - Update documentation
