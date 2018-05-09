# ![nf-core/vipr](docs/images/vipr_logo.png)

[![Build Status](https://travis-ci.org/nf-core/vipr.svg?branch=master)](https://travis-ci.org/nf-core/vipr)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.6-brightgreen.svg)](https://www.nextflow.io/)
[![Gitter](https://img.shields.io/badge/gitter-%20join%20chat%20%E2%86%92-4fb99a.svg)](https://gitter.im/nf-core/Lobby)
[![Docker Container available](https://img.shields.io/docker/automated/nfcore/vipr.svg)](https://hub.docker.com/r/nfcore/vipr/)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)


### Introduction

**nf-core/vipr** is a bioinformatics best-practice analysis pipeline
used for viral amplicon analysis and intrahost / low-frequency variant
calling.

The pipeline uses [Nextflow](https://www.nextflow.io), a
bioinformatics workflow tool. It pre-processes raw data from FastQ
inputs, aligns the reads and performs extensive quality-control on the
results.


### Pipeline Steps

| Step                                            | Main program/s          |
|-------------------------------------------------|-------------------------|
| Trimming and combining of read-pairs per sample | Skewer, FastQC          |
| Decontamination                                 | BWA                     |
| Metagenomics classifcation / Sample purity      | Kraken                  |
| Assembly to contigs                             | BBtools' Tadpole        |
| Assembly polishing (two processes)              | Mummer's Nucmer, LoFreq |
| Mapping to assembly                             | BWA, LoFreq             |
| Low frequency variant calling                   | LoFreq                  |
| Coverage and variant AF plots (two processes)   | Bedtools, Plotly        |



### Documentation

The nf-core/vipr pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation and configuration](docs/installation.md)
2. [Running the pipeline](docs/usage.md)
3. [Output and how to interpret the results](docs/output.md)


### Credits

These scripts were originally written for use at the [Genome Institute of Singapore](http://a-star.edu.sg/gis)

ViPR had three predecessors:
- ViPR1: FIXME
- ViPR2: FIXME
- ViPR3: FIXME

Supporting publications for older versions:
- FIXME
- FIXME


* Main author:
  * Andreas Wilm ([@ewels](https://github.com/andreas-wilm/))
