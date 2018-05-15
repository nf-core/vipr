# nf-core/vipr Usage

Please refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for generic Nextflow options, like `-resume` etc.

## General Nextflow info

Nextflow handles job submissions on compute environments, and
supervises running the jobs. Thus the Nextflow process must run until
the pipeline is finished. We recommend that you put the process
running in the background through `screen` / `tmux` or similar
tool. Alternatively you can run nextflow within a cluster job
submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines
memory. We recommend adding the following line to your environment
(typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run nf-core/vipr -params-file params.yaml -profile docker
```

This will launch the pipeline with the `docker` configuration profile
using input paramters as defined in `params.yaml`. See below for more
information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the
pipeline code from GitHub and stores it as a cached version. When
running the pipeline after this, it will always use the cached version
if available - even if the pipeline has been updated since. To make
sure that you're running the latest version of the pipeline, make sure
that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/vipr
```

### Reproducibility

It's a good idea to specify a pipeline version when running the
pipeline on your data. This ensures that a specific version of the
pipeline code and software are used when you run your pipeline. If you
keep using the same tag, you'll be running the same version of the
pipeline, even if there have been changes to the code since.

First, go to the
[nf-core/vipr releases page](https://github.com/nf-core/vipr/releases)
and find the latest version number - numeric only (e.g. `1.0`). Then
specify this when running the pipeline with `-r` (one hyphen) -
eg. `-r 1.0`.

This version number will be logged in reports when you run the
pipeline, so that you'll know what you used when you look back in the
future.


## Main Arguments

### `-params-file`

You can modify program behaviour and specify input files in a yaml
configuration file. An example is given in `example_params.yaml`.

Please note: this is currently the only way to specify read input.
The corresponding entries in params.yaml looks as follows:

```nextflow
params {
  samples:
    sample-name-1:
      readunits:
        readunit-1:
          fq1: path-to-R1.fastq.gz
          fq2: path-to-R2.fastq.gz
        ...
        readunit-n:
          fq1: path-to-R1.fastq.gz
          fq2: path-to-R2.fastq.gz
    ...
    sample-name-n:
      readunits:
        ..
}
```

So you can specify multiple samples and each samples can contain multiple fastq pairs (AKA readunits)


### `-profile`
Use this parameter to choose a configuration profile. Each profile is designed for a different compute environment - follow the links below to see instructions for running on that system. Available profiles are:

* `standard`
    * The default profile, used if `-profile` is not specified at all. Runs locally and expects all software to be installed and available on the `PATH`.
    * This profile is mainly designed to be used as a starting point for other configurations and is inherited by most of the other profiles.
* `none`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config profile (not recommended).

See `nextflow.config` for more available profiles

### `--skip-kraken`

Skips the optional [Kraken](https://ccb.jhu.edu/software/kraken/) metagenomics classifaction of your reads.



## Reference Genomes

The pipeline requires you to specify close reference for samples
(`--ref-fasta`, e.g. Zika) and a fasta reference for decontamination
(`--cont-fasta` e.g. human) If you in addition also use Kraken (see
above) you will need to specify the path to your Kraken database (`--kraken-db`).

You can add the above parameters to your params-file (see above). Then entries looks as follows
```nextflow
params {
  ref_fasta: "data/ref/DENV2-NC_001474.2.fa"
  cont_fasta: "/mnt/projects/rpd/genomes/human_g1k_v37/human_g1k_v37.fasta"
  kraken_db: "/mnt/projects/rpd/genomes/kraken/minikraken_20171019_8GB/"
  skip_kraken: false
  }
}
```

An example is given in `example_params.yaml`.

### `--outdir`
The output directory where the results will be saved.

## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples.
