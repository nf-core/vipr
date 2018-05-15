# Troubleshooting

## Pipeline run completed but output files seem to be missig

Under certain circumstance the assembled contigs cannot be oriented or
aligned against the given reference. In such cases, you will only get
the contigs (`*_contigs.fa`) but the downstream analysis stops,
i.e. no polished assembly is created or variant calls are made. This
can for example happen if you picked the wrong reference for your
sample. For troubleshooting, take the contigs and BLAST them against
NR and also have a look at the Kraken classifcation output
`*_kraken.report` to see whether their is a mismatch between the
contig classification and the user provided reference.


## Getting help

If you have an issue with running the pipeline then feel free to
contact us.  Have look at the
[issue tracker for our repo](https://github.com/nf-core/vipr/issues). Maybe
someone has already had the same problem?

If you have problems that are related to Nextflow and not our pipeline
then check out the
[Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or
the [google group](https://groups.google.com/forum/#!forum/nextflow).
