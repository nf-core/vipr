**Under Construction**

Tests are implemented in Snakemake (which thus needs to be installed).

Run with e.g.: `snakemake --rerun-incomplete --timestamp --printshellcmds --latency-wait 60 --max-jobs-per-second 1 --max-status-checks-per-second 0.1 --cores 3`

This will run Nextflow on multiple data-sets (see profile setting there) and subsequently test the output
