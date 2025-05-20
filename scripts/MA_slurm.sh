#!/bin/bash
#SBATCH -p emergence
#SBATCH -o %x.%N.%j.out
#SBATCH -e %x.%N.%j.err
#SBATCH --cpus-per-task=1
#SBATCH --job-name=BacWork

snakemake -s workflow.smk  -k \
--use-conda --latency-wait 30 --resources load=100 \
--cluster "sbatch --cpus-per-task={threads} --job-name={rule}-{wildcards} -p prod" \
--jobs 8  --conda-prefix /global/bio/anaconda3/ \
--conda-frontend mamba 


snakemake -s workflow.smk  -k \
--use-conda --latency-wait 30 --resources load=100 \
--cluster "sbatch --cpus-per-task={threads} --job-name={rule}-{wildcards} -p prod" \
--jobs 8  --conda-prefix /global/bio/anaconda3/bacwork/ \
--conda-frontend mamba 