#!/bin/bash
#SBATCH -p bioinfo
#SBATCH -o %x.%N.%j.out
#SBATCH -e %x.%N.%j.err
#SBATCH --cpus-per-task=1
#SBATCH --job-name=BacWork

mkdir -p logs

/nfs/conda/envs/snakemake_7.25.0/bin/snakemake -s workflow.smk  -k \
--use-conda --latency-wait 60 --resources load=100 \
--cluster "sbatch --cpus-per-task={threads} --job-name={rule}-{wildcards} -p $1 --output=logs/%j-{rule}.out" \
--jobs 16  --conda-prefix /nfs/kebab \
--conda-frontend mamba 