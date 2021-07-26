#!/bin/bash
snakemake -s bulkHiCpipeline.smk --cluster 'sbatch --cpus-per-task={threads} -t 7-00:00 -J HiC!' --jobs 188 --resources nodes=188 --rerun-incomplete
