#!/bin/bash
#SBATCH -J make_flat_file       ## Name of job
#SBATCH -A b1042               ## Allocation
#SBATCH -p genomics                ## Queue
#SBATCH -t 04:00:00             ## Walltime/duration of the job
#SBATCH --mem-per-cpu=32G               ## Memory per node in GB needed for a job.
#SBATCH --cpus-per-task=1           ## Number of processors for your task

module load R/4.0.0

Rscript 20210323_MakeFlatFile.R
Rscript flat_file_gene.R
Rscript flat_file_consequence_scoring.R
Rscript add_divergent_flag.R
