#!/bin/bash
#SBATCH -J ncsq224       ## Name of job
#SBATCH -A b1042               ## Allocation
#SBATCH -p genomics                ## Queue
#SBATCH -t 02:00:00             ## Walltime/duration of the job
#SBATCH --mem-per-cpu=32G               ## Memory per node in GB needed for a job.
#SBATCH --cpus-per-task=1            ## Number of processors for your task

module purge all

conda activate bcftools1_12

bcftools csq --ncsq 224 -f c_elegans.PRJNA13758.WS276.genomic.fa -g c_elegans.PRJNA13758.WS276.annotations.NEW_CSQ.gff3  -p a -Ov -o WI.20210121.hard-filter.isotype.re_an224.vcf.gz WI.20210121.hard-filter.isotype.pre_an.vcf.gz
