#!/bin/bash
#SBATCH -A b1042               # Allocation
#SBATCH --partition=genomics                # Queue
#SBATCH -t 01:00:00             # Walltime/duration of the job
#SBATCH --mem-per-cpu=8G               ## Memory per node in GB needed for a job.
#SBATCH --cpus-per-task=1           ## Number of processors for your task

module purge all

module load bcftools/1.10.1 

bcftools view -e 'INFO/BCSQ="."' -Ou /projects/b1059/projects/Ryan/csq/GitreanWI.20210121.hard-filter.isotype.bcsq.vcf.gz | bcftools query -i 'GT ="alt"' -f'%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE:%TBCSQ{*}=]\n' > /projects/b1059/projects/Sophie/csq/re-an_strain_BCSQ.tsv
