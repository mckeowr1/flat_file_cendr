#!/bin/bash
#SBATCH -A b1042               # Allocation
#SBATCH --partition=genomics                # Queue
#SBATCH -t 01:00:00             # Walltime/duration of the job
#SBATCH --mem-per-cpu=8G               ## Memory per node in GB needed for a job.
#SBATCH --cpus-per-task=1           ## Number of processors for your task

module purge all

conda activate bcftools1_12

bcftools view -e 'INFO/BCSQ="."' -Ou /projects/b1059/analysis/WI-20210121/isotype_only/WI.20210121.hard-filter.isotype.bcsq.20210401.vcf.gz  | bcftools query -i 'GT ="alt"' -f'%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE:%TBCSQ{*}=]\n' > /projects/b1059/projects/Ryan/csq/flat_file/WI.20210121.hard-filter.isotype.bcsq.20210401_strain_bcsq.tsv
