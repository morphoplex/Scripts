#!/bin/sh
#SBATCH --job-name=trinotate
#SBATCH --account=uio
#SBATCH --output=slurm-%j.base
##SBATCH --partition=hugemem
#SBATCH --time=10:00:00
##SBATCH --ntasks=4
#SBATCH --mem-per-cpu=8G
module purge
module load perlmodules

~/Trinotate-3.0.0/Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls
