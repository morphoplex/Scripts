#!/bin/sh
#SBATCH --job-name=sratools
#SBATCH --account=uio
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --output=slurm-%j.base

## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit

/usit/abel/u1/jonbra/sratoolkit.2.5.5-ubuntu64/bin/fastq-dump --split-3 -O /usit/abel/u1/jonbra/Sarctica_RNA_revision/mapping --gzip --accession SRR350209

