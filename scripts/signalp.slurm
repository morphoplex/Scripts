#!/bin/sh
#SBATCH --job-name=signalP
#SBATCH --account=uio
#SBATCH --output=slurm-%j.base
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=8G

~/signalp-4.1/signalp merged_fixed_no_gene_name_200nt_150nt_unique_tcons.fasta.transdecoder.pep
