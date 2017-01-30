#!/bin/sh
#SBATCH --job-name=hmm
#SBATCH --account=uio
#SBATCH --output=slurm-%j.base
##SBATCH --partition=hugemem
#SBATCH --time=100:00:00
##SBATCH --ntasks=4
#SBATCH --mem-per-cpu=16G

module purge
module load hmmer

hmmscan --cpu 8 --domtblout TrinotatePFAM.out /work/databases/bio/trinotate/Pfam-A.hmm ../trinity/Trinity.longest.isoform.cut.fasta.transdecoder_dir/longest_orfs.pep > pfam.log

