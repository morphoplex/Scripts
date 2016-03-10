#!/bin/sh
#SBATCH --job-name=hmmer
#SBATCH --account=uio
#SBATCH --output=slurm-%j.base
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=8G

module load hmmer

hmmscan --tblout /projects/researchers/researchers01/jonbra/MULTICELL/Andreas/potential_lncRNA/second_try/hmmer_PfamB.out -E 0.0001 /usit/abel/u1/jonbra/Pfam-B.hmm /projects/researchers/researchers01/jonbra/MULTICELL/Andreas/potential_lncRNA/second_try/merged_fixed_no_gene_name_200nt_150nt_unique_tcons.fasta.transdecoder.pep
