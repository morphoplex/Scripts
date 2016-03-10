#!/bin/sh
#SBATCH --job-name=blastx
#SBATCH --account=uio
#SBATCH --output=slurm-%j.base
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=16


module load blast+/2.2.29

blastx -query new_assembled_transcripts_potential_lncRNA.fa -evalue 0.0001 -db /work/databases/bio/ncbi-data13042015/refseq_protein -outfmt 5 -out new_assembled_transcripts_potential_lncRNA_blastx_refseq_protein_1E-4.xml -max_target_seqs 1 -num_threads $SLURM_CPUS_PER_TASK

