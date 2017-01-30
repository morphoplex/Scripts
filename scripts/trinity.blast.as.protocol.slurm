#!/bin/sh
#SBATCH --job-name=Blastx
#SBATCH --account=uio
#SBATCH --output=slurm-%j.base
##SBATCH --partition=hugemem
#SBATCH --time=100:00:00
##SBATCH --ntasks=4
#SBATCH --mem-per-cpu=6G

module purge
module load blast+/2.2.29

blastx -query ../trinity/Trinity.longest.isoform.cut.fasta -outfmt 6 -evalue 1e-20 -out blastx.protocol.outfmt6 -db /work/databases/bio/trinotate/uniprot_sprot.trinotate.pep -num_threads 8 -max_target_seqs 1

