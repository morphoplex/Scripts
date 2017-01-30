#!/bin/sh
#SBATCH --job-name=Blastx
#SBATCH --account=uio
#SBATCH --output=slurm-%j.base
##SBATCH --partition=hugemem
#SBATCH --time=100:00:00
##SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3G

module purge
module load blast+/2.2.29

blastx -query ../trinity/Trinity.longest.isoform.cut.fasta -outfmt 6 -evalue 1e-20 -out blastx.outfmt6 -db /work/databases/bio/ncbi-data/swissprot -num_descriptions 1 -num_alignments 1 -num_threads $SLURM_CPUS_PER_TASK

