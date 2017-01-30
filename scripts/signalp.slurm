#!/bin/sh
#SBATCH --job-name=signalP
#SBATCH --account=uio
#SBATCH --output=slurm-%j.base
##SBATCH --partition=hugemem
#SBATCH --time=100:00:00
##SBATCH --ntasks=4
#SBATCH --mem-per-cpu=16G

module purge

perl ~/signalp-4.1/signalp -f short -n signalp.out ../trinity/Trinity.longest.isoform.cut.fasta.transdecoder_dir/longest_orfs.pep
