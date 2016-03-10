#!/bin/sh
#SBATCH --job-name=cpc
#SBATCH --account=uio
#SBATCH --output=slurm-%j.base
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=32G

module load blast
bash ~/cpc-0.9-r2/bin/run_predict.sh merged_fixed_no_gene_name_200nt_150nt_unique_tcons.fasta cpc_out /projects/researchers/researchers01/jonbra/MULTICELL/Andreas/potential_lncRNA/second_try evidence_out
