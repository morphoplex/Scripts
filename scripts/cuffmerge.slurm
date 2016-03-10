#!/bin/sh
#SBATCH --job-name=tophat
#SBATCH --account=uio
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --output=slurm-%j.base

## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit

module load cufflinks

cuffmerge -p 8 -g /projects/researchers/researchers01/jonbra/Mnemiopsis/genome_transcriptome/ML2.2.nogene.gff3 -s /projects/researchers/researchers01/jonbra/Mnemiopsis/bowtie2_index/Mlgenome.fa -o /projects/researchers/researchers01/jonbra/Mnemiopsis/assemblies/MergedAssemblies/ /projects/researchers/researchers01/jonbra/Mnemiopsis/assemblies/assemblies.txt 
