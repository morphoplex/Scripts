#!/bin/sh
#SBATCH --job-name=fastqc
#SBATCH --account=uio
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --output=slurm-%j.base

## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit

module load fastqc

for i in ../SequenceData/RawSeqs/*.fastq.gz
do
        fastqc -o ../Analyses/FastQC/ $i
done

