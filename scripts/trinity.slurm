#!/bin/sh
#SBATCH --job-name=trinity
#SBATCH --account=uio
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=61G
#SBATCH --output=slurm-%j.base

## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit

module load trinityrnaseq

Trinity --seqType fq --left ../QC/Trimmomatic/trimmed_6_1_GTTTCG_L006_R1_001.fastq.gz,../QC/Trimmomatic/trimmed_7_1_GAGTGG_L006_R1_001.fastq.gz,../QC/Trimmomatic/trimmed_6_2_CGTACG_L006_R1_001.fastq.gz,../QC/Trimmomatic/trimmed_5_1_GTGGCC_L006_R1_001.fastq.gz,../QC/Trimmomatic/trimmed_4_2_GGCTAC_L006_R1_001.fastq.gz,../QC/Trimmomatic/trimmed_4_1_TAGCTT_L006_R1_001.fastq.gz,../QC/Trimmomatic/trimmed_1_2_ATCACG_L006_R1_001.fastq.gz,../QC/Trimmomatic/trimmed_2_1_TTAGGC_L006_R1_001.fastq.gz,../QC/Trimmomatic/trimmed_2_2_ACTTGA_L006_R1_001.fastq.gz,../QC/Trimmomatic/trimmed_3_1_GATCAG_L006_R1_001.fastq.gz --right ../QC/Trimmomatic/trimmed_6_1_GTTTCG_L006_R2_001.fastq.gz,../QC/Trimmomatic/trimmed_7_1_GAGTGG_L006_R2_001.fastq.gz,../QC/Trimmomatic/trimmed_6_2_CGTACG_L006_R2_001.fastq.gz,../QC/Trimmomatic/trimmed_5_1_GTGGCC_L006_R2_001.fastq.gz,../QC/Trimmomatic/trimmed_4_2_GGCTAC_L006_R2_001.fastq.gz,../QC/Trimmomatic/trimmed_4_1_TAGCTT_L006_R2_001.fastq.gz,../QC/Trimmomatic/trimmed_1_2_ATCACG_L006_R2_001.fastq.gz,../QC/Trimmomatic/trimmed_2_1_TTAGGC_L006_R2_001.fastq.gz,../QC/Trimmomatic/trimmed_2_2_ACTTGA_L006_R2_001.fastq.gz,../QC/Trimmomatic/trimmed_3_1_GATCAG_L006_R2_001.fastq.gz --SS_lib_type RF --max_memory 61G --CPU 8 --output ../trinity/ --normalize_reads

