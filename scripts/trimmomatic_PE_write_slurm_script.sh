#!/bin/bash
array=("7_1_GAGTGG_L006_R1_001.fastq.gz" "6_2_CGTACG_L006_R1_001.fastq.gz" "6_1_GTTTCG_L006_R1_001.fastq.gz" "5_1_GTGGCC_L006_R1_001.fastq.gz" "4_2_GGCTAC_L006_R1_001.fastq.gz" "4_1_TAGCTT_L006_R1_001.fastq.gz" "3_1_GATCAG_L006_R1_001.fastq.gz" "2_2_ACTTGA_L006_R1_001.fastq.gz" "2_1_TTAGGC_L006_R1_001.fastq.gz" "1_2_ATCACG_L006_R1_001.fastq.gz")
for SAMPLE_ID_R1 in "${array[@]}"
do
	SAMPLE_ID_R2=${SAMPLE_ID_R1%%R1_001.fastq.gz}"R2_001.fastq.gz"
#set a prefix to make understanding what has been done to the file easier
        PREFIX=trimmed_
        PREFIX_SE=SE_

        #create the new file name (e.g. trimmedgly7a.fq.gz)
        NEWTRIMFILE_R1=${PREFIX}$SAMPLE_ID_R1
        NEWTRIMFILE_R2=${PREFIX}$SAMPLE_ID_R2
        SE_TRIMFILE_R1=${PREFIX_SE}$SAMPLE_ID_R1
        SE_TRIMFILE_R2=${PREFIX_SE}$SAMPLE_ID_R2
                
cat > /projects/researchers/researchers01/jonbra/Beroe/scripts/trimmomatic_${SAMPLE_ID_R1}.slurm <<EOF
#!/bin/sh
#SBATCH --job-name=trim
#SBATCH --account=uio
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --output=slurm-%j.base
## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit
#do the trimming
java -jar ~/Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 -threads 8 ../rawdata/$SAMPLE_ID_R1 ../rawdata/$SAMPLE_ID_R2 ../QC/Trimmomatic/$NEWTRIMFILE_R1 ../QC/Trimmomatic/$SE_TRIMFILE_R1 ../QC/Trimmomatic/$NEWTRIMFILE_R2 ../QC/Trimmomatic/$SE_TRIMFILE_R2 ILLUMINACLIP:/usit/abel/u1/jonbra/Trimmomatic-0.35/adapters/TruSeq3-SE.fa:2:28:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:4:28 MINLEN:36
EOF

sbatch /projects/researchers/researchers01/jonbra/Beroe/scripts/trimmomatic_${SAMPLE_ID_R1}.slurm

done

