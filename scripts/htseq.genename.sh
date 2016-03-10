#!/bin/bash
array=("1-1-1-1_S1_L00X" "1-2-1-3_S2_L00X" "1-3-2-1_S3_L00X" "1-4-2-3_S4_L00X" "1-5-3-1_S5_L00X" "1-6-3-3_S6_L00X" "1-7-4-1_S7_L00X" "1-8-4-3_S8_L00X" "1-4B_S7" "2-5B_S10" "3-6B_S13" "4-7B_S16" "5-8B_S19" "6-10B_S22" "7-13B_S1" "8-15B_S4" "9-17B_S8" "10-21B_S11" "11-12M_S14" "12-13M_S17" "13-14M_S20" "15-17M_S23" "16-3T_S2" "17-4T_S5" "18-7T_S9" "19-8T_S12" "20-9T_S15" "21-10T_S18" "22-11T_S21" "23-12T_S24" "24-13T_S3" "25-17T_S6")
for SAMPLE_ID in "${array[@]}"
do
cat > /projects/researchers/researchers01/jonbra/Mnemiopsis/scripts/htseq_${SAMPLE_ID}.slurm <<EOF
#!/bin/sh
#SBATCH --job-name=htseq
#SBATCH --account=uio
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --output=slurm-%j.base

## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit

module load python2
module load samtools

samtools sort -O bam -T /projects/researchers/researchers01/jonbra/Mnemiopsis/tophat/${SAMPLE_ID}/${SAMPLE_ID} -n /projects/researchers/researchers01/jonbra/Mnemiopsis/tophat/${SAMPLE_ID}/accepted_hits.bam -o /projects/researchers/researchers01/jonbra/Mnemiopsis/tophat/${SAMPLE_ID}/accepted_hits_sorted.bam

python -m HTSeq.scripts.count -f bam -r name -s reverse -t exon -i gene_name /projects/researchers/researchers01/jonbra/Mnemiopsis/tophat/${SAMPLE_ID}/accepted_hits_sorted.bam /projects/researchers/researchers01/jonbra/Mnemiopsis/assemblies/MergedAssemblies/merged_fixed_genename_only.gtf > /projects/researchers/researchers01/jonbra/Mnemiopsis/assemblies/HTSeq_counts_genename/HTSeq_${SAMPLE_ID}.txt

rm /projects/researchers/researchers01/jonbra/Mnemiopsis/tophat/${SAMPLE_ID}/accepted_hits_sorted.bam

EOF

sbatch /projects/researchers/researchers01/jonbra/Mnemiopsis/scripts/htseq_${SAMPLE_ID}.slurm

done

