#!/bin/bash
for f1 in *R1_001.fastq.gz
do
	f2=${f1%%R1_001.fastq.gz}"R2_001.fastq.gz"

        #set a prefix to make understanding what has been done to the file easier
        PREFIX=trimmed_
        PREFIX_SE=SE_

        #create the new file name (e.g. trimmedgly7a.fq.gz)
        NEWTRIMFILE_R1=${PREFIX}$f1
        NEWTRIMFILE_R2=${PREFIX}$f2
        SE_TRIMFILE_R1=${PREFIX_SE}$f1
        SE_TRIMFILE_R2=${PREFIX_SE}$f2

        #do the trimming
	java -jar ~/Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 -threads 8 $f1 $f2 ../QC/Trimmomatic/$NEWTRIMFILE_R1 ../QC/Trimmomatic/$SE_TRIMFILE_R1 ../QC/Trimmomatic/$NEWTRIMFILE_R2 ../QC/Trimmomatic/$SE_TRIMFILE_R2 ILLUMINACLIP:/usit/abel/u1/jonbra/Trimmomatic-0.35/adapters/TruSeq3-SE.fa:2:28:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:4:28 MINLEN:36

done

