#!/bin/bash
#$ -cwd
for i in *R1*
do
R2=${i//R1_001.fastq.gz/R2_001.fastq.gz}
FOLDER=${i//_S*_L004_R1_001.fastq.gz}
LOG=${i//_S*_L004_R1_001.fastq.gz/.log}

kallisto quant -i Mus_musculus.GRCm39.cdna.all.index -o $FOLDER -t 48 --genomebam --gtf Mus_musculus.GRCm39.105.gtf.gz $i $R2 &> $LOG

done
      
