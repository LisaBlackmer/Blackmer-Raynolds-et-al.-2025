#!/bin/bash


#SBATCH -N 1
#SBATCH --cpus-per-task 40
#SBATCH -o /home/pchopr4/logs/cellbender.o-%j
#SBATCH -e /home/pchopr4/logs/cellbender.e-%j
#SBATCH --time=72:00:00
#SBATCH --array=1-3

## ref: https://davetang.org/muse/2018/08/09/getting-started-with-cell-ranger/
##  derived from xmo.sh
## sbatch array.xxx.sh
## created by Pankaj Chopra


source activate CellBender

echo "cellbender version:"
cellbender -v

### Project: 21173FL-04 (Lisa - mouse; 14 May, 2024) - 15 samples (set --array=1-15)
#DDIR=/home/pchopr4/data/scrna/ssloan/lisa.blackmer-raynolds.march24/fastq
#SFILE=${DDIR}/sampleids.txt

## lisa: Re-do 3 problematic samples using long fastq reads (i.e. not trimmed by florida U.
DDIR=/home/pchopr4/data/scrna/ssloan/lisa.blackmer-raynolds.march24/fastq
SFILE=${DDIR}/sampleids.3.txt

MYID=${SLURM_ARRAY_TASK_ID}
#MYFASTQ=`sed '"${MYID}"q;d' ${SFILE}`
MYFASTQ=$(sed "${MYID}q;d" "${SFILE}")

#MYFASTQ="21173FL-04-03"
#FASTQDIR=${DDIR}/${MYFASTQ}/
#FASTQDIR=${DDIR}/${MYFASTQ}.paired.trimmed/


#FILTER="none"

echo "${SFILE}"
echo "${MYID}"
echo "${MYFASTQ}"

cd ${DDIR}
ls -lh ${MYFASTQ}

cellranger --version
cellbender -v

#exit 0

### 1. Run CellRanger
#cellranger count --id=${MYFASTQ}_out \
#	--transcriptome=/home/pchopr4/data/references/refdata-gex-mm10-2020-A \
#	--fastqs=${MYFASTQ} \
#	--expect-cells=10000

DDIR2=${DDIR}/${MYFASTQ}_out/outs

### 2. Run cellbender on output from CellRanger above
cd ${DDIR2}
rm *cellbender*
echo pwd


## For cuda/GPU submit job with 'sbatch -p gpu -w gpu01 cellbender.sh'
## not working (cuda is NVIDIA while skynet GPU is intel) If using GPU: sbatch -p gpu -w gpu02 cellbinder.sh (single node)
#cellbender remove-background --cuda --fpr 0.05 0.1 --input raw_feature_bc_matrix.h5 --output ${DDIR2}/${MYFASTQ}.cellbender.h5


cellbender remove-background --input raw_feature_bc_matrix.h5 \
	--expected-cells 5500 --total-droplets-included 12000 \
	--output ${DDIR2}/${MYFASTQ}.cellbender.h5


