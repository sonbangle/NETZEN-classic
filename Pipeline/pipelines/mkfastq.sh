#!/bin/bash
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=1000G
#SBATCH --output=mkfastq_HNJ2MDSX3_%j.log
#SBATCH --error=mkfastq_HNJ2MDSX3_%j.log
#SBATCH --qos=dtran
#SBATCH --job-name=HNJ2MDSX3

## Convert BCL files to fastq. Assumes Illumina run dir is current directory. This is for 5' Dual Index Library. It will affect the base mask

RUNID=fastq_HNJ2MDSX3  
#the path of illumina bcl run folder
RUNDIR=/orange/dtran/sequencing/NS2738/HNJ2MDSX3/	 	
#output directory
OUTDIR=.
#Base Mask
BM=Y26N*,I10N*,I10N*,Y151
#Sample Sheet
SSHEET=SampleSheet-HNJ2MDSX3.csv


module purge
module load cellranger/7.0.0
module list





echo "Running mkfastq with following params:"

echo "RUNDIR: ${RUNDIR}, RUNID:${RUNID}, samplesheet:${SSHEET}, base mask:${BM}"

cellranger mkfastq \
  --run="$RUNDIR" \
  --csv=$SSHEET \
  --id=$RUNID \
  --localcores=128 \
  --localmem=900 \
  --use-bases-mask=$BM \
  --output-dir=$OUTDIR \



