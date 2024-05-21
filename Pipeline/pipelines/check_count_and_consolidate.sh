#!/bin/sh
CONDA=check_count_consolidate
source ~/.bashrc
conda activate ${CONDA}




Usage="""
This script check and consolidate counts after star alignment step.\n
Usage:\n
check_count_and_consolidate.sh count_dir organism outdir.\n\n
counts dir : COUNTS directory; \n
organism: species of sample (human, mouse)\n
outdir: output directory\n
"""
if [ "$#" -eq 0 ]; then
echo -e $Usage
echo "Please enter the input"
exit
fi
main_dir=`pwd`
COUNTS=$1
organism=$2
outdir=$3

echo "organism : ${organism}"

# Checking for number of file downloaded,
cd ${COUNTS}
ls -lh  ./*/count*/gene*/*count.txt > final_count_file.txt
grep unmapped  */*/Log.final.out > unmapped.log
# total number of sample having count:
wc -l final_count_file.txt > number_of_samples_with_counts.txt


cd $main_dir
#module load 
echo "pipe post processing starting "
Rscript $SOURCE/pipe_post_processing_052918.R --cmd=consolidate --indir=${COUNTS} --outdir=${outdir} --organism=$organism

echo "done pipe post processing"
#module load python
submit_cluster $SOURCE/multiqc.sh ${COUNTS}



