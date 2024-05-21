
#module load multiqc pandoc
source ~/.bashrc
CONDA=multiqc
conda activate ${CONDA}

COUNTS=$1

for d in ${COUNTS}/*/ 
do echo ${d%*/}
while true; do
missing_file=0
#check if finished rseqc
check_files=( "read_distribution.csv" "geneBodyCoverage.txt" "eRPKM.xls" "junctionSaturation_plot.pdf" "bam_stat.csv" "clipping_profile.xls" "tin.csv" )
    for check_file in "${check_files[@]}"; do
        if [ $(find ${d} -name "*${check_file}"|wc -l) -gt 0 ]; then
        printf .
        else
        echo ${check_file} does not exist in the folder ${d} 
        missing_file=1
        break
        fi
    done
    echo ""
if  [ ${missing_file} -eq 1 ]; then
sleep 60
else
break
fi

done
done

echo all required files exist
multiqc -d  ${COUNTS}
multiqc  -d  --pdf ${COUNTS}
