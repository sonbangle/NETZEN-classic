#!/bin/bash
export DISPLAY=:0.0

Usage="""Usage: $(basename $0) sampleID  organism outdir ncpu remove_bam fastq1 fastq2 ; organism can be mouse or human \n
sampleID: sampleID \n
organism: species of sample \n
outdir:  output directory for this sample
ncpu: number of cpu for alignment
remove_bam: removing bam file after alignment. Take value 0 : not removing, 1: removing

fastq1: fastq file
fastq2: fastq file
"""
if [ "$#" -eq 0 ]; then
echo -e $Usage
echo "Please enter the input"
exit

fi



#initial values

do_BedGraph=1
do_RseQC=1
do_qualimap=1
remove_bam=$5
CONDA=bulkRNAseq

source ~/.bashrc

conda activate ${CONDA}

n_end=2
#organism: name of species for reference genome, can be mouse or human

organism=$2

echo "MY_BIN:${MY_BIN}"

if [ "$organism" == "mouse" ] ; then
echo "mouse"
genomeDir=${MY_BIN}/STAR/GRCm38_genecodevM11
sjdbGTFfile=${MY_BIN}/STAR/GRCm38_genecodevM11/gencode.vM11.primary_assembly.annotation.gtf
elif [ "$organism" == "human" ] ; then
echo "human"
genomeDir=${MY_BIN}/STAR/human_Genecode28/Genomedir
sjdbGTFfile=${MY_BIN}/STAR/human_Genecode28/gencode.v28.primary_assembly.annotation.gtf
fi 

echo "genomeDir:${genomeDir}"
module purge
export PATH=$PATH:$MY_BIN


let "n_end = $# - 5"
fastq1=$(readlink -f $6)
echo $fastq1
if [ "$n_end" -eq 1 ]; then
fastq_files=($fastq1)
else
fastq2=$(readlink -f $7)
echo $fastq2
fastq_files=($fastq1 $fastq2)
fi


outdir=$3
sampleID=$1
mkdir $outdir
cd $outdir
mkdir $sampleID
cd $sampleID
ncpu=$4


#quality control ==================================
echo "beginning fasqc before trimming at $(date "+%Y-%m-%d %H:%M:%S")"
module load gcc/11.3.0 fastqc/0.11.9
mkdir BEFORE_TRIM

echo "ncpu: ${ncpu}"

if [ $(find ./BEFORE_TRIM -name "*_fastqc.html"|wc -l)  -gt 0 ] ; then 
    echo " .fastqc.html exists "
else
echo "./BEFORE_TRIM/${fastq_files[0]}_fastqc.html does not exists"
echo current path: $(pwd)
fastqc -t $ncpu  -o  ./BEFORE_TRIM ${fastq_files[0]}
if [ $n_end -eq 2 ] ; then
fastqc -t $ncpu  -o  ./BEFORE_TRIM ${fastq_files[1]}
fi
fi
echo "done fasqc before trimming at $(date "+%Y-%m-%d %H:%M:%S")"

#trimming =======================================================q
echo "beginning  trimming at $(date "+%Y-%m-%d %H:%M:%S")"
#module load trimmomatic/0.39
echo "number of read end: $n_end" 
if [ $n_end -eq 2 ]
then
if [ -f "${sampleID}_forward_paired.fq.gz" ]; then
    echo "${sampleID}_forward_paired.fq.gz exits "
else

trimmomatic PE -threads $ncpu  ${fastq_files[0]}  ${fastq_files[1]}  ${sampleID}_forward_paired.fq.gz        ${sampleID}_forward_unpaired.fq.gz ${sampleID}_reverse_paired.fq.gz    ${sampleID}_reverse_unpaired.fq.gz \
    ILLUMINACLIP:/spack/2206/apps/linux-centos7-x86_64_v3/gcc-11.3.0/trimmomatic-0.39-nughsr4/share/adapters/TruSeq3-PE-2.fa:2:40:15   LEADING:3   TRAILING:3   SLIDINGWINDOW:4:15 MINLEN:20
 #This will perform the following:

          #Remove adapters
          #Remove leading low quality or N bases (below quality 3)
          #Remove trailing low quality or N bases (below quality 3)
          #Scan the read with a 4-base wide sliding window, cutting when the average quality #per base drops below 15
          # Drop reads below the 20 bases long
rm ${sampleID}_forward_unpaired.fq.gz ${sampleID}_reverse_unpaired.fq.gz
fi
elif [ $n_end -eq 1 ]
then
echo "single end"
if [ -f "${sampleID}_forward.fq.gz" ]; then
    echo "${sampleID}_forward.fq.gz exits"
else
trimmomatic SE  -threads $ncpu   ${fastq_files[0]}   ${sampleID}_forward.fq.gz    ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-SE.fa:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
fi
fi

echo "done  trimming at $(date "+%Y-%m-%d %H:%M:%S")"
#rm ${fastq_files[0]}  ${fastq_files[1]} 

#quality control after trimming ==================================
mkdir AFTER_TRIM
echo "beginning fasqc after trimming at $(date "+%Y-%m-%d %H:%M:%S")"
if [ $(find ./AFTER_TRIM -name "*_fastqc.html"|wc -l)  -gt 0 ] ; then 
echo " AFTER_TRIM/fastqc file exits"
else
if [ $n_end -eq 2 ] ; then
fastqc -t $ncpu -o ./AFTER_TRIM *_paired.fq.gz 
echo 1
elif [ $n_end -eq 1 ]
then
echo "single end"

fastqc  -t $ncpu -o ./AFTER_TRIM ${sampleID}_forward.fq.gz
fi
fi
echo "done fasqc after trimming at $(date "+%Y-%m-%d %H:%M:%S")"



#STAR aligment ============================================================================

echo "beginning STAR alignment at $(date "+%Y-%m-%d %H:%M:%S")"
#module load star/2.7.6a
ulimit -n 65535
echo "Current path: $(pwd)"
if [ -f "Aligned.sortedByCoord.out.bam.bai" ]; then
echo " Aligned.sortedByCoord.out.bam.bai from STAR already exits, skipping alignment"
else
if [ $n_end -eq 2 ]
then
STAR \
--runThreadN $ncpu \
--genomeLoad NoSharedMemory \
--genomeDir $genomeDir \
--sjdbGTFfile $sjdbGTFfile --sjdbOverhang 100 \
--readFilesIn ${sampleID}_forward_paired.fq.gz  ${sampleID}_reverse_paired.fq.gz  \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate Unsorted \
--outBAMsortingThreadN  24  \
--outFilterMultimapNmax 20 \
--outReadsUnmapped Fastx \
--alignMatesGapMax 4000 \
--limitBAMsortRAM 90000000000 \
--chimSegmentMin 20 \
--outSAMstrandField intronMotif \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--outFilterMatchNmin 0 \
--quantMode GeneCounts 


elif [ $n_end -eq 1 ]
then
STAR \
--runThreadN $ncpu \
--genomeLoad NoSharedMemory \
--genomeDir $genomeDir \
--sjdbGTFfile $sjdbGTFfile --sjdbOverhang 100 \
--readFilesIn ${sampleID}_forward.fq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate Unsorted \
--outBAMsortingThreadN  24 \
--outFilterMultimapNmax 20 \
--outReadsUnmapped Fastx \
--alignMatesGapMax 40000 \
--limitBAMsortRAM 90000000000 \
--chimSegmentMin 20 \
--outSAMstrandField intronMotif \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--outFilterMatchNmin 0 \
--quantMode GeneCounts  

 
fi

mkdir ../star_readcount
cut -f 1,2 ReadsPerGene.out.tab > ../star_readcount/$1.ReadsPerGene.csv


echo "end of STAR alignment at $(date "+%Y-%m-%d %H:%M:%S")"

if [ ${do_BedGraph} -gt 0 ]; then

# STAR --generate BedGraph file bg for visualization on IGV browser ======================================================================
echo "generating BedGraph  file at $(date "+%Y-%m-%d %H:%M:%S")"
STAR \
--runMode inputAlignmentsFromBAM \
--runThreadN $ncpu \
--inputBAMfile Aligned.sortedByCoord.out.bam \
--outWigType bedGraph  --outWigStrand Unstranded
fi

#generate BAI sorted and indexed for IGV viewer:
echo "generating BAM indexed file at $(date "+%Y-%m-%d %H:%M:%S")"
#module load samtools/1.17
samtools index Aligned.sortedByCoord.out.bam

fi


# QC after alignment by RseQC ===============================================================================================
#module load R
#module load python
echo before RseqQC: organism: $organism
if [ ${do_RseQC} -gt 0 ]; then
if [ -f "${sampleID}.junctionSaturation_plot.pdf" ]; then
echo "RSQC already done, junctionSaturation_plot.pdf exist"
else

echo "beginning RseQC at $(date "+%Y-%m-%d %H:%M:%S")"
#module load rseqc
if [ $organism == "mouse" ] ; then
bed_file=${MY_BIN}/STAR/GRCm38_genecodevM11/GRCm38_EnsemblGenes.bed
elif [ $organism == "human" ] ; then
bed_file=${MY_BIN}/STAR/human_Genecode28/hg38_RefSeq.bed
fi
echo bed_file: $bedfile

# submit_cluster --modules=rseqc,R --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf read_distribution.py -r ${bed_file} -i Aligned.sortedByCoord.out.bam > read_distribution.csv
# submit_cluster --modules=rseqc,R --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf geneBody_coverage.py -r ${bed_file} -i Aligned.sortedByCoord.out.bam -o $1
# submit_cluster --modules=rseqc,R --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf RPKM_saturation.py -r ${bed_file} -i Aligned.sortedByCoord.out.bam  -o $1
# submit_cluster --modules=rseqc,R --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf junction_annotation.py -r ${bed_file} -i Aligned.sortedByCoord.out.bam -o $1
# submit_cluster --modules=rseqc,R --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf junction_saturation.py -r ${bed_file} -i Aligned.sortedByCoord.out.bam -o $1
# submit_cluster --modules=rseqc,R --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf bam_stat.py  -i Aligned.sortedByCoord.out.bam > bam_stat.csv
# submit_cluster --modules=rseqc,R --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf clipping_profile.py -i Aligned.sortedByCoord.out.bam -s "PE" -o $1
# submit_cluster --modules=rseqc,R --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf tin.py -r ${bed_file} -i Aligned.sortedByCoord.out.bam > tin.csv


submit_cluster --conda=${CONDA} --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf read_distribution.py -r ${bed_file} -i Aligned.sortedByCoord.out.bam > read_distribution.csv
submit_cluster --conda=${CONDA} --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf geneBody_coverage.py -r ${bed_file} -i Aligned.sortedByCoord.out.bam -o $1
submit_cluster --conda=${CONDA} --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf RPKM_saturation.py -r ${bed_file} -i Aligned.sortedByCoord.out.bam  -o $1
submit_cluster --conda=${CONDA} --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf junction_annotation.py -r ${bed_file} -i Aligned.sortedByCoord.out.bam -o $1
submit_cluster --conda=${CONDA} --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf junction_saturation.py -r ${bed_file} -i Aligned.sortedByCoord.out.bam -o $1
submit_cluster --conda=${CONDA} --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf bam_stat.py  -i Aligned.sortedByCoord.out.bam > bam_stat.csv
submit_cluster --conda=${CONDA} --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf clipping_profile.py -i Aligned.sortedByCoord.out.bam -s "PE" -o $1
submit_cluster --conda=${CONDA} --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf tin.py -r ${bed_file} -i Aligned.sortedByCoord.out.bam > tin.csv



if [ ${n_end} -eq 2 ] ; then
#submit_cluster --modules=rseqc,R --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf inner_distance.py -r ${bed_file} -i Aligned.sortedByCoord.out.bam -o $1
submit_cluster --conda=${CONDA} --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf inner_distance.py -r ${bed_file} -i Aligned.sortedByCoord.out.bam -o $1
fi 
echo "done RseqQC at $(date "+%Y-%m-%d %H:%M:%S")"
fi
fi

#qualimap
echo " Checking qualimap"
if [ ${do_qualimap} -gt 0 ]; then
mkdir qualimap_reports
echo "making qualimap"
if [ $(find ./qualimap_reports -name  "*.html"|wc -l) -gt 0 ]; then
echo "qualimap is already done "
else
#module load qualimap
#qualimap bamqc -bam Aligned.sortedByCoord.out.bam -c -outdir qualimap_reports --java-mem-size=20G
#submit_cluster --modules=qualimap,R --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf qualimap rnaseq -bam Aligned.sortedByCoord.out.bam -gtf $sjdbGTFfile -outdir qualimap_reports --java-mem-size=85G
#submit_cluster --conda=${CONDA} --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf unset DISPLAY; qualimap rnaseq -bam Aligned.sortedByCoord.out.bam -gtf $sjdbGTFfile -outdir qualimap_reports --java-mem-size=85G --java_options="-Djava.awt.headless=true"

submit_cluster --conda=${CONDA} --cluster_config=$SOURCE/pipelines/config_templates/cluster_configs/cluster_config_90Gb.conf qualimap  rnaseq -bam Aligned.sortedByCoord.out.bam -gtf $sjdbGTFfile -outdir qualimap_reports --java-mem-size=85G 
fi
fi
echo "doing feature counts"
#Feature Counts to output FPKM, RPKM
mkdir ../count_cpm_rpkm_tpm

if [ $organism == "human" ] ; then
gtf=${MY_BIN}/STAR/human_Genecode28/gencode.v28.primary_assembly.annotation.gtf
elif [ $organism == "mouse" ] ; then
gtf=${MY_BIN}/STAR/GRCm38_genecodevM11/gencode.vM11.primary_assembly.annotation.gtf
fi

if [ $(find ../count_cpm_rpkm_tpm/transcript_level/ -name  "*.txt"|wc -l) -eq 4 ]; then
echo " already calculated count, cpm, rpkm, tpm"
else

module purge
# Loading the old R version as Rsubread package only works with old R version.
#module load R/3.6  

#Rscript ${MY_BIN}/get_normalized_RNAseq_data.R --cmd=run --bam=Aligned.sortedByCoord.out.bam --GTF_annotation=$gtf --outdir=../count_cpm_rpkm_tpm --sampleID=$1
Rscript ${SOURCE}/get_normalized_RNAseq_data.R --cmd=run --bam=Aligned.sortedByCoord.out.bam --GTF_annotation=$gtf --outdir=../count_cpm_rpkm_tpm --sampleID=$1 --n_end=$n_end
fi

# Remove unwanted file to clean up space
if [ "${remove_bam}" = "1" ]; then
echo removing bam file
rm -rf Aligned.out.bam	Aligned.sortedByCoord.out.bam	Aligned.sortedByCoord.out.bam.bai Chimeric.out.junction Chimeric.out.sam Chimeric.out.sam  *.fq.gz  *.fastq.gz _STARgenome
else
echo keeping bam files
fi
