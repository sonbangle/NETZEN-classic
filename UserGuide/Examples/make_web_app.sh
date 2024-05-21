#!/bin/sh
#SBATCH --job-name=make_web_app # Job name
#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --cpus-per-task=1
#SBATCH --mem=12Gb
#SBATCH --time=96:00:00
#SBATCH --qos=dtran-b
#SBATCH --output=cluster_logs/logs/1_generic_job_%j.out

dataset_name=ATOX_KO
data_folder=ATOX1_KO_WebApp/NETZEN_WebApp/data/ATOX1_KO/
outdir=ATOX1_KO_WebApp


echo started:`date` 
module purge 
source /etc/profile.d/modules.sh
module load python
$NETZEN/make_Web_app $dataset_name  $data_folder $outdir

echo terminated:`date` 
