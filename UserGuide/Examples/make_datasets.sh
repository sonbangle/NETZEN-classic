outdir=.
datasets=$(ls  -d ../NETZEN_analysis/cell_types/Cluster_*/NETZEN_WebApp/data/Cluster*)
echo -e 'name\tlocation' > ${outdir}/datasets.csv
rm -f ${outdir}/datasets.html


for dataset in ${datasets}
do
ln -s ${dataset} data/
dataset_name=$(basename ${dataset})
echo ${dataset_name}
echo -e ${dataset_name}'\t'${dataset_name} >> ${outdir}/datasets.csv
done

$NETZEN/json_export  ${outdir}/datasets.csv  ${outdir}/datasets.html
echo end exporting dataset
