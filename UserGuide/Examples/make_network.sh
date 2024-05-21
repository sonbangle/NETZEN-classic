#/bin/bash
net_file=neurons_net.csv

cp NETWORK/genereprun/consolidated_tpm_table_translated.cyto.csv ${net_file}.tmp
cut ${net_file}.tmp -f 1,2,3 > $net_file
sed -i '1d'  $net_file
rm ${net_file}.tmp
