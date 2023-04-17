#!/bin/bash
echo "Start quantification of filtered data"

for fn in filtered_data/filtered_S{1..12};
do
samp=${fn}
echo "Processing sample ${samp}"
read_1=${samp}_R1.fastq.gz
read_2=${samp}_R2.fastq.gz

echo "read 1: ${read_1}"
echo "read 2: ${read_2}"
salmon quant -i mouse_index_GRCm39 -l A \
         -1 ${read_1} \
         -2 ${read_2} \
         -p 8 --validateMappings --gcBias -o quants_filtered_GRCm39/${samp}_quant_filtered
done
echo "End" 
#read -p "Press enter to continue";
