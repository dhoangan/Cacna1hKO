#!/bin/bash
echo "Start fastp for reporting and filtering"
i=1
for fn in data/CH-US-HAD-0{01..12}_S;
do
samp=${fn}${i}
echo "Processing sample ${samp}"
read_1=${samp}_L001_R1_001.fastq.gz
read_2=${samp}_L001_R2_001.fastq.gz

echo "read 1: ${read_1}"
echo "read 2: ${read_2}"
fastp -g -h fastp_report/fastp_report${i}.html -j fastp_report/fastp_report${i}.json\
         -i ${read_1} -o filtered_data/filtered_S${i}_R1.fastq.gz \
         -I ${read_2} -O filtered_data/filtered_S${i}_R2.fastq.gz\
		 
((i=i+1))
done
echo "End" 
#read -p "Press enter to continue";
