#!/usr/bin/env bash

DIR='/home/a1645424/hpcfs/analysis/shannon'
FQ='/home/a1645424/hpcfs/analysis/shannon/results/qc/fastp'

IDS=$(cat ${DIR}/data/sample-lists/*.txt)
for i in $IDS; do
    if [ ! -f ${FQ}/${i}.fastq.gz ]; then
        echo "$i" >> 'need-processing.txt' 
    fi
done
