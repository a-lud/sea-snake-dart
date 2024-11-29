#!/usr/bin/env bash

RESSTORE='/uofaresstor/sanders_lab/sequencing-datasets/radseq'
FQDIR="/hpcfs/users/${USER}/analysis/shannon/data/fastq"

mkdir -p "${FQDIR}"

cd "${FQDIR}" || exit 1

# Read 'sample-sheet.csv' line by line
while IFS=',' read -r -a LINE; do
  if [[ ! -f "${RESSTORE}/DaRT-${LINE[0]}/${LINE[1]}.FASTQ.gz" ]]; then
    echo "File for ${LINE[]} does not exist in ${LINE[0]}"
    exit 1
  fi
  echo -e "FROM: ${RESSTORE}/DaRT-${LINE[0]}/${LINE[1]}.FASTQ.gz\tTO: ${LINE[2]}.FASTQ.gz" 
  rsync -aP "${RESSTORE}/DaRT-${LINE[0]}/${LINE[1]}.FASTQ.gz" .
  rename -v "${LINE[1]}" "${LINE[2]}" "${LINE[1]}.FASTQ.gz"
done < <(tail -n +2 '/hpcfs/users/a1645424/analysis/shannon/data/240524-sample-linkage.csv')

