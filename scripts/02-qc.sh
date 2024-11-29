#!/usr/bin/env bash
#SBATCH --job-name=qc
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -a 1-497%25
#SBATCH --time=02:00:00
#SBATCH --mem=60GB
#SBATCH -o ./joblogs/qc/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

DIR='/home/a1645424/hpcfs/analysis/shannon'
FQDIR="${DIR}/data/fastq"
OUTDIR="${DIR}/results/qc"
MULTIQC="${OUTDIR}/multiqc"

# Databases
UNIVEC='/home/a1645424/al-biohub/database/univec/univec.fasta'
KDB='/home/a1645424/al-biohub/database/k2_standard_20210517'

FQ=$(find "${FQDIR}" -type f -name '*.gz' | tr '\n' ' ' | cut -d' ' -f "${SLURM_ARRAY_TASK_ID}")
BN=$(basename "${FQ%%.*}")
SEQ_RUN=$(grep --no-filename "${BN}" ${DIR}/data/sample-sheets/*.csv | head -n 1 | cut -d',' -f 1)

echo "SAMPLE: ${BN} RUN: ${SEQ_RUN}"
mkdir -p "${OUTDIR}" "${MULTIQC}" "${OUTDIR}/fastp" "${OUTDIR}/kraken2" "${OUTDIR}/bbduk"

source "/home/a1645424/hpcfs/micromamba/etc/profile.d/micromamba.sh"

# Newer DaRT runs remove barcodes and adapter content
micromamba activate bbmap
if [[ "${SEQ_RUN}" =~ ("DNote23-8392"|"DNote23-8556"|"DNote23-8773"|"DNote24-9763") ]]; then
  # DON'T use hard filter to trim adapter (already removed)
  bbduk.sh \
    --in="${FQ}" \
    --out="${OUTDIR}/bbduk/${BN}.fastq.gz" \
    --ref="${UNIVEC}" \
    --ktrim=r \
    --mink=11 \
    --threads="${SLURM_CPUS_PER_TASK}" 2> "${MULTIQC}/${BN}.bbduk.log"
else
  # Get total length of barcode and cut site overhang
  B9l=$(grep "${BN}" ${DIR}/data/sample-sheets/*.csv | head -n 1 | cut -d',' -f 4)
  LENGTH=$(echo -n "${B9l}" | wc -c)

  bbduk.sh \
    --in="${FQ}" \
    --out="${OUTDIR}/bbduk/${BN}.fastq.gz" \
    --ref="${UNIVEC}" \
    --ktrim=r \
    --mink=11 \
    --ftl="${LENGTH}" \
    --threads="${SLURM_CPUS_PER_TASK}" 2> "${MULTIQC}/${BN}.bbduk.log"
fi
micromamba deactivate

# Kraken
micromamba activate kraken2
kraken2 \
  --db "${KDB}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --gzip-compressed \
  --output '-' \
  --unclassified-out "${OUTDIR}/kraken2/${BN}-unclassified.fastq" \
  --report "${MULTIQC}/${BN}.report" \
  "${OUTDIR}/bbduk/${BN}.fastq.gz"
pigz -p "${SLURM_CPUS_PER_TASK}" "${OUTDIR}/kraken2/${BN}-unclassified.fastq"

micromamba activate fastp
fastp \
  -i "${OUTDIR}/kraken2/${BN}-unclassified.fastq.gz" \
  -o "${OUTDIR}/fastp/${BN}.fastq.gz" \
  --average_qual 10 \
  --length_required 25 \
  --thread "${SLURM_CPUS_PER_TASK}" \
  --json "${MULTIQC}/${BN}.fastp.json"
micromamba deactivate

if [[ -f "${OUTDIR}/bbduk/${BN}.fastq.gz" ]]; then rm -v "${OUTDIR}/bbduk/${BN}.fastq.gz"; fi
if [[ -f "${OUTDIR}/kraken2/${BN}-unclassified.fastq.gz" ]]; then rm -v "${OUTDIR}/kraken2/${BN}-unclassified.fastq.gz"; fi
