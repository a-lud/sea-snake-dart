#!/usr/bin/env bash
#SBATCH --job-name=pixy-allsite-vcfs
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --ntasks-per-core=2
#SBATCH -a 3
#SBATCH --time=20:00:00
#SBATCH --mem=4GB
#SBATCH -o ./joblogs/allsite-vcfs/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

DIR='/home/a1645424/hpcfs/analysis/shannon'

VCF=$(find "${DIR}/results/diversity-statistics" -name '*all_loci.vcf.gz' | tr '\n' ' ' | cut -d ' ' -f "${SLURM_ARRAY_TASK_ID}")
BN=$(basename "${VCF}" -stringent-all_loci.vcf.gz)
OUT="${DIR}/results/diversity-statistics/pixy"

mkdir -p "${OUT}"

source '/hpcfs/users/a1645424/micromamba/etc/profile.d/micromamba.sh'
micromamba activate pixy

echo "${BN}"
pixy \
    --stats pi \
    --vcf "${VCF}" \
    --populations "${DIR}/data/popmaps/${BN}-popmap.tsv" \
    --window_size 50000 \
    --n_cores "${SLURM_CPUS_PER_TASK}" \
    --output_prefix "${BN}" \
    --output_folder "${OUT}"

micromamba deactivate
