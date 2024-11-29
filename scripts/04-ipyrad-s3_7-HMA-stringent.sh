#!/usr/bin/env bash
#SBATCH --job-name=HMA-stringent-s3_7
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --time=4:00:00
#SBATCH --mem=60GB
#SBATCH -o ./joblogs/ipyrad-s37/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

# Variables
DIR='/home/a1645424/hpcfs/analysis/shannon'
cd "${DIR}/results/ipyrad" || exit 1

source "/home/a1645424/hpcfs/micromamba/etc/profile.d/micromamba.sh"
micromamba activate ipyrad

ipyrad \
    -s 34567 \
    -p 'params-HMA-stringent.txt' \
    -c "${SLURM_CPUS_PER_TASK}"

micromamba deactivate
