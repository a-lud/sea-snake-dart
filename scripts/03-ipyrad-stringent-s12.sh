#!/usr/bin/env bash
#SBATCH --job-name=all_samples_stringent-s12
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --time=4:00:00
#SBATCH --mem=20GB
#SBATCH -o ./joblogs/ipyrad-s12/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

# Variables
DIR='/home/a1645424/hpcfs/analysis/shannon'
cd "${DIR}/results/ipyrad" || exit 1

source "/home/a1645424/hpcfs/micromamba/etc/profile.d/micromamba.sh"
micromamba activate ipyrad

# Use the '--ipcluster' argument to connect to the cluster we started above
ipyrad \
    -s 12 \
    -p 'params-all_samples_stringent-s12.txt' \
    -c "${SLURM_CPUS_PER_TASK}"

micromamba deactivate
