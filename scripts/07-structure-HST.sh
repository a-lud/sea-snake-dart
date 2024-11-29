#!/usr/bin/env bash
#SBATCH --job-name=structure-HST
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --ntasks-per-core=2
#SBATCH --time=10:00:00
#SBATCH --mem=20GB
#SBATCH -o ./joblogs/structure/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

DIR='/home/a1645424/hpcfs/analysis/shannon/results/ipyrad/HST-stringent_outfiles'
mkdir -p "${DIR}/structure"

# Copy required data and files to working directory
cp structure-HST.py "${DIR}/structure"
cp "${DIR}/HST-stringent.highQ.filtered.LD50k.snps.hdf5" "${DIR}/structure"
cp /home/a1645424/hpcfs/analysis/shannon/data/popmaps/HST-popmap.txt "${DIR}/structure"

# Change into working directory
cd "${DIR}/structure" || exit 1

source '/hpcfs/users/a1645424/micromamba/etc/profile.d/micromamba.sh'
micromamba activate ipyrad

ipcluster start -n "${SLURM_CPUS_PER_TASK}" --cluster-id="ipyrad_HST" --daemonize

sleep 60

python structure-HST.py

ipcluster stop

micromamba deactivate
