#!/usr/bin/env bash
#SBATCH --job-name=structure-HST-all
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 38
#SBATCH --ntasks-per-core=1
#SBATCH --time=48:00:00
#SBATCH --mem=20GB
#SBATCH -o ./joblogs/structure/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

SCRIPT_DIR=$(pwd)

HDF5="${SCRIPT_DIR}/../results/variant-filtering/HST-ld-maf.snps.hdf5"
POPMAP="${SCRIPT_DIR}/../data/popmaps/HST-clean-popmap.tsv"
OUTDIR="${SCRIPT_DIR}/../results/population-structure/structure"

# Change into working directory
cd "${OUTDIR}" || exit 1

source "/hpcfs/users/$USER/micromamba/etc/profile.d/micromamba.sh"
micromamba activate ipyrad

ipcluster start -n "${SLURM_CPUS_PER_TASK}" --cluster-id="ipyrad_HST" --daemonize

sleep 60

python ${SCRIPT_DIR}/ipa-structure.py 'HST_all' "${HDF5}" "${POPMAP}" ipyrad_HST -k 2 10

ipcluster stop --cluster-id="ipyrad_HST"

micromamba deactivate
