#!/usr/bin/env bash
#SBATCH --job-name=structure-HMA-refined
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

HDF5='/home/a1645424/hpcfs/analysis/shannon/results/ipyrad/HMA-stringent_outfiles'
WD='/home/a1645424/hpcfs/analysis/shannon/results/population-structure/structure'
mkdir -p "${WD}"

# Copy required data and files to working directory
cp structure-HMA.py "${WD}"
cp "${HDF5}/HMA-stringent.highQ.filtered.LD50k.snps.hdf5" "${WD}"
cp /home/a1645424/hpcfs/analysis/shannon/data/popmaps/HMA-popmap.txt "${WD}"
cp /home/a1645424/hpcfs/analysis/shannon/data/popmaps/dropped-samples.tsv "${WD}"
cp "${WD}/../HMA-ignore.txt" "${WD}"

# Change into working directory
cd "${WD}" || exit 1

source '/hpcfs/users/a1645424/micromamba/etc/profile.d/micromamba.sh'
micromamba activate ipyrad

ipcluster start -n "${SLURM_CPUS_PER_TASK}" --cluster-id="ipyrad_HMA" --daemonize

sleep 60

python structure-HMA.py

ipcluster stop
ipcluster clean

micromamba deactivate
