#!/usr/bin/env bash
#SBATCH --job-name=structure-ALA-refined
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 40
#SBATCH --ntasks-per-core=2
#SBATCH --time=25:00:00
#SBATCH --mem=20GB
#SBATCH -o ./joblogs/structure/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

HDF5='/home/a1645424/hpcfs/analysis/shannon/results/ipyrad/ALA-stringent_outfiles'
WD='/home/a1645424/hpcfs/analysis/shannon/results/population-structure/structure'
mkdir -p "${WD}"

# Copy required data and files to working directory
cp structure-ALA.py "${WD}"
cp "${HDF5}/ALA-stringent.highQ.filtered.LD50k.snps.hdf5" "${WD}"
cp /home/a1645424/hpcfs/analysis/shannon/data/popmaps/ALA-popmap.txt "${WD}"
cp /home/a1645424/hpcfs/analysis/shannon/data/popmaps/dropped-samples.tsv "${WD}"
cp "${WD}/../ALA-ignore.txt" "${WD}"

# Change into working directory
cd "${WD}" || exit 1

source '/hpcfs/users/a1645424/micromamba/etc/profile.d/micromamba.sh'
micromamba activate ipyrad

ipcluster start -n "${SLURM_CPUS_PER_TASK}" --cluster-id="ipyrad_ALA" --daemonize

sleep 60

python structure-ALA.py

ipcluster stop
ipcluster clean

micromamba deactivate
