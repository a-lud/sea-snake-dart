#!/usr/bin/env bash
#SBATCH --job-name=write-allsite-vcfs
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --ntasks-per-core=1
#SBATCH --time=02:00:00
#SBATCH --mem=4GB
#SBATCH -o ./joblogs/allsite-vcfs/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

DIR='/home/a1645424/hpcfs/analysis/shannon'
SCRIPT="${DIR}/scripts"
SAMPLES="${DIR}/data/sample-lists"
OUTDIR="${DIR}/results/diversity-statistics"

mkdir -p ${OUTDIR}

source '/hpcfs/users/a1645424/micromamba/etc/profile.d/micromamba.sh'
micromamba activate ipyrad

cd ${OUTDIR} || exit 1
echo "[Processing] ALA"
"${SCRIPT}/loci2vcf.py" \
    "${DIR}/results/ipyrad/ALA-stringent_outfiles/ALA-stringent.loci" \
    "${DIR}/results/ipyrad/ALA-stringent_outfiles/ALA-stringent-all_loci.vcf" \
    "${SAMPLES}/ALA-samples.txt" \
    True

cp "${DIR}/results/ipyrad/ALA-stringent_outfiles/ALA-stringent-all_loci.vcf" .
bgzip 'ALA-stringent-all_loci.vcf'
tabix 'ALA-stringent-all_loci.vcf.gz'

echo "[Processing] HMA"
"${SCRIPT}/loci2vcf.py" \
    "${DIR}/results/ipyrad/HMA-stringent_outfiles/HMA-stringent.loci" \
    "${DIR}/results/ipyrad/HMA-stringent_outfiles/HMA-stringent-all_loci.vcf" \
    "${SAMPLES}/HMA-samples.txt" \
    True

cp "${DIR}/results/ipyrad/HMA-stringent_outfiles/HMA-stringent-all_loci.vcf" .
bgzip 'HMA-stringent-all_loci.vcf'
tabix 'HMA-stringent-all_loci.vcf.gz'

echo "[Processing] HST"
"${SCRIPT}/loci2vcf.py" \
    "${DIR}/results/ipyrad/HST-stringent_outfiles/HST-stringent.loci" \
    "${DIR}/results/ipyrad/HST-stringent_outfiles/HST-stringent-all_loci.vcf" \
    "${SAMPLES}/HST-samples.txt" \
    True

cp "${DIR}/results/ipyrad/HST-stringent_outfiles/HST-stringent-all_loci.vcf" .
bgzip 'HST-stringent-all_loci.vcf'
tabix 'HST-stringent-all_loci.vcf.gz'

micromamba deactivate
