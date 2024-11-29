#!/usr/bin/env bash

# Set a sample to missing (./.) if their depth at a variant site is <5 coverage
# Add in 'tags' - specifically the fraction of missing data per population
# Remove any site where more than 10% of samples have missing data
#   - Applied across the WHOLE dataset - not by population

DIR='/home/a1645424/hpcfs/analysis/shannon'
VCF=$(find "${DIR}/results/ipyrad" -name '*-stringent.vcf.gz')

for V in $VCF; do 
    BN=$(basename "${V}" -stringent.vcf.gz)

    echo "${BN}"

    OUT="${DIR}/results/ipyrad/${BN}-stringent_outfiles"
    POP="${DIR}/data/popmaps/${BN}-popmap.tsv"

    bcftools view "${V}" | 
        bcftools +setGT -- -t q -n . -i 'FMT/DP<5' |
        bcftools +fill-tags -- -t F_MISSING -S "${POP}" |
        bcftools filter -i 'F_MISSING<=0.1'| 
        bgzip > ${OUT}/${BN}-stringent.highQ.filtered.vcf.gz
    
    tabix ${OUT}/${BN}-stringent.highQ.filtered.vcf.gz
done

