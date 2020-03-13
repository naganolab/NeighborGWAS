#######################################
# shell script to prepare input files #
#######################################

# 23-Feb-2020
# This shell script was written for 
# Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana
# which was co-authored by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano
# correspondence: sato.yasuhiro.36c@kyoto-u.jp

# subset SNP matrix from the full imputed data.
python subsetSNP.py

# cut off the SNP data with MAF and LD
Rscript reshapeSNP.R
rm sub_snps.csv
rm positions.csv

# prepare neighbor covariates
Rscript prepNeighborsJPN.R

