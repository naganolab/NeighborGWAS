#######################################
# control run_lmm.R from the terminal #
#######################################

# 16-Oct-2019
# This shell script was written for 
# Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana
# which was co-authored by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano
# correspondence: sato.yasuhiro.36c@kyoto-u.jp

# split whole GWAS into multiple single-core iterations

Rscript run_lmm.R 1 60000&
Rscript run_lmm.R 60001 120000&
Rscript run_lmm.R 120001 180000&
Rscript run_lmm.R 180001 240000&
sleep 20m
Rscript run_lmm.R 240001 300000&
Rscript run_lmm.R 300001 360000&
Rscript run_lmm.R 360001 420000&
Rscript run_lmm.R 420001 480000&
sleep 20m
Rscript run_lmm.R 480001 540000&
Rscript run_lmm.R 540001 600000&
Rscript run_lmm.R 600001 660000&
Rscript run_lmm.R 660001 720000&
sleep 20m
Rscript run_lmm.R 720001 780000&
Rscript run_lmm.R 780001 840000&
Rscript run_lmm.R 840001 900000&
Rscript run_lmm.R 900001 960000&
sleep 20m
Rscript run_lmm.R 960001 1020000&
Rscript run_lmm.R 1020001 1080000&
Rscript run_lmm.R 1080001 1140000&
Rscript run_lmm.R 1140001 1200000&
Rscript run_lmm.R 1200001 1242127&
wait

# GWAS output file is made at "./lmekinGWAS/output/XXXX".
Rscript ./output/lmm_output.R
cd ./output
rm lmmGastonJPN2017ScoreMAF5_*
