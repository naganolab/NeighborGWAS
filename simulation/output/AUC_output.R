######################################################
# R script to combine multiple output files into one #
######################################################

# 16-Oct-2019
# This R script was written for 
# Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana
# which was co-authored by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano
# correspondence: sato.yasuhiro.36c@kyoto-u.jp

file = "AUCs_"
n_file = 72

out = data.frame()
for(i in 1:n_file) {
	f = read.csv(paste0(file,i,".csv"),header=T)
	out = rbind(out, f)
}

write.csv(out,paste0(file, ".csv"))
