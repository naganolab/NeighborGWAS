######################################################
# R script to combine multiple output files into one #
######################################################

# 16-Oct-2019
# This R script was written for 
# Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana
# which was co-authored by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano
# correspondence: sato.yasuhiro.36c@kyoto-u.jp

position = read.csv("../positionsMAF5LD80.csv",header=T)
position = position[,-1]

file = "lmmGastonJPN2017ScoreMAF5"
per_pos = 60000
n_file = 21
path = paste0("./output/", file, "_", 1+per_pos*c(0:(n_file-1)), ".csv")

p_list_lm = data.frame()
for(i in 1:n_file) {
	f = read.csv(path[i],header=T)
	p_list_lm = rbind(p_list_lm, f)
}

p_list = data.frame(position,p_list_lm)
colnames(p_list) = c("Chr", "Position", "MAF", "beta_self", "beta_nei", "P_self", "P_nei")
write.csv(p_list,paste0("./output/", file, ".csv"))

