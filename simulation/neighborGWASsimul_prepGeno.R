######################################
# Preparing genotypes for simulation #
######################################

# 23-Feb-2020
# This R script was written for 
# Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana
# which was co-authored by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano
# correspondence: sato.yasuhiro.36c@kyoto-u.jp

library(Matrix)

# load data
g = read.csv("call_method_75_TAIR9_250k.csv",header=T, skip=1) #data downloaded from http://bergelson.uchicago.edu/?page_id=790
line_names = colnames(g)[-c(1:2)]
info = g[,1:2]
g = g[,-c(1:2)]

# convert rare alleles into -1 degit
snp_degit = function(vec) {
  if(table(vec)[1]>table(vec)[2]) {
    vec[which(vec==levels(factor(vec))[1])] = 0
    vec[which(vec==levels(factor(vec))[2])] = 1
  } else {
    vec[which(vec==levels(factor(vec))[1])] = 1
    vec[which(vec==levels(factor(vec))[2])] = 0
  }
  return(as.numeric(vec))
}

g_bin = apply(g,1,snp_degit)
g_bin = t(g_bin)
colnames(g_bin) = line_names

rm(g)
gc();gc()

# calc. minor allele frequency
maf = apply(g_bin,1,sum)/ncol(g_bin)

tg = which((info$Chromosome==1)|(info$Chromosome==2))
g_bin = g_bin[tg,][maf[tg]>0.1,]

save(g_bin,file="Ara250kRegMap_chr12_MAF10.RData",compress=T)


