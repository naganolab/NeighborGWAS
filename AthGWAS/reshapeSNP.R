################################################################################
# Cuting off SNP data with MAF and LD, and then preparing SNPs for gl1 mutants #
################################################################################

# 16-Oct-2019
# This R script was written for 
# Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana
# which was co-authored by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano
# correspondence: sato.yasuhiro.36c@kyoto-u.jp

# (1) load data and cut off with MAF and LD
# load SNP matrix, col = gwasID, row = position
d_geno = read.csv("sub_snps.csv",header=T)
d_geno = d_geno[,-1]
gwasID = read.csv("gwasIDlist.csv",header=T)
colnames(d_geno) = paste0("X",gwasID$GWASid)

# read a position file
all_pos = read.csv("positions.csv",header=T)
all_pos = all_pos[,2]

n_loci = nrow(d_geno)
n_acc = ncol(d_geno)
a_freq = apply(d_geno,1,sum)/n_acc

# cut off at MAF > 0.05
n_sub_snps = sum(a_freq>0.05&a_freq<0.95)
sub_pos = which(a_freq>0.05&a_freq<0.95)
d_geno_maf = d_geno[sub_pos,]

# cut off at r2 > 0.8 of LD between adjacent loci
cor_LD = function(i) {
  r = cor(as.numeric(d_geno_maf[i,]),as.numeric(d_geno_maf[i+1,]))
  return(r)
}
cor_pos = mapply(cor_LD,1:(n_sub_snps-1))

sub_pos = sub_pos[cor_pos^2<0.80]
d_geno_maf = d_geno[sub_pos,]

# cut off the position
position = all_pos[sub_pos]
n_sub_snps = nrow(d_geno_maf)

chr_which = which((position[-length(position)] - position[-1])>0)

# (2) consider single-gene mutants
# insert Col(gl1-2) SNPs
n_pos = sum(position[chr_which[2]:chr_which[3]]<10362188)

position[(chr_which[2]+n_pos)] # just before gl1-2

position = c(position[1:(chr_which[2]+n_pos)], 10362188, position[(chr_which[2]+n_pos+1):n_sub_snps])
d_geno_maf = rbind(d_geno_maf[1:(chr_which[2]+n_pos),],rep(0,196),d_geno_maf[(chr_which[2]+n_pos+1):n_sub_snps,])

gl1_2 = d_geno_maf[,"X6909"] # Col-0
gl1_2[chr_which[2]+n_pos+1] = 1

d_geno_maf = cbind(d_geno_maf, gl1_2)

# insert Ler(gl1-1) SNPs
n_pos1 = sum(position[chr_which[2]:chr_which[3]]<10361426)
n_pos2 = sum(position[chr_which[2]:chr_which[3]]<10364564)

position[(chr_which[2]+n_pos1+1):(chr_which[2]+n_pos2+1)] # whole GL1 region
gl1_1 = d_geno_maf[,"X6932"] #Ler-1
gl1_1[(chr_which[2]+n_pos1+1):(chr_which[2]+n_pos2+1)] = 1

d_geno_maf = cbind(d_geno_maf, gl1_1)

# assign gwasID of 1 and 2 to gl1-2 and gl1-1, respectively.
colnames(d_geno_maf) = c(colnames(d_geno), "X2", "X1")

rm(d_geno)
gc();gc()

# (3) export the input genotype data
# export SNP matrix
geno_d = d_geno_maf
save(geno_d, file="sub_snpMAF5LD80.RData", compress = T)

# export the position file with the choromosome no.
chr_which = which((position[-length(position)] - position[-1])>0)
chr = c(rep(1, chr_which[1]), rep(2, (chr_which[2]-chr_which[1])), rep(3, (chr_which[3]-chr_which[2])), rep(4, (chr_which[4]-chr_which[3])), rep(5, (nrow(d_geno_maf)-chr_which[4])))
n_acc = ncol(d_geno_maf)
maf = apply(d_geno_maf,1,sum)/n_acc
posMAF5 = cbind(chr, position, maf)
write.csv(posMAF5, "positionsMAF5LD80.csv")
