#####################################
### Preparing neighbor covariates ###
#####################################

# 16-Oct-2019
# This R script was written for 
# Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana
# which was co-authored by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano
# correspondence: sato.yasuhiro.36c@kyoto-u.jp

# (1) set functions 
# set a fild map
set_field = function(block_id) {
  field = matrix(rep(NA,10*40),10,40)
  for(i in 2*c(1:20)) {
    field[seq(2,10,2),i] = seq((i-1)*5+1,i*5,1)+(200*(block_id-1))
  }
  for(i in 2*c(1:20)-1) {
    field[seq(1,9,2),i] = seq((i-1)*5+1,i*5,1)+(200*(block_id-1))
  }
  return(field)
}

# function for calculating neighbor covariates
nei_coval = function(id,s=1) {
  s_seq = sort(seq(-s,s,by=1),decreasing = T)
  place = which(field==id,arr.ind=TRUE)
  n = rep(0,n_marker)
  k = 0
  for(i in s_seq) {
    for(j in s_seq) {
      nei_id = try(field[place[1]+i,place[2]+j],silent=T)
      if((class(nei_id)!="try-error")) { # skip if calling below the field
        if(length(nei_id)==1) { # skip if calling above the field
          if(is.na(nei_id)==FALSE) { # skip if a neighbor is absent
            if(nei_id!=id) { # excl. a focal plant itself
              n = n + snp_list_degit[,nei_id]
              k = k + 1
            }
          }
        }
      }
    }
  }
  coval = snp_list_degit[,id]*n/k
  print(id); print(k)
  return(coval)
}

# (2) load data

# load SNP data
load(file="sub_snpMAF5LD80.RData")
geno_d[geno_d==0]=-1 #replace 0 into -1

# load a position file
position = read.csv("./positionsMAF5LD80.csv",header=T)
position = position[,-1]

# load phenotype with gwas ID
pheno_d = read.csv("./phenotype.csv",header=T)
pheno_d$gwasID = paste0("X",pheno_d$gwasID)

n_marker = nrow(geno_d)
n_plants = nrow(pheno_d)


# (2) subsetting the data and calculating neighbor covariates
# make snp data for individuals
snp_list_degit = matrix(rep(0,n_marker*n_plants),n_marker,n_plants)
for(i in 1:n_plants) {
  col_id = pheno_d$gwasID[i]
  snp_vec = try(geno_d[,col_id])
  if(class(snp_vec)!="try-error") { snp_list_degit[,i] = snp_vec } else { print(i) }
}

# calc. neighbor covariates for block 1-8
geno_nei = c()
for(i in 1:8) {
  ph = subset(pheno_d, Block==i)
  field = set_field(i)
  g_nei = mapply(nei_coval,ph$IndivID,s=1)
  geno_nei = cbind(geno_nei, g_nei)
}

# reshape pheno data
pheno_d = subset(pheno_d, gwasID!="X7329")
n_marker = nrow(geno_d)
n_plants = nrow(pheno_d)

# setup snp data for individuals by excluding NAs
snp_list_degit = matrix(rep(NA,n_marker*n_plants),n_marker,n_plants)
for(i in 1:n_plants) {
  col_id = pheno_d$gwasID[i]
  snp_vec = try(geno_d[,col_id])
  if(class(snp_vec)!="try-error") { snp_list_degit[,i] = snp_vec }
}

nei_list_degit = geno_nei[,pheno_d$IndivID]

# release memory
rm(geno_nei)
gc();gc()

pheno = pheno_d
self = snp_list_degit
nei = nei_list_degit

# release memory
rm(snp_list_degit)
rm(nei_list_degit)
gc();gc()

# export input data for lmm
save(self, file = "lmmJPN2017selfMAF5.RData",compress=T)
save(nei, file = "lmmJPN2017neiMAF5.RData",compress=T)
save(pheno, file = "lmmJPN2017phenoMAF5.RData")
