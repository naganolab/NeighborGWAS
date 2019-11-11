####################################################
### Variation partitioning and association tests ###
####################################################

# 16-Oct-2019
# This R script was written for 
# Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana
# which was co-authored by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano
# correspondence: sato.yasuhiro.36c@kyoto-u.jp

# (1) load input data and arguments
# run GWAS from args1-th to args2-th SNPs
args1 = commandArgs(trailingOnly=TRUE)[1]
args2 = commandArgs(trailingOnly=TRUE)[2]

# load libraties
library(Matrix); library(gaston)
setThreadOptions(numThreads = 1)

# load data
load("../lmmJPN2017selfMAF5.RData")
load("../lmmJPN2017neiMAF5.RData")
load("../lmmJPN2017phenoMAF5.RData")

# (2) variation partitioning using lmm.aireml()
selfCross = crossprod(self)
neiCross = crossprod(nei)

selfK = ((nrow(self)-1)+selfCross)/(2*(nrow(self)-1))
neiK = neiCross/(nrow(nei)-1)

# approximate the kinship matrices as a near postive-definite matrix
selfPDK = as.matrix(nearPD(selfK,maxit=10^6)$mat)
neiPDK = as.matrix(nearPD(neiK,maxit=10^6)$mat)

Y = scale(pheno$Score)
X = model.matrix(~pheno$Bolting+scale(pheno$InitLeafLen)+factor(pheno$Block))
res = lmm.aireml(Y = Y, X = X, K=list(selfPDK,neiPDK), verbose = T)
tau = res$tau

Vs = res$tau[1]/sum(res$tau,res$sigma2)
paste0("Vs = ", Vs)
Vn = res$tau[2]/sum(res$tau,res$sigma2)
paste0("Vn = ", Vn)

# (3) association test using lmm.diago()
eiKs = eigen(as.matrix(selfK)) # self kinship matrix

kin_mat = tau[1]*selfPDK + tau[2]*neiPDK
eiKsn = eigen(as.matrix(kin_mat)) # mixed kinship matrix based on tau[1] and tau[2]

self = self[c(args1:args2),]
nei = nei[c(args1:args2),]

gc();gc()

# estimate a likelihood of the null model
res_null = lmm.diago(Y = Y, X = X, eigenK = eiKs, method="newton", verbose = F)
LL_null = lmm.diago.profile.likelihood(tau=res_null$tau, s2=res_null$sigma2, Y = Y, X = X, eigenK = eiKs)[1,1]

res_null2 = lmm.diago(Y = Y, X = X, eigenK = eiKsn, method="newton", verbose = F)

# run ordinal GWAS and then test neighbor effects with LL-tests
p_list_lm = c()
for(i in 1:nrow(self)) { print(sum(as.numeric(args1)+as.numeric(i)))
  
  X = model.matrix(~self[i,]+pheno$Bolting+scale(pheno$InitLeafLen)+factor(pheno$Block))
  LL_self = lmm.diago.profile.likelihood(tau=res_null$tau, s2=res_null$sigma2, Y = Y, X = X, eigenK = eiKs)[1,1]

  X = model.matrix(~self[i,]+pheno$Bolting+scale(pheno$InitLeafLen)+factor(pheno$Block))
  LL_self2 = lmm.diago.profile.likelihood(tau=res_null2$tau, s2=res_null2$sigma2, Y = Y, X = X, eigenK = eiKsn)[1,1]
  
  X = model.matrix(~self[i,]+nei[i,]+pheno$Bolting+scale(pheno$InitLeafLen)+factor(pheno$Block))
  LL_nei = lmm.diago.profile.likelihood(tau=res_null2$tau, s2=res_null2$sigma2, Y = Y, X = X, eigenK = eiKsn)[1,1]
  
  res = lmm.diago(Y = Y, X = X, eigenK = eiKsn, verbose = F)
  
  p_self = pchisq(-2*(LL_null-LL_self),1,lower.tail = F)
  p_nei = pchisq(-2*(LL_self2-LL_nei),1,lower.tail = F)
  
  p_list_lm = rbind(p_list_lm, c(res$BLUP_beta[2:3],p_self,p_nei))
}

# (4) export results
colnames(p_list_lm) = c("beta_self", "beta_nei", "p_self", "p_nei")
write.csv(p_list_lm, paste0("./output/lmmGastonJPN2017ScoreMAF5_", args1, ".csv"), row.names = F)

