####################################################
### Variation partitioning and association tests ###
####################################################

# 23-Feb-2020
# This R script was written for 
# Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana
# which was co-authored by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano
# correspondence: sato.yasuhiro.36c@kyoto-u.jp

# (1) load input data and arguments
# run GWAS from args1-th to args2-th SNPs
args1 = commandArgs(trailingOnly=TRUE)[1]
args2 = commandArgs(trailingOnly=TRUE)[2]

library(Matrix); library(gaston)
setThreadOptions(numThreads = 1)

# load data
load("../lmmJPN2017selfMAF5.RData")
load("../lmmJPN2017neiMAF5.RData")
int = self*nei # set X for testing beta_12
self = self[c(args1:args2),]
nei = nei[c(args1:args2),]
int = int[c(args1:args2),]

load("../lmmJPN2017phenoMAF5.RData")
gc();gc()

# load kinship matrices
load("selfPDK.RData")
load("neiPDK.RData")
load("intPDK.RData")
load("eiKs.RData")

# (2) variation partitioning using lmm.aireml()
Y = scale(pheno$Score)
X = model.matrix(~pheno$Bolting+scale(pheno$InitLeafLen)+factor(pheno$Block)+pheno$edge)
res = lmm.aireml(Y = Y, X = X, K=list(selfPDK,neiPDK), verbose = T)
tau = res$tau

Vs = res$tau[1]/sum(res$tau,res$sigma2)
Vn = res$tau[2]/sum(res$tau,res$sigma2)
paste0("Vs = ", Vs)
paste0("Vn = ", Vn)

kin_mat = tau[1]*selfPDK + tau[2]*neiPDK
eiKsn = eigen(as.matrix(kin_mat)) #mixed kinship matrix based on tau[1] and tau[2]

X = model.matrix(~pheno$Bolting+scale(pheno$InitLeafLen)+factor(pheno$Block)+pheno$edge)
res = lmm.aireml(Y = Y, X = X, K=list(selfPDK,neiPDK,intPDK), verbose = T)
tau = res$tau

Vs = res$tau[1]/sum(res$tau,res$sigma2)
paste0("Vs = ", Vs)
Vn = res$tau[2]/sum(res$tau,res$sigma2)
paste0("Vn = ", Vn)
Vsxn = res$tau[3]/sum(res$tau,res$sigma2)
paste0("Vsxn = ", Vsxn)

# (3) association test using lmm.diago()
kin_mat = tau[1]*selfPDK + tau[2]*neiPDK + tau[3]*intPDK
eiKsxn = eigen(as.matrix(kin_mat)) #mixed kinship matrix based on tau[1] and tau[2]

# estimate likelihood of the null model
res_null = lmm.diago(Y = Y, X = X, eigenK = eiKs, method="newton", verbose = F)
LL_null = lmm.diago.profile.likelihood(tau=res_null$tau, s2=res_null$sigma2, Y = Y, X = X, eigenK = eiKs)[1,1]

res_null2 = lmm.diago(Y = Y, X = X, eigenK = eiKsn, method="newton", verbose = F)

# run normal GWAS and then test neighbor effects with LL-tests
p_list_lm = c()
for(i in 1:nrow(self)) { print(sum(as.numeric(args1)+as.numeric(i)))
  
  X = model.matrix(~self[i,]+pheno$Bolting+scale(pheno$InitLeafLen)+factor(pheno$Block)+pheno$edge)
  LL_self = lmm.diago.profile.likelihood(tau=res_null$tau, s2=res_null$sigma2, Y = Y, X = X, eigenK = eiKs)[1,1]

  X = model.matrix(~self[i,]+pheno$Bolting+scale(pheno$InitLeafLen)+factor(pheno$Block)+pheno$edge)
  LL_self2 = lmm.diago.profile.likelihood(tau=res_null2$tau, s2=res_null2$sigma2, Y = Y, X = X, eigenK = eiKsn)[1,1]
  
  X = model.matrix(~self[i,]+nei[i,]+pheno$Bolting+scale(pheno$InitLeafLen)+factor(pheno$Block)+pheno$edge)
  LL_nei = lmm.diago.profile.likelihood(tau=res_null2$tau, s2=res_null2$sigma2, Y = Y, X = X, eigenK = eiKsn)[1,1]
  
  X = model.matrix(~self[i,]+nei[i,]+pheno$Bolting+scale(pheno$InitLeafLen)+factor(pheno$Block)+pheno$edge)
  LL_nei2 = lmm.diago.profile.likelihood(tau=res_null2$tau, s2=res_null2$sigma2, Y = Y, X = X, eigenK = eiKsxn)[1,1]

  X = model.matrix(~self[i,]+nei[i,]+int[i,]+pheno$Bolting+scale(pheno$InitLeafLen)+factor(pheno$Block)+pheno$edge)
  LL_int = lmm.diago.profile.likelihood(tau=res_null2$tau, s2=res_null2$sigma2, Y = Y, X = X, eigenK = eiKsxn)[1,1]

  res = lmm.diago(Y = Y, X = X, eigenK = eiKsn, verbose = F)
  
  p_self = pchisq(-2*(LL_null-LL_self),1,lower.tail = F)
  p_nei = pchisq(-2*(LL_self2-LL_nei),1,lower.tail = F)
  p_int = pchisq(-2*(LL_nei2-LL_int),1,lower.tail = F)
  
  p_list_lm = rbind(p_list_lm, c(res$BLUP_beta[2:3],p_self,p_nei,p_int))
}
colnames(p_list_lm) = c("beta_self", "beta_nei", "p_self", "p_nei", "p_int")
write.csv(p_list_lm, paste0("./output/lmmGastonJPN2017ScoreMAF5_", args1, ".csv"), row.names = F)

