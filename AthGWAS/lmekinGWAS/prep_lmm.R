######################################
#prepare kinship matrices before GWAS#
######################################

# 23-Feb-2020
# This R script was written for 
# Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana
# which was co-authored by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano
# correspondence: sato.yasuhiro.36c@kyoto-u.jp


library(Matrix); library(gaston)
setThreadOptions(numThreads = 1)

# load data
load("../lmmJPN2017selfMAF5.RData")
load("../lmmJPN2017neiMAF5.RData")
load("../lmmJPN2017phenoMAF5.RData")
int = self*nei # set X for testing beta_12

# for variation partitioning with lmm.aireml
selfCross = crossprod(self)
neiCross = crossprod(nei)
intCross = crossprod(int)

selfK = ((nrow(self)-1)+selfCross)/(2*(nrow(self)-1))
neiK = neiCross/(nrow(nei)-1)
intK = intCross/(nrow(int)-1)

# approximate the kinship matrices as a near postive-definite matrix
selfPDK = as.matrix(nearPD(selfK,maxit=10^6)$mat)
neiPDK = as.matrix(nearPD(neiK,maxit=10^6)$mat)
intPDK = as.matrix(nearPD(intK,maxit=10^6)$mat)

save(selfPDK,file="selfPDK.RData")
save(neiPDK,file="neiPDK.RData")
save(intPDK,file="intPDK.RData")

eiKs = eigen(as.matrix(selfK))
save(eiKs,file="eiKs.RData")

############################
