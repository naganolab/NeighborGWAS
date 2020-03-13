############################################
### Simulating the power of neighborGWAS ###
############################################

# 23-Feb-2020
# This R script was written for 
# Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana
# which was co-authored by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano
# correspondence: sato.yasuhiro.36c@kyoto-u.jp


# (1) load libraries
library(Matrix)
library(gaston)
setThreadOptions(numThreads = 1)
library(pROC)


# (2) set functions developed for this simulation

# 'qtl_pheno_simu' simulates phenotypic values.
# Arguments
# g : n x 1 genotype vector to design major additive genetic effect.
# sigma : a vector of explanatory variable for beta_2 (non-genic effects)
# eigen_K.. : products of function 'eigen' with covariance matricies that were used as explanatory variables for the phenotype.
# tau_ratio : ratio for contribution of engenK_self,  eigenK_nei and eigenK_sxn to the phenotype.
# pveB : proportion of variance explained by genetic effect designed by 'g' vector.
# pveM : proportion of variance explained by all genetic effects (i.e., g and eigenK_..)
# print.PVE : print the realized PVEs. Sometimes, the calculations do not work well and this is useful to check the result
qtl_pheno_simu = function(g, sigma, eigenK_self, eigenK_nei, eigenK_sxn, tau_ratio = c(1,1,1), pveB, pveM, print.PVE = FALSE)
{
  beta_g = g / sd(g)
  beta_sigma = sigma / sd(sigma)
  Maj_eff = beta_g + beta_sigma
  
  omega_self = eigenK_self$vectors %*% rnorm(nrow(eigenK_self$vectors), sd = sqrt(tau_ratio[1]*eigenK_self$values))
  omega_nei = eigenK_nei$vectors %*% rnorm(nrow(eigenK_nei$vectors), sd = sqrt(tau_ratio[2]*eigenK_nei$values))
  omega_sxn = eigenK_sxn$vectors %*% rnorm(nrow(eigenK_sxn$vectors), sd = sqrt(tau_ratio[3]*eigenK_sxn$values))	
  omega_all = omega_self + omega_nei + omega_sxn
  omega_all = omega_all / sd(omega_all)
  adj.omega = function(adj, Maj_eff, omega_all, pveB, pveM)	{
    (pveB/pveM - (var(Maj_eff) / (var(Maj_eff + (omega_all / adj)))))^2
  }#adj.omega()
  adj = optimize(adj.omega , interval=c(0.0000001, 10000000), Maj_eff, omega_all, pveB, pveM)$minimum
  adj.omega_all = omega_all / adj
  Model_effect = Maj_eff + adj.omega_all
  
  epsilon = rnorm(nrow(g), sd = 1)
  pveM.adj = function(adj, Model_effect, epsilon, pveM)	{
    (pveM - (var(Model_effect) / (var(Model_effect + (epsilon / adj)))))^2
  }#pveM.adj()
  adj = optimize(pveM.adj, interval=c(0.0000001, 10000000),  Model_effect, epsilon, pveM)$minimum
  adj.epsilon = epsilon / adj
  
  y = Model_effect + adj.epsilon
  
  realized_pveB = var(Maj_eff) / var(y)
  realized_pveM = var(Model_effect) / var(y)
  
  if (print.PVE) {
    print(c(paste("pveB =",realized_pveB),paste("pve.beta_g =",as.double(var(beta_g)/var(y))),paste("pve.beta_sigma =",as.double(var(beta_sigma)/var(y))),paste("cor.g_x_sigma =",cor(beta_g,beta_sigma))))
  }
  return(list(y=y, beta_g=beta_g, beta_sigma=beta_sigma, omega_slef=omega_self, omega_nei=omega_nei, omega_sxn=omega_sxn, epsilon=epsilon, res_pveB = rep(realized_pveB,length(y)), res_pveM = rep(realized_pveM,length(y))))	
}#qtl_pheno_simu()

# 'real_geno_simu' simulates power of Neighbor GWAS
# Arguments
# s & a : parameters for distant decay of gene effects
# n_causal : number of causal genes
# rect : experimental plot deign (rect * rect square matrix)
# following argmuments are used for function 'qtl_pheno_simu'
# tau_ratio : ratio for contribution of engenK_self,  eigenK_nei and eigenK_sxn to the phenotype.
# pveB : proportion of variance explained by genetic effect designed by 'g' vector.
# pveM : proportion of variance explained by all genetic effects (i.e., g and eigenK_..)
real_geno_simu = function(s_max, a, n_causal, tau_ratio=c(1,1,1), pveB, pveM, rect)
{
  
  perm = sample(colnames(g_bin),rect*rect,replace = FALSE)
  g_self = g_bin[,perm]
  
  field = matrix(rep(NA,rect*rect*4),rect*2,rect*2)
  for(i in 2*c(1:rect)) {
    field[seq(2,rect*2,2),i] = seq((i-1)*rect+1,i*rect,1)
  }
  for(i in 2*c(1:rect)-1) {
    field[seq(1,rect*2,2),i] = seq((i-1)*rect+1,i*rect,1)
  }
  
  field = field[,1:36]
  
  p = nrow(g_self)
  n = ncol(g_self)
  maf = apply(g_self,1,sum)/n
  
  g_self[g_self==0] = -1
  K = crossprod(g_self)
  K = (p/2+K/2)/p
  
  true = sample(c(1:p),n_causal,replace=F)
  
  b_int = rep(0,p)
  
  both_true = sample(true,round(n_causal*0.3))
  int_true = sample(both_true,round(length(both_true)*0.5))
  b_int[int_true] = runif(length(int_true),0.5,2)*sample(c(-1,1),length(int_true),replace=T)
  
  main_true = true[is.element(true,both_true)==FALSE]
  mainTF = sample(c(TRUE,FALSE),length(main_true),replace=TRUE)
  
  b_self = rep(0,p)
  self_true = c(both_true, main_true[mainTF==TRUE])
  b_self[self_true] = runif(length(self_true),0.5,2)*sample(c(-1,1),length(self_true),replace=T)
  
  b_nei = rep(0,p)
  nei_true = c(both_true, main_true[mainTF==FALSE])
  b_nei[nei_true] = runif(length(nei_true),0.5,2)*sample(c(-1,1),length(nei_true),replace=T)
  true = sample(c(1:p),n_causal,replace=F)
  
  nei_coval = function(id,s) {
    s_seq = sort(seq(-s,s,by=1),decreasing = T)
    place = which(field==id,arr.ind=TRUE)
    n = rep(0,p)
    k = 0
    for(i in s_seq) {
      for(j in s_seq) {
        nei_id = try(field[place[1]+i,place[2]+j],silent=T)
        if((class(nei_id)!="try-error")) {
          if(length(nei_id)==1) {
            if(is.na(nei_id)==FALSE) {
              if(nei_id!=id) {
                n = n + g_self[,nei_id]
                k = k + 1
              }
            }
          }
        }
      }
    }
    coval = g_self[,id]*n/k
    print(id); print(k)
    return(coval)
  }
  
  w = function(a,s) {
    w = exp(-a*(s-1))
    return(w)
  }
  
  nei_covalS = function(id,s) {
    s_seq = sort(seq(-s,s,by=1),decreasing = T)
    place = which(field==id,arr.ind=TRUE)
    n = rep(0,p)
    k = 0
    for(i in s_seq) {
      for(j in s_seq) {
        nei_id = try(field[place[1]+i,place[2]+j],silent=T)
        if((class(nei_id)!="try-error")) {
          if(length(nei_id)==1) {
            if(is.na(nei_id)==FALSE) {
              if(nei_id!=id) {
                d = max(abs(i),abs(j))
                n = n + w(a=a,s=d)*g_self[,nei_id]
                k = k + 1
              }
            }
          }
        }
      }
    }
    coval = g_self[,id]*n/k
    print(id); print(k)
    return(coval)
  }
  
  g_neiS = mapply(nei_covalS,1:n,s=s_max)
  g_neiS = t(g_neiS)
  g_neiS = (g_neiS-mean(g_neiS))/sd(g_neiS)
  
  g_self = t(g_self)
  g_self = (g_self-mean(g_self))/sd(g_self)
  g_selfxnei = (g_self*g_neiS)
  
  K_self = tcrossprod(g_self)/(p-1)
  K_nei = tcrossprod(g_neiS)/(p-1)
  K_sxn = tcrossprod(g_selfxnei)/(p-1)
  
  eigenK_self=eigen(K_self); eigenK_nei=eigen(K_nei); eigenK_sxn = eigen(K_sxn)
  
  sigma = g_neiS%*%b_nei + g_selfxnei%*%b_int
  pheno = qtl_pheno_simu(g=g_self%*%b_self, sigma=sigma, eigenK_self=eigenK_self, eigenK_nei=eigenK_nei, eigenK_sxn=eigenK_sxn, tau_ratio=tau_ratio, pveB=pveB, pveM=pveM)
  
  while(((round(pheno$res_pveB[1],2)==pveB)&(round(pheno$res_pveM[1],2)==pveM))==FALSE) {
    pheno = qtl_pheno_simu(g=g_self%*%b_self, sigma=sigma,eigenK_self=eigenK_self, eigenK_nei=eigenK_nei, eigenK_sxn=eigenK_sxn, tau_ratio=tau_ratio, pveB=pveB, pveM=pveM)
  }
  
  Y = pheno$y
  
  # s = 0 without neighbor effects
  res = lmm.aireml(Y=Y,K=list(K_self),verbose=T)
  tau = c(res$sigma2,res$tau,0)
  
  eiK = eigen(K_self)
  res_null = lmm.diago(Y=Y, eigenK=eiK, verbose = F)
  LL_null = lmm.diago.profile.likelihood(tau=res_null$tau, s2=res_null$sigma2, Y=Y, eigenK=eiK)[1,1] 
  
  res_null_lm = lm(Y~1)
  LL_null_lm = logLik(res_null_lm)
  
  p_list = c()
  for(i in 1:p) { print(i)
    LL =  lmm.diago.profile.likelihood(tau=res_null$tau, s2=res_null$sigma2, Y=Y, X=model.matrix(~g_self[,i]), eigenK=eiK)[1,1]
    p_LL = pchisq(-2*(LL_null-LL),1,lower.tail = F)

    LL_lm = logLik(lm(Y~g_self[,i]))
    p_LL_lm = pchisq(-2*(LL_null_lm-LL_lm),1,lower.tail = F)
    
    p_list = rbind(p_list, c(p_LL,p_LL_lm))
  }
  
  p_self = -log10(p_list[,1])
  p_self_lm = -log10(p_list[,2])
  ans_self = round(b_self,0)
  ans_self[ans_self!=0] = 1

  ans = data.frame(p_self,p_self_lm,ans_self)
  roc_self = roc(ans_self~p_self,data=ans)
  roc_self_lm = roc(ans_self~p_self_lm,data=ans)

  auc = roc_self$auc
  auclm = roc_self_lm$auc
  
  resList = c(0, a, n_causal, tau_ratio[1], tau_ratio[2], tau_ratio[3], pveB, pveM, tau, auc, auclm)
  # save(ans,file=paste0("ans_s",0,"_n",n_causal,"_pveB",pveB,"_alpha",alpha,".RData"),compress = T)

  #s = 1 to s_max
  for(s in 1:s_max) {
    g_self = g_bin[,perm]
    g_self[g_self==0] = -1
    
    g_nei = mapply(nei_coval,1:n,s=s)
    g_nei = t(g_nei)
    g_nei = (g_nei-mean(g_nei))/sd(g_nei)
    
    K_nei = tcrossprod(g_nei)/(p-1)
    
    g_self = t(g_self)
    g_self = (g_self-mean(g_self))/sd(g_self)
    
    res = lmm.aireml(Y=Y,K=list(K_self,K_nei),verbose=T)
    tau = c(res$sigma2,res$tau)
    
    eiK = eigen(res$tau[1]*K_self+res$tau[2]*K_nei)
    
    p_list = c()
    for(i in 1:p) { print(i)
      LL_self =  lmm.diago.profile.likelihood(tau=res_null$tau, s2=res_null$sigma2, Y=Y, X=model.matrix(~g_self[,i]), eigenK=eiK)[1,1]
      LL_nei =  lmm.diago.profile.likelihood(tau=res_null$tau, s2=res_null$sigma2, Y=Y, X=model.matrix(~g_self[,i]+g_nei[,i]), eigenK=eiK)[1,1]
      p_nei = pchisq(-2*(LL_self-LL_nei),1,lower.tail = F)
      
      LL_self_lm = logLik(lm(Y~g_self[,i]))
      LL_nei_lm = logLik(lm(Y~g_self[,i]+g_nei[,i]))
      p_nei_lm = pchisq(-2*(LL_self_lm-LL_nei_lm),1,lower.tail = F)
      
      p_list = rbind(p_list, c(p_nei,p_nei_lm))
    }
    
    p_nei = -log10(p_list[,1])
    p_nei_lm = -log10(p_list[,2])
    ans_nei = round(b_nei,0)
    ans_nei[ans_nei!=0] = 1
    
    ans = data.frame(p_nei,p_nei_lm,ans_nei)
    roc_nei = roc(ans_nei~p_nei,data=ans)
    roc_nei_lm = roc(ans_nei~p_nei_lm,data=ans)
    
    auc = roc_nei$auc
    auclm = roc_nei_lm$auc
    
    resList = rbind(resList, c(s, a, n_causal, tau_ratio[1], tau_ratio[2], tau_ratio[3], pveB, pveM, tau, auc, auclm))
	# save(ans,file=paste0("ans_s",s,"_n",n_causal,"_pveB",pveB,"_alpha",alpha,".RData"),compress = T)
  }
  return(resList)
}#real_geno_simu()


# (3) load arguents from terminal command

a = commandArgs(trailingOnly = TRUE)[1]
n_causal = commandArgs(trailingOnly = TRUE)[2]
tau1 = commandArgs(trailingOnly = TRUE)[3]
tau2 = commandArgs(trailingOnly = TRUE)[4]
tau3 = commandArgs(trailingOnly = TRUE)[5]
pveB = commandArgs(trailingOnly = TRUE)[6]
pveM = commandArgs(trailingOnly = TRUE)[7]
iter = commandArgs(trailingOnly = TRUE)[8]


# (4) load Arabidopsis RegMap data (chr1 & 2 only)

load("./Ara250kRegMap_chr12_MAF10.RData")


# (5) set parameters for simulation

s_max = 3
rect = 36
tau_ratio = c(as.numeric(tau1), 
            as.numeric(tau2), 
            as.numeric(tau3))


# (6) execute simulation

resList = c()
for(i in 1:10) {
  res = real_geno_simu(s_max = s_max,
                         a = as.numeric(a), 
                         n_causal = as.numeric(n_causal), 
                         tau_ratio = tau_ratio, 
                         pveB = as.numeric(pveB), 
                         pveM = as.numeric(pveM),
                         rect = rect)
  resList = rbind(resList, res)
}#i


# (7) return results
rownames(resList) = c()
colnames(resList) = c("s",
                      "alpha", 
                      "n_causal", 
                      "tau1", 
                      "tau2", 
                      "tau3", 
                      "pveB", 
                      "pveM",
                      "sigma2",
                      "tau_s",
                      "tau_n",
                      "AUC", 
                      "AUC_lm")
write.csv(resList, file = paste0("./output/AUCs_", iter,".csv"), row.names = FALSE)

