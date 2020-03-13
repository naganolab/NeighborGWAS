##############################
# Examples of 2D Ising model #
##############################

# 23-Feb-2020
# This R script was written for 
# Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana
# which was co-authored by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano
# correspondence: sato.yasuhiro.36c@kyoto-u.jp

# set coefficients
J = 0.28 #= beta_2
h = 0.14 #= beta_1

# set a field
rect1 = 50 #rect1 = 10
rect2 = 50 #rect2 = 40 for 10 x 40 grid for Figure 5
N = rect1*rect2
field = matrix(c(1:N),nrow=rect1,ncol=rect2)
Xi = sample(c(-1,1),N,replace = T)

# function to count neighbors
nei_states = function(id,s) {
  s_seq = sort(seq(-s,s,by=1),decreasing = T)
  place = which(field==id,arr.ind=TRUE)
  n = 0
  for(i in s_seq) {
    for(j in s_seq) {
      nei_id = try(field[place[1]+i,place[2]+j],silent=T)
      if((class(nei_id)!="try-error")) { # skip if calling below the field
        if(length(nei_id)==1) { # skip if calling above the field
          if(nei_id!=id) { # excl. a focal plant itself
            n = n + Xi[nei_id]
          }
        }
      }
    }
  }
  coval = Xi[id]*n
  print(id); print(n)
  return(coval)
}

xixj = mapply(nei_states,1:N,s=1)
Ei_t0 = -(J*xixj+h*Xi)

# MCMC by Gibbs sampling
for(j in 1:1000) {
  Xi_t0 = Xi
  xixj_t0 = xixj
  
  perm = sample(1:N,N)
  for(i in perm) {
    Xi[i] = 2*sample(0:1,1)-1
    xixj[i] = nei_states(i,s=1)
    Ei = -(J*xixj[i]+h*Xi[i]) # set the energy negative to minimize itself
    
    # Metropolis algorithm
    if(Ei_t0[i]<Ei) {
      Xi[i] = Xi[i]; xixj[i] = xixj[i]; Ei_t0[i] = Ei
    } else if(runif(1,0,1)<exp(Ei-Ei_t0[i])) {
      Xi[i] = Xi[i]; xixj[i] = xixj[i]; Ei_t0[i] = Ei
    } else {
      Xi[i] = Xi_t0[i]; xixj[i] = xixj_t0[i]
    }
  }
  E = sum(+J*(xixj)+h*Xi)
  image(matrix(Xi,rect1,rect2),main=paste(j,E),col=c("black","white"))
}

svg(file=paste0("IsingJ",J,"_h",h,".svg"),width = 10, height = 4)
image(t(matrix(Xi,rect1,rect2)),main=paste(j,E),col=c("black","white"))
dev.off()
