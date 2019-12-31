#################################
# GWAS Manhattan plot & QQ plot #
#################################

# 16-Oct-2019
# This R script was written for 
# Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana
# which was co-authored by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano
# correspondence: sato.yasuhiro.36c@kyoto-u.jp

# set a file path
folder_path = "../AthGWAS/lmekinGWAS/output/"
file_name = "lmmGastonJPN2017ScoreMAF5"
file_path = paste0(folder_path,file_name,".csv")

# load file
gwas_out = read.csv(file_path, header=T)
position = gwas_out[,1:4]
colnames(position) = c("X","chr","pos","maf")
position$maf[position$maf>0.5] = 1-position$maf[position$maf>0.5]

p_adj = 0.05/(nrow(gwas_out))
chr_rep = cumsum(table(position$chr))
cols2 = c(rgb(1,0,0, 2*position$maf[1:chr_rep[1]]), rgb(0,1,0, 2*position$maf[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*position$maf[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*position$maf[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*position$maf[(chr_rep[4]+1):(chr_rep[5])]))

# (1) depicting Manhattan plots
png(paste0(file_name,".png"), width=7, height=8, res=600, units="in")
par(mfcol=c(2,1))

plot(-log(gwas_out$P_self,10), xlab="position", ylab="-log10(P_self)", col=cols2, pch=16, las=1)
abline(h = -log(p_adj,10),lty=2,col=grey(0.5,0.5)) #Bonferroni correction
plot(-log(gwas_out$P_nei,10), xlab="position", ylab="-log10(P_nei)", col=cols2, pch=16, las=1)
abline(h = -log(p_adj,10),lty=2,col=grey(0.5,0.5)) #Bonferroni correction

dev.off()

# (2) depicting QQ-plots
png(paste0(file_name,"_QQ.png"), width=7, height=8, res=600, units="in")
par(mfcol=c(2,1))

o = -log(sort(gwas_out$P_self,decreasing=F),10)
e = -log(ppoints(length(gwas_out$P_self)),10)
plot(e, o, pch=16, las=1, main=file_name, ylab="obs. -log10(P_self)", xlab="exp. -log10(P)", col=cols2[order(gwas_out$P_self,decreasing=F)])
abline(0,1,lty=2)

o = -log(sort(gwas_out$P_nei,decreasing=F),10)
e = -log(ppoints(length(gwas_out$P_nei)),10)
plot(e, o, pch=16, las=1, main=file_name, ylab="obs. -log10(P_nei)", xlab="exp. -log10(P)", col=cols2[order(gwas_out$P_nei,decreasing=F)])
abline(0,1,lty=2)

dev.off()

