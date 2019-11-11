##############################
### GO enrichment analysis ###
##############################

# 16-Oct-2019
# This R script was written for 
# Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana
# which was co-authored by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano
# correspondence: sato.yasuhiro.36c@kyoto-u.jp

### load a library and functions
library(GO.db)

# multiple fisher test --------------------------
ng.mft = function(
  cgt, #output from ng.MakeContigGOidTable, [,1]:"locus", [,2]:"GOid"
  gn.test, #contig names for test
  alternative="greater"
){
    
  gid.u = unique(cgt[,"GOid"])
  
  ft.in = matrix(0, nrow=length(gid.u), ncol=9)
  colnames(ft.in) = c("xtt", "xft", "xtf", "xff", "xnt", "xnf", "xtn", "xfn", "xnn")
  rownames(ft.in) = gid.u
  
  #               gn.test
  #             TRUE FALSE
  #Group  TRUE   xtt   xft   xnt
  #      FALSE   xtf   xff   xnf
  #              xtn   xfn   xnn
  
  ft.in[,"xnn"] = length(unique(cgt[, "locus"]))
  
  gn.pp.gid = table(cgt[, "GOid"])
  ft.in[names(gn.pp.gid), "xnt"] = gn.pp.gid
  ft.in[,"xnf"] = ft.in[,"xnn"] - ft.in[,"xnt"]
  
  ft.in[,"xtn"] = length(intersect(gn.test, unique(cgt[, "locus"])))
  ft.in[,"xfn"] = ft.in[,"xnn"] - ft.in[,"xtn"]
  
  gsea.test = cgt[is.element(cgt[,"locus"], gn.test), ]
  gn.test.gid = table(gsea.test[, "GOid"])
  ft.in[names(gn.test.gid), "xtt"] = gn.test.gid
  
  ft.in[,"xtf"] = ft.in[,"xtn"] - ft.in[,"xtt"]
  ft.in[,"xft"] = ft.in[,"xnt"] - ft.in[,"xtt"]
  ft.in[,"xff"] = ft.in[,"xnf"] - ft.in[,"xtf"]
    
  #Fisher's exact test
  fr = rep(1, nrow(ft.in))
  dt = rep(1, nrow(ft.in))
  for(i in 1:nrow(ft.in)){
    start = Sys.time()
    if(ft.in[i,"xtn"] > 1 && ft.in[i,"xnt"] > 1){ 
      contable = matrix(ft.in[i, 1:4], ncol=2)
      tmp = fisher.test(contable, alternative = alternative)
      fr[i] = tmp$p.value
    } else {
    }
    end = Sys.time()
    dt[i] = end - start
  }
  
  out = cbind(fr, ft.in, dt)
  colnames(out) = c("p.value", colnames(ft.in), "time")
  rownames(out) = rownames(ft.in)
    
  return(out)
  
}


# get GO terms -----------------------------
ng.GetGOTerms = function(GOid){
  
  out = NULL
  for(i in GOid){
    tmp = try(get(i, GOTERM), silent=TRUE)
    if(class(tmp)=="try-error"){
      out = c(out, "NA")
    } else {
      out = c(out, Term(tmp))
    }
  }
  return(out)
}


# prep. GO test output table ------------------------------
ng.prepGOtestOutTable = function(r, alpha=0.01){
  
  adp = r[,"p.value"] #if adjust p-value, use "p.adjust(r[,"p.value"], method="BH")" instead
  
  tmp.id = rownames(r)[adp < alpha]
  tmp.adp = adp[adp < alpha]
  tmp.description = ng.GetGOTerms(tmp.id)
  tmp.xnn = r[adp < alpha, c("xtt", "xtn", "xnt", "xnn")]
  
  out = cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  colnames(out) = c("p-value", "ID", "Description",  "A & B", "A",	"B", "U")
  
  return(out)
}

##########


# load a set of AGI codes from clipboard
gl = read.table("clipboard",header=F)
gl = gl[,1]

# load a list of AGI and GO. 
# The input file available at https://github.com/naganolab/AthRNAseq2016Zurich_Sato_et_al
load("./180907_GO/ulg.TAIR_180907") 

# Fisher tests with FDR
fisher.res = ng.mft(ulg, gl)
GO.list = ng.prepGOtestOutTable(fisher.res[(fisher.res[,"xtt"]>20)&(fisher.res[,"xtt"]<200),],alpha=0.05)
GO.list[order(as.numeric(GO.list[,1])),]
