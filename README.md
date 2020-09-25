# NeighborGWAS
R source codes, accession list, and phenotype data used in "Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana" by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano. R package version is available at https://cran.r-project.org/package=rNeighborGWAS.  

## simulation
neighbor GWAS of simulated phenotypes
1. Ising2Dexample.R  
R source code to show examples of spatial arrangements in a 2D Ising model. This works independently from 2.-4.  
2. neighborGWASsimul_prepGeno.R  
R source code to subset the RegMap data  
3. neighborGWAS250kSimulation.R  
R source code to run simulations   
4. neighborGWAS250kSimulation.sh  
A shell script to control simulations from terminal. Current strategy is simultanous runnning with a single core  
#### /output
A folder to store output files  
1. AUC_output.R  
R source code to join multiple output files into one .csv file  

## AthGWAS
Application to A. thaliana herbivory data  
1. subsetSNP.R  
R source code to subset the full imputed SNP matrix  
2. reshapeSNP.R  
R source code to cut-off the SNP data with MAF and LD  
3. prepNeighbors.R  
R source code to reshape and prepare neighbor covariates 
4. gwasID.csv  
.csv list of genotype ID used in this study  
5. phenotype.csv  
.csv file including the leaf damage score with some covariates  
6. prep_AthNeighborGWAS.sh  
A shell script to prepare input files before running lmm  

### /lmekinGWAS
linear mixed model with a kinship structure using the gaston package in R  
1. prep_lmm.R  
R source code to prepare kinship matrices before running GWAS  
2. run_lmm.R  
R source code to estimate variance components and then to perform association tests  
3. run_lmm.sh  
A shell script to control lmm iterations from terminal. Current strategy is simultanous runnning with a single core  
#### /output
A folder to store output files  
1. lmm_output.R  
R source code to combine multiple output files into one .csv file  

## postGWAS
postGWAS analyses
1. Manhattan_lmm.R  
R source code to depict Manhattan plots and QQ-plots  
2. GeneList_lmm.py  
Python code to make a list of candidate genes from GWAS output and .gff file  
