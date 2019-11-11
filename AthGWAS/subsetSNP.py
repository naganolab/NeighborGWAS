####################################
# Subsetting the full imputed SNPs #
####################################

# 16-Oct-2019
# This Python script was written for 
# Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on Arabidopsis thaliana
# which was co-authored by Yasuhiro Sato, Eiji Yamamoto, Kentaro K. Shimizu and Atsushi J. Nagano
# correspondence: sato.yasuhiro.36c@kyoto-u.jp

# (1) load libraries
import h5py
import pandas as pd
import numpy as np
import csv

# (2) load data
# The input .hdf5 file was provided by Uemit Seren at GMI Molecular Plant Biology
# As recommended, the input file was handled using Python
input_file = "all_chromosomes_binary.hdf5"
h5file = h5py.File(input_file,"r")

snps = h5file["snps"].value
accs = h5file["accessions"].value
pos = h5file["positions"].value

acID = []
for i in accs:
    acID.append(int(i))

acID = np.array(acID)
gwasID = pd.read_csv("gwasIDlist.csv") # load gwasIDs

acc_list = []
for i in gwasID["GWASid"]:
    place = np.where(acID == i) # search a focal accession
    acc_list.append(int(place[0]))
    
# (2) export the subset data
sub_snps = snps[:,acc_list]
del snps
sub_snps = pd.DataFrame(sub_snps)
sub_snps.to_csv("sub_snps.csv") # export SNPs

pos = pd.DataFrame(pos)
pos.to_csv("positions.csv") # export positions
