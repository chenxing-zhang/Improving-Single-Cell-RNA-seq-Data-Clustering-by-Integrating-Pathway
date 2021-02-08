# SNFtool
install.packages("SNFtool")

install.packages("BiocManager")
# GSEABase 
library("BiocManager")
BiocManager::install("GSEABase")

# AUCelln
BiocManager::install("AUCell")

# fastcluster
install.packages("fastcluster")
# ARI, NMI
install.packages("aricode") 

# dbscan
install.packages("dbscan")

# wordspace
install.packages("wordspace")

# cor 
BiocManager::install("impute")
BiocManager::install("preprocessCore")
BiocManager::install("GO.db") # de novo pathway
install.packages("WGCNA")

# SingleCellExperiment
BiocManager::install("SingleCellExperiment")

# SC3
BiocManager::install("SC3")

# Seurat
install.packages("Seurat")

# cidr
install.packages("devtools")
library("devtools")
install.packages("ade4")
install.packages("clusterCrit")
install.packages("minpack.lm")
install.packages("NbClust")
install.packages("RcppParallel")
devtools::install_github("VCCRI/CIDR")

# pcaReduce
# 'pcaMethods', 'mnormt', 'mclust
BiocManager::install("pcaMethods")
install.packages("mnormt")
install.packages("mclust")
# "pcaReduce_1.0.tar.gz" can be downloaded from https://github.com/JustinaZ/pcaReduce

# SOUP
# 'limSolve', 'quantreg', 'PMA', 'ggpubr', 'ape', 'princurve'
install.packages("limSolve")
install.packages("quantreg")
install.packages("PMA")
install.packages("ggpubr")
install.packages("ape")
install.packages("princurve")
library("devtools")
devtools::install_github("lingxuez/SOUP")

# SNN-Cliq
install.packages("reticulate")


