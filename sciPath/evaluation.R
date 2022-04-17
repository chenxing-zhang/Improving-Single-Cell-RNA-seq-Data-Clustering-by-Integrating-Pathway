# evaluation
library(aricode) # ARI, NMI
  # for MES
library(Seurat) # FindAllMarkers, find markers
library(plyr) # ddply
library(clusterProfiler) # GSEA
library(DOSE) # GSEA_internal

# evaluation, accuracy, ARI
evaluation_ARI <- function(clust_results,label_int){
  return(ARI(clust_results,label_int))
}

# evaluation, accuracy, NMI
evaluation_NMI <- function(clust_results,label_int){
  return(NMI(clust_results,label_int))
}

# evaluation, accuracy, MES
# thre: the top ? gene are used to as marker gene,  c(25,50,100)
evaluation_MES <- function(mat_gene,clust_results,label_int,thre){
  # identify marker based on cell label
  pbmc_small <- CreateSeuratObject(counts=mat_gene)
  b = factor(x=label_int)
  names(b) = colnames(mat_gene)
  pbmc_small@active.ident = b
  markerTable_label = FindAllMarkers(pbmc_small,logfc.threshold=0.0, return.thresh=1) # all marker
  
  # identify marker based on clustering result
  pbmc_small <- CreateSeuratObject(counts=mat_gene)
  b = factor(x=clust_results)
  names(b) = colnames(mat_gene)
  pbmc_small@active.ident = b
  markerTable_cluster = FindAllMarkers(pbmc_small,logfc.threshold=0.0, return.thresh=1) # all marker
  
  # MES
  cellNameList = unique(labels)
  eScore_table = data.frame()
  
  for(i in seq(1:length(cellNameList))){
    cellName = cellNameList[i]
    markerGenes = markerTable_label[(markerTable_label['cluster']==i)&(markerTable_label['p_val_adj']<=0.05)&(abs(markerTable_label['avg_logFC'])>1),'gene']
    markerGenes = as.character(markerGenes)
    if(length(markerGenes)>thre){markerGenes=markerGenes[1:thre]}
    markerGenesSet = data.frame(cellName = rep(cellName,length(markerGenes)),
                                markerGenes = markerGenes)
    
    # GSEA, markerGenes to markerTable_cluster
    tempSNF_table = data.frame()
    for(clustName_SNF in unique(markerTable_cluster[['cluster']])){
      DEclust_SNF_table = markerTable_cluster[markerTable_cluster['cluster']==clustName_SNF,]
      if(dim(DEclust_SNF_table[DEclust_SNF_table[['gene']]%in%markerGenes,])[1]==0){next}
      DEclust_SNF_geneList = -log(DEclust_SNF_table[['p_val']]+10e-100)
      names(DEclust_SNF_geneList) = DEclust_SNF_table[['gene']]
      b = GSEA(DEclust_SNF_geneList, TERM2GENE = markerGenesSet, minGSSize = 1, pvalueCutoff = 1,verbose='FALSE')
      if(dim(data.frame(b))[1]==0){next}
      eScore_SNF = b@result$enrichmentScore
      esPvalue_SNF = b@result$p.adjust
      temp_table = data.frame(scName=scName,Method=method,ARI_i=ARI_i,type='gene_pathway',cellName=cellName,clustName=clustName_SNF,eScore=eScore_SNF,espvalue=esPvalue_SNF)
      tempSNF_table = rbind(tempSNF_table,temp_table)
    }
    tempSNFMax_table = ddply(.data = tempSNF_table, .variables = .(cellName), subset, eScore == max(eScore))
    eScore_table = rbind(eScore_table,tempSNFMax_table)
  }
  MES = mean(eScore[['eScore']])
  return(MES)
}

# evaluation, robustness, AUC
# accuracy_List, A vector of accuracy indicators, arranged in ascending or descending order of noise proprotion
evaluation_AUC<- function(accuracy_List){
  len = length(accuracy_List)
  AUC = (2 * sum(accuracy_List) - accuracy_List[1] - accuracy_List[len]) / (2.0 * len)
  return(AUC)
}
