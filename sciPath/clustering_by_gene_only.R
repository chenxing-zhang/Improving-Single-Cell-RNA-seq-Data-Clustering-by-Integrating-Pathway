library(SNFtool) # SNF;spectralClustering
library(GSEABase) # getGmt, load pathway information
library(AUCell) # AUCell, pathway scoring method 
library(SingleCellExperiment)

# clustering method
library(stats) # kmeans
library(fastcluster) # fastcluster::hclust
library(dbscan) # dbscan
library(wordspace) # dist.matrix, fast distance calculation function
library(SC3) # SC3
library(Seurat)# Seurat
library(cidr) # CIDR
library(pcaReduce) # pcaReduce
library(SOUP) # SOUP
source("..//SOUP_ori.R") 
library(reticulate) # for python SNN-Cliq
py_config() # config python
source_python("..//SNN-Cliq.py") # SNN-Cliq


# load cell label
load_label <- function(path,name){
  label = read.table(paste(path,name,sep='\\'),sep='\t')
  return(as.vector(label$'cell_type1'))
}

# load single cell expresion matrix
load_matrix <- function(path, name){
  expr_matrix = read.table(paste(path,name,sep='\\'),sep='\t')
  expr_matrix = as.matrix(expr_matrix)
  return(expr_matrix)
}

clustering_by_gene<- function(mat_gene,name,k){
  dis_gene = dist.matrix(t(mat_gene),method = "euclidean",as.dist = TRUE)
  m_list = c('kmeans','hierarchical','spectral','DBSCAN',
             'SC3','Seurat','CIDR','pcaReduce','SOUP','SNN-Cliq')
  for(m in m_list){
    switch(which(m_list == m),
           # 1.keams
           {c = lapply(seq(1:5),kmeans,x=dis_gene,centers=k); 
           clust_results = lapply(c,FUN=function(x) x$cluster)},
           
           # 2.hierarchical
           {gcl <- fastcluster::hclust(dis_gene, method = 'ward.D');
           clust_results <- cutree(gcl, k)},
           
           # 3.spectral
           {K = ceiling(dim(mat_gene)[2]/10)
           sim_gene = affinityMatrix(as.matrix(dis_gene), K, 0.5);
           clust_results = spectralClustering(sim_gene, k)},
           
           # 4.DBSCAN
           {b = sort(kNNdist(dis_gene,5)); # dbscan::kNNdistplot(as.dist(1/w_gene), k =  5)
           eps_vec = b[ceiling(seq(0.5,0.9,by=0.1)*length(b))];
           d = lapply(eps_vec,dbscan,x=dis_gene);
           clust_results = lapply(d,FUN=function(x) x$cluster)},
           
           # 5.SC3
           {sce <- SingleCellExperiment(assays = list(counts = mat_gene,logcounts = mat_gene),colData = colnames(mat_gene))
           rowData(sce)$feature_symbol <- rownames(sce)
           a <- sc3(sce, ks = k, biology = FALSE, gene_filter = FALSE,kmeans_iter_max = 100) #SC3
           clust_results = as.numeric(as.vector(a@colData[[paste('sc3',as.character(k),'clusters',sep='_')]]))
           gc()},
           
           # 6. Seurat
           {pbmc <- CreateSeuratObject(counts=mat_gene)
           pbmc <- FindVariableFeatures(object = pbmc)
           pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
           pbmc <- ScaleData(object = pbmc, vars.to.regress = "percent.mt")
           pbmc <- RunPCA(object = pbmc)
           pbmc <- FindNeighbors(object = pbmc, reduction = "pca")
           pbmc <- FindClusters(object = pbmc)
           clust_results= as.numeric(as.vector(pbmc$seurat_clusters))},
           
           # 7. CIDR
           {sData <- new("scData", tags = mat_gene, tagType = 'CPM')
           sData@sampleSize <- ncol(sData@tags)
           sData@librarySizes <- rep(1e+06, sData@sampleSize)
           sData@nData <- mat_gene
           sData <- determineDropoutCandidates(sData)
           sData <- wThreshold(sData,cutoff=0.5, plotTornado=FALSE)
           sData <- scDissim(sData)
           sData <- scPCA(sData,plotPC = FALSE)
           sData <- nPC(sData)
           sData <- scCluster(sData,nCluster =k)
           clust_results= as.numeric(sData@clusters)},
           
           # 8. pcaReduce
           {Output_S <- PCAreduce(t(mat_gene), nbt=5, q=k-1, method='S')
           clust_results = lapply(seq(1:5),FUN=function(n) as.numeric(Output_S[[n]][,1]))},
           
           # 9. SOUP
           {mat_gene_t = t(mat_gene)
           spca.out = SPCAselect(mat_gene_t, type="log")
           spca.genes = spca.out$select.genes
           log.select.expr = mat_gene_t[, colnames(mat_gene_t) %in% spca.genes]
           k_soup = min(nrow(log.select.expr)-1, ncol(log.select.expr)-1, k)
           soup.out = try(SOUP(log.select.expr, Ks=k_soup, type="log"))
           if('try-error' %in% class(soup.out)){soup.out  = SOUP(log.select.expr,type="log")}
           soup.labels = soup.out$major.labels[[1]]
           clust_results = as.numeric(soup.labels)},
		   
		   # 10. SNN-Cliq
		   {clust_results = snn_cliq(as.matrix(dis_gene), 0.5, 0.7)}
    )
  }
  return(clust_results)
}


# main function
# cName: Name of clustering method
#     'kmeans','hierarchical','spectral','DBSCAN',
#     'SC3','Seurat','CIDR','pcaReduce','SOUP','SNN-Cliq'
# scName: Name of singel cell dataset
#     'yan', 'biase'
cName = 'DBSCAN'
scName= 'yan'
labelPath = '..\\Demo_data\\label'
scPath = '..\\Demo_data\\matrix'
paPath = "..\\Demo_data\\pathway"

main<-function(cName, paName, scName,s,
               labelPath, paPath, scPath){
  
  # load single cell data
  mat_gene = load_matrix(scPath, scName)
  
  # load cell label
  label = load_label(labelPath, paste(scName,'.rds_label',sep=''))
  label_int = as.numeric(as.vector(factor(label,levels=unique(label),labels=seq(1:length(unique(label))))))
  k = length(unique(label_int)) # real k
  
  # clustering
  clust_results = clustering_by_gene(mat_gene,scName,k)
  print(clust_results)
}






