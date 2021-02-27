# sciPath
Improving Single-Cell RNA-seq Clustering by Integrating Pathways

## Descriptions
We designed a framework (sciPath) to study the accuracy and robustness of existing single-cell clustering by integrating pathways, including 10 state-of-art single-cell clustering methods and 4 pathway databases, pathway integration method and a complete set of validation indicators.

## Preparations
**1. Datasets**  
Demo datasets are saved in `".//Demo_data"`, including scRNA-swq matrix `(".//Demo_data//matrix")`, pathway `(".//Demo_data//pathway")` and cell labels `(".//Demo_data//label")`.

**2. Package installation**  
The installation code of the offline package and the online package are saved in `".//package//package_install.R"`.
 
## Codes
**1. clustering_by_gene_only.R**  
Single cell clustering that only consider gene-level information, including 'kmeans','hierarchical','spectral','DBSCAN','SC3','Seurat','CIDR','pcaReduce','SOUP','SNN-Cliq'. 
  
**2. clustering_by_integreting_pathway.R**  
Single cell clustering integrating pathway-level information, including 10 clustering methods, pathway score method (AUCell, function `"pathway_scoring"`) and integration method (SNF, function `"integrating_pathway"`).
  
**3. data_simulation_with_noise.R**  
The generation of three different expression profiles with noise, including "randomly set to 0" (function `"simulation_0"`), "gaussian noise" (function `"simulation_gaussian"`) and "randomly amplify" (function `"simulation_gaussian"`)

**4. evaluation.R**  
The evaluation of accuracy and robustness of clustering methods, including ARI(function `"evaluation_ARI"`), NMI(function `"evaluation_NMI"`), MES(function `"evaluation_MES"`) and AUC(function `"evaluation_AUC"`).

## Reference
Zhang C, Gao L, Wang B, Gao Y. Improving Single-Cell RNA-seq Clustering by Integrating Pathway. Briefings in Bioinformatics. (Revision)
