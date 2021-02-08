# data simulation with noise

# simulation, randomly set 0
# p: the proprotion of random selection, c(0.05,0.10,0.15,0.20)
# t: simulation times, c(1,2,3)
simulation_0 <- function(mat_gene,p,t){
  p = as.numeric(p)
  t = as.numeric(t)
  set.seed(p*100+t)
  index = matrix(sample(c(TRUE,FALSE), length(mat_gene), p = c(1,1-p), replace=TRUE), nrow=nrow(mat_gene))
  newMat = mat_gene
  colnames(newMat) = colnames(mat_gene)
  rownames(newMat) = rownames(mat_gene)
  rm(mat_gene,gMean,gSd);gc()
  newMat[index==TRUE]=0
  return(newMat)
}

# simulation, gaussian
# p: noise proprotion, c(0.05,0.10,0.15,0.20)
# t: simulation times, c(1,2,3)
simulation_gaussian <- function(mat_gene,p,t){
  p = as.numeric(p)
  t = as.numeric(t)
  set.seed(p*100+t)
  cNum = length(colnames(mat_gene))
  gMean = apply(mat_gene,1,mean)
  gSd = apply(mat_gene,1,sd)
  noise = t(sapply(1:length(gMean),function(i) rnorm(cNum,mean=gMean[i],sd=gSd[i])))
  newMat = mat_gene*(1-p) + noise*(p)
  colnames(newMat) = colnames(mat_gene)
  rownames(newMat) = rownames(mat_gene)
  rm(mat_gene,gMean,gSd);rm(noise);gc()
  newMat[newMat < 0] <- 0
  return(newMat)
}

# simulation, amplify
# p: the proprotion of random selection,  0.05 or 0.01
# f: amplify factor, c(5,10,50,100)
# t: simulation times, c(1,2,3)
simulation_amplify <- function(mat_gene,p,f,t){
  p = as.numeric(p)
  f = as.numeric(f)
  t = as.numeric(t)
  set.seed(p*1000+f*10+t)
  index = matrix(sample(c(TRUE,FALSE), length(mat_gene), p = c(p,1-p), replace=TRUE), nrow=nrow(mat_gene_ori))
  newMat = mat_gene
  colnames(newMat) = colnames(mat_gene)
  rownames(newMat) = rownames(mat_gene)
  rm(mat_gene,gMean,gSd);gc()
  newMat[index==TRUE] = newMat[index==TRUE]*f
  return(newMat)
}