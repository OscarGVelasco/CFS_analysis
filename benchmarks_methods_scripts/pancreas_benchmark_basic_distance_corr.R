#####################
##
## https://github.com/MarioniLab/MouseGastrulationData/blob/devel/vignettes/MouseGastrulationData.Rmd
## Paper: https://www.nature.com/articles/s41587-023-01766-z#data-availability

# BiocManager::install("MouseGastrulationData")

# install.packages("igraph")
# BiocManager::install("Rcpp")
# BiocManager::install("dqrng")
# BiocManager::install(pkgs = c("RcppAnnoy","RcppProgress","FNN"))
# install.packages("Matrix-1.6.4 ", type = "source")
# BiocManager::install("irlba")
# install.packages("https://cran.r-project.org/src/contrib/Archive/uwot/uwot_0.1.16.tar.gz")
# install.packages("Seurat")
# install.packages("leidenAlg")

PATH <- "/home/o872o/o872o/cluster_similarity_sc/review/"
setwd(PATH)
library(dplyr)
library(scran)
library(bluster)
library(Seurat)

source("/home/o872o/o872o/cluster_similarity_sc/review/initialice.R")

# Load the Pancreas single-cell objects
fname <- load("/home/o872o/o872o/cluster_similarity_sc/rdata/Pancreas_sc_objects_list_processed.RData");fname
# Select common cell types for comparison
cellstypes <- names(table(sc_list_input[[1]]@colData$cell.type))[names(table(sc_list_input[[1]]@colData$cell.type)) %in% names(table(sc_list_input[[2]]@colData$cell.type))]

seu1 <- sc_list_input[[1]]
seu1 <- Seurat::as.Seurat(seu1[,seu1$cell.type %in% cellstypes],data = NULL)

seu2 <- sc_list_input[[2]]
seu2 <- Seurat::as.Seurat(seu2[,seu2$cell.type %in% cellstypes],data = NULL)

######################################
######################################
pancreas_list <- list(seu1, seu2)
samples = c("pancreas_large", "pancreas_small")
names(pancreas_list) <- samples
# Get variable features:
vfeat <- lapply(pancreas_list, function(x){Seurat::VariableFeatures(Seurat::FindVariableFeatures(x, nfeatures=4000))})
hvgs <- intersect(vfeat[[1]], vfeat[[2]])[1:1500]

k = 8
labels = "cell.type"
nGenes_all = c(25, 50, 100, 250, 500, 1000, 1500)
nGenes_all <- rep(nGenes_all, each=6)
names(nGenes_all) <- paste0("Sim_", seq_len(length(nGenes_all)))

labels = "cell.type"

resFile = paste0(PATH,"pancreas.benchmark.results.CORRS.res2.pc15.1.5kfeats.seed.04022025.Rds")

if (!file.exists(resFile)) {
  res = NULL
} else {
  res = readRDS(resFile)
}

BiocParallel::register(BPPARAM =  BiocParallel::MulticoreParam(workers = 16))
library(Seurat)
embeddings_names=""
######### SETTING SEED FOR REPRODUCIBILITY
set.seed(04022025)
######### SETTING SEED FOR REPRODUCIBILITY

for (sim in names(nGenes_all)) {
  
  nGenes = nGenes_all[sim]
  
  for (querySample in unique(samples)) {
    
    for (referenceSample in setdiff(samples,querySample)) {
      
      print(nGenes)
      print(paste0("Query: ", querySample))
      print(paste0("Reference: ", referenceSample))
      
      # check has this already been done?
      if (any(res[,"Sim"] == sim & 
              res[,"querySample"] == querySample & 
              res[,"referenceSample"] == referenceSample)) next
      
      # select HVGs from the referenceSCE and subset
      #decomp <- modelGeneVar(referenceSCE)
      #decomp <- Seurat::VariableFeatures(atlas_seu)[1:1500]
      #hvgs <- rownames(decomp)[decomp$mean > 0.05 & decomp$bio > 0]
      #hvgs <- decomp
      # genes = sample(hvgs, size = nGenes)
      # referenceSCE = atlas_seu[,atlas_seu$sample %in% referenceSample]
      genes = sample(hvgs)[seq_len(nGenes)]
      referenceSCE = pancreas_list[[referenceSample]]
      referenceSCE <- Seurat::NormalizeData(referenceSCE)
      referenceSCE <- Seurat::ScaleData(referenceSCE)
      #referenceSCE <- RunPCA(referenceSCE, features = VariableFeatures(object = referenceSCE), npcs = 30)
      referenceSCE <- referenceSCE[genes,]
      # select random genes from among the HVGs for the query data
      querySCE = pancreas_list[[querySample]]
      #querySCE = atlas_seu[, atlas_seu$sample %in% querySample]
      querySCE <- Seurat::NormalizeData(querySCE)
      querySCE <- Seurat::ScaleData(querySCE)
      querySCE <- Seurat::FindVariableFeatures(querySCE)
      querySCE <- RunPCA(querySCE, features = VariableFeatures(object = querySCE), npcs = 30)
      querySCE <- RunUMAP(querySCE, features = VariableFeatures(object = querySCE))
      querySCE <- FindNeighbors(querySCE, dims = seq(15), k.param = 20)
      querySCE <- FindClusters(querySCE, resolution = 2, method = 4)
      
      querySCE <- querySCE[genes,]
      # Compute base accuracy by cluster 
      querySCE@meta.data$assign = ""
      for (cls in unique(querySCE@meta.data$seurat_clusters)){
        ccount <- table(querySCE@meta.data$cell.type[querySCE@meta.data$seurat_clusters == cls])
        querySCE@meta.data$assign[querySCE@meta.data$seurat_clusters == cls] <- names(which.max(ccount))
      }
      sum(querySCE@meta.data$assign == querySCE@meta.data$cell.type) / nrow(querySCE@meta.data)
      
      # data type factor
      type = factor(ifelse(c(colnames(referenceSCE), colnames(querySCE)) %in% colnames(referenceSCE),
                           "Reference cells", "Query cells"), levels = c("Reference cells", "Query cells"))
      names(type) <- c(colnames(referenceSCE), colnames(querySCE))
      
      # input files as lists
      SCE_list = list(Reference = referenceSCE, Query = querySCE)
      assayNames = list(Reference = assayNameReference, Query = assayNameQuery)
      Idents(SCE_list[[1]]) <- SCE_list[[1]]$cell.type
      Idents(SCE_list[[2]]) <- SCE_list[[2]]$seurat_clusters
      sim.results <- ClusterFoldSimilarity::clusterFoldSimilarity(scList = SCE_list, sampleNames = names(SCE_list), parallel = T, nSubsampling = 25)
      
      sim.community <- ClusterFoldSimilarity::findCommunitiesSimmilarity(sim.results)
      
      tmp1 <- sim.results[sim.results$datasetL == "Query",c("clusterL","clusterR")]
      tmp2 <- querySCE@meta.data[,c("seurat_clusters", "cell.type")]
      tmp2$assign = ""
      
      for(i in 1:nrow(tmp1)){
        tmp2$assign[tmp2$seurat_clusters == as.numeric(tmp1[i,1])] = tmp1[i,2]
      }
      
      acc_clusterFoldSim = sum(tmp2$cell.type == tmp2$assign)/nrow(tmp2);acc_clusterFoldSim


      ################## Corr Pearson
      cor1 <- Seurat::AverageExpression(object = referenceSCE, features = genes, group.by = "cell.type")$originalexp
      cor2 <- Seurat::AverageExpression(object = querySCE, features = genes, group.by = "seurat_clusters")$originalexp
      colnames(cor2) <- 0:(ncol(cor2)-1)
      tmp2 <- querySCE@meta.data[,c("seurat_clusters", "cell.type")]
      tmp2$assign = ""
      for(i in 1:ncol(cor2)){
        cell_target <- names(sort(apply(cor1, 2, function(ccr)cor(cor2[,i], ccr, method = "pearson")), decreasing = T)[1])
        clust_target <- colnames(cor2)[i]
        tmp2$assign[tmp2$seurat_clusters == clust_target] = cell_target
      }      
      acc_pearson_cor = sum(tmp2$cell.type == tmp2$assign)/nrow(tmp2);acc_pearson_cor


      ################## Cosine Sim
      tmp2$assign = ""
      for(i in 1:ncol(cor2)){
        cell_target <- names(sort(apply(cor1, 2, function(ccr)lsa::cosine(x = cor2[,i],y =  ccr)), decreasing = T)[1])
        clust_target <- colnames(cor2)[i]
        tmp2$assign[tmp2$seurat_clusters == clust_target] = cell_target
      }      
      acc_cosine_cor = sum(tmp2$cell.type == tmp2$assign)/nrow(tmp2);acc_cosine_cor

      ################## Cosine Sim
      tmp2$assign = ""
      for(i in 1:ncol(cor2)){
        cell_target <- names(sort(apply(cor1, 2, function(ccr)dist(x = t(cbind(cor2[,i],ccr)), method = "minkowski", p = 1.5)), decreasing = F)[1])
        clust_target <- colnames(cor2)[i]
        tmp2$assign[tmp2$seurat_clusters == clust_target] = cell_target
      }      
      acc_minkowski_dist = sum(tmp2$cell.type == tmp2$assign)/nrow(tmp2);acc_cosine_cor
      
      
      
      embeddings_accuracy <- c(ClusterFoldSimil=acc_clusterFoldSim, pearson=acc_pearson_cor, cosine=acc_cosine_cor, minkowski=acc_minkowski_dist)
      embeddings_names_corrected <- c("ClusterFoldSimilarity" , "pearson", "cosine", "minkowski")
      res = rbind(res,
                  data.frame(
                    genes = rep(length(genes), length(embeddings_names_corrected)),
                    type = embeddings_names_corrected,
                    Accuracy = unlist(embeddings_accuracy),
                    querySample = rep(querySample, length(embeddings_names_corrected)),
                    referenceSample = rep(referenceSample, length(embeddings_names_corrected)),
                    Sim = sim
                  ))
      
      print(res)
      
      # save the results as we go:
      saveRDS(res, file = resFile)
    }
  }
}
