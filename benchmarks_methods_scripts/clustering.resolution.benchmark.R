PATH <- "/cluster_similarity_sc/"
setwd(PATH)
library(dplyr)
source("./initialice.R")

fname <- load(paste0(PATH,"rdata/Pancreas_sc_objects_list_processed.RData"));fname

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

nGenes=1500
querySample="pancreas_small"
referenceSample="pancreas_large"
results_reso=NULL

vfeat <- lapply(pancreas_list, function(x){Seurat::VariableFeatures(Seurat::FindVariableFeatures(x, nfeatures=5000))})
hvgs <- intersect(vfeat[[1]], vfeat[[2]])[1:1500]

genes = sample(hvgs)[seq_len(nGenes)]
referenceSCE = pancreas_list[[referenceSample]]
referenceSCE <- Seurat::NormalizeData(referenceSCE)
referenceSCE <- Seurat::ScaleData(referenceSCE)
referenceSCE <- Seurat::FindVariableFeatures(referenceSCE)
referenceSCE <- RunPCA(referenceSCE, features = VariableFeatures(object = referenceSCE), npcs = 30)
referenceSCE <- referenceSCE[hvgs,]
# select random genes from among the HVGs for the query data
querySCE = pancreas_list[[querySample]]
querySCE <- Seurat::NormalizeData(querySCE)
querySCE <- Seurat::FindVariableFeatures(querySCE)
querySCE <- Seurat::ScaleData(querySCE)
querySCE <- RunPCA(querySCE, features = VariableFeatures(object = querySCE), npcs = 30)

gc(reset=T, full=T)

labels = "cell.type"
assayNameReference = "logcounts"
assayNameQuery = "logcounts"

for(resol in c(0.6, 1, 1.4, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.4)){
  for(pc_test in c(10, 15, 20, 25, 28, 30)){
    querySCE <- FindNeighbors(querySCE, dims = seq(pc_test), k.param = 30)
    querySCE <- FindClusters(querySCE, resolution = resol, method = 4)
    querySCE <- querySCE[genes,]
    nclusts <- length(unique(querySCE$seurat_clusters))
    # Compute base accuracy by cluster 
    querySCE@meta.data$assign = ""
    for (cls in unique(querySCE@meta.data$seurat_clusters)){
      ccount <- table(querySCE@meta.data$cell.type[querySCE@meta.data$seurat_clusters == cls])
      querySCE@meta.data$assign[querySCE@meta.data$seurat_clusters == cls] <- names(which.max(ccount))
    }
    base_accu = sum(querySCE@meta.data$assign == querySCE@meta.data$cell.type) / nrow(querySCE@meta.data)
    
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
    results_reso <- rbind(results_reso, data.frame(resolution=resol, acc=acc_clusterFoldSim, 
                                                   base_acc=base_accu, nclusters=nclusts, pca=pc_test, ngenes=nGenes))
  }
  }  

saveRDS(results_reso, 
        file = paste0(PATH, "clustering.resolution.and.pca.benchmark.pancreas.multi.1500.genes.Rds"))
