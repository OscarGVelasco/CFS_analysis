PATH <- "/home/o872o/o872o/cluster_similarity_sc/review/"
setwd(PATH)
library(dplyr)
# devtools::install_github("MarioniLab/StabMap")
source("/home/o872o/o872o/cluster_similarity_sc/review/initialice.R")

# Load pancreas data sets
fname <- load("/home/o872o/o872o/cluster_similarity_sc/rdata/Pancreas_sc_objects_list_processed.RData");fname

# Find and subset common cell types
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
acc_df_all = NULL
nGenes_all = c(25, 50, 100, 250, 500, 1000, 1500)
nGenes_all <- rep(nGenes_all, each=6)
names(nGenes_all) <- paste0("Sim_", seq_len(length(nGenes_all)))

labels = "cell.type"
assayNameReference = "logcounts"
assayNameQuery = "logcounts"
doUMAP = FALSE

resFile = paste0(PATH,"mouse.pancreas.benchmark.results.res.2.pca.30.seed.04022025.Rds")

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
gc(reset = T, full = T)

for (sim in names(nGenes_all)) {
  
  nGenes = nGenes_all[sim]
  # As defined in the original benchmark paper
  if (nGenes <= 100) {
    nPCs = 20
    k.param = 30
  } else {
    nPCs = 40
    k.param = 20
  }
  
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
      # hvgs - defined at the beginning
      genes = sample(hvgs)[seq_len(nGenes)]
      referenceSCE = pancreas_list[[referenceSample]]
      referenceSCE <- Seurat::NormalizeData(referenceSCE)
      referenceSCE <- Seurat::ScaleData(referenceSCE)
      referenceSCE <- Seurat::FindVariableFeatures(referenceSCE)
      referenceSCE <- RunPCA(referenceSCE, features = VariableFeatures(object = referenceSCE), npcs = 15)
      referenceSCE <- referenceSCE[hvgs,]
      # select random genes from among the HVGs for the query data
      querySCE = pancreas_list[[querySample]]
      querySCE <- Seurat::NormalizeData(querySCE)
      querySCE <- Seurat::ScaleData(querySCE)
      querySCE <- Seurat::FindVariableFeatures(querySCE)
      querySCE <- RunPCA(querySCE, features = VariableFeatures(object = querySCE), npcs = 30)        
      querySCE <- FindNeighbors(querySCE, dims = seq(15), k.param = 30) 
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
      
      #sim.community <- ClusterFoldSimilarity::findCommunitiesSimmilarity(sim.results)
      
      tmp1 <- sim.results[sim.results$datasetL == "Query",c("clusterL","clusterR")]
      tmp2 <- querySCE@meta.data[,c("seurat_clusters", "cell.type")]
      tmp2$assign = ""
      
      for(i in 1:nrow(tmp1)){
        tmp2$assign[tmp2$seurat_clusters == as.numeric(tmp1[i,1])] = tmp1[i,2]
      }
      
      acc_clusterFoldSim = sum(tmp2$cell.type == tmp2$assign)/nrow(tmp2);acc_clusterFoldSim

      # Seurat Integration
      message("COMPUTING: Seurat")
      Seurat_embedding = Seurat_wrap(SCE_list, method="cca")
      dim(Seurat_embedding)
      embeddings_names <- c(Seurat="Seurat_embedding")
      
      SCE_list[[1]] <- as.SingleCellExperiment(SCE_list[[1]])
      SCE_list[[2]] <- as.SingleCellExperiment(SCE_list[[2]])
      assay_list = mapply(assay, SCE_list, assayNames)
      counts_list = mapply(assay, SCE_list, "counts")
      labels_list = list(Reference = setNames(referenceSCE$cell.type, colnames(referenceSCE)))
      
      # perform PCA
      message("COMPUTING: PCA emb")
      PC_embedding = mapPCA(lapply(assay_list, as.matrix), nPCs = nPCs)
      dim(PC_embedding)
      colnames(PC_embedding) <- paste0("PC_", seq_len(ncol(PC_embedding)))
      embeddings_names <- c(embeddings_names, PCA="PC_embedding")
      
      library(StabMap)
      StabMap_embedding = stabMap(lapply(assay_list, as.matrix),
                                  reference_list = "Reference",
                                  ncomponentsReference = nPCs,
                                  projectAll = TRUE,
                                  scale.center = FALSE,
                                  scale.scale = FALSE,
                                  plot = FALSE
      )
      dim(StabMap_embedding)
      embeddings_names <- c(embeddings_names, StabMap="StabMap_embedding")

      MultiMAP_embedding = MultiMAP_wrap(assay_list, verbose = TRUE)
      dim(MultiMAP_embedding)
      embeddings_names <- c(embeddings_names, MultiMap="MultiMAP_embedding")
      
      nobatch = c("MultiMAP_embedding","Seurat_embedding")#, "UINMF_embedding")
      noumap = c("MultiMAP_embedding")
      
      library(Seurat)
      referenceSCE <- as.SingleCellExperiment(referenceSCE)
      querySCE <- as.SingleCellExperiment(querySCE)
      
      cnames = unlist(lapply(assay_list, colnames))
      batchFactor_ref <- ifelse(cnames %in% colnames(referenceSCE),"Reference", "Query")
      batchFactor <- batchFactor_ref
      names(batchFactor) <- cnames

      sapply(embeddings_names, function(nm) {
        print(nm)
        # can include optional flags for certain embeddings to 
        # not be passed through a batch correction stage
        # e.g. from MultiMAP
        if (nm %in% nobatch) {
          # i.e. do nothing
          print(paste0("no batch correction for ", nm))
          assign(paste0(nm, "_corrected"), get(nm), envir = .GlobalEnv)
        } else {
          assign(paste0(nm, "_corrected_none"), get(nm),
                 envir = .GlobalEnv)
          assign(paste0(nm, "_corrected_MNN"), reducedMNN_batchFactor(get(nm), batchFactor),
                 envir = .GlobalEnv)
          assign(paste0(nm, "_corrected_Harmony"), Harmony_batchFactor(get(nm), batchFactor),
                 envir = .GlobalEnv)
        }
      }, simplify = FALSE)
      
      PC_embedding_corrected_Harmony <- t(PC_embedding_corrected_Harmony)
      embeddings_names_corrected = c(nobatch, apply(expand.grid(setdiff(embeddings_names, nobatch),
                                                                paste0("_corrected_", c("none", "MNN", "Harmony", "Seurat"))), 1, paste0, collapse = ""))
      
      # predict cell type labels using knn with k = 5
      referenceLabels = colData(referenceSCE)[,labels]
      names(referenceLabels) = colnames(referenceSCE)
      
      queryLabels = colData(querySCE)[,labels]
      names(queryLabels) = colnames(querySCE)
      
      embeddings_names_corrected <- embeddings_names_corrected[!embeddings_names_corrected %in% 
                                                                 c("StabMap_embedding_corrected_Harmony",
                                                                   "PC_embedding_corrected_Seurat","StabMap_embedding_corrected_Seurat")]
      
      # calculate accuracy of the query cells
      embeddings_accuracy = sapply(embeddings_names_corrected, function(nm) {
        if (is.null(get(nm))) return(NA)
        
        print(nm)
        
        # only use cells with labels
        data_all = get(paste0(nm))
        
        labels_train = referenceLabels[!is.na(referenceLabels) & names(referenceLabels) %in% rownames(data_all)]
        
        knn_out = embeddingKNN(data_all,
                               labels_train,
                               type = "uniform_fixed",
                               k_values = 5)
        
        acc = mean(isEqual(knn_out[names(queryLabels),"predicted_labels"], queryLabels), na.rm = TRUE)
        
        print(acc)
        
        return(acc)
      }, simplify = FALSE)
      
      print(embeddings_accuracy)
      embeddings_accuracy <- c(ClusterFoldSimil=acc_clusterFoldSim, embeddings_accuracy)
      embeddings_names_corrected <- c("ClusterFoldSimilarity" , embeddings_names_corrected)
      res = rbind(res,
                  data.frame(
                    genes = rep(length(genes), length(embeddings_names_corrected)),
                    type = embeddings_names_corrected,
                    Accuracy = unlist(embeddings_accuracy),
                    #Accuracy_resub = unlist(embeddings_resub_accuracy),
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

