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

# sc.gastro <- MouseGastrulationData::EmbryoAtlasData()
####### SAVE POINT ########
# save(sc.gastro, file = paste0(PATH,"gastrulation.mouse.data.RData"))
####### SAVE POINT ########

####### LOAD POINT ########
fname <- load(file = paste0(PATH,"gastrulation.mouse.data.RData")); fname
####### LOAD POINT ########

pheno.e8 <- colData(sc.gastro) %>% as.data.frame() %>% dplyr::filter(stage == "E8.5")
samples = sc.gastro$sample[sc.gastro$stage == "E8.5"]

atlas = MouseGastrulationData::EmbryoAtlasData(type = "processed", samples = unique(samples))
atlas <- logNormCounts(atlas)

# Remove small cell populations as in the original benchmark from Stabmap
atlas <- atlas[, !atlas$celltype %in% c("Parietal endoderm", "PGC",  "Visceral endoderm", "Rostral neurectoderm", "Notochord")]

gc(reset = T,full = T)
counts(atlas) <- as(counts(atlas, withDimnames = FALSE), "CsparseMatrix")
logcounts(atlas) <- as(logcounts(atlas, withDimnames = FALSE), "CsparseMatrix")

atlas_seu <- Seurat::as.Seurat(atlas)
atlas_seu <- atlas_seu[,!is.na(atlas_seu$celltype)]
atlas_seu <- FindVariableFeatures(atlas_seu, selection.method = "vst", nfeatures = 4000)
atlas_seu <- ScaleData(atlas_seu, features = VariableFeatures(atlas_seu))
atlas_seu <- RunPCA(atlas_seu, features = VariableFeatures(object = atlas_seu))
atlas_seu <- RunUMAP(atlas_seu, features = VariableFeatures(object = atlas_seu))
atlas_seu <- FindNeighbors(atlas_seu, dims = seq(25))
atlas_seu <- FindClusters(atlas_seu, resolution = 4)

atlas_seu_8_5 <- atlas_seu
# ####### SAVE POINT ########
#save(atlas_seu_8_5, file = paste0(PATH,"gastrulation.mouse.8.5.data.Seurat.atlas.1500.features.RData"))
# ####### SAVE POINT ########

atlas_seu_8_5=NULL
atlas=NULL
gc(reset = T, full = T)

####### LOAD POINT ########
# fname <- load(file = paste0(PATH,"gastrulation.mouse.8.5.data.Seurat.atlas.1500.features.RData"));fname
# atlas_seu <- atlas_seu_8_0
####### LOAD POINT ########

samples = unique(atlas_seu$sample)
k = 8
acc_df_all = NULL
nGenes_all = c(25, 50, 100, 250, 500, 1000, 1500)
nGenes_all = rep(nGenes_all, each=2)
names(nGenes_all) <- paste0("Sim_", seq_len(length(nGenes_all)))

labels = "celltype"
assayNameReference = "logcounts"
assayNameQuery = "logcounts"
doUMAP = FALSE

### ALGORITHMS TO USE:
# devtools::install_github("MarioniLab/StabMap")
###

resFile = paste0(PATH,"mouse.gastrulation.benchmark.results.res1.5.pc15.2.5kfeats.seed.29072024.Rds")

if (!file.exists(resFile)) {
  res = NULL
} else {
  res = readRDS(resFile)
}

BiocParallel::register(BPPARAM =  BiocParallel::MulticoreParam(workers = 16))
library(Seurat)
embeddings_names=""

######### SETTING SEED FOR REPRODUCIBILITY
set.seed(29072024)
######### SETTING SEED FOR REPRODUCIBILITY

for (sim in names(nGenes_all)) {
  
  nGenes = nGenes_all[sim]
  
  if (nGenes <= 100) {
    nPCs = 20
    k.param = 30
  } else {
    nPCs = 15
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
      hvgs <- Seurat::VariableFeatures(atlas_seu)[1:2000]
      genes = sample(hvgs, size = nGenes)
      referenceSCE = atlas_seu[,atlas_seu$sample %in% referenceSample]
      referenceSCE <- Seurat::NormalizeData(referenceSCE)
      referenceSCE <- Seurat::ScaleData(referenceSCE)
      referenceSCE <- referenceSCE[hvgs,]
      # select random genes from among the HVGs for the query data
      querySCE = atlas_seu[, atlas_seu$sample %in% querySample]
      querySCE <- Seurat::NormalizeData(querySCE)
      querySCE <- Seurat::ScaleData(querySCE)
      querySCE <- RunPCA(querySCE, features = VariableFeatures(object = querySCE), npcs = nPCs)
      querySCE <- RunUMAP(querySCE, features = VariableFeatures(object = querySCE))
      querySCE <- FindNeighbors(querySCE, dims = seq(nPCs), k.param = k.param)
      querySCE <- FindClusters(querySCE, resolution = 1.5, method = 4)
      querySCE <- querySCE[genes,]
      # Compute base accuracy by cluster 
      querySCE@meta.data$assign = ""
      for (cls in unique(querySCE@meta.data$seurat_clusters)){
        ccount <- table(querySCE@meta.data$celltype[querySCE@meta.data$seurat_clusters == cls])
        querySCE@meta.data$assign[querySCE@meta.data$seurat_clusters == cls] <- names(which.max(ccount))
      }
      sum(querySCE@meta.data$assign == querySCE@meta.data$celltype) / nrow(querySCE@meta.data)
      
      # data type factor
      type = factor(ifelse(c(colnames(referenceSCE), colnames(querySCE)) %in% colnames(referenceSCE),
                           "Reference cells", "Query cells"), levels = c("Reference cells", "Query cells"))
      names(type) <- c(colnames(referenceSCE), colnames(querySCE))
      
      # input files as lists
      SCE_list = list(Reference = referenceSCE, Query = querySCE)
      assayNames = list(Reference = assayNameReference, Query = assayNameQuery)
      Idents(SCE_list[[1]]) <- SCE_list[[1]]$celltype
      Idents(SCE_list[[2]]) <- SCE_list[[2]]$seurat_clusters
      sim.results <- ClusterFoldSimilarity::clusterFoldSimilarity(scList = SCE_list, sampleNames = names(SCE_list), parallel = T, nSubsampling = 20)
      
      tmp1 <- sim.results[sim.results$datasetL == "Query",c("clusterL","clusterR")]
      tmp2 <- querySCE@meta.data[,c("seurat_clusters", "celltype")]
      tmp2$assign = ""
      
      for(i in 1:nrow(tmp1)){
        tmp2$assign[tmp2$seurat_clusters == as.numeric(tmp1[i,1])] = tmp1[i,2]
      }
      
      acc_clusterFoldSim = sum(tmp2$celltype == tmp2$assign)/nrow(tmp2);acc_clusterFoldSim

      # Seurat Integration
      Seurat_embedding = Seurat_wrap(SCE_list, method="cca")
      dim(Seurat_embedding)
      embeddings_names <- c(Seurat="Seurat_embedding")

      SCE_list[[1]] <- as.SingleCellExperiment(SCE_list[[1]])
      SCE_list[[2]] <- as.SingleCellExperiment(SCE_list[[2]])
      assay_list = mapply(assay, SCE_list, assayNames)
      counts_list = mapply(assay, SCE_list, "counts")
      labels_list = list(Reference = setNames(referenceSCE$celltype, colnames(referenceSCE)))
      
      # perform PCA
      PC_embedding = mapPCA(lapply(assay_list, as.matrix), nPCs = nPCs)
      dim(PC_embedding)
      colnames(PC_embedding) <- paste0("PC_", seq_len(ncol(PC_embedding)))
      embeddings_names <- c(embeddings_names, PCA="PC_embedding")
      
      # 
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
      batchFactor = as.character(interaction(batchFactor_ref, atlas_seu@meta.data[cnames, "sample"]))
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
      
      # # calculate the resubstitution accuracy
      embeddings_resub_accuracy = sapply(embeddings_names_corrected, function(nm) {
        if (is.null(get(nm))) return(NA)
        print(nm)

        # only use cells with labels
        data_all = get(paste0(nm))
        labels_train = referenceLabels[!is.na(referenceLabels) & names(referenceLabels) %in% rownames(data_all)]

        knn_out = embeddingKNN(data_all,
                               labels_train,
                               type = "uniform_fixed",
                               k_values = 5)

        acc = mean(isEqual(knn_out[names(labels_train),"resubstituted_labels"],
                           labels_train), na.rm = TRUE)
        return(acc)
      }, simplify = FALSE)
      
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
      
      # save the results on the go:
      saveRDS(res, file = resFile)
    }
  }
}




####### Plot
res.t <- res
res.t$type <- factor(res.t$type)
res.t$genes <- factor(res.t$genes)
res.t$referenceSample <- factor(res.t$referenceSample)

ggplot(res.t, aes(x=genes, y=Accuracy, fill = type, color=type)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("n. of genes") +
  #ylim(0.25,1) +
  ggplot2::scale_y_continuous(breaks =  seq(0.30, 1, by=0.1), labels = as.character(seq(0.30, 1, by=0.1))) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(length(levels(res.t$type)), "Pastel1")) +# Pastel2
  scale_color_manual(values = RColorBrewer::brewer.pal(length(levels(res.t$type)), "Pastel1")) # Pastel2

# check specific sims
#res.t %>% dplyr::filter(genes==1500 & type == "ClusterFoldSimilarity") %>% summarise(mean(Accuracy))
#res.t %>% dplyr::filter(genes==1500 & type == "ClusterFoldSimilarity") %>% pull(Accuracy) %>% summary()
# Barplot with means
means.all <- res.t %>% group_by(type, genes) %>% summarise(accuracy=mean(Accuracy), sd=sd(Accuracy)/sqrt(length(Accuracy)))
pdf(file = "/Users/oscargonzalezvelasco/Desktop/ClusterFoldSimilarity/review/barplot.means.gastrulation.acc.pdf", width = 10,height = 2.2)
g <- ggplot(means.all, aes(x=genes, y=accuracy, fill = type)) +
  geom_bar(stat="identity", position = "dodge") +
  theme_minimal() +
  xlab("n. of genes") +
  #ylim(0.25,1) +
  ggplot2::geom_errorbar( aes(x=genes, ymin=accuracy-sd, ymax=accuracy+sd), 
                          width=0.4, colour="grey", alpha=0.9, position=position_dodge(.9), linewidth=0.5) +
  ggplot2::scale_y_continuous(breaks =  seq(0.30, 1, by=0.1), labels = as.character(seq(0.30, 1, by=0.1))) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(length(levels(res.t$type)), "Pastel1")) # Pastel2
print(g)
dev.off()
