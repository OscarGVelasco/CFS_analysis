#!/usr/bin/env Rscript
### BENCHMARKING PARTY **
library(Seurat)
library(SeuratObject)
library(ClusterFoldSimilarity)
library(Matrix)

PATH <- "SET_PATH"

# data dimensions
nfeatures <- 5000
ncells_all <- c(5000, 25000, 50000, 100000, 150000, 200000, 250000)

#nfeatures <- 2000; ncells <- 4000
ncpus <- 17
for (ncells in ncells_all){
  # single-cell 1
  counts <- Matrix::Matrix(rpois(n=nfeatures * ncells, lambda=1), nfeatures, sparse = TRUE)
  gc()
  Matrix::nnzero(counts)
  rownames(counts) <- paste0("gene",seq(nfeatures))
  colnames(counts) <- paste0("cell",seq(ncells))
  # Cluster Definition
  nClusters <- 20
  colData <- data.frame(cluster=sample(paste("Cluster",seq(nClusters), sep = ""),size = ncells,replace = TRUE),
                        row.names=paste0("cell",seq(ncells)))
  seu1 <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = colData)
  Idents(object = seu1) <- "cluster"
  gc()
  
  # single-cell 2
  counts <- Matrix::Matrix(rpois(n=nfeatures * ncells, lambda=1), nfeatures)
  gc()
  rownames(counts) <- paste0("gene",seq(nfeatures))
  colnames(counts) <- paste0("cell",seq(ncells))
  # Cluster Definition
  nClusters <- 20
  colData <- data.frame(cluster=sample(paste("Cluster",seq(nClusters), sep = ""),size = ncells,replace = TRUE),
                        row.names=paste0("cell",seq(ncells)))
  seu2 <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = colData)
  Idents(object = seu2) <- "cluster"
  rm(counts)
  # Create a list with the unprocessed single-cell datasets
  singlecellObjectList <- list(seu1, seu2)
  rm(seu1)
  rm(seu2)
  gc()
  #i=26
  singlecellObjectList[[3]] <- singlecellObjectList[[1]]
  # # Test with 1 CPU
  # time_tmp <- system.time(
  #   similarityTable <- clusterFoldSimilarity(scList=singlecellObjectList, 
  #                                            sampleNames = c("sc1","sc2","sc3"), 
  #                                            nSubsampling = 30, 
  #                                            parallel = FALSE)
  # )
  # save(time_tmp, file = file.path(PATH, paste0("time.3.datasets.CPU.1.RData")))
  # Setting the number of CPUs with BiocParallel:
  message(paste("Trying with size: ", ncells))
  time <- data.frame()
  BiocParallel::register(BPPARAM =  BiocParallel::MulticoreParam(workers = ncpus))
  for (j in 1:2){
    message(paste("- Round ", j))
    gc1 <- gc(reset = TRUE)
    time_tmp <- system.time(
      similarityTable <- clusterFoldSimilarity(scList=singlecellObjectList, 
                                               sampleNames = c("sc1","sc2","sc3"), 
                                               nSubsampling = 30, 
                                               parallel = TRUE)
    )
    gc2 <- gc()
    time <- rbind(time, c(time_tmp, size=ncells, cpus=ncpus))
    colnames(time) <- c(names(time_tmp), "size", "cpus")
    #cat(sprintf("mem: %.1fMb.\n", sum(gc2[,6] - gc1[,2])))
    #print(paste("Time Measured:", time))
    #inception.list[[i]] <- time
    save(time, file = file.path(PATH, paste0("time.3.datasets.SIZE.",ncells,".",j,".using.17.CPUs.RData")))
  }
}
# save(inception.list, file = paste0(PATH, "time.benchmarking.results.RData"))

#head(similarityTable)
# load(file = paste0(PATH, "time.benchmarking.results.RData")
# inception.list <- cbind.data.frame(do.call("rbind", inception.list),cpus=seq(1,56,by=2))
# library(ggplot2)
# ggplot(inception.list, aes(x=cpus,y=elapsed/60)) +
#   geom_point(col="blue", size=3) +
#   theme_minimal() +
#   geom_hline(yintercept=(min(inception.list$elapsed)/60)-0.2, linetype="dashed",
#              color = "red", size=0.5) +
#   ggtitle("Computing time for 400k cells 10 cluster per dataset")

