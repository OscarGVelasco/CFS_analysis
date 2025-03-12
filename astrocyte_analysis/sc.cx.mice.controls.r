# Load libraries ---------------------------------------------------------------------
library(SingleCellExperiment)
library(scater)
library(gridExtra)
library(Seurat)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(SingleR)
# Set folder to load and save data/plots
MASTERPATH <- "/SET_FOLDER/"
setwd(MASTERPATH)
PATH="cellRanger_data/"

############################################################
############# Loading and merging data:

############ Batch 1
# Sample 26
sc_26 <- Read10X(data.dir = paste0(MASTERPATH,PATH,"26_sc/"))
pheno26 <- data.frame(stringsAsFactors = F,row.names = colnames(sc_26),mouse = rep("26_sc",ncol(sc_26)),status= rep("Control",ncol(sc_26)))
sc_26 <- Seurat::CreateSeuratObject(counts = sc_26,meta.data = pheno26)
# Sample 36
sc_36 <- Read10X(data.dir = paste0(MASTERPATH,PATH,"36_sc/"))
pheno36 <- data.frame(stringsAsFactors = F,row.names = colnames(sc_36),mouse = rep("36_sc",ncol(sc_36)),status= rep("Control",ncol(sc_36)))
sc_36 <- Seurat::CreateSeuratObject(counts = sc_36,meta.data = pheno36)
# Sample 699
sc_699 <- Read10X(data.dir = paste0(MASTERPATH,PATH,"699_sc/"))
pheno699 <- data.frame(stringsAsFactors = F,row.names = colnames(sc_699),mouse = rep("699_sc",ncol(sc_699)),status= rep("Control",ncol(sc_699)))
sc_699 <- Seurat::CreateSeuratObject(counts = sc_699,meta.data = pheno699)
# Cortex
# Sample 26
cx_26 <- Read10X(data.dir = paste0(MASTERPATH,PATH,"26_cx/"))
pheno26 <- data.frame(stringsAsFactors = F,row.names = colnames(cx_26),mouse = rep("26_ctx",ncol(cx_26)),status= rep("Control",ncol(cx_26)))
cx_26 <- Seurat::CreateSeuratObject(counts = cx_26,meta.data = pheno26)
# Sample 36
cx_36 <- Read10X(data.dir = paste0(MASTERPATH,PATH,"36_cx/"))
pheno36 <- data.frame(stringsAsFactors = F,row.names = colnames(cx_36),mouse = rep("36_ctx",ncol(cx_36)),status= rep("Control",ncol(cx_36)))
cx_36 <- Seurat::CreateSeuratObject(counts = cx_36,meta.data = pheno36)
# Sample 699
cx_699 <- Read10X(data.dir = paste0(MASTERPATH,PATH,"699_cx/"))
pheno699 <- data.frame(stringsAsFactors = F,row.names = colnames(cx_699),mouse = rep("699_ctx",ncol(cx_699)),status= rep("Control",ncol(cx_699)))
cx_699 <- Seurat::CreateSeuratObject(counts = cx_699,meta.data = pheno699)

# # Visualization of ALL cells merged + including Astrocyte groups detected in analysis with ClusterFoldSimilarity (tmp_meta)
all_merged <- merge(sc_36,list(sc_26,sc_699,cx_36,cx_26,cx_699),)
all_merged = JoinLayers(all_merged)
all_merged[["percent.mt"]] <- PercentageFeatureSet(all_merged, pattern = "^mt")
all_merged <- subset(all_merged, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 1.5)
all_merged <- NormalizeData(all_merged)
all_merged <- FindVariableFeatures(all_merged, selection.method = "vst", nfeatures = 2000)
all_merged <- ScaleData(all_merged,features = VariableFeatures(all_merged))
all_merged <- RunPCA(all_merged, features = VariableFeatures(object = all_merged))
all_merged <- RunUMAP(all_merged, features = VariableFeatures(object = all_merged))
all_merged@meta.data$tissue <- unlist(lapply(strsplit(all_merged@meta.data$mouse, split = "_"),function(cell)cell[2]))
all_merged@meta.data$cell_type <- "other"
all_merged@meta.data[rownames(tmp_meta),]$cell_type <- "Astrocytes"
all_merged@meta.data$mouse_id <- unlist(lapply(strsplit(all_merged@meta.data$mouse, "_"), function(mos)mos[1]))
DimPlot(all_merged, group.by = "mouse_id") + ggtitle("Mouse ID all cells")
DimPlot(all_merged, group.by = "mouse") + ggtitle("Mouse ID all cells")
DimPlot(all_merged, group.by = "tissue", cols=scales::alpha(c("#D5BADB","#7EB6D9"))) + ggtitle("Tissue type all cells")
DimPlot(all_merged, group.by = "cell_type", cols=scales::alpha(c("#F08080","#D9E8F5", 0.6))) + ggtitle("Cell type all cells")
all_merged@meta.data$astrocyte_group <- "NA"
all_merged@meta.data[rownames(tmp_meta),]$astrocyte_group <- tmp_meta$group
DimPlot(all_merged, group.by = "astrocyte_group", cols=scales::alpha(c("#F08080","#F2D377","#7EB6D9","#D9E8F5", 0.6))) + ggtitle("Astrocyte groups")
# # Visualization END


all(rownames(all_merged@meta.data) %in% rownames(tmp_meta))
#all_merged <- single.process(x = all_merged, resolution = 0.8, dims = 16)
g_merged <- DimPlot(all_merged, group.by = "mouse")
ggsave(plot = ggsave, filename="yourfile.pdf",device = cairo_pdf, width = 297, height = 210, units = "mm") # A4 landscape

tmp_meta

sc_list <- list()
sc_list2 <- list()
sc_list[[1]] <- sc_36
sc_list[[2]] <- sc_26
sc_list[[3]] <- sc_699
sc_list[[4]] <- cx_36
sc_list[[5]] <- cx_26
sc_list[[6]] <- cx_699

sum(unlist(lapply(sc_list,ncol)))
mean(unlist(lapply(sc_list,ncol)))
mean(unlist(lapply(sc_list[1:3],ncol))) # Spinal Cord
mean(unlist(lapply(sc_list[4:6],ncol))) # Motor Cortex
#
sd(unlist(lapply(sc_list,ncol)))
sd(unlist(lapply(sc_list[1:3],ncol))) # Spinal Cord
sd(unlist(lapply(sc_list[4:6],ncol))) # Motor Cortex

ref <- celldex::MouseRNAseqData(ensembl = F)

single.process <- function(x,resolution,dims){
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^mt")
  x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 1.5)
  x <- NormalizeData(x)
  pred <- SingleR(test=as.SingleCellExperiment(x), ref=ref, labels=ref$label.main,check.missing = T)
  pred[is.na(pred$pruned.labels),]$pruned.labels <- "unknown"
  x@meta.data <- cbind.data.frame(x@meta.data,cell.type = pred[rownames(x@meta.data),]$pruned.labels)
  tmp <- x@meta.data %>% dplyr::filter(cell.type %in% "Astrocytes") %>% rownames()
  x <- x[,tmp]
  #
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x,features = VariableFeatures(x))
  x <- RunPCA(x, features = VariableFeatures(object = x))
  x <- RunUMAP(x, features = VariableFeatures(object = x))
  x <- JackStraw(x, num.replicate = 100, dims = 25)
  x <- ScoreJackStraw(x, dims = 1:25)
  print(JackStrawPlot(x, dims = 1:25))
  print(ElbowPlot(x, ndims=25))
  x <- FindNeighbors(x, dims = seq(dims))
  x <- FindClusters(x, resolution = resolution)
}
DimPlot(x,group.by = "cell.type") +
DimPlot(x, cells.highlight = sc_list2[[5]]@meta.data %>% filter(group %in% "group 2") %>% rownames())


sc_list2[[1]] <- single.process(x = sc_list[[1]], resolution = 0.8, dims = 8)
Seurat::DimPlot(sc_list2[[1]])
sc_list2[[2]] <- single.process(x = sc_list[[2]], resolution = 0.8, dims = 8)
Seurat::DimPlot(sc_list2[[2]])
sc_list2[[3]] <- single.process(x = sc_list[[3]], resolution = 1.6, dims = 6)
Seurat::DimPlot(sc_list2[[3]])
sc_list2[[4]] <- single.process(x = sc_list[[4]], resolution = 0.6, dims = 6)
Seurat::DimPlot(sc_list2[[4]])
sc_list2[[5]] <- single.process(x = sc_list[[5]], resolution = 0.4, dims = 5)
Seurat::DimPlot(sc_list2[[5]])
sc_list2[[6]] <- single.process(x = sc_list[[6]], resolution = 0.5, dims = 5)
Seurat::DimPlot(sc_list2[[6]])

names(sc_list2) <- c("sc.36","sc.26","sc.699","cx.36","cx.26","cx.699")

######################### - Save Point - ##########################
# save(sc_list2, file = paste0(MASTERPATH, "mice.controls.data.cortex.spinalcord.processed.6.samples.RData"))
######################### - Save Point - ##########################
fileName <- load(file = paste0(MASTERPATH, "mice.controls.data.cortex.spinalcord.processed.6.samples.RData"));fileName
sprintf(fileName)

table(sc_list2[[1]]@meta.data$seurat_clusters) *1/3
table(sc_list2[[2]]@meta.data$seurat_clusters) *1/3
table(sc_list2[[3]]@meta.data$seurat_clusters) *1/3
table(sc_list2[[4]]@meta.data$seurat_clusters) *1/3
table(sc_list2[[5]]@meta.data$seurat_clusters) *1/3
table(sc_list2[[6]]@meta.data$seurat_clusters) *1/3

lapply(sc_list2,dim)

n.of.astros <- unlist(lapply(sc_list2,ncol))
mean(n.of.astros)
sd(n.of.astros)

# Clusters:
(Seurat::DimPlot(sc_list2[[1]]) + ggtitle(names(sc_list2)[1])) +
  (Seurat::DimPlot(sc_list2[[2]])  + ggtitle(names(sc_list2)[2])) +
  (Seurat::DimPlot(sc_list2[[3]]) + ggtitle(names(sc_list2)[3])) +
  (Seurat::DimPlot(sc_list2[[4]]) + ggtitle(names(sc_list2)[4])) +
  (Seurat::DimPlot(sc_list2[[5]]) + ggtitle(names(sc_list2)[5])) +
  (Seurat::DimPlot(sc_list2[[6]]) + ggtitle(names(sc_list2)[6]))

(Seurat::DimPlot(sc_list2[[1]]) + ggtitle(names(sc_list2)[1])) +
  (Seurat::DimPlot(sc_list2[[2]])  + ggtitle(names(sc_list2)[2])) +
  (Seurat::DimPlot(sc_list2[[4]]) + ggtitle(names(sc_list2)[4])) +
  (Seurat::DimPlot(sc_list2[[5]]) + ggtitle(names(sc_list2)[5]))


set.m <- c(1,2,3,4,5,6)
variable.features <- Reduce(union,lapply(sc_list2[set.m],function(x){Seurat::VariableFeatures(x)}));length(variable.features)
sc_list3 <- lapply(sc_list2[set.m], function(x)return(x[variable.features[1:6000],]))

# devtools::install_github(repo = "https://github.com/OscarGVelasco/ClusterFoldSimilarity")
# Set number of CPUs to use:
BiocParallel::register(BPPARAM =  BiocParallel::MulticoreParam(workers = 16))
simil <- ClusterFoldSimilarity::clusterFoldSimilarity(scList = sc_list2, sampleNames = names(sc_list2),topN = 1, nSubsampling = 18, parallel = T)
comms <- findCommunitiesSimmilarity(similarityTable = simil, leidenResolution = 3)
comms <- ClusterFoldSimilarity::findCommunitiesSimmilarity(similarityTable = simil, leidenResolution = 3)

simil3 <- ClusterFoldSimilarity::clusterFoldSimilarity(scList = sc_list3, sampleNames = names(sc_list2),topN = 1, nSubsampling = 20, parallel = T)
comms3 <- ClusterFoldSimilarity::findCommunitiesSimmilarity(similarityTable = simil3, leidenResolution = 3)


simil.all <- ClusterFoldSimilarity::clusterFoldSimilarity(sceList = sc_list2[c(2,5)], sampleNames = names(sc_list2)[c(2,5)],topN = Inf, nSubsampling = 22)

simil <- ClusterFoldSimilarity::clusterFoldSimilarity(sceList = sc_list3, sampleNames = names(sc_list3),topN = 1, nSubsampling = 18)

simil.all <- ClusterFoldSimilarity::clusterFoldSimilarity(sceList = sc_list3,sampleNames = names(sc_list3),topN = Inf)


ClusterFoldSimilarity::plotClustersGraph(simil)
simil.top10 <- ClusterFoldSimilarity::clusterFoldSimilarity(sceList = sc_list3,sampleNames = names(sc_list3),topN = 1,topNFeatures = 20, nSubsampling = 18)
simil.top10 <- ClusterFoldSimilarity::clusterFoldSimilarity(scList = sc_list2,sampleNames = names(sc_list2),topN = 1,topNFeatures = 20, nSubsampling = 16)

####### Testing canonical astrocyte markers
astro.markers <- c("Atp1b2", "Apoe", "Slc1a3","Slc1a2", "Gpc5", "Gria1", "S100b", "Aldh1l1", "Aqp4", "Gfap", "Slc7a10","Glud1","Aldoc","Sox9","Glul")#,"Glt1")
astro.markers <- c("Frzb", "Slc1a3","Dab1")
astro.markers <- c("Gria1","Gria4")
# Canonical ones...
astro.markers <- c("Atp1b2", "Apoe", "Slc1a3","Slc1a2", "Gpc5", "Aldh1l1", "Aqp4", "Gfap","Glud1","Glul","Dbx2")
neuronal.progenitor.markers <- c("Ascl1","Gfap","Dcx","Calb2", "Sox2", "Pcna") # "Mki67" = 0 exprs
schwann <- c( "Gap43", "S100") #"Mpz", "Ncam",
# Nkx3-1 = 0
extra.astro.markers <- c("Vim","Nkx6-1")
astro.markers <- neuronal.progenitor.markers
astro.markers <- schwann

g.1 <- Seurat::VlnPlot(sc_list2[[1]],features = astro.markers,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  theme(text = element_text(size = 18), axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(), axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[1]])
g.2 <- Seurat::VlnPlot(sc_list2[[2]],features = astro.markers,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  theme(text = element_text(size = 18), axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[2]])
g.3 <- Seurat::VlnPlot(sc_list2[[3]],features = astro.markers,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  theme(text = element_text(size = 18), axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[3]])
g.4 <- Seurat::VlnPlot(sc_list2[[4]],features = astro.markers,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  theme(text = element_text(size = 18),axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[4]])
g.5 <- Seurat::VlnPlot(sc_list2[[5]],features = astro.markers,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  theme(text = element_text(size = 18), axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[5]])
g.6 <- Seurat::VlnPlot(sc_list2[[6]],features = astro.markers,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  theme(text = element_text(size = 18), axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[6]])

g.1 + g.2 + g.3 + g.4 + g.5 + g.6

multiViolines(geneList = astro.markers, sc_list = sc_list2, highlight = "group 1")

#### # Function to highlight specific clusters:
clust.highlight <- function(clust,obj){
  true.color <- "#3C7DA6"
  other.color <- "#D9E8F5"
  clust.highlight <- clust
  color.cl <- rep(other.color,length(levels(Idents(obj))))
  color.cl[clust.highlight+1] <- true.color
  return(color.cl)
}
sc_list3 <- sc_list2
# Group definition by similarity graph groups:
sc_list3 <- lapply(sc_list3,function(x){
  x@meta.data$group <- "no group"
  x@meta.data$group2 <- NA
                   return(x)})
# Group 1
sc_list3[[1]]@meta.data[sc_list3[[1]]@meta.data$seurat_clusters %in% c(5),]$group <- "group 1"
sc_list3[[2]]@meta.data[sc_list3[[2]]@meta.data$seurat_clusters %in% c(4),]$group <- "group 1"
sc_list3[[3]]@meta.data[sc_list3[[3]]@meta.data$seurat_clusters %in% c(5),]$group <- "group 1"
sc_list3[[4]]@meta.data[sc_list3[[4]]@meta.data$seurat_clusters %in% c(3),]$group <- "group 1"
sc_list3[[5]]@meta.data[sc_list3[[5]]@meta.data$seurat_clusters %in% c(5),]$group <- "group 1"
sc_list3[[6]]@meta.data[sc_list3[[6]]@meta.data$seurat_clusters %in% c(6),]$group <- "group 1"

###### Check populations on merged data:
library(Seurat)
DimPlot(object = all_merged,cells.highlight = unlist(lapply(1:6,function(x){
  return(paste(sc_list3[[x]]@meta.data %>% filter(group=="group 1") %>% rownames(),x,sep = "_"))
})), cols.highlight = "#F2D377") + ggtitle("Group of interest 1")
FeaturePlot(all_merged,features = "Gfap")
DimPlot(object = all_merged,cells.highlight = unlist(lapply(1:6,function(x){
  return(paste(sc_list3[[x]]@meta.data %>% filter(group=="group 2") %>% rownames(),x,sep = "_"))
})), cols.highlight = "#F2D377") + ggtitle("Group of interest 2")

(DimPlot(object = all_merged,cells.highlight = unlist(lapply(1:6,function(x){
  return(paste(sc_list3[[x]]@meta.data %>% filter(group2=="group 3 high") %>% rownames(),x,sep = "_"))
})), cols.highlight = "#F2D377") + ggtitle("Group 3 low exprs.")) +
(DimPlot(object = all_merged,cells.highlight = unlist(lapply(1:6,function(x){
  return(paste(sc_list3[[x]]@meta.data %>% filter(group2=="group 3 med") %>% rownames(),x,sep = "_"))
})), cols.highlight = "#F28D35") + ggtitle("Group 3 medium exprs.")) +
(DimPlot(object = all_merged,cells.highlight = unlist(lapply(1:6,function(x){
  return(paste(sc_list3[[x]]@meta.data %>% filter(group2=="group 3 high") %>% rownames(),x,sep = "_"))
})), cols.highlight = "#D94D1A") + ggtitle("Group 3 high exprs.")
)## check completed

# Group 2
sc_list3[[1]]@meta.data[sc_list3[[1]]@meta.data$seurat_clusters %in% c(2,3),]$group <- "group 2"
sc_list3[[2]]@meta.data[sc_list3[[2]]@meta.data$seurat_clusters %in% c(2),]$group <- "group 2"
sc_list3[[3]]@meta.data[sc_list3[[3]]@meta.data$seurat_clusters %in% c(3),]$group <- "group 2"
sc_list3[[4]]@meta.data[sc_list3[[4]]@meta.data$seurat_clusters %in% c(2),]$group <- "group 2"
sc_list3[[5]]@meta.data[sc_list3[[5]]@meta.data$seurat_clusters %in% c(4),]$group <- "group 2"
sc_list3[[6]]@meta.data[sc_list3[[6]]@meta.data$seurat_clusters %in% c(4,3),]$group <- "group 2"



# Group 3.1 - Low Exprs
sc_list3[[2]]@meta.data[sc_list3[[2]]@meta.data$seurat_clusters %in% c(0),]$group2 <- "group 3 low"
sc_list3[[3]]@meta.data[sc_list3[[3]]@meta.data$seurat_clusters %in% c(2),]$group2 <- "group 3 low"
sc_list3[[5]]@meta.data[sc_list3[[5]]@meta.data$seurat_clusters %in% c(2,3),]$group2 <- "group 3 low"
sc_list3[[6]]@meta.data[sc_list3[[6]]@meta.data$seurat_clusters %in% c(5),]$group2 <- "group 3 low"

# Group 3.2 - Med. Exprs
sc_list3[[1]]@meta.data[sc_list3[[1]]@meta.data$seurat_clusters %in% c(4),]$group2 <- "group 3 med"
sc_list3[[2]]@meta.data[sc_list3[[2]]@meta.data$seurat_clusters %in% c(3),]$group2 <- "group 3 med"
sc_list3[[3]]@meta.data[sc_list3[[3]]@meta.data$seurat_clusters %in% c(0),]$group2 <- "group 3 med"
sc_list3[[4]]@meta.data[sc_list3[[4]]@meta.data$seurat_clusters %in% c(1),]$group2 <- "group 3 med"
sc_list3[[5]]@meta.data[sc_list3[[5]]@meta.data$seurat_clusters %in% c(0),]$group2 <- "group 3 med"
sc_list3[[6]]@meta.data[sc_list3[[6]]@meta.data$seurat_clusters %in% c(2),]$group2 <- "group 3 med"

# Group 3.2 - high Exprs
sc_list3[[1]]@meta.data[sc_list3[[1]]@meta.data$seurat_clusters %in% c(0,1),]$group2 <- "group 3 high"
sc_list3[[2]]@meta.data[sc_list3[[2]]@meta.data$seurat_clusters %in% c(1),]$group2 <- "group 3 high"
sc_list3[[3]]@meta.data[sc_list3[[3]]@meta.data$seurat_clusters %in% c(1,4),]$group2 <- "group 3 high"
sc_list3[[4]]@meta.data[sc_list3[[4]]@meta.data$seurat_clusters %in% c(0),]$group2 <- "group 3 high"
sc_list3[[5]]@meta.data[sc_list3[[5]]@meta.data$seurat_clusters %in% c(1),]$group2 <- "group 3 high"
sc_list3[[6]]@meta.data[sc_list3[[6]]@meta.data$seurat_clusters %in% c(0,1),]$group2 <- "group 3 high"

########################
## Group of interest 1 *
library(dplyr)
markers.group.1 <- list()
markers.group.1[[1]] <- simil.top10 %>% filter(datasetL == "cx.26" & clusterL == 5) %>% pull(topFeatureConserved)
markers.group.1[[2]] <- simil.top10 %>% filter(datasetL == "sc.26" & clusterL == 4) %>% pull(topFeatureConserved)
markers.group.1[[3]] <- simil.top10 %>% filter(datasetL == "cx.36" & clusterL == 3) %>% pull(topFeatureConserved)
markers.group.1[[4]] <- simil.top10 %>% filter(datasetL == "sc.36" & clusterL == 5) %>% pull(topFeatureConserved)
markers.group.1[[5]] <- simil.top10 %>% filter(datasetL == "cx.699" & clusterL == 6) %>% pull(topFeatureConserved)
markers.group.1[[6]] <- simil.top10 %>% filter(datasetL == "cx.699" & clusterL == 5) %>% pull(topFeatureConserved)
markers.group.1[[7]] <- simil.top10 %>% filter(datasetL == "sc.699" & clusterL == 5) %>% pull(topFeatureConserved)

markers.group.1.common <- Reduce(intersect, lapply(markers.group.1,function(x)return(x)))
markers.group.1.common <- unique(unlist(markers.group.1))

sort(table(unlist(markers.group.1)),decreasing = T)

markers.group.1.common <- names(table(unlist(markers.group.1)))[table(unlist(markers.group.1)) > 5]
markers.group.1.selected.high <- c("Dock10", "Edil3", "Elmo1", "Frmd5", "Nkain2", "Pde4b", "Pex5l", "Plcl1", "Prr5l", "Slc24a2", "St18", "Tmeff2", "Zfp536","Id3", "Tmem117")
mouse_to_human_genes(markers.group.1.selected.high)
markers.group.1.selected.low <- c("Mertk", "Prex2", "Rgs6", "Rmst","Itga6","Aldh1l1")
markers.group.1.common <- c(markers.group.1.selected.high,markers.group.1.selected.low)

markers.group.1.common <- markers.group.1.selected.high
markers.group.1.common <- markers.group.1.selected.low

g.1 <- Seurat::VlnPlot(sc_list2[[1]],features = markers.group.1.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(5),sc_list2[[1]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[1]])
g.2 <- Seurat::VlnPlot(sc_list2[[2]],features = markers.group.1.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(4),sc_list2[[2]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[2]])
g.3 <- Seurat::VlnPlot(sc_list2[[3]],features = markers.group.1.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(5),sc_list2[[3]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[3]])
g.4 <- Seurat::VlnPlot(sc_list2[[4]],features = markers.group.1.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(3),sc_list2[[4]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[4]])
g.5 <- Seurat::VlnPlot(sc_list2[[5]],features = markers.group.1.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(5),sc_list2[[5]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[5]])
g.6 <- Seurat::VlnPlot(sc_list2[[6]],features = markers.group.1.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(6),sc_list2[[6]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[6]])

g.1 + g.2 + g.3 + g.4 + g.5 + g.6

########################
## Group of interest 2 *
markers.group.2 <- list()
markers.group.2[[1]] <- simil.top10 %>% filter(datasetL == "cx.26" & clusterL == 4) %>% pull(topFeatureConserved)
markers.group.2[[2]] <- simil.top10 %>% filter(datasetL == "sc.26" & clusterL == 2) %>% pull(topFeatureConserved)
markers.group.2[[3]] <- simil.top10 %>% filter(datasetL == "cx.36" & clusterL == 2) %>% pull(topFeatureConserved)
markers.group.2[[4]] <- simil.top10 %>% filter(datasetL == "sc.36" & clusterL == 2) %>% pull(topFeatureConserved)
#markers.group.2[[8]] <- simil.top10 %>% filter(datasetL == "sc.36" & clusterL == 3) %>% pull(topFeatureConserved)
markers.group.2[[5]] <- simil.top10 %>% filter(datasetL == "cx.699" & clusterL == 4) %>% pull(topFeatureConserved)
#markers.group.2[[7]] <- simil.top10 %>% filter(datasetL == "cx.699" & clusterL == 3) %>% pull(topFeatureConserved)
markers.group.2[[6]] <- simil.top10 %>% filter(datasetL == "sc.699" & clusterL == 3) %>% pull(topFeatureConserved)
#
sort(table(unlist(markers.group.2)),decreasing = T)

markers.group.2.common <- Reduce(intersect, lapply(markers.group.2,function(x)return(x)))
markers.group.2.common <- names(table(unlist(markers.group.2)))[table(unlist(markers.group.2)) > 6]
markers.group.2.selected.high <- c("Ablim2", "Aqp4", "Gfap", "Grik2", "Kcnj3", "Prune2", "Robo2", "Slc38a1", "Sorbs2")
markers.group.2.selected.low <- c("Brinp3", "Gm12239", "Gm20713",	"Gpc5", "Gria2", "Kcnd2", "Kirrel3",	"Lsamp", "Mdga2","Nhsl1",
                                  "Nrxn1","Pde7b","Pitpnc1","Plcb1",	"Rgs7","Rora","Trpm3")

markers.group.2.common <- c(markers.group.2.selected.high,markers.group.2.selected.low)
markers.group.2.common <- c(markers.group.2.selected.high)
markers.group.2.common <- c(markers.group.2.selected.low)

g.1 <- Seurat::VlnPlot(sc_list2[[1]],features = markers.group.2.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(2),sc_list2[[1]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[1]])
g.2 <- Seurat::VlnPlot(sc_list2[[2]],features = markers.group.2.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(2,sc_list2[[2]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[2]])
g.3 <- Seurat::VlnPlot(sc_list2[[3]],features = markers.group.2.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(3,sc_list2[[3]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[3]])
g.4 <- Seurat::VlnPlot(sc_list2[[4]],features = markers.group.2.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(2,sc_list2[[4]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[4]])
g.5 <- Seurat::VlnPlot(sc_list2[[5]],features = markers.group.2.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(4,sc_list2[[5]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[5]])
g.6 <- Seurat::VlnPlot(sc_list2[[6]],features = markers.group.2.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(4),sc_list2[[6]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[6]])

g.1 + g.2 + g.3 + g.4 + g.5 + g.6


########################
## Group of interest 3 *
markers.group.3 <- list()
sc_list3[[1]]@meta.data[sc_list3[[1]]@meta.data$seurat_clusters %in% c(0,1,4),]$group <- "group 3"
sc_list3[[2]]@meta.data[sc_list3[[2]]@meta.data$seurat_clusters %in% c(0,1,3),]$group <- "group 3"
sc_list3[[3]]@meta.data[sc_list3[[3]]@meta.data$seurat_clusters %in% c(0,1,2),]$group <- "group 3"
sc_list3[[4]]@meta.data[sc_list3[[4]]@meta.data$seurat_clusters %in% c(0,1),]$group <- "group 3"
sc_list3[[5]]@meta.data[sc_list3[[5]]@meta.data$seurat_clusters %in% c(0,1,2,3),]$group <- "group 3"
sc_list3[[6]]@meta.data[sc_list3[[6]]@meta.data$seurat_clusters %in% c(0,1,2),]$group <- "group 3"

markers.group.3[[1]] <- simil.top10 %>% filter(datasetL == "cx.26" & clusterL == 0) %>% pull(topFeatureConserved)
markers.group.3[[2]] <- simil.top10 %>% filter(datasetL == "cx.26" & clusterL == 1) %>% pull(topFeatureConserved)
markers.group.3[[18]] <- simil.top10 %>% filter(datasetL == "cx.26" & clusterL == 2) %>% pull(topFeatureConserved)
markers.group.3[[3]] <- simil.top10 %>% filter(datasetL == "cx.26" & clusterL == 3) %>% pull(topFeatureConserved)
markers.group.3[[4]] <- simil.top10 %>% filter(datasetL == "sc.26" & clusterL == 0) %>% pull(topFeatureConserved)
markers.group.3[[5]] <- simil.top10 %>% filter(datasetL == "sc.26" & clusterL == 1) %>% pull(topFeatureConserved)
markers.group.3[[6]] <- simil.top10 %>% filter(datasetL == "sc.26" & clusterL == 3) %>% pull(topFeatureConserved)
markers.group.3[[7]] <- simil.top10 %>% filter(datasetL == "cx.36" & clusterL == 0) %>% pull(topFeatureConserved)
markers.group.3[[8]] <- simil.top10 %>% filter(datasetL == "cx.36" & clusterL == 1) %>% pull(topFeatureConserved)
markers.group.3[[9]] <- simil.top10 %>% filter(datasetL == "sc.36" & clusterL == 0) %>% pull(topFeatureConserved)
markers.group.3[[10]] <- simil.top10 %>% filter(datasetL == "sc.36" & clusterL == 1) %>% pull(topFeatureConserved)
markers.group.3[[11]] <- simil.top10 %>% filter(datasetL == "sc.36" & clusterL == 4) %>% pull(topFeatureConserved)
markers.group.3[[12]] <- simil.top10 %>% filter(datasetL == "cx.699" & clusterL == 0) %>% pull(topFeatureConserved)
markers.group.3[[13]] <- simil.top10 %>% filter(datasetL == "cx.699" & clusterL == 1) %>% pull(topFeatureConserved)
markers.group.3[[14]] <- simil.top10 %>% filter(datasetL == "cx.699" & clusterL == 2) %>% pull(topFeatureConserved)
markers.group.3[[15]] <- simil.top10 %>% filter(datasetL == "sc.699" & clusterL == 0) %>% pull(topFeatureConserved)
markers.group.3[[16]] <- simil.top10 %>% filter(datasetL == "sc.699" & clusterL == 1) %>% pull(topFeatureConserved)
markers.group.3[[17]] <- simil.top10 %>% filter(datasetL == "sc.699" & clusterL == 2) %>% pull(topFeatureConserved)

#

sort(table(unlist(markers.group.3)),decreasing = T)

markers.group.3.common <- Reduce(intersect, lapply(markers.group.3,function(x)return(x)))
markers.group.3.common <- names(table(unlist(markers.group.3)))[table(unlist(markers.group.3)) > 20]

g.1 <- Seurat::VlnPlot(sc_list2[[1]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(0,1,4),sc_list2[[1]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[1]])
g.2 <- Seurat::VlnPlot(sc_list2[[2]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(0,1,3),sc_list2[[2]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[2]])
g.3 <- Seurat::VlnPlot(sc_list2[[3]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(0,1,2),sc_list2[[3]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[3]])
g.4 <- Seurat::VlnPlot(sc_list2[[4]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(0,1),sc_list2[[4]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[4]])
g.5 <- Seurat::VlnPlot(sc_list2[[5]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(0,1,2,3),sc_list2[[5]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[5]])
g.6 <- Seurat::VlnPlot(sc_list2[[6]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(0,1,2),sc_list2[[6]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[6]])

g.1 + g.2 + g.3 + g.4 + g.5 + g.6
#
violin.per.clusters <- function(obj,clust,gene){
  cells.grp <- obj@meta.data %>% filter(seurat_clusters %in% clust) %>% rownames()
  exprs.plot.data <- cbind.data.frame(exprs=Seurat::GetAssayData(object = obj,slot = "data")[gene,cells.grp],
                   cluster= factor(obj@meta.data[cells.grp,"seurat_clusters"]),
                   group=factor(obj@meta.data[cells.grp,"group2"],levels = c("group 3 low","group 3 med","group 3 high"),ordered = T))
  g <- ggplot(exprs.plot.data, aes(x=reorder(cluster, exprs, mean),y=exprs,fill=group)) +
    geom_violin(alpha=0.6, draw_quantiles =T) +
    scale_fill_manual(values = c("group 3 high" = "#D94D1A",
                                  "group 3 med"="#FDD49E",
                                  "group 3 low"="#7EB6D9")) +
    stat_summary(fun = "mean",
                 geom = "crossbar", 
                 width = 0.2,
                 colour = "black") +
    theme_minimal() + 
    xlab("")
  return(g)
}
#############################################
# PLOT of group 3 by gene expression Sub-group
#       LOW Exprs.      >>.      Med. Exprs.        >>.      High Exprs.
gene <- "Gpc5"
gene <- "Gria2"
gene <- "Vav3"
gene <- "Eps8"
p1 <- violin.per.clusters(obj = sc_list3[[1]],c(0,1,4),gene = gene)
p2 <- violin.per.clusters(obj = sc_list3[[2]],c(0,1,3),gene = gene)
p3 <- violin.per.clusters(obj = sc_list3[[3]],c(0,1,2,4),gene = gene)
p4 <- violin.per.clusters(obj = sc_list3[[4]],c(0,1),gene = gene)
p5 <- violin.per.clusters(obj = sc_list3[[5]],c(0,1,2,3),gene = gene)
p6 <- violin.per.clusters(obj = sc_list3[[6]],c(0,1,2,5),gene = gene)
#
legend_b <- cowplot::get_legend(p2 + theme(legend.position="bottom"))
# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
prow <- cowplot::plot_grid( p1 + theme(legend.position="none"),
                   p2 + theme(legend.position="none"),
                   p3 + theme(legend.position="none"),
                   p4 + theme(legend.position="none"),
                   p5 + theme(legend.position="none"),
                   p6 + theme(legend.position="none"),
                   hjust = -1,
                   nrow = 2
)
p <- cowplot::plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2), labels = paste("Gene",gene));p
################################## End of plot group 3 sub-groups

################################## Check the same genes but in all groups
g.1 <- Seurat::VlnPlot(sc_list3[[1]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(0,1,4),sc_list3[[1]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[1]])
g.2 <- Seurat::VlnPlot(sc_list3[[2]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(0,1,3),sc_list3[[2]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[2]])
g.3 <- Seurat::VlnPlot(sc_list3[[3]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(0,1,2,4),sc_list3[[3]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[3]])
g.4 <- Seurat::VlnPlot(sc_list3[[4]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(0,1),sc_list3[[4]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[4]])
g.5 <- Seurat::VlnPlot(sc_list3[[5]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(0,1,2,3),sc_list3[[5]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[5]])
g.6 <- Seurat::VlnPlot(sc_list3[[6]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(0,1,2,5),sc_list3[[6]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[6]])
g.1 + g.2 + g.3 + g.4 + g.5 + g.6
##################################

markers.group.3.common <- Reduce(intersect, lapply(markers.group.3,function(x)return(x)))
markers.group.3.common <- names(table(unlist(markers.group.3)))[table(unlist(markers.group.3)) > 5]
markers.group.3.common <- markers.group.3.common[-9]


g.1 <- Seurat::VlnPlot(sc_list3[[1]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(0,sc_list3[[1]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[1]])

g.2 <- Seurat::VlnPlot(sc_list3[[2]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(1,sc_list3[[2]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[2]])

g.3 <- Seurat::VlnPlot(sc_list3[[3]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(4,1),sc_list3[[3]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[3]])

g.4 <- Seurat::VlnPlot(sc_list3[[4]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(0),sc_list3[[4]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[4]])

g.5 <- Seurat::VlnPlot(sc_list3[[5]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(1),sc_list3[[5]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[5]])

g.6 <- Seurat::VlnPlot(sc_list3[[6]],features = markers.group.3.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(0),sc_list3[[6]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[6]])

g.1 + g.2 + g.3 + g.4 + g.5 + g.6

##########
######
###
#

sc_list_group3 <- lapply(sc_list3,function(x)subset(x,subset = group == "group 3"))

sc_list_group3 <- lapply(sc_list2,function(x)subset(x,subset = group == "group 3"))
variable.features3 <- Reduce(union,lapply(sc_list_group3,function(x){Seurat::VariableFeatures(x)}));length(variable.features3)
sc_list_group3 <- lapply(sc_list_group3, function(x)return(x[variable.features3,]))

simil.group.3 <- ClusterFoldSimilarity::clusterFoldSimilarity(scList = sc_list_group3,sampleNames = names(sc_list_group3), topN = 1,
                                                              nSubsampling = 17)
findCommunitiesSimmilarity(similarityTable = simil.group.3, leidenResolution = 1)

ClusterFoldSimilarity::plotClustersGraph(simil.group.3)
simil.group.3.top.10 <- ClusterFoldSimilarity::clusterFoldSimilarity(sceList = sc_list_group3,sampleNames = names(sc_list_group3), topN = 15,
                                                              nSubsampling = 17)

markers.group.3.sub <- list()
markers.group.3.sub[[1]] <- simil.group.3.top.10 %>% filter(datasetL == "cx.26" & clusterL == 3) %>% pull(topFeatureConserved)
markers.group.3.sub[[2]] <- simil.group.3.top.10 %>% filter(datasetL == "sc.26" & clusterL == 0) %>% pull(topFeatureConserved)
markers.group.3.sub[[3]] <- simil.group.3.top.10 %>% filter(datasetL == "sc.26" & clusterL == 3) %>% pull(topFeatureConserved)
markers.group.3.sub[[4]] <- simil.group.3.top.10 %>% filter(datasetL == "cx.36" & clusterL == 1) %>% pull(topFeatureConserved)
markers.group.3.sub[[5]] <- simil.group.3.top.10 %>% filter(datasetL == "sc.36" & clusterL == 4) %>% pull(topFeatureConserved)
markers.group.3.sub[[6]] <- simil.group.3.top.10 %>% filter(datasetL == "cx.699" & clusterL == 2) %>% pull(topFeatureConserved)
markers.group.3.sub[[7]] <- simil.group.3.top.10 %>% filter(datasetL == "sc.699" & clusterL == 4) %>% pull(topFeatureConserved)

sort(table(unlist(markers.group.3.sub)),decreasing = T)
markers.group.3.sub.common <- names(table(unlist(markers.group.3.sub)))[table(unlist(markers.group.3.sub)) > 2]

g.1 <- Seurat::VlnPlot(sc_list3[[1]],features = markers.group.3.sub.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(4,sc_list3[[1]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[1]])

g.2 <- Seurat::VlnPlot(sc_list3[[2]],features = markers.group.3.sub.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(3,sc_list3[[2]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[2]])

g.3 <- Seurat::VlnPlot(sc_list3[[3]],features = markers.group.3.sub.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(4),sc_list3[[3]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[3]])

g.4 <- Seurat::VlnPlot(sc_list3[[4]],features = markers.group.3.sub.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(1),sc_list3[[4]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[4]])

g.5 <- Seurat::VlnPlot(sc_list3[[5]],features = markers.group.3.sub.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(3),sc_list3[[5]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[5]])

g.6 <- Seurat::VlnPlot(sc_list3[[6]],features = markers.group.3.sub.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=clust.highlight(c(2),sc_list3[[6]])) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list3)[[6]])

g.1 + g.2 + g.3 + g.4 + g.5 + g.6


# Group definition by similarity graph groups:
sc_list_group3 <- lapply(sc_list_group3,function(x){
  x@meta.data$group2 <- "group 3 high"
  return(x)})
# Group 3.2 - Low. Exprs
sc_list_group3[[1]]@meta.data[sc_list_group3[[1]]@meta.data$seurat_clusters %in% c(4),]$group2 <- "group 3 low"
sc_list_group3[[2]]@meta.data[sc_list_group3[[2]]@meta.data$seurat_clusters %in% c(0,3),]$group2 <- "group 3 low"
sc_list_group3[[3]]@meta.data[sc_list_group3[[3]]@meta.data$seurat_clusters %in% c(4),]$group2 <- "group 3 low"
sc_list_group3[[4]]@meta.data[sc_list_group3[[4]]@meta.data$seurat_clusters %in% c(1),]$group2 <- "group 3 low"
sc_list_group3[[5]]@meta.data[sc_list_group3[[5]]@meta.data$seurat_clusters %in% c(3),]$group2 <- "group 3 low"
sc_list_group3[[6]]@meta.data[sc_list_group3[[6]]@meta.data$seurat_clusters %in% c(2),]$group2 <- "group 3 low"

#
violin.per.clusters <- function(obj,gene){
  cells.grp <- obj@meta.data %>% rownames()
  exprs.plot.data <- cbind.data.frame(exprs=Seurat::GetAssayData(object = obj,slot = "data")[gene,cells.grp],
                                      cluster= factor(obj@meta.data[cells.grp,"seurat_clusters"]),
                                      group=factor(obj@meta.data[cells.grp,"group2"],levels = c("group 3 low","group 3 high"),ordered = T))
  g <- ggplot(exprs.plot.data, aes(x=reorder(cluster, exprs, mean),y=exprs,fill=group)) +
    geom_violin(alpha=0.6, draw_quantiles =T) +
    scale_fill_manual(values = c("group 3 high" = "#D94D1A",
                                 "group 3 low"="#7EB6D9")) +
    stat_summary(fun = "mean",
                 geom = "crossbar", 
                 width = 0.2,
                 colour = "black") +
    #geom_point(alpha = 0.5, position= position_jitter(width = 0.5)) + #position_jitter(height=.5, width=.5)) +
    theme_minimal() + 
    #theme(legend.position = "none") + 
    xlab("")
  return(g)
}
#############################################
# PLOT of group 3 by gene expression Sub-group
#       LOW Exprs.      >>.      Med. Exprs.        >>.      High Exprs.
gene <- "Gpc5" # Works nice
gene <- "Gria2"
gene <- "Vav3"
gene <- "Eps8"
gene <- markers.group.3.common[8]
p1 <- violin.per.clusters(obj = sc_list_group3[[1]],gene = gene)
p2 <- violin.per.clusters(obj = sc_list_group3[[2]],gene = gene)
p3 <- violin.per.clusters(obj = sc_list_group3[[3]],gene = gene)
p4 <- violin.per.clusters(obj = sc_list_group3[[4]],gene = gene)
p5 <- violin.per.clusters(obj = sc_list_group3[[5]],gene = gene)
p6 <- violin.per.clusters(obj = sc_list_group3[[6]],gene = gene)
#
#
legend_b <- cowplot::get_legend(p2 + theme(legend.position="bottom"))
# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
prow <- cowplot::plot_grid( p1 + theme(legend.position="none"),
                            p2 + theme(legend.position="none"),
                            p3 + theme(legend.position="none"),
                            p4 + theme(legend.position="none"),
                            p5 + theme(legend.position="none"),
                            p6 + theme(legend.position="none"),
                            hjust = -1,
                            nrow = 2
)
p <- cowplot::plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2), labels = paste("Gene",gene));p


