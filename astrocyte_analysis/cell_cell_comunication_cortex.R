# Load libraries ---------------------------------------------------------------------
library(SingleCellExperiment)
#library(scater)
# library(ggfortify)
library(gridExtra)
#library(scran)
library(Seurat)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(SingleR)
MASTERPATH <- "/omics/groups/OE0436/internal/o872o/KIF5A_analysis/"
setwd(MASTERPATH)
PATH="cellRanger_data/"

############################################################
############# Loading and merging data:

pheno.data <- cbind.data.frame(mouse = c("25",     "27",    "36",     "37",   "26",     "28",   "70",   "699"),
                               status = c("KIF5A","Control","Control","KIF5A","Control","KIF5A","KIF5A","Control"))
rownames(pheno.data) <- pheno.data$mouse
############ Batch 1
# Sample 26
sc_26 <- Read10X(data.dir = paste0(MASTERPATH,PATH,"26_sc/"))
pheno26 <- data.frame(stringsAsFactors = F,row.names = colnames(sc_26),mouse = rep("26_sc",ncol(sc_26)),status= rep("Control",ncol(sc_26)))
sc_26 <- Seurat::CreateSeuratObject(counts = sc_26,meta.data = pheno26)
# Sample 36
sc_36 <- Read10X(data.dir = paste0(MASTERPATH,PATH,"36/"))
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

tmp <- all_merged@meta.data %>% dplyr::filter(tissue %in% "ctx") %>% rownames()
all_merged <- all_merged[,tmp]

all_merged <- FindNeighbors(all_merged, dims = seq(1:25))
all_merged <- FindClusters(all_merged, resolution = 0.5)
g <- DimPlot(all_merged, group.by = "RNA_snn_res.0.5", label = T)

######### SELECTING NEURONS
# NEURONAL CELLS
neuronal_markers <- c(
  "Slc17a7", # excitatory neurons
  "Camk2a", #excitatory neurons
  "Gad1",# inhibitory neurons 
  "Gad2", # inhibitory neurons
  "Rbfox3",
  "Dlg4",
  "Slc12a5",
  "Sgip1"
)
g + DotPlot(all_merged, features = neuronal_markers, group.by = "RNA_snn_res.0.5")+ scale_color_viridis_c()
###### SELECTING NEURONS END
all_merged@meta.data$cell.type <- "Other"
all_merged@meta.data[rownames(tmp_meta),]$cell.type <- tmp_meta$group

neuron_id <- all_merged@meta.data %>% dplyr::filter(RNA_snn_res.0.5 %in% c(2,3,4,5,7,8,10,11,13,16,17,20)) %>% rownames()
all_merged@meta.data[neuron_id,]$cell.type <- all_merged@meta.data[neuron_id,]$RNA_snn_res.0.5
save(all_merged, file = "/omics/groups/OE0436/internal/o872o/cluster_similarity_sc/rdata/all_merged_clustered.RData")
DimPlot(all_merged, group.by = "cell.type", label = T)
#
interest_id <- all_merged@meta.data %>% dplyr::filter(cell.type %in% c("group 1", "group 2", "group 3",
                                                                       "2","3","4","5","7","8","10","11","13","16","17","20")) %>% rownames()
all_merged_neurons_astro <- all_merged[, interest_id]
DimPlot(all_merged_neurons_astro, group.by = "cell.type", label = T)
save(all_merged_neurons_astro, file = "/omics/groups/OE0436/internal/o872o/cluster_similarity_sc/rdata/all_merged_neurons_astro_only.RData")


##################################################################
###################### CELL - CELL COMMUNICATION ANALYSIS
devtools::install_github("jinworks/CellChat")
library(CellChat)
# Extract the data from Seurat:
data.input <- as.matrix(GetAssayData(all_merged_neurons_astro, assay = "RNA", slot = "data")) # normalized data matrix
labels <- all_merged_neurons_astro@meta.data[colnames(data.input),"cell.type",drop=FALSE]
labels$cell.type[labels$cell.type %in% c("11", "17", "3", "4", "5", "8")] <- paste("neuron cluster",labels$cell.type[labels$cell.type %in% c("11", "17", "3", "4", "5", "8")])
labels$cell.type[labels$cell.type %in% c("group 1")] <- "astrocyte group 1"
labels$cell.type[labels$cell.type %in% c("group 2")] <- "astrocyte group 2"
labels$cell.type[labels$cell.type %in% c("group 3")] <- "astrocyte group 3"

#meta <- data.frame(labels = labels, row.names = rownames(all_merged_neurons_astro@meta.data)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = labels, group.by = "cell.type")
cellchat <- addMeta(cellchat, meta = labels, meta.name = "cell.type")
cellchat <- setIdent(cellchat, ident.use = "cell.type") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels


groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse #CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
#
# set the used database in the object
cellchat@DB <- CellChatDB
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multiprocess", workers = 4) # do parallel
future::multisession("multiprocess", workers = 16)
#
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#
cellchat <- computeCommunProb(cellchat, population.size = TRUE)


# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
df.netP %>% dplyr::filter(source %in% c("group 1","group 2","group 3") & pval < 0.005) %>% group_by(source, pathway_name) %>% dplyr::count(pathway_name)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. 
# Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


############### Advance Visualisation
cellchat@netP$pathways
# Group 1
# NRG # in group 1 !
# MAG
# MPZ
# CDH
# SEMA6
# CNTN # in group 3
#
# Group 1 & 2 only
# BMP
pathways.show <- c("NRG") 
pathways.show <- c("SEMA6")
pathways.show <- c("CDH")
pathways.show <- c("MAG") #NOP
pathways.show <- c("MPZ")
pathways.show <- c("CNTN") #AG3
pathways.show <- c("NCAM") #AG3
pathways.show <- c("NEGR") #AG3
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1:3) # a numeric vector. 
vertex.receiver = seq(1:9) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
#> 

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = T,enriched.only = T)
LR.show <- pairLR.CXCL$geneLR # show one ligand-receptor pair

Seurat::VlnPlot(all_merged_neurons_astro, features = c("Sema6a","Plxna2", "Plxna4"), group.by = "cell.type", stack = T,flip = T, fill.by = "ident",pt.size = 0.6) + 
  scale_color_viridis_c() +
  theme(text = element_text(size = 18), axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(), axis.title.x = element_blank(),legend.position = "none") +
  scale_x_discrete(labels=c("8" = "neuron cluster 8", 
                            "5" = "neuron cluster 5",
                            "4" = "neuron cluster 4",
                            "3" = "neuron cluster 3",
                            "17" = "neuron cluster 17",
                            "11" = "neuron cluster 11",
                            "group 3" = "astrocyte group 3",
                            "group 2" = "astrocyte group 2",
                            "group 1" = "astrocyte group 1")) +
  ggtitle("SEMA 6 - ligand receptor")

# Hierarchy plot
vertex.receiver = seq(1:9) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
netVisual_bubble(cellchat, sources.use = c(7,8,9), targets.use = c(1:6), remove.isolate = FALSE)


# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = c(1:3), targets.use = c(4:9), signaling = c("NRG"), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(1,2,3), targets.use = c(4:9), signaling = c("SEMA6"), remove.isolate = FALSE, angle.x = 45)


