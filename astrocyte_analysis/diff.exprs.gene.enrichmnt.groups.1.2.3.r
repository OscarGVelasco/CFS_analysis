library(gprofiler2)

######################### - Load Point - ##########################
fileName <- load(file = paste0("/omics/groups/OE0436/internal/o872o/cluster_similarity_sc/","mice.Kif5a.controls.data.cortex.spinalcord.processed.6.samples.RData"))
sprintf(fileName)

### Merged Analysis

sc_list2[[1]]@meta.data$tissue <- "spinal cord"
sc_list2[[2]]@meta.data$tissue <- "spinal cord"
sc_list2[[3]]@meta.data$tissue <- "spinal cord"
sc_list2[[4]]@meta.data$tissue <- "cortex"
sc_list2[[5]]@meta.data$tissue <- "cortex"
sc_list2[[6]]@meta.data$tissue <- "cortex"
# Group definition by similarity graph groups:
sc_list2 <- lapply(sc_list2,function(x){
  x@meta.data$group <- "group 3"
  #x@meta.data$group2 <- NA
  return(x)})
# Group 1
sc_list2[[1]]@meta.data[sc_list2[[1]]@meta.data$seurat_clusters %in% c(5),]$group <- "group 1"
sc_list2[[2]]@meta.data[sc_list2[[2]]@meta.data$seurat_clusters %in% c(4),]$group <- "group 1"
sc_list2[[3]]@meta.data[sc_list2[[3]]@meta.data$seurat_clusters %in% c(5),]$group <- "group 1"
sc_list2[[4]]@meta.data[sc_list2[[4]]@meta.data$seurat_clusters %in% c(3),]$group <- "group 1"
sc_list2[[5]]@meta.data[sc_list2[[5]]@meta.data$seurat_clusters %in% c(5),]$group <- "group 1"
sc_list2[[6]]@meta.data[sc_list2[[6]]@meta.data$seurat_clusters %in% c(6),]$group <- "group 1"
# Group 2
sc_list2[[1]]@meta.data[sc_list2[[1]]@meta.data$seurat_clusters %in% c(2),]$group <- "group 2"
sc_list2[[2]]@meta.data[sc_list2[[2]]@meta.data$seurat_clusters %in% c(2),]$group <- "group 2"
sc_list2[[3]]@meta.data[sc_list2[[3]]@meta.data$seurat_clusters %in% c(3),]$group <- "group 2"
sc_list2[[4]]@meta.data[sc_list2[[4]]@meta.data$seurat_clusters %in% c(2),]$group <- "group 2"
sc_list2[[5]]@meta.data[sc_list2[[5]]@meta.data$seurat_clusters %in% c(4),]$group <- "group 2"
sc_list2[[6]]@meta.data[sc_list2[[6]]@meta.data$seurat_clusters %in% c(4),]$group <- "group 2"

sc_list2[[3]]@meta.data[sc_list2[[3]]@meta.data$seurat_clusters %in% c(4),]$group <- "unknown"
sc_list2[[6]]@meta.data[sc_list2[[6]]@meta.data$seurat_clusters %in% c(5,3),]$group <- "unknown"

astrocyte_seurat_objects_mc_sc <- sc_list2
######################### - Save Point - ##########################
save(astrocyte_seurat_objects_mc_sc, file = paste0("/omics/groups/OE0436/internal/o872o/cluster_similarity_sc/rdata/","mice.Kif5a.controls.data.cortex.spinalcord.processed.6.samples.astrocytes.labeled.RData"))
######################### - Save Point - ##########################
astrocyte_seurat_objects_mc_sc <- NULL

cell.groups.percent <- cbind.data.frame(as.data.frame(table(sc_list2[[1]]@meta.data$group)),sample="sc.36")
cell.groups.percent <- rbind.data.frame(cell.groups.percent, cbind.data.frame(as.data.frame(table(sc_list2[[2]]@meta.data$group)),sample="sc.26"))
cell.groups.percent <- rbind.data.frame(cell.groups.percent, cbind.data.frame(as.data.frame(table(sc_list2[[3]]@meta.data$group)),sample="sc.699"))
cell.groups.percent <- rbind.data.frame(cell.groups.percent, cbind.data.frame(as.data.frame(table(sc_list2[[4]]@meta.data$group)),sample="cx.36"))
cell.groups.percent <- rbind.data.frame(cell.groups.percent, cbind.data.frame(as.data.frame(table(sc_list2[[5]]@meta.data$group)),sample="cx.26"))
cell.groups.percent <- rbind.data.frame(cell.groups.percent, cbind.data.frame(as.data.frame(table(sc_list2[[6]]@meta.data$group)),sample="cx.699"))

colores = c("#301934","#86608E","#D5BADB","#3C7DA6","#7EB6D9","#D9E8F5","#4AA147","#92C791",
            "#DBECDA","#D94D1A","#F28D35","#F2D377")
percent_all <- ggplot(cell.groups.percent,aes(x=sample,y=Freq,fill=Var1)) +
  geom_bar(position = position_fill(), stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 12),  
  axis.title.x = element_text(size = 12),
  axis.title.y = element_text(size = 12), 
  legend.title = element_text(size=12), #change legend title font size
  legend.text = element_text(size=12)) +
  scale_fill_manual(values=colores) +
  ylab("percentage") +
  ggtitle("Astrocyte groups percentages");percent_all


tmp_meta <- all_merged@meta.data
library(Seurat)
all_merged <- merge(sc_list2[[1]],c(sc_list2[[2]],sc_list2[[3]],sc_list2[[4]],sc_list2[[5]],sc_list2[[6]]),)
all_merged = JoinLayers(all_merged)
all_merged <- ScaleData(all_merged,features = VariableFeatures(all_merged))
all_merged <- FindVariableFeatures(all_merged, selection.method = "vst", nfeatures = 3000)
all_merged <- RunPCA(all_merged, features = VariableFeatures(object = all_merged))
all_merged <- RunUMAP(all_merged, features = VariableFeatures(object = all_merged))
#all_merged <- single.process(x = all_merged, resolution = 0.8, dims = 16)
DimPlot(all_merged,group.by = "mouse") + ggtitle("Mouse ID Astrocytes")
DimPlot(all_merged,group.by = "tissue", pt.size = 0.8, cols = c("#D5BADB","#7EB6D9")) + ggtitle("Tissue type Astrocytes")
DimPlot(all_merged,group.by = "group", cols = c("#F08080","#F2D377","#7EB6D9","#301934")) + ggtitle("Astrocyte groups")

# Save subset of group 3 for Monocle
spinalc_cortex_controls_group3_merged <- subset(all_merged, subset = group == "group 3")
save(spinalc_cortex_controls_group3_merged, file = "/omics/groups/OE0436/internal/o872o/cluster_similarity_sc/spinalc_cortex_controls_group3_merged.Rdata")

sc_list2.group3 <- lapply(sc_list2,function(x){
  x <- subset(x, subset = group == "group 3") 
  return(x)
  })
variable.features <- Reduce(union,lapply(sc_list2.group3,function(x){Seurat::VariableFeatures(x)}));length(variable.features)
sc_list2.group3 <- lapply(sc_list2.group3, function(x)return(x[variable.features,]))


ClusterFoldSimilarity::clusterFoldSimilarity(sceList = sc_list2.group3[1:3], sampleNames = names(sc_list2.group3)[1:3],
                                             nSubsampling = 20, topN = 1)
sc_list2.group3 <- sc_list2.group3[1:3]
sc_list2.group3[[1]]@meta.data$group2 <- ""
sc_list2.group3[[2]]@meta.data$group2 <- ""
sc_list2.group3[[3]]@meta.data$group2 <- ""

# 1
sc_list2.group3[[1]]@meta.data[sc_list2.group3[[1]]@meta.data$seurat_clusters %in% c(0),]$group2 <- "group 1"
sc_list2.group3[[2]]@meta.data[sc_list2.group3[[2]]@meta.data$seurat_clusters %in% c(1),]$group2 <- "group 1"
sc_list2.group3[[3]]@meta.data[sc_list2.group3[[3]]@meta.data$seurat_clusters %in% c(1),]$group2 <- "group 1"
# 2
sc_list2.group3[[1]]@meta.data[sc_list2.group3[[1]]@meta.data$seurat_clusters %in% c(4,1),]$group2 <- "group 2"
sc_list2.group3[[2]]@meta.data[sc_list2.group3[[2]]@meta.data$seurat_clusters %in% c(3),]$group2 <- "group 2"
sc_list2.group3[[3]]@meta.data[sc_list2.group3[[3]]@meta.data$seurat_clusters %in% c(0),]$group2 <- "group 2"
# 3
sc_list2.group3[[1]]@meta.data[sc_list2.group3[[1]]@meta.data$seurat_clusters %in% c(3),]$group2 <- "group 3"
sc_list2.group3[[2]]@meta.data[sc_list2.group3[[2]]@meta.data$seurat_clusters %in% c(0),]$group2 <- "group 3"
sc_list2.group3[[3]]@meta.data[sc_list2.group3[[3]]@meta.data$seurat_clusters %in% c(2),]$group2 <- "group 3"

sc_list2.group3 <- merge(sc_list2.group3[[1]],c(sc_list2.group3[[2]],sc_list2.group3[[3]]))
sc_list2.group3 <- ScaleData(sc_list2.group3,features = VariableFeatures(sc_list2.group3))
sc_list2.group3 <- FindVariableFeatures(sc_list2.group3, selection.method = "vst", nfeatures = 3000)
sc_list2.group3 <- RunPCA(sc_list2.group3, features = VariableFeatures(object = sc_list2.group3))
sc_list2.group3 <- RunUMAP(sc_list2.group3, features = VariableFeatures(object = sc_list2.group3))
DimPlot(sc_list2.group3, group.by = "group2")

save(sc_list2.group3, file = "/omics/groups/OE0436/internal/o872o/cluster_similarity_sc/spinalcord.only.list2.group3.merged.Rdata")

####### Transcription Factor List
tf.mice <- read.table(file = "/omics/groups/OE0436/internal/o872o/cluster_similarity_sc/Mus_musculus_TF.txt",header = T,sep = "\t")

## Spinal cord vs Cortex
astro.Control.sc.vs.ctx <- FindMarkers(SetIdent(all_merged, value = "tissue"), ident.1 = "spinal cord",ident.2 = "cortex",
                                       test.use = "negbinom", latent.vars = "group", densify=TRUE)
gostres.pos <- gost(query = c(astro.Control.sc.vs.ctx %>% filter(p_val_adj < 0.01 & avg_log2FC > 0) %>% rownames()), 
                organism = "mmusculus", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostres.neg <- gost(query = c(astro.Control.sc.vs.ctx %>% filter(p_val_adj < 0.01 & avg_log2FC < 0) %>% rownames()), 
                    organism = "mmusculus", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)

pos <- gostres.pos$result %>% arrange(p_value) %>% pull(term_name) %>% head(20)
neg <- gostres.neg$result %>% arrange(p_value) %>% pull(term_name) %>% head(20)
setdiff(pos,neg) # Upregulated in Spinal Cord
setdiff(neg,pos) # Upregulated in Cortex
##

### Group 1
group1.astro.Control.sc.vs.ctx <- FindMarkers(SetIdent(subset(all_merged,subset = group == "group 1"),value = "tissue"),ident.1 = "spinal cord",ident.2 = "cortex")

group1.astro.Control.sc.vs.ctx.POS <- group1.astro.Control.sc.vs.ctx %>% filter(p_val_adj < 0.005 & avg_log2FC > 0)
multiViolines(geneList = rownames(group1.astro.Control.sc.vs.ctx.POS)[1:5],sc_list = sc_list2,highlight = "group 1")

group1.astro.Control.sc.vs.ctx.NEG <- group1.astro.Control.sc.vs.ctx %>% filter(p_val_adj < 0.005 & avg_log2FC < 0)
sum(group1.astro.Control.sc.vs.ctx$p_val_adj < 0.01)
gostres <- gost(query = c(group1.astro.Control.sc.vs.ctx %>% filter(p_val_adj < 0.01 & avg_log2FC > 0) %>% rownames()), 
                organism = "mmusculus", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
View(gostres$result)


##################
# group 1 vs ALL
# change the current plan to access parallelization
plan("multiprocess", workers = 8)
plan()

group1.astro.Control.1.vs.all <- FindMarkers(SetIdent(all_merged,value = "group"),ident.1 = "group 1",
                                           test.use = "negbinom", latent.vars = "tissue", densify=TRUE)

group1.astro.Control.sc.vs.ctx.POS

# # Watch out!!
# group1.astro.Control.1.vs.all <- FindMarkers(SetIdent(all_merged,value = "group"),ident.1 = "group 1",
#                                              test.use = "negbinom", latent.vars = c("tissue","mouse"),
#                                              min.pct = -Inf, logfc.threshold = -Inf, min.cells.feature = 1, min.cells.group = 1,
#                                              densify=TRUE)

genedes <- mouse_genes_description(rownames(group1.astro.Control.1.vs.all))
group1.astro.Control.1.vs.all$description <- genedes[rownames(group1.astro.Control.1.vs.all),]$description

multiViolines(geneList = intersect(rownames(group1.astro.Control.sc.vs.ctx.POS),c(group1.astro.Control.1.vs.all %>% filter(p_val_adj < 0.01) %>% rownames())),
              sc_list = sc_list2,highlight = "group 1")

c("Cdk4", "Sirt2", "Dab1") %in% (group1.astro.Control.1.vs.all %>% filter(p_val_adj < 0.005) %>% rownames())
write.table(x = group1.astro.Control.1.vs.all %>% filter(p_val_adj < 0.005), sep = "\t",
            file = "/omics/groups/OE0436/internal/o872o/cluster_similarity_sc/table_results/group1.astro.Control.1.vs.all.significant.txt")

gostres.pos <- gprofiler2::gost(query = c(group1.astro.Control.1.vs.all %>% filter(p_val_adj < 0.005 & avg_log2FC > 0) %>% rownames()), 
                    organism = "mmusculus", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostres.neg <- gprofiler2::gost(query = c(group1.astro.Control.1.vs.all %>% filter(p_val_adj < 0.005 & avg_log2FC < 0) %>% rownames()), 
                    organism = "mmusculus", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)
pos <- gostres.pos$result %>% arrange(p_value) %>% pull(term_name) %>% head(20)
neg <- gostres.neg$result %>% arrange(p_value) %>% pull(term_name) %>% head(20)
setdiff(pos,neg) # Upregulated in Group 1
setdiff(neg,pos) # Upregulated in Group 2
gostplot(gostres.pos)
df <- gostres.pos$result %>% filter(p_value < 0.005 & source != "TF") %>% arrange(p_value) %>% head(10)
df$term_name <- factor(df$term_name, levels = df$term_name)
hallmarks.plot.gost <- ggplot(df, aes(x=-log10(p_value),y=term_name)) +
  geom_point(aes(size=intersection_size,fill = p_value),pch=21) +
  theme_minimal() +
  scale_fill_gradientn(colours=rev(c(#"#3C7DA6", # Dark Blue
    #"#D9E8F5", # Light Blue
    "#F2D377",
    "#F2D377", # Yellow
    "#F2D377", # Yellow
    "#F2D377", # Yellow
    "#F2D377", # Yellow
    "#F28D35", # Ligh red
    "#F28D35",
    "#F28D35", # Ligh red
    "#D94D1A"))
  ) +
  theme(legend.position = "none",
  axis.text.y = element_text(size = 12),  
  axis.title.x = element_text(size = 12),
  axis.title.y = element_blank()) +
  ggtitle("Gene Ontology enrichment analysis on Astrocyte group 1 \n upregulated genes. Spinal cord and cortex.")
hallmarks.plot.gost


neurogenesis.genes <- unlist(strsplit(gostres.pos$result %>% filter(term_name == "neurogenesis") %>% pull(intersection), split = ","))
group1.astro.Control.1.vs.all[neurogenesis.genes,]
multiViolines(geneList = neurogenesis.genes[1:15],sc_list = sc_list2, highlight = "group 1")

indx <- unlist(strsplit(gostres.pos$result[which(gostres.pos$result$term_name == "neurogenesis"),]$intersection, split = ","))
indx2 <- unlist(strsplit(gostres.pos$result[which(gostres.pos$result$term_name == pos[1]),]$intersection, split = ","))
indx3 <- unlist(strsplit(gostres.pos$result[which(gostres.pos$result$term_name == pos[4]),]$intersection, split = ","))
intersect(indx,indx3)

# TRANSCRIPTION FACTORS
tfs.1 <- group1.astro.Control.1.vs.all %>% filter(p_val_adj < 0.005 & row.names(group1.astro.Control.1.vs.all) %in% tf.mice$Symbol) %>% rownames()
multiViolines(geneList = tfs.1,sc_list = sc_list2, highlight = "group 1")

intersect(
group1.astro.Control.1.vs.all %>% filter(p_val_adj < 0.005 & row.names(group1.astro.Control.1.vs.all) %in% tf.mice$Symbol) %>% rownames()
,
group2.astro.Control.2.vs.all %>% filter(p_val_adj < 0.005 & row.names(group2.astro.Control.2.vs.all) %in% tf.mice$Symbol) %>% rownames()
)
gostres.pos <- gost(query = group1.astro.Control.1.vs.all %>% 
                      filter(p_val_adj < 0.005 & row.names(group1.astro.Control.1.vs.all) %in% tf.mice$Symbol & avg_log2FC > 0) %>% rownames(), 
                    organism = "mmusculus", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)
######

group1.astro.Control.1.vs.all[indx,]
multiViolines(indx[1:15], sc_list2)

                                                                                                           
which(indx %in% markers.group.1.common)

group1.astro.Control.1.vs.2 <- FindMarkers(SetIdent(all_merged,value = "group"),ident.1 = "group 1", ident.2 = "group 2",
                                             test.use = "negbinom", latent.vars = "tissue")
                                             # min.pct = -Inf, logfc.threshold = -Inf, min.cells.feature = 1, min.cells.group = 1,
                                             # densify=TRUE)

write.table(x = group1.astro.Control.1.vs.2 %>% filter(p_val_adj < 0.005), sep = "\t",
            file = "/omics/groups/OE0436/internal/o872o/cluster_similarity_sc/table_results/group1and2.astro.Control.1.vs.2.significant.txt")

group1.astro.Control.1.vs.2[c("Unc13c","Agt","Ascl1"),]
group1.astro.Control.1.vs.all[c("Unc13c","Agt","Ascl1"),]
gostres.pos <- gost(query = c(group1.astro.Control.1.vs.2 %>% filter(p_val_adj < 0.005 & avg_log2FC > 0) %>% rownames()), 
                    organism = "mmusculus", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostres.neg <- gost(query = c(group1.astro.Control.1.vs.2 %>% filter(p_val_adj < 0.005 & avg_log2FC < 0) %>% rownames()), 
                    organism = "mmusculus", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)
pos <- gostres.pos$result %>% filter(p_value < 0.005) %>% arrange(p_value) %>% pull(term_name) %>% head(20)
neg <- gostres.neg$result %>% filter(p_value < 0.01) %>% arrange(p_value) %>% pull(term_name) %>% head(20)
setdiff(pos,neg) # Upregulated in Group 1
setdiff(neg,pos) # Upregulated in Group 2


library(fgsea)
gmt.file <- "enrichment_analysis/mh.all.v2023.1.Mm.symbols.gmt"
pathways <- gmtPathways(gmt.file)
str(head(pathways));length(pathways)
# Use reactome pathways??
# Obtaining entrez IDs
ensembl = useMart(biomart = "ensembl", dataset="mmusculus_gene_ensembl")
genes_hgn = getBM(attributes=c("external_gene_name",'entrezgene_id'),
                  filters = 'external_gene_name',
                  values = rownames(group1.astro.Control.1.vs.all),
                  mart = ensembl)
gseaDat <- cbind.data.frame(group1.astro.Control.1.vs.all[genes_hgn$external_gene_name,],genes_hgn)
gseaDat <- gseaDat[!is.na(gseaDat$entrezgene_id),]
#
ranks <- sign(gseaDat$avg_log2FC) * -log10(gseaDat$p_val)
#ranks <- sign(group1.astro.Control.1.vs.all$avg_log2FC) * -log10(group1.astro.Control.1.vs.all$p_val)
names(ranks) <- gseaDat$entrezgene_id
#names(ranks) <- rownames(gseaDat)
ranks <- ranks[!(is.infinite(ranks)|is.na(ranks))]
head(ranks)
barplot(sort(ranks, decreasing = T))
ranks <- sort(ranks, decreasing = T)
plot(density(ranks))
#
library(fgsea)
pathways <- reactomePathways(names(ranks));length(pathways)
fgseaRes.reactome <- fgsea(pathways, ranks, minSize=2, maxSize = 500,scoreType = "std")
fgseaResPlot.reactome <- fgseaRes.reactome[fgseaRes.reactome$padj<0.1,];dim(fgseaResPlot.reactome)

fgseaResPlot.reactome <- fgseaResPlot.reactome[order(fgseaResPlot.reactome$NES,decreasing = T),]
fgseaResPlot.reactome$pathway <- factor(fgseaResPlot.reactome$pathway,levels = fgseaResPlot.reactome$pathway)
hallmarks.plot.reactome <- ggplot(fgseaResPlot.reactome, aes(x=NES,y=pathway)) +
  geom_point(aes(size=size,fill = padj),pch=21) +
  theme_minimal() +
  scale_fill_gradientn(colours=rev(c(#"#3C7DA6", # Dark Blue
    #"#D9E8F5", # Light Blue
    "#F2D377",
    "#F2D377", # Yellow
    "#F2D377", # Yellow
    "#F2D377", # Yellow
    "#F2D377", # Yellow
    "#F28D35", # Ligh red
    "#F28D35",
    "#F28D35", # Ligh red
    "#D94D1A"))
  ) +
  ggtitle("GSEA analysis Hallmark Pathways on Cortex Neurons \n Reactome Pathways")
hallmarks.plot.reactome


# Gene annotation

group1.astro.Control.1.vs.all[markers.group.1.selected.high,]
markers.group.1.selected.high %in% (group1.astro.Control.1.vs.all %>% filter(p_val_adj < 0.01) %>% rownames())
markers.group.1.selected.low %in% (group1.astro.Control.1.vs.all %>% filter(p_val_adj < 0.01) %>% rownames())
group2.astro.Control.2.vs.all["Ntf3",]

install.packages("gprofiler2")
library(gprofiler2)
gostres.pos <- gost(query = c(group1.astro.Control.1.vs.all %>% filter(p_val_adj < 0.01 & avg_log2FC > 0) %>% rownames()), 
                organism = "mmusculus", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostres.neg <- gost(query = c(group1.astro.Control.1.vs.all %>% filter(p_val_adj < 0.01 & avg_log2FC < 0) %>% rownames()), 
                    organism = "mmusculus", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)

pos <- gostres.pos$result %>% arrange(p_value) %>% pull(term_name) %>% head(20)
neg <- gostres.neg$result %>% arrange(p_value) %>% pull(term_name) %>% head(20)
setdiff(pos,neg) # Upregulated in Group 1
setdiff(neg,pos) # Upregulated in Group 2

##################
# group 2 vs group 1
group2.astro.Control.2.vs.1 <- FindMarkers(SetIdent(all_merged,value = "group"),ident.1 = "group 2",ident.2 = "group 1",
                                             test.use = "negbinom", latent.vars = "tissue")
gostres.pos <- gost(query = c(group2.astro.Control.2.vs.1 %>% filter(p_val_adj < 0.005 & avg_log2FC > 0) %>% rownames()), 
                organism = "mmusculus", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostres.neg <- gost(query = c(group2.astro.Control.2.vs.1 %>% filter(p_val_adj < 0.005 & avg_log2FC < 0) %>% rownames()), 
                    organism = "mmusculus", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "fdr", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)
pos <- gostres.pos$result %>% filter(p_value < 0.005) %>% arrange(p_value) %>% pull(term_name) %>% head(20)
neg <- gostres.neg$result %>% filter(p_value < 0.005) %>% arrange(p_value) %>% pull(term_name) %>% head(20)
setdiff(pos,neg) # Upregulated in Group 2
setdiff(neg,pos) # Upregulated in other groups
intersect(pos,neg)
##
# Group 2 VS ALL
group2.astro.Control.2.vs.all <- FindMarkers(SetIdent(all_merged,value = "group"),ident.1 = "group 2",
                                             test.use = "negbinom", latent.vars = "tissue")

group2.astro.Control.2.vs.all %>% 
  filter(p_val_adj < 0.005 & row.names(group2.astro.Control.2.vs.all) %in% tf.mice$Symbol & avg_log2FC >0) %>% rownames()
gostres.TF <- gost(query = c(group2.astro.Control.2.vs.all %>% 
                               filter(p_val_adj < 0.005 & row.names(group2.astro.Control.2.vs.all) %in% tf.mice$Symbol & avg_log2FC >0) %>% rownames()), 
                    organism = "mmusculus", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)

genedes <- mouse_genes_description(rownames(group2.astro.Control.2.vs.all))
group2.astro.Control.2.vs.all$description <- genedes[rownames(group2.astro.Control.2.vs.all),]$description

write.table(x = group2.astro.Control.2.vs.all %>% filter(p_val_adj < 0.005), sep = "\t",
            file = "/omics/groups/OE0436/internal/o872o/cluster_similarity_sc/table_results/group2.astro.Control.2.vs.all.significant.txt")

gostres.pos <- gost(query = c(group2.astro.Control.2.vs.all %>% filter(p_val_adj < 0.005 & avg_log2FC > 0) %>% rownames()), 
                    organism = "mmusculus", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostres.neg <- gost(query = c(group2.astro.Control.2.vs.all %>% filter(p_val_adj < 0.005 & avg_log2FC < 0) %>% rownames()), 
                    organism = "mmusculus", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)
pos <- gostres.pos$result %>% filter(p_value < 0.005) %>% arrange(p_value)f %>% pull(term_name) %>% head(20)
neg <- gostres.neg$result %>% filter(p_value < 0.005) %>% arrange(p_value) %>% pull(term_name) %>% head(20)
setdiff(pos,neg) # Upregulated in Group 2
setdiff(neg,pos) # Upregulated in other groups

###################
# group 3 vs ALL
group3.astro.Control.3.vs.all <- FindMarkers(SetIdent(all_merged,value = "group"),ident.1 = "group 3",
                                             test.use = "negbinom", latent.vars = "tissue")

genedes <- mouse_genes_description(rownames(group3.astro.Control.3.vs.all))
group3.astro.Control.3.vs.all$description <- genedes[rownames(group3.astro.Control.3.vs.all),]$description

write.table(x = group3.astro.Control.3.vs.all %>% filter(p_val_adj < 0.005), sep = "\t",
            file = "/omics/groups/OE0436/internal/o872o/cluster_similarity_sc/table_results/group3.astro.Control.3.vs.all.significant.txt")

group3.astro.Control.3.vs.all %>% 
  filter(p_val_adj < 0.005 & row.names(group3.astro.Control.3.vs.all) %in% tf.mice$Symbol & avg_log2FC >0) %>% rownames()
gostres.TF <- gost(query = c(group3.astro.Control.3.vs.all %>% 
                               filter(p_val_adj < 0.005 & row.names(group3.astro.Control.3.vs.all) %in% tf.mice$Symbol & avg_log2FC >0) %>% rownames()), 
                   organism = "mmusculus", ordered_query = FALSE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                   measure_underrepresentation = FALSE, evcodes = TRUE, 
                   user_threshold = 0.05, correction_method = "g_SCS", 
                   domain_scope = "annotated", custom_bg = NULL, 
                   numeric_ns = "", sources = NULL, as_short_link = FALSE)

multiViolines(sc_list = sc_list2, 
              geneList = group3.astro.Control.3.vs.all %>% filter(p_val_adj < 0.005 & avg_log2FC > 0) %>% 
                arrange(desc(avg_log2FC)) %>% rownames() %>% head(10))

gostres.pos <- gost(query = c(group3.astro.Control.3.vs.all %>% filter(p_val_adj < 0.005 & avg_log2FC > 0) %>% rownames()), 
                organism = "mmusculus", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostres.neg <- gost(query = c(group3.astro.Control.3.vs.all %>% filter(p_val_adj < 0.005 & avg_log2FC < 0) %>% rownames()), 
                    organism = "mmusculus", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "fdr", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)

pos <- gostres.pos$result %>% filter(p_value < 0.005) %>% arrange(p_value) %>% pull(term_name) %>% head(20)
neg <- gostres.neg$result %>% filter(p_value < 0.005) %>% arrange(p_value) %>% pull(term_name) %>% head(20)
setdiff(pos,neg) # Upregulated in Group 3
setdiff(neg,pos) # Upregulated in the other groups

