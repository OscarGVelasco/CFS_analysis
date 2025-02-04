

markers.group.1.selected.high <- c("Dock10", "Edil3", "Elmo1", "Frmd5", "Nkain2", "Pde4b", "Pex5l", "Plcl1", "Prr5l", "Slc24a2", "St18", "Tmeff2", "Zfp536","Id3", "Tmem117")
markers.group.1.selected.low <- c("Mertk", "Prex2", "Rgs6", "Rmst","Itga6","Aldh1l1")
markers.group.1.common <- c(markers.group.1.selected.high,markers.group.1.selected.low)
#
markers.group.2.selected.high <- c("Ablim2", "Aqp4", "Gfap", "Grik2", "Kcnj3", "Prune2", "Robo2", "Slc38a1", "Sorbs2")
markers.group.2.selected.low <- c("Brinp3", "Gm12239", "Gm20713",	"Gpc5", "Gria2", "Kcnd2", "Kirrel3",	"Lsamp", "Mdga2","Nhsl1",
                                  "Nrxn1","Pde7b","Pitpnc1","Plcb1",	"Rgs7","Rora","Trpm3")
markers.group.2.common <- c(markers.group.2.selected.high,markers.group.2.selected.low)
#

markers.group.1.common
markers.group.2.common
markers.group.3.common


AG1 <- setdiff(setdiff(markers.group.1.common, markers.group.2.common), markers.group.3.common)
AG1[AG1 %in% markers.group.1.selected.high]
AG1[AG1 %in% markers.group.1.selected.low]

AG2 <- setdiff(setdiff(markers.group.2.common, markers.group.1.common), markers.group.3.common)

AG3 <- setdiff(setdiff(markers.group.3.common, markers.group.1.common), markers.group.2.common)

markers.common <- c(AG1, AG2, AG3)

markers.common <- c("Dock10", "Edil3", "Elmo1", "Frmd5", #AG1 selected
"Ablim2", "Aqp4", "Gfap", "Sorbs2", # AG2 selected
"Slc1a2", "Slc4a4","Pcdh7","Ntm")

RColorBrewer::brewer.pal(6, "Pastel1")
agnac="#D9E8F5"
agnac="#FFFFFF"
ag1c="#CCEBC5"
ag2c="#FDCDAC"
ag3c="#DECBE4"

update_geom_defaults("violin", aes(linewidth = 0.1))

g.1 <- Seurat::VlnPlot(sc_list2[[1]],features = markers.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=c(ag3c,ag3c,ag2c,agnac,ag3c,ag1c)) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[1]])
#
g.2 <- Seurat::VlnPlot(sc_list2[[2]],features = markers.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=c(ag3c,ag3c,ag2c,ag3c,ag1c)) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[2]])
#
g.3 <- Seurat::VlnPlot(sc_list2[[3]],features = markers.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=c(ag3c,ag3c,ag3c,ag2c,agnac,ag1c)) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[3]])
#
g.4 <- Seurat::VlnPlot(sc_list2[[4]],features = markers.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=c(ag3c,ag3c,ag2c,ag1c)) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[4]])
#
g.5 <- Seurat::VlnPlot(sc_list2[[5]],features = markers.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=c(ag3c,ag3c,ag3c,ag3c,ag2c,ag1c)) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[5]])
#
g.6 <- Seurat::VlnPlot(sc_list2[[6]],features = markers.common,stack = T,flip = T,fill.by = "ident",pt.size = 0.6) +
  scale_fill_manual(values=c(ag3c,ag3c,ag3c,agnac,ag2c,agnac,ag1c)) +
  theme(axis.title.y = element_blank(), axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y =element_blank(),axis.title.x = element_blank(),legend.position = "none") +
  ggtitle(names(sc_list2)[[6]])
# PLOT
g.1 + g.2 + g.3 + g.4 + g.5 + g.6


