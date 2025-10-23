library(Rmagic)
library(Seurat)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(MAST)
library(DESeq2)
library(EnhancedVolcano)
library(limma)
library(scales)
library(metR)
library(ggpubr)
library(rstatix)
library(svglite)
library(viridis)

# Load in the Data
Vyom_color_Umap <- CustomPalette(low = "purple", high = "yellow", mid = "black", k = 100)

contr4.data <- Read10X(data.dir ="~/data/intestine/SB04/Control/filtered_feature_bc_matrix/")
contr4 <- CreateSeuratObject(counts = contr4.data, project = "Control4")
contr4

AA4.data <- Read10X(data.dir = "~/data/intestine/SB04/AA/filtered_feature_bc_matrix/")
AA4 <- CreateSeuratObject(counts = AA4.data, project = "AA4")
AA4

Data4 <- merge(contr4, y = c(AA4), add.cell.ids = c("Control", "AA"), project = "SB04")
Data4
Data4[["percent.mt"]] <- PercentageFeatureSet(Data4, pattern = "mt-")
VlnPlot(Data4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Data4 <- subset(Data4, subset = nCount_RNA > 10000 & nCount_RNA < 100000 & nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 15)
Data4
VlnPlot(organoid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'type', pt.size = 0)
AA5.data <- Read10X(data.dir = "~/data/intestine/SB05/AA/")
AA5 <- CreateSeuratObject(counts = AA5.data, project = "AA5")
AA5

Control5.data <- Read10X(data.dir = "~/data/intestine/SB05/Control/")
C5 <- CreateSeuratObject(counts = Control5.data, project = "C5")
C5

P5.data <- Read10X(data.dir = "~/data/intestine/SB05/Pge2/")
P5 <- CreateSeuratObject(counts = P5.data, project = "P5")
P5

Data5 <- merge(C5, y = c(AA5, P5), add.cell.ids = c("Control", "AA", "PGE2"), project = "SB05")
Data5
Data5[["percent.mt"]] <- PercentageFeatureSet(Data5, pattern = "mt-")

VlnPlot(Data5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Data5 <- subset(Data5, subset = nCount_RNA > 3000 & nCount_RNA < 35000 & nFeature_RNA > 1500 & nFeature_RNA < 7000 & percent.mt < 15)
Data5

Data <- merge(Data4, y = c(Data5), add.cell.ids = c("SB04", "SB05"), project = "intest")
Data
# Begin Preprocessing and Filtering
Data[["percent.mt"]] <- PercentageFeatureSet(Data, pattern = "mt-")
VlnPlot(Data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
head(Data@meta.data, 5)
plot1 <- FeatureScatter(Data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
Data

# Batch correct data
Data.list <- SplitObject(Data, split.by = "ident")
Data.list <- Data.list[c("AA4", "AA5", "C5", "Control4", "P5")]
for (i in 1:length(Data.list)) {
  Data.list[[i]] <- SCTransform(Data.list[[i]], verbose = FALSE)
}

# Normilization
#select highly variable genes 
Data.features <- SelectIntegrationFeatures(object.list = Data.list, nfeatures = 3000)
options (future.globals.maxSize = 4000 * 1024^10)
Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                verbose = FALSE)
Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                       anchor.features = Data.features, verbose = FALSE)
organoid <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                                 verbose = TRUE)

# Visulization and Clustering
organoid <- RunPCA(organoid, verbose = FALSE)

organoid <- RunUMAP(organoid, dims = 1:10)

levels(factor(organoid@meta.data$orig.ident))
Idents(organoid) <- organoid$orig.ident
organoid[["Type"]] <- Idents(organoid)
new.cluster.ids <- c("AA", "AA", "Control", "Control", "Pge2")
names(new.cluster.ids) <- levels(organoid)
organoid <- RenameIdents(organoid, new.cluster.ids)

organoid[["Type"]] <- Idents(organoid)
head(Idents(organoid))
plots <- DimPlot(organoid, group.by = c("Type"))
plots 

organoid <- FindNeighbors(organoid, dims = 1:10)
organoid <- FindClusters(organoid, resolution = .3)#.3
DimPlot(organoid, reduction = "umap", label = TRUE)

stem.markers <- c("Lgr5" ,"Ascl2" ,"Slc12a2" ,"Axin2" ,"Olfm4" ,"Gkn3")
TA.markers <- c("Tubb5" ,"Hmgb2" ,"Stmn1" ,"H2afz1" ,"Tuba1b" ,"Hmgb11" ,"Hmgn22" ,"Ptma1" ,"Kiaa0101" ,"Tk1" ,"Cenpw" ,"Tyms" ,"Ranbp11" ,"Idh21" ,"Ran1" ,"Dtymk" ,"Nucks11" ,"Dut1" ,"Cks1b")
paneth.markers <- c("Lyz1" ,"Defa17" ,"Defa22" ,"Defa24" ,"Ang4")
Enterocyte.markers <- c("Alpi" ,"Apoa1" ,"Apoa4" ,"Fabp1" ,"Adh6a")
EP.markers <- c("Ccnb1" ,"Cdc20" ,"Cenpa" ,"Cdkn3" ,"Ccnb2" ,"Cdc25c" ,"Kif22" ,"Ube2c" ,"Sapcd2" ,"Rbp7" ,"Aurka" ,"Ccna2" ,"Cdkn2d" ,"Kif23" ,"Nek2" ,"Slc16a1")
Enteroendocrine.markers <- c("Chga" ,"Chgb" ,"Tac1" ,"Tph1" ,"Neurog3" ,"Gch1" ,"Slc39a2")
Goblet.markers <- c("Muc2" ,"Clca3" ,"Tff3" ,"Agr2" ,"Spink4" ,"Fcgbp" ,"Zg16" ,"Ccl9" ,"Atoh1")
Tuft.markers <- c("Trpm5" ,"Gfi1b" ,"Il25" ,"Alox5ap" ,"Lrmp" ,"Rgs13" ,"Ltc4s" ,"Adh1")
library(readxl)
Roulis_gene_lists_for_metagenes <- read_excel("Desktop/Roulis_gene_lists_for_metagenes.xlsx")
revival.markers_roulis <- c(Roulis_gene_lists_for_metagenes$Gene)

organoid <- AddModuleScore(object = organoid, features = list(pathways$LGR5_stem_cell_signature), name = 'Stem')
organoid <- AddModuleScore(object = organoid, features = list(revival.markers), name = 'revival.markers', assay =  'RNA')
organoid <- AddModuleScore(object = organoid, features = list(TA.markers), name = 'TA.markers')
organoid <- AddModuleScore(object = organoid, features = list(paneth.markers), name = 'paneth.markers')
organoid <- AddModuleScore(object = organoid, features = list(Enterocyte.markers), name = 'Enterocyte.markers')
organoid <- AddModuleScore(object = organoid, features = list(EP.markers), name = 'EP.markers')
organoid <- AddModuleScore(object = organoid, features = list(Enteroendocrine.markers), name = 'Enteroendocrine.markers')
organoid <- AddModuleScore(object = organoid, features = list(Goblet.markers), name = 'Goblet.markers')
organoid <- AddModuleScore(object = organoid, features = list(Tuft.markers), name = 'Tuft.markers')
plasma <- viridis(10, direction = 1, option = "C")
FeaturePlot(object = organoid, features = 'stem1', pt.size = .001)  + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'Organoid_Stem.pdf', width=2, height=2, units="in")

FeaturePlot(object = organoid, features = 'Goblet.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'Organoid_Goblet.pdf', width=2, height=2, units="in")

FeaturePlot(object = organoid, features = 'paneth.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'Organoid_paneth.pdf', width=2, height=2, units="in")

FeaturePlot(object = organoid, features = 'Enterocyte.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'Organoid_Enterocyte.pdf', width=2, height=2, units="in")

FeaturePlot(object = organoid, features = 'Enteroendocrine.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'Organoid_Enteroendocrine.pdf', width=2, height=2, units="in")

FeaturePlot(object = organoid, features = 'revival.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'Organoid_RevStem.pdf', width=2, height=2, units="in")

DotPlot(organoid, group.by = "Cell_Type", features = c('stem1' ,'revival.markers1' ,'reserve.markers1' ,'TA.markers1' ,'paneth.markers1' ,'Enterocyte.markers1' ,'EP.markers1' ,'Enteroendocrine.markers1' ,'Goblet.markers1' ,'Tuft.markers1'))

Cell_Type_Sig_Score <- data.frame(Stem=organoid$stem1, Revival=organoid$revival.markers1, Reserve=organoid$reserve.markers1, Transit_Amplifying =organoid$TA.markers1, Paneth=organoid$paneth.markers1, Enterocyte=organoid$Enterocyte.markers1, EP=organoid$EP.markers1 ,Enteroendocrine= organoid$Enteroendocrine.markers1 ,Goblet=organoid$Goblet.markers1 ,Tuft=organoid$Tuft.markers1)

# module score distribution
modulescores <- Cell_Type_Sig_Score %>%
  rownames_to_column(var="id") %>%
  pivot_longer(-id, names_to="celltype", values_to="score")


p <- ggplot(modulescores)
#p <- ggplot(onescore)
p + geom_point(aes(x=fct_inorder(id), y=sort(score))) +
  facet_wrap(~celltype) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p
Idents(organoid) <- organoid$Type

write.csv(z_RevSC_cluster)
z_Stem_cluster<-WhichCells(object = organoid, expression = stem1 > .5)
z_RevSC_cluster<-WhichCells(object = organoid , expression = revival.markers1 > 11) #idents = c('AA' ,'Pge2')  ,
z_reserve_cluster<-WhichCells(object = organoid , expression = reserve.markers1 > 1)
z_TA_cluster<-WhichCells(object = organoid, expression = TA.markers1 > 3.1)
z_paneth_cluster<-WhichCells(object = organoid, expression = paneth.markers1 > 13.5)
z_Enterocyte_cluster<-WhichCells(object = organoid, expression = Enterocyte.markers1 > 1)
z_EP_cluster<-WhichCells(object = organoid, expression = EP.markers1 > 3.7)
z_Enteroendocrine_cluster<-WhichCells(object = organoid, expression = Enteroendocrine.markers1 > 2)
z_Goblet_cluster<-WhichCells(object = organoid, expression = Goblet.markers1 > 12.9)
z_Tuft_cluster<-WhichCells(object = organoid, expression = Tuft.markers1 > 14)

table(organoid$Cell_type)

organoid <- FindClusters(organoid, resolution = .3)#.3
DimPlot(organoid, reduction = "umap", label = TRUE)

DimPlot(organoid, label=T, group.by="Cell_Type", cells.highlight= list(z_Tuft_cluster, z_Stem_cluster, z_RevSC_cluster, z_reserve_cluster, z_TA_cluster, z_paneth_cluster, z_Enterocyte_cluster, z_EP_cluster, z_Enteroendocrine_cluster, z_Goblet_cluster  ),  cols.highlight = c("darkblue", "darkred", 'green' ,'lightgreen', 'darkgreen', 'yellow', 'pink', 'purple', 'orange' ,'lightblue') ,cols= "grey")
z_RevSC_cluster1<-WhichCells(object = organoid ,idents = c('Pge2'), expression = revival.markers1 > 5)
z_RevSC_cluster2<-WhichCells(object = organoid ,idents = c('AA'), expression = revival.markers1 > 5)
z_RevSC_cluster3<-WhichCells(object = organoid ,idents = c('Control'), expression = revival.markers1 > 5)

DimPlot(organoid, label=T, cells.highlight= list(z_RevSC_cluster),  cols.highlight = c("darkblue") ,cols= "grey")

#Annotate Clusters
Cluster_5 <- subset(organoid, idents = 5)
Cluster_5
Cluster_5 <- RunPCA(Cluster_5, verbose = FALSE)

Cluster_5 <- RunUMAP(Cluster_5, dims = 1:10)

Cluster_5 <- FindNeighbors(Cluster_5, dims = 1:10)
Cluster_5 <- FindClusters(Cluster_5, resolution = .3)
DimPlot(Cluster_5, reduction = "umap", label = TRUE)


Cluster_5 <- AddModuleScore(object = Cluster_5, features = list(Goblet.markers), name = 'Goblet.markers')
FeaturePlot(object = Cluster_5, features = 'Goblet.markers1')
Cluster_5 <- AddModuleScore(object = Cluster_5, features = list(Tuft.markers), name = 'Tuft.markers')
FeaturePlot(object = Cluster_5, features = 'Tuft.markers1')
Cluster_5 <- AddModuleScore(object = Cluster_5, features = list(paneth.markers), name = 'paneth.markers')
FeaturePlot(object = Cluster_5, features = 'paneth.markers1')
new.cluster.ids <- c("Goblet" ,"Goblet" ,"Secretory Progenitor" ,"Paneth" ,"Paneth" ,"Tuft")
Cluster_5[["Cell_Type"]] <- Idents(Cluster_5)
names(new.cluster.ids) <- levels(Cluster_5)
Cluster_5 <- RenameIdents(Cluster_5, new.cluster.ids)
Cluster_5[["Cell_Type"]] <- Idents(Cluster_5)
DimPlot(Cluster_5, reduction = "umap", label = TRUE)


#pre annotate clusters
new.cluster.ids <- c("Stem 1" ,"Transit Amplifying" ,"Enterocyte Progenitor" ,"Stem 1" ,"Enterocyte", "" ,"Stem 1", "Enteroendocrine", "Transit Amplifying", "Stem Like Progenitor")
organoid[["Cell_type"]] <- Idents(organoid)
names(new.cluster.ids) <- levels(organoid)
organoid <- RenameIdents(organoid, new.cluster.ids)
organoid[["Cell_type"]] <- Idents(organoid)

# Generate a new column called Cell_type in the metadata
organoid$Cell_Type <- as.character(Idents(organoid))

# Change the information of cells containing sub-cluster information z_reserve_cluster
organoid$Cell_Type[WhichCells(Cluster_5)] <- paste(Idents(Cluster_5))
Idents(organoid) <- organoid$Cell_type
Idents(organoid) <- organoid$integrated_snn_res.0.3
#6,8
Stem_2_x <- WhichCells(object = organoid ,idents = c('6'), z_RevSC_cluster)

Stem3_x <- WhichCells(object = organoid ,idents = c('8'), z_RevSC_cluster)

Stem2_cells <- c(Stem_2_x, z_EP_cluster)
Stem3_cells <- c(Stem3_x, z_Tuft_cluster)

Idents(organoid) <- organoid$Cell_Type
organoid$Cell_Type[Stem2_cells] <- paste('Stem 2')
organoid$Cell_Type[Stem3_cells] <- paste('Stem 3')

table(organoid$Cell_Type, organoid$Type)
#organoid$Cell_Type[WhichCells(object = organoid ,idents = c('AA' ,'Pge2') , expression = reserve.markers1 > 1.5)] <- paste('Stem 3')

DimPlot(organoid, group.by = "Cell_Type", label = FALSE, pt.size=1, label.size = 3)
Idents(organoid) <- organoid$Cell_Type
organoid$Cell_type <- organoid$Cell_Type
my_levels <- c("Stem 1", "Stem 2", "Stem 3", "Stem Like Progenitor", "Secretory Progenitor", "Transit Amplifying", "Enterocyte Progenitor" , "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Paneth" )
organoid$Cell_type <- factor(x = organoid$Cell_Type, levels = my_levels)
DimPlot(organoid, group.by = "Cell_type", label = FALSE, pt.size=1, label.size = 2)
Idents(organoid) <- organoid$Type
my_levels <- c('Control' ,'AA' ,'Pge2')
organoid$type <- factor(x = organoid$Type, levels = my_levels)

#Cluster Markers
Idents(organoid) <- organoid$Cell_type

markers_stem <- FindAllMarkers(organoid, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
markers_stem %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

write.csv(markers_stem,'organoid_Sigs_Per_Clust.csv')

#Gene level analysis
Idents(Gene_de) <- Gene_de$type
Gene_de <- subset(organoid,  idents =c('Control', 'AA'))
Gene_de <- NormalizeData(object = Gene_de ,normalization.method = "LogNormalize", assay = "RNA")
# impute for violin plots
my_levels <- c('Control' ,'AA')
Gene_de$type <- factor(x = Gene_de$Type, levels = my_levels)
DefaultAssay(Gene_de) <- "RNA"
Gene_de <- magic(Gene_de)

#pseudotime analysis
library(SeuratWrappers)
library(monocle3)
library(Signac)
library(org.Mm.eg.db)
gene_symbol <- as.list(org.Mm.egSYMBOL)
Idents(Gene_de) <- Gene_de$Cell_type
pseudo_stem <- WhichCells(object = Gene_de, idents = 'Stem 1')


cds <- as.cell_data_set(Gene_de)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(Gene_de[["RNA"]])
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE)
cds <- order_cells(cds, reduction_method = "UMAP")

saveRDS(cds, file = "Organoid_pseudotime.rds")
plot_cells(cds = cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, cell_size = 2)
pseudo_umap <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = FALSE, label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, trajectory_graph_color = "#A8A8A8", graph_label_size=1.5,  cell_size = 2) 
pseudo_umap + scale_fill_gradientn(colors = annColors$Pseudotime)  + geom_polygon(data = hulls2, aes(x = x, y = y, group=CellType), fill=NA, color="black", size=0.3, alpha = 0.3) + geom_shadowtext(data = clusterMedian, aes(x = UMAP_1,y= UMAP_2, label = CellType, group = as.factor(CellType)), size = 3.5, bg.colour="black") + labs(fill = "Pseudotime")
ggsave(file="Gene_de_Pseudotime_Umap_no_trajectory.svg", width=10, height=10)
plot_cells(cds,
           color_cells_by = "Cell_type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

Gene_de <- AddMetaData(
  object = Gene_de,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)

