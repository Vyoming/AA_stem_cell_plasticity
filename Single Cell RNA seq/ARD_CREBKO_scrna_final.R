#ARD vs ARDKO single cell RNA sequencing
library(Rmagic)
library(plyr)
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
library(reshape)
library(harmony)
library(nichenetr)
library(RColorBrewer)
library(Libra)

#load data
ARD1 <- Read10X_h5("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/ARD_CREBKO_Multi/RNA/count/to_local/Beyaz_OE02_7703_Ctrl_WT/outs/filtered_feature_bc_matrix.h5")
ARD2 <- Read10X(data.dir ="/Users/vyom/data/ARD_KO_updated/AG01_1759_arasco/outs/filtered_feature_bc_matrix/")
ARDKO1 <- Read10X_h5("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/ARD_CREBKO_Multi/RNA/count/to_local/Beyaz_OE02_7700_ARD_Creb1KO/outs/filtered_feature_bc_matrix.h5")
ARDKO2 <- Read10X_h5("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/OE06_creb_multi/RNA/count/Beyaz_OE06_7962_Creb1KO/outs/filtered_feature_bc_matrix.h5")
ARDKO3 <- Read10X_h5("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/OE06_creb_multi/RNA/count/Beyaz_OE06_7963_Creb1KO/outs/filtered_feature_bc_matrix.h5")

ARD1 <- CreateSeuratObject(ARD1, project = "ARD1")
ARD2 <- CreateSeuratObject(ARD2, project = "ARD2")

ARDKO1 <- CreateSeuratObject(ARDKO1, project = "ARDKO1")
ARDKO2 <- CreateSeuratObject(ARDKO2, project = "ARDKO2")
ARDKO3 <- CreateSeuratObject(ARDKO3, project = "ARDKO3")

d <- merge(ARD1, y = c(ARD2, ARDKO1,ARDKO2, ARDKO3), add.cell.ids = c("ARD1", "ARD2", "ARDKO1","ARDKO2", "ARDKO3"), project = "ARD_CREBKO")
d

#filtering and qc
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "mt-")
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0 ,ncol = 3)
d <- subset(d, subset = nCount_RNA > 500 & nCount_RNA < 50000 & nFeature_RNA > 2000 & nFeature_RNA < 7500 & percent.mt < 15)
d
plot1 <- FeatureScatter(d, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(d, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = 'orig.ident',pt.size = 0 ,ncol = 3)


d$orig.ident
Data.list <- SplitObject(d, split.by = "ident")
Data.list <- Data.list[c("ARD1", "ARD2", "ARDKO1","ARDKO2", "ARDKO3")]
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
ARDKO_obj <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                           verbose = TRUE)
# Visulization and Clustering
ARDKO_obj <- RunPCA(ARDKO_obj)
VizDimLoadings(ARDKO_obj, dims = 1:12, reduction = "pca")

#determine number of PCs to use
DimPlot(ARDKO_obj, reduction = "pca")
ElbowPlot(ARDKO_obj, ndims = 50, reduction = "pca")


ARDKO_obj <- RunUMAP(ARDKO_obj, dims = 1:35, n.epochs = 500, n.neighbors = 10)

ARDKO_obj <- FindNeighbors(ARDKO_obj, dims = 1:35)

ARDKO_obj <- FindClusters(ARDKO_obj, resolution = 1)
DimPlot(ARDKO_obj, reduction = "umap", label = TRUE)
DimPlot(ARDKO_obj, reduction = "umap", label = TRUE, group.by = 'Treatment')
DefaultAssay(ARDKO_obj) <- 'RNA'
FeaturePlot(ARDKO_obj, features = c('Ptprc')) + scale_color_gradientn(colors = rev(brewer.pal(11,'Spectral')))
ARDKO_obj$nCount_SCT

markers <- FindAllMarkers(ARDKO_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = .5)
library(readxl)
write.csv(markers,'Gene_de_Sigs_Per_Clust.csv')

new.cluster.ids <- c("Ent.", "Stem 1", "EP", "Ent.", "TA", "Ent.", "Gob.", "TA", "EP", "Ent.", "T", "Ent.", "TA", "TA", "Ent.", "Stem 2", "Ent.", "Ent.", "T", "Ent.", "End.", "EP", "End.", "Tuft", "End.", "End.", "Stem 1", "Mac.", "Gob.", "Pan.")

ARDKO_obj[["Cell_Type"]] <- Idents(ARDKO_obj)
names(new.cluster.ids) <- levels(ARDKO_obj)
ARDKO_obj <- RenameIdents(ARDKO_obj, new.cluster.ids)
ARDKO_obj[["Cell_Type"]] <- Idents(ARDKO_obj)
DimPlot(ARDKO_obj, reduction = "umap", group.by= 'Cell_Type')
ARDKO_obj@meta.data$Cell_Type <- factor(ARDKO_obj@meta.data$Cell_Type, levels = c("Stem 1", "Stem 2", "TA", "EP", "Ent.", "End.", "Gob.",'Pan.',"Tuft",'T', 'Mac.')) 


#MAGIC Imputation
ARDKO_obj@assays$SCT <- NULL
use_python("/Users/vyom/miniconda3/bin/python")
DefaultAssay(ARDKO_obj) <- 'RNA'
ARDKO_obj <- NormalizeData(ARDKO_obj, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(ARDKO_obj)
ARDKO_obj <- ScaleData(ARDKO_obj, features = all.genes)
ARDKO_obj[["RNA"]] <- as(object = ARDKO_obj[["RNA"]], Class = "Assay")
ARDKO_obj <- magic(ARDKO_obj)


#pseudotime:
#pseudotime analysis
library(SeuratWrappers)
library(monocle3)
#install.packages("Signac")
library(Signac)
library(org.Mm.eg.db)
gene_symbol <- as.list(org.Mm.egSYMBOL)
Idents(ARDKO_obj) <- ARDKO_obj$Cell_Type

cds <- as.cell_data_set(ARDKO_obj)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(ARDKO_obj[["RNA"]])
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE,   close_loop = FALSE)
cds <- order_cells(cds, reduction_method = "UMAP")


plot_cells(cds = cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, cell_size = 2)
pseudo_umap <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = FALSE, label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, trajectory_graph_color = "#A8A8A8", graph_label_size=1.5,  cell_size = 2) 
pseudo_umap + scale_fill_gradientn(colors = annColors$Pseudotime)  + geom_polygon(data = hulls2, aes(x = x, y = y, group=CellType), fill=NA, color="black", size=0.3, alpha = 0.3) + geom_shadowtext(data = clusterMedian, aes(x = UMAP_1,y= UMAP_2, label = CellType, group = as.factor(CellType)), size = 3.5, bg.colour="black") + labs(fill = "Pseudotime")
ggsave(file="ARDKO_obj_Pseudotime_Umap_no_trajectory.svg", width=10, height=10)
plot_cells(cds,
           color_cells_by = "Cell_Type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

ARDKO_obj <- AddMetaData(
  object = ARDKO_obj,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)
ARDKO_obj$Pseudotime
Idents(ARDKO_obj) <- ARDKO_obj$Cell_Type
FeaturePlot(ARDKO_obj, c("Pseudotime"), pt.size = .0001, label = TRUE, repel = TRUE) +scale_fill_gradientn(colors = rev(brewer.pal(11,'Spectral'))) + scale_color_gradientn(colors =  rev(brewer.pal(11,'Spectral')))
ggsave(file = 'ARDKO_obj_pseudotime_umap.pdf', width=4, height=4, units="in")

#differential expression
Idents(ARDKO_obj) <- ARDKO_obj$Treatment
DE_all <- FindMarkers(ARDKO_obj, ident.1 = "ARD", ident.2 = "ARDKO", test.use = "MAST", logfc.threshold = .01, min.pct = .01, assay = 'RNA')
write.csv(DE_all, paste0('All_ARDvCREBKO_scrna_DE.csv')) 

#housekeeping
#saveRDS(ARDKO_obj, file = "/Users/vyom/data/Seurat_Objects/may23_AA_V_ARD_crebko_scRNA.rds")
#ARDKO_obj <- readRDS('/Users/vyom/CSHL Dropbox Team Dropbox/Vyom Shah/Seurat_Objects/may23_AA_V_ARD_crebko_scRNA.rds', refhook = NULL)

