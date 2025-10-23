#CREBKO_IR_scrna - Vyom Shah
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
library(Nebulosa)

# Load in the Data
WT_1 <- Read10X_h5("/Users/vyom/data/crebIR_scrna/Beyaz_OE09_SB1365_Villin_CreERT2neg_Creb_WT/filtered_feature_bc_matrix.h5")
KO_1 <- Read10X_h5("/Users/vyom/data/crebIR_scrna/Beyaz_OE09_SB1366_Villin_CreERT2pos_Creb_KO/filtered_feature_bc_matrix.h5")

WT_2 <- Read10X_h5("/Users/vyom/data/crebIR_scrna/Beyaz_OE09_SB1367_Villin_CreERT2neg_Creb_WT/filtered_feature_bc_matrix.h5")
KO_2 <- Read10X_h5("/Users/vyom/data/crebIR_scrna/Beyaz_OE09_SB1368_Villin-CreERT2pos_Creb_KO/filtered_feature_bc_matrix.h5")


WT_1 <- CreateSeuratObject(WT_1, project = "WT_1")
KO_1 <- CreateSeuratObject(KO_1, project = "KO_1")

WT_2 <- CreateSeuratObject(WT_2, project = "WT_2")
KO_2 <- CreateSeuratObject(KO_2, project = "KO_2")


d <- merge(WT_1, y = c(KO_1, WT_2, KO_2 ), add.cell.ids = c("WT_1", 'KO_1', 'WT_2', 'KO_2'), project = "Immune")
d

#QC and Filtering
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "mt-")
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0 ,ncol = 3)
d <- subset(d, subset = nCount_RNA > 500 & nCount_RNA < 12500 & nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 15)
d

Data.list <- SplitObject(d, split.by = "ident")
Data.list <- Data.list[c("WT_1", 'KO_1', 'WT_2', 'KO_2')]
for (i in 1:length(Data.list)) {
  
  Data.list[[i]] <- SCTransform(Data.list[[i]], verbose = FALSE)
}

# Normilization
#select highly variable genes 
Data.features <- SelectIntegrationFeatures(object.list = Data.list, nfeatures = 2000)
Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                verbose = FALSE)
Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                       anchor.features = Data.features, verbose = FALSE)
KOIR_obj <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                          verbose = TRUE)

# Visulization and Clustering
KOIR_obj <- RunPCA(KOIR_obj)
VizDimLoadings(KOIR_obj, dims = 1:2, reduction = "pca")

DimPlot(KOIR_obj, reduction = "pca")
ElbowPlot(KOIR_obj, ndims = 50, reduction = "pca")

levels(factor(KOIR_obj@meta.data$orig.ident))
Idents(KOIR_obj) <- KOIR_obj$orig.ident
KOIR_obj[["Treatment"]] <- Idents(KOIR_obj)
new.cluster.ids <- c("WT",  "KO",  "WT", "KO")
names(new.cluster.ids) <- levels(KOIR_obj)
KOIR_obj <- RenameIdents(KOIR_obj, new.cluster.ids)

KOIR_obj[["Treatment"]] <- Idents(KOIR_obj)
head(Idents(KOIR_obj))
DimPlot(KOIR_obj, group.by = c("Treatment"))

KOIR_obj@meta.data$Treatment <- factor(KOIR_obj@meta.data$Treatment, levels = c("WT",  "KO")) 

DefaultAssay(KOIR_obj) <- 'integrated'
KOIR_obj <- RunUMAP(KOIR_obj, dims = 1:25)

KOIR_obj <- FindNeighbors(KOIR_obj, dims = 1:25)

KOIR_obj <- FindClusters(KOIR_obj, resolution = 1)
DimPlot(KOIR_obj, reduction = "umap", label = TRUE)
DimPlot(KOIR_obj, reduction = "umap",group.by = 'Treatment', label = TRUE)

markers <- FindAllMarkers(KOIR_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = .5)
library(readxl)
write.csv(markers,'Gene_de_Sigs_Per_Clust_crebKo_IR.csv')

my_vector <- c( "Imm.", "Imm.", "Imm.", "Imm.", "Imm.", "Imm.", "Imm.", "Imm.", "Imm.", "Imm.", "Imm.", "Imm.",  "Imm.","Ent.", "Imm.", "Imm.", "Sec.", "Imm.", "Imm.", "Sec.", "Ent.", "Ent.", "Clu+/Atoh1+ Stem", "End.", "Stem 2")

unique(my_vector)
Idents(KOIR_obj) <- KOIR_obj$seurat_clusters
KOIR_obj[["Cell_Type"]] <- Idents(KOIR_obj)
names(my_vector) <- levels(KOIR_obj)
KOIR_obj <- RenameIdents(KOIR_obj, my_vector)
KOIR_obj[["Cell_Type"]] <- Idents(KOIR_obj)
DimPlot(KOIR_obj)

my_levels <- c('Clu+/Atoh1+ Stem', 'Stem 2', 'Ent.', 'End.', 'Sec.', 'Imm.')
KOIR_obj$Cell_Type <- factor(KOIR_obj$Cell_Type, levels = my_levels)

KOIR_obj <- subset(KOIR_obj,  idents = c('Clu+/Atoh1+ Stem', 'Stem 2', 'Ent.', 'End.', 'Sec.'))
my_levels <- c('Clu+/Atoh1+ Stem', 'Stem 2', 'Ent.', 'End.', 'Sec.')
KOIR_obj$Cell_Type <- factor(KOIR_obj$Cell_Type, levels = my_levels)

KOIR_obj <- RunUMAP(KOIR_obj, dims = 1:25)
DimPlot(KOIR_obj, reduction = "umap", label = TRUE, group.by = 'Cell_Type')

#normalize and scale all
KOIR_obj@assays$SCT <- NULL
DefaultAssay(KOIR_obj) <- 'RNA'
KOIR_obj <- NormalizeData(KOIR_obj, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(KOIR_obj)

KOIR_obj[["RNA"]] <- as(object = KOIR_obj[["RNA"]], Class = "Assay")
use_python("/Users/vyom/miniconda3/bin/python")
KOIR_obj <- magic(KOIR_obj)

KOIR_obj <- ScaleData(KOIR_obj, features = all.genes)


DefaultAssay(KOIR_obj) <- 'RNA'
#differenetial expression analysis
KOIR_obj <- PrepSCTFindMarkers(KOIR_obj)
Idents(KOIR_obj) <- KOIR_obj$Treatment
DE_all <- FindMarkers(KOIR_obj, ident.1 = "KO", ident.2 = "WT", test.use = "MAST", logfc.threshold = .1, min.pct = .01, assay = 'RNA')
write.csv(DE_all, paste0('All_CREBKO_IR_scrna_DE.csv')) 
if (max(-log10(DE_all$p_val_adj)) < 321){
  ylim_num <- max(-log10(DE_all$p_val_adj))} else {
    ylim_num <- 320}
EnhancedVolcano(DE_all, lab = rownames(DE_all), x = 'avg_log2FC', y = 'p_val_adj', title = 'KO vs WT',
                pCutoff = .05, FCcutoff = 0.25, xlim = c(min(DE_all$avg_log2FC)-.05,max(DE_all$avg_log2FC)+.05),ylim = c(0,ylim_num+5),
                subtitle = 'ALL CELLS', gridlines.minor = FALSE, gridlines.major = FALSE)
ggsave(file = paste0('All_DE_peak_volcano_SCRNA.pdf'), width=6, height=6, units="in")

min(DE_all$p_val_adj)

Cell_Types <- c(levels(as.factor(KOIR_obj$Cell_Type)))
Cell_Types <- Cell_Types[2:5]
Idents(KOIR_obj) <- KOIR_obj$Cell_Type
for(i in Cell_Types){
  subset_cell <- subset(KOIR_obj,  idents = i)
  Idents(subset_cell) <- subset_cell$Treatment
  DE_subset <- FindMarkers(subset_cell, ident.1 = "KO", ident.2 = "WT",  test.use = "MAST",recorrect_umi=FALSE, logfc.threshold = .1, min.pct = .01, assay = 'RNA')
  write.csv(DE_subset, paste0(i,'_CREBKO_IR_scrna_DE.csv')) 
  if (max(-log10(DE_subset$p_val_adj)) < 321){
    ylim_num <- max(-log10(DE_subset$p_val_adj))} else {
      ylim_num <- 320}
  EnhancedVolcano(DE_subset, lab = rownames(DE_subset), x = 'avg_log2FC', y = 'p_val_adj', title = 'KO vs WT',
                  pCutoff = .05, FCcutoff = 0.1, xlim = c(min(DE_subset$avg_log2FC)-.05,max(DE_subset$avg_log2FC)+.05),ylim = c(0,ylim_num+5),
                  subtitle = i, gridlines.minor = FALSE, gridlines.major = FALSE)
  ggsave(file = paste0(i,'_DE_peak_volcano_SCRNA.pdf'), width=6, height=6, units="in")
}


Cell_Types <- c(levels(as.factor(KOIR_obj$Cell_Type)))
Cell_Types <- c('All', "clu_atoh1_stem",Cell_Types)
i = 'All'
for (i in Cell_Types) {
  res_format <- read.csv(paste0('~/',i, '_CREBKO_IR_scrna_DE.csv'))
  res_format <- as.data.frame(res_format[complete.cases(res_format),])
  dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023", "MSigDB_Hallmark_2020")
  enriched <- enrichr(unique(res_format[(res_format$p_val_adj < .05 & res_format$avg_log2FC < -0.01 ),]$X), dbs)
  Enriched1 <- rbind(enriched[[1]],enriched[[2]],enriched[[3]],enriched[[4]] )
  Enriched_filter <- Enriched1[Enriched1$Adjusted.P.value < .05,]
  plotEnrich(Enriched_filter, showTerms = 100, numChar = 100, y = "Count", orderBy = "Adjusted.P.value")
  write.csv(Enriched_filter, paste0(i,'_enrichR_Down_gsea.csv'))
  
  dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023", "MSigDB_Hallmark_2020")
  enriched <- enrichr(unique(res_format[(res_format$p_val_adj < .05 & res_format$avg_log2FC > 0.01 ),]$X), dbs)
  Enriched1 <- rbind(enriched[[1]],enriched[[2]],enriched[[3]],enriched[[4]] )
  Enriched_filter <- Enriched1[Enriched1$Adjusted.P.value < .05,]
  plotEnrich(Enriched_filter, showTerms = 100, numChar = 100, y = "Count", orderBy = "Adjusted.P.value")
  write.csv(Enriched_filter, paste0(i,'_enrichR_Up_gsea.csv'))
}

#pseudotime analysis
library(SeuratWrappers)
library(monocle3)
#install.packages("Signac")
library(Signac)
library(org.Mm.eg.db)
gene_symbol <- as.list(org.Mm.egSYMBOL)
Idents(KOIR_obj) <- KOIR_obj$Cell_Type

cds <- as.cell_data_set(KOIR_obj)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(KOIR_obj[["SCT"]])
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE,   close_loop = FALSE)
cds <- order_cells(cds, reduction_method = "UMAP")


plot_cells(cds = cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, cell_size = 2)
pseudo_umap <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = FALSE, label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, trajectory_graph_color = "#A8A8A8", graph_label_size=1.5,  cell_size = 2) 
pseudo_umap + scale_fill_gradientn(colors = annColors$Pseudotime)  + geom_polygon(data = hulls2, aes(x = x, y = y, group=CellType), fill=NA, color="black", size=0.3, alpha = 0.3) + geom_shadowtext(data = clusterMedian, aes(x = UMAP_1,y= UMAP_2, label = CellType, group = as.factor(CellType)), size = 3.5, bg.colour="black") + labs(fill = "Pseudotime")
ggsave(file="S100KO_obj_Pseudotime_Umap_no_trajectory.pdf", width=10, height=10)
plot_cells(cds,
           color_cells_by = "Cell_Type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

KOIR_obj <- AddMetaData(
  object = KOIR_obj,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)
KOIR_obj$Pseudotime
Idents(KOIR_obj) <- KOIR_obj$Cell_Type
FeaturePlot(KOIR_obj, c("Pseudotime"), pt.size = .0001, label = TRUE, repel = TRUE) +scale_fill_gradientn(colors = rev(brewer.pal(11,'Spectral'))) + scale_color_gradientn(colors =  rev(brewer.pal(11,'Spectral')))
ggsave(file = 'S100KO_obj_pseudotime_umap.pdf', width=4, height=4, units="in")

#housekeeping
#saveRDS(KOIR_obj, file = "/Users/vyom/CSHL Dropbox Team Dropbox/Vyom Shah/Seurat_Objects/CREBKO_IR_scrna_sep16.rds")
#KOIR_obj <- readRDS('/Users/vyom/CSHL Dropbox Team Dropbox/Vyom Shah/Seurat_Objects/CREBKO_IR_scrna_sep16.rds', refhook = NULL)


