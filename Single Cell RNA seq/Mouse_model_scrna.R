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
library(readxl)

AA_10x <- Read10X_h5(file = "/Users/vyom/data/SB10_SB11_results/Individual_Samples/SB10/Beyaz_10_10xgex_araasco/filtered_feature_bc_matrix.h5")

Control_10X <- Read10X_h5(file = "/Users/vyom/data/SB10_SB11_results/Individual_Samples/SB10/Beyaz_10_10xgex_control/filtered_feature_bc_matrix.h5")

AA_11x <- Read10X_h5(file = "/Users/vyom/data/SB10_SB11_results/Individual_Samples/SB11/Beyaz_11_10xgex_arasco/filtered_feature_bc_matrix.h5")

Control_11X <- Read10X_h5(file = "/Users/vyom/data/SB10_SB11_results/Individual_Samples/SB11/Beyaz_11_10xgex_control/filtered_feature_bc_matrix.h5")

AA_10x <- CreateSeuratObject(AA_10x, project = "AA_10x")
Control_10X <- CreateSeuratObject(Control_10X, project = "Control_10X")
AA_11x <- CreateSeuratObject(AA_11x, project = "AA_11x")
Control_11X <- CreateSeuratObject(Control_11X, project = "Control_11X")

d <- merge(Control_10X, y = c(AA_10x,Control_11X, AA_11x ), add.cell.ids = c("Control_10", "Arasco_10","Control_11", "Arasco_11"), project = "Arasco")
d
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "mt-")
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0 ,ncol = 3)
d <- subset(d, subset = nCount_RNA > 2000 & nCount_RNA < 40000 & nFeature_RNA > 1000 & nFeature_RNA < 9000 & percent.mt < 15)
d
d$orig.ident
Data.list <- SplitObject(d, split.by = "ident")
Data.list <- Data.list[c("Control_10X", "AA_10x","Control_11X", "AA_11x")]
for (i in 1:length(Data.list)) {
  Data.list[[i]] <- SCTransform(Data.list[[i]], verbose = FALSE, ncells=Inf)
}

# Normilization
#select highly variable genes 
Data.features <- SelectIntegrationFeatures(object.list = Data.list, nfeatures = 3000)
options (future.globals.maxSize = 4000 * 1024^5)
Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                verbose = FALSE)
Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                       anchor.features = Data.features, verbose = FALSE)

arasco_obj <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                                 verbose = TRUE)

# Visulization and Clustering
arasco_obj <- RunPCA(arasco_obj, verbose = FALSE)
arasco_obj <- RunUMAP(arasco_obj, dims = 1:10)

arasco_obj[["Type"]] <- Idents(arasco_obj)
new.cluster.ids <- c("Arasco", "Arasco", "Control", "Control")
names(new.cluster.ids) <- levels(arasco_obj)
arasco_obj <- RenameIdents(arasco_obj, new.cluster.ids)

arasco_obj[["Type"]] <- Idents(arasco_obj)
head(Idents(arasco_obj))
DimPlot(arasco_obj, group.by = c("Cell_Type"))

arasco_obj <- FindNeighbors(arasco_obj, dims = 1:10)
arasco_obj <- FindClusters(arasco_obj, resolution = .6)#.3
DimPlot(arasco_obj, reduction = "umap", label = TRUE)

#Annotation decisions
stem.markers <- c("Lgr5" ,"Ascl2" ,"Slc12a2" ,"Axin2" ,"Olfm4" ,"Gkn3")
TA.markers <- c("Tubb5" ,"Hmgb2" ,"Stmn1" ,"H2afz1" ,"Tuba1b" ,"Hmgb11" ,"Hmgn22" ,"Ptma1" ,"Kiaa0101" ,"Tk1" ,"Cenpw" ,"Tyms" ,"Ranbp11" ,"Idh21" ,"Ran1" ,"Dtymk" ,"Nucks11" ,"Dut1" ,"Cks1b")
paneth.markers <- c("Lyz1" ,"Defa17" ,"Defa22" ,"Defa24" ,"Ang4")
Enterocyte.markers <- c("Alpi" ,"Apoa1" ,"Apoa4" ,"Fabp1" ,"Adh6a")
EP.markers <- c("Ccnb1" ,"Cdc20" ,"Cenpa" ,"Cdkn3" ,"Ccnb2" ,"Cdc25c" ,"Kif22" ,"Ube2c" ,"Sapcd2" ,"Rbp7" ,"Aurka" ,"Ccna2" ,"Cdkn2d" ,"Kif23" ,"Nek2" ,"Slc16a1")
Enteroendocrine.markers <- c("Chga" ,"Chgb" ,"Tac1" ,"Tph1" ,"Neurog3" ,"Gch1" ,"Slc39a2")
Goblet.markers <- c("Muc2" ,"Clca3" ,"Tff3" ,"Agr2" ,"Spink4" ,"Fcgbp" ,"Zg16" ,"Ccl9" ,"Atoh1")
Tuft.markers <- c("Trpm5" ,"Gfi1b" ,"Il25" ,"Alox5ap" ,"Lrmp" ,"Rgs13" ,"Ltc4s" ,"Adh1")

Roulis_gene_lists_for_metagenes <- read_excel("Desktop/Roulis_gene_lists_for_metagenes.xlsx")
revival.markers <- c('S100a6','Ly6a','Clu','Anxa3','Areg')

arasco.obj <- AddModuleScore(object = arasco.obj, features = list(pathways$LGR5_stem_cell_signature), name = 'Stem')
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(revival.markers), name = 'revival.markers')
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(TA.markers), name = 'TA.markers')
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(paneth.markers), name = 'paneth.markers')
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(Enterocyte.markers), name = 'Enterocyte.markers')
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(EP.markers), name = 'EP.markers')
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(Enteroendocrine.markers), name = 'Enteroendocrine.markers')
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(Goblet.markers), name = 'Goblet.markers')
arasco.obj <- AddModuleScore(object = arasco.obj, features = list(Tuft.markers), name = 'Tuft.markers')

FeaturePlot(object = arasco.obj, features = 'stem1', pt.size = .001) + scale_color_viridis(option = 'B') + 
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'arasco_Stem.pdf', width=2, height=2, units="in")

FeaturePlot(object = arasco.obj, features = 'Goblet.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'arasco_Goblet.pdf', width=2, height=2, units="in")

FeaturePlot(object = arasco.obj, features = 'paneth.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'arasco_paneth.pdf', width=2, height=2, units="in")

FeaturePlot(object = arasco.obj, features = 'Enterocyte.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'arasco_Enterocyte.pdf', width=2, height=2, units="in")

FeaturePlot(object = arasco.obj, features = 'Enteroendocrine.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'arasco_Enteroendocrine.pdf', width=2, height=2, units="in")

FeaturePlot(object = arasco.obj, features = 'revival.markers1', pt.size = .001) + scale_color_viridis(option = 'B') +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.key.size = unit(.0, "cm"), legend.text=element_text(size=0), legend.title = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'arasco_RevStem.pdf', width=2, height=2, units="in")

DotPlot(arasco.obj, features = c('stem1','revival.markers1','reserve.markers1','TA.markers1','paneth.markers1','Enterocyte.markers1','EP.markers1','Enteroendocrine.markers1','Goblet.markers1','Tuft.markers1'))

Cell_Type_Sig_Score <- data.frame(Stem=arasco_obj$stem1, Revival=arasco_obj$revival.markers1, Reserve=arasco_obj$reserve.markers1, Transit_Amplifying =arasco_obj$TA.markers1, Paneth=arasco_obj$paneth.markers1, Enterocyte=arasco_obj$Enterocyte.markers1, EP=arasco_obj$EP.markers1,Enteroendocrine= arasco_obj$Enteroendocrine.markers1,Goblet=arasco_obj$Goblet.markers1,Tuft=arasco_obj$Tuft.markers1)

# module score distribution
modulescores <- Cell_Type_Sig_Score %>%
  rownames_to_column(var="id") %>%
  pivot_longer(-id, names_to="celltype", values_to="score")


p <- ggplot(modulescores)
#p <- ggplot(onescore)
p + geom_point(aes(x=fct_inorder(id), y=sort(score))) +
  facet_wrap(~celltype) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

Idents(arasco_obj) <- arasco_obj$Type

z_Stem_cluster<-WhichCells(object = arasco_obj, expression = stem1 > 1)
z_RevSC_cluster<-WhichCells(object = arasco_obj,expression = revival.markers1 > 1, )
z_reserve_cluster<-WhichCells(object = arasco_obj, expression = reserve.markers1 > 1)
z_TA_cluster<-WhichCells(object = arasco_obj, expression = TA.markers1 > .5)
z_paneth_cluster<-WhichCells(object = arasco_obj, expression = paneth.markers1 > 5)
z_Enterocyte_cluster<-WhichCells(object = arasco_obj, expression = Enterocyte.markers1 > 1)
z_EP_cluster<-WhichCells(object = arasco_obj, expression = EP.markers1 > .5)
z_Enteroendocrine_cluster<-WhichCells(object = arasco_obj, expression = Enteroendocrine.markers1 > 11)
z_Goblet_cluster<-WhichCells(object = arasco_obj, expression = Goblet.markers1 > 1)
z_Tuft_cluster<-WhichCells(object = arasco_obj, expression = Tuft.markers1 > 1)

table(arasco_obj$Cell_type)

arasco_obj <- FindClusters(arasco_obj, resolution = .6)#.3
DimPlot(arasco_obj, reduction = "umap", label = TRUE)

DimPlot(arasco_obj, label=T,group.by = 'Cell_type', cells.highlight= list( z_Stem_cluster, z_RevSC_cluster,  z_TA_cluster, z_paneth_cluster, z_Enterocyte_cluster, z_EP_cluster, z_Enteroendocrine_cluster, z_Goblet_cluster, z_Tuft_cluster),  cols.highlight = c("darkblue", "darkred", 'green','lightgreen', 'darkgreen', 'yellow', 'pink', 'purple', 'orange','lightblue'),cols= "grey")
DimPlot(arasco_obj, label=T,group.by = 'Cell_type', cells.highlight= list(z_paneth_cluster),  cols.highlight = c("darkblue", "darkred"),cols= "grey")
VlnPlot(arasco_obj, features = c('paneth.markers1','Goblet.markers1'), pt.size = 0)

z_RevSC_clusterx<-WhichCells(object = arasco_obj ,expression = revival.markers1 > 1)
x_rev <- WhichCells(object = arasco_obj,idents = c('x'), z_RevSC_clusterx )
write.csv(x_rev)
#______________________________Clustering________________________________________#
new.cluster.ids <- c('Enterocyte (Proximal)', 'Stem 1', 'Goblet', 'Enterocyte Progenitor', 'Enterocyte (Proximal)', 'Enterocyte (Proximal)', 'Stem 1', 'Transit Amplifying','Enterocyte Progenitor','Transit Amplifying', 'Stem 1', 'Enterocyte (Proximal)', 'Enterocyte (Distal)', 'Enteroendocrine', 'Enterocyte (Distal)', 'Enterocyte Progenitor', 'Tuft', 'CD8+ intraepithelial cells', 'Macrophage', 'Goblet', 'Enterocyte (Proximal)', 'x', 'Paneth')
arasco_obj[["Cell_type"]] <- Idents(arasco_obj)
names(new.cluster.ids) <- levels(arasco_obj)
arasco_obj <- RenameIdents(arasco_obj, new.cluster.ids)
arasco_obj[["Cell_type"]] <- Idents(arasco_obj)
DimPlot(arasco_obj, reduction = "umap", group.by= 'Cell_type')

Idents(arasco_obj) <- arasco_obj$Cell_type
arasco_obj$Cell_Type <- as.character(Idents(arasco_obj))
stem2_cells1 <- WhichCells(object = arasco_obj,idents = c('x'), z_RevSC_cluster )
stem2_cells <- c(stem2_cells1, z_Enteroendocrine_cluster)
arasco_obj$Cell_Type[stem2_cells] <- paste('Stem 2')
DimPlot(arasco_obj, reduction = "umap", group.by= 'Cell_Type')
arasco_obj$Cell_type


levels(factor(arasco_obj@meta.data$Type))
arasco_obj@meta.data$type <- factor(arasco_obj@meta.data$Type, levels = c("Control", "Arasco")) 
Idents(arasco_obj) <- arasco_obj$Cell_Type
prop.table(x = table(Idents(arasco_obj), arasco_obj$type), margin = 2)
tail(arasco_obj$type)
arasco.obj <- subset(arasco_obj, idents = c('Stem 1','Stem 2', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft'))
DimPlot(arasco.obj)
my_levels <- c('Stem 1','Stem 2', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft')
arasco.obj$Cell_type <- factor(x = arasco.obj$Cell_Type, levels = my_levels)
Idents(arasco.obj) <- arasco.obj$Cell_type
Proportion_Arasco<- prop.table(x = table(arasco.obj$Cell_type, arasco.obj$type), margin = 2)
DimPlot(arasco.obj, group.by = "Cell_type", label = FALSE, pt.size=1, label.size = 3)
arasco.obj@meta.data$type <- factor(arasco.obj@meta.data$Type, levels = c("Control", "Arasco")) 

# imputation for violin plots
DefaultAssay(arasco.obj) <- "RNA"
arasco.obj <- NormalizeData(object = arasco.obj,normalization.method = "LogNormalize", assay = "RNA")
arasco.obj <- magic(arasco.obj)

# Pseudotime Analysis
library(SeuratWrappers)
library(monocle3)
library(Signac)
library(org.Mm.eg.db)
gene_symbol <- as.list(org.Mm.egSYMBOL)
Idents(arasco.obj) <- arasco.obj$Cell_type
pseudo_stem <- WhichCells(object = arasco.obj, idents = 'Stem 1')

arasco.obj1 <- arasco.obj
cds <- as.cell_data_set(arasco.obj1)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(arasco.obj1[["RNA"]])
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE)

cds <- order_cells(cds, reduction_method = "UMAP")

cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)

#visualize Seurat object using the pseudotime variable
arasco.obj <- AddMetaData(
  object = arasco.obj,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)