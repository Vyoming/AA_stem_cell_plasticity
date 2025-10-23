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

#Load in the data
ARD3 <- Read10X_h5("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/YAPKO_AA_ovy_scrna/Beyaz_OE07_6812_YAP_WT_ARD/filtered_feature_bc_matrix.h5")
ARD1 <- Read10X_h5("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/ARD_YapKO/RNA/count/Beyaz_OE05_6800_ARD_WT/outs/filtered_feature_bc_matrix.h5")
ARD2 <- Read10X(data.dir ="/Users/vyom/data/ARD_KO_updated/AG01_1759_arasco/outs/filtered_feature_bc_matrix")
ARDKO1 <- Read10X_h5("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/YAPKO_AA_ovy_scrna/Beyaz_OE07_6813_YAP_KO_ARD/filtered_feature_bc_matrix.h5")
ARDKO2 <- Read10X_h5("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/YAPKO_AA_ovy_scrna/Beyaz_OE07_6406_YAP_KO_ARD/filtered_feature_bc_matrix.h5")
ARDKO3 <- Read10X_h5("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/ARD_YapKO/RNA/count/Beyaz_OE05_6801_ARD_YAP_KO/outs/filtered_feature_bc_matrix.h5")

ARD1 <- CreateSeuratObject(ARD1, project = "ARD1")
ARD2 <- CreateSeuratObject(ARD2, project = "ARD2")
ARD3 <- CreateSeuratObject(ARD3, project = "ARD3")

ARDKO1 <- CreateSeuratObject(ARDKO1, project = "ARDKO1")
ARDKO2 <- CreateSeuratObject(ARDKO2, project = "ARDKO2")
ARDKO3 <- CreateSeuratObject(ARDKO3, project = "ARDKO3")

#filtering and QC
ARD1 <- subset(ARD1, subset = nCount_RNA > 250 & nCount_RNA < 80000 & nFeature_RNA > 2000 & nFeature_RNA < 80000)
ARD2 <- subset(ARD2, subset = nCount_RNA > 250 & nCount_RNA < 80000 & nFeature_RNA > 2000 & nFeature_RNA < 80000)
ARD3 <- subset(ARD3, subset = nCount_RNA > 250 & nCount_RNA < 80000 & nFeature_RNA > 2000 & nFeature_RNA < 80000)
ARDKO1 <- subset(ARDKO1, subset = nCount_RNA > 250 & nCount_RNA < 80000 & nFeature_RNA > 2000 & nFeature_RNA < 80000)
ARDKO2 <- subset(ARDKO2, subset = nCount_RNA > 250 & nCount_RNA < 80000 & nFeature_RNA > 3250 & nFeature_RNA < 80000)
ARDKO3 <- subset(ARDKO3, subset = nCount_RNA > 250 & nCount_RNA < 80000 & nFeature_RNA > 1000 & nFeature_RNA < 80000)

d <- merge(ARD1, y = c(ARD2, ARD3, ARDKO1, ARDKO2, ARDKO3), add.cell.ids = c('ARD1', 'ARD2', 'ARD3', 'ARDKO1', 'ARDKO2', 'ARDKO3'), project = "ARD_YapKO")
d
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "mt-")
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0 ,ncol = 3)
d <- subset(d, subset =  percent.mt < 15)
d
plot1 <- FeatureScatter(d, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(d, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = 'orig.ident',pt.size = 0 ,ncol = 3)


d$orig.ident
Data.list <- SplitObject(d, split.by = "ident")
Data.list <- Data.list[c('ARD1', 'ARD2', 'ARD3', 'ARDKO1', 'ARDKO2', 'ARDKO3')]
for (i in 1:length(Data.list)) {
  
  Data.list[[i]] <- SCTransform(Data.list[[i]], verbose = FALSE)
}

# Normilization
#select highly variable genes 
Data.features <- SelectIntegrationFeatures(object.list = Data.list, nfeatures = 2500)
options (future.globals.maxSize = 4000 * 1024^20)
Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                verbose = FALSE)
Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                       anchor.features = Data.features, verbose = FALSE)
rm(Data.list, d)
AYKO_obj <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                          verbose = TRUE)
rm(Data.anchors)
# Visulization and Clustering
AYKO_obj <- RunPCA(AYKO_obj)
VizDimLoadings(AYKO_obj, dims = 1:12, reduction = "pca")

#determine num of PCs for clustering
DimPlot(AYKO_obj, reduction = "pca")
ElbowPlot(AYKO_obj, ndims = 40, reduction = "pca")

AYKO_obj <- RunUMAP(AYKO_obj, dims = 1:40,n.neighbors = 25, n.epochs = 250)

AYKO_obj <- FindNeighbors(AYKO_obj, dims = 1:40,  n.trees = 200)

AYKO_obj <- FindClusters(AYKO_obj, resolution = 1.25)
DimPlot(AYKO_obj, reduction = "umap", label = TRUE)

markers <- FindAllMarkers(AYKO_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = .5)
library(readxl)
write.csv(markers,'Gene_de_Sigs_Per_Clust.csv')
DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3",'S100a6','Ly6a','Clu','Anxa3','Areg',"Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a",'Tmigd1', 'Fabp6', 'Slc51b', 'Slc51a', "Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Lyz1","Defa17","Defa24","Ang4","Pou2f3","Avil","Tuba1a","Adh1",'Epcam','Krt19','Ptprc','Cd3e','Cd8a','Cd4','Cd68','Ccl5','Itgam','Itgax','Ighm') 

DotPlot(AYKO_obj, features = DotPlot_Sig, assay = 'SCT') + labs(y= "Cell Treatment", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1)) +
  theme(text = element_text(size=5), axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())

new.cluster.ids <- c("EP", "TA", "TA", "Gob.", "Ent.", "Ent.", "Stem 1", "Stem 1", "Gob.", "Ent.", "Ent.", "TA", "Stem 1", "Ent.", "Ent.", "Ent.", "Stem 1", "Stem 1", "Tuft", "Ent.", "End.", "Ent.", "T", "EP", "Gob.", "Gob.", "End.", "T", "SP", "Ent.", "SP", "Stem 2", "Stem 1", "End.", "Myeloid", "Gob.", "Gob.")

AYKO_obj[["Cell_Type"]] <- Idents(AYKO_obj)
names(new.cluster.ids) <- levels(AYKO_obj)
AYKO_obj <- RenameIdents(AYKO_obj, new.cluster.ids)
AYKO_obj[["Cell_Type"]] <- Idents(AYKO_obj)
DimPlot(AYKO_obj, reduction = "umap", group.by= 'Cell_Type')

#magic imputation
use_python("/Users/vyom/miniconda3/bin/python")
AYKO_obj <- magic(AYKO_obj, assay = 'SCT')

#differential expression
AYKO_obj <- PrepSCTFindMarkers(AYKO_obj)
Idents(AYKO_obj) <- AYKO_obj$Treatment
DE_all <- FindMarkers(AYKO_obj, ident.1 = "ARD", ident.2 = "ARDKO", test.use = "MAST", logfc.threshold = .01, min.pct = .01, assay = 'SCT')
DE_all<- DE_all[!grepl("Rik", rownames(DE_all)), ]
DE_all<- DE_all[!grepl("mt-", rownames(DE_all)), ]
DE_all<- DE_all[!grepl("Rpl", rownames(DE_all)), ]
DE_all<- DE_all[!grepl("Rps", rownames(DE_all)), ]
DE_all<- DE_all[!grepl("Atp", rownames(DE_all)), ]
DE_all<- DE_all[!grepl("AY0", rownames(DE_all)), ]
DE_all<- DE_all[!grepl("Gm", rownames(DE_all)), ]
write.csv(DE_all, paste0('All_ARDvYAPKO_scrna_DE.csv')) 

#pseudotime:
#pseudotime analysis
library(SeuratWrappers)
library(monocle3)
#install.packages("Signac")
library(Signac)
library(org.Mm.eg.db)
gene_symbol <- as.list(org.Mm.egSYMBOL)
Idents(AYKO_obj) <- AYKO_obj$Cell_Type

cds <- as.cell_data_set(AYKO_obj)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(AYKO_obj[["SCT"]])
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

AYKO_obj <- AddMetaData(
  object = AYKO_obj,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)
AYKO_obj$Pseudotime
Idents(AYKO_obj) <- AYKO_obj$Cell_Type
FeaturePlot(AYKO_obj, c("Pseudotime"), pt.size = .0001, label = TRUE, repel = TRUE) +scale_fill_gradientn(colors = rev(brewer.pal(11,'Spectral'))) + scale_color_gradientn(colors =  rev(brewer.pal(11,'Spectral')))
ggsave(file = 'ARDKO_obj_pseudotime_umap.pdf', width=4, height=4, units="in")

#create scores for each value in signals across peaks and reads 
Fetal <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "Fetal_Spheroid")$Gene_Name
Granuloma <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "Granuloma_Induced")$Gene_Name
Radiation <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "Radiation_induced")$Gene_Name
ECM_Induced <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "ECM_Induced")$Gene_Name
ASCL_Targets <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "ASCL_Targets")$Gene_Name
Regeneration <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "Regeneration_Induced")$Gene_Name
Homeostatic <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "Homeostatic_Signature")$Gene_Name
AA_induced_invivo <- c('Hspa1b',	'Gm2000',	'Defa24',	'Banp',	'Lockd',	'Dpm3',	'Rps27rt',	'Snrpg',	'mt-Nd3',	'Cenpw',	'Slirp',	'Pet100',	'Tomm7',	'Tmem256',	'Cops9',	'Atp5k',	'Snhg6',	'Smim4',	'Mrpl52',	'Ndufb1-ps',	'Epb41l4aos',	'Rpl38',	'Pin4',	'Snhg8',	'Snrpf',	'Ndufa2',	'Psmg4',	'Atp5mpl',	'Romo1',	'Rps29',	'Rps21',	'Rps28',	'Anapc13',	'Atp5md',	'Lsm7',	'Tmem258',	'Uqcc2',	'Plp2',	'S100a6',	'Phgdh',	'Ndufa3',	'Eif3j1',	'Nop10',	'Sec61g',	'Aspm',	'Atox1',	'Polr2l',	'Ndufb2',	'Ccdc167',	'Fkbp5',	'Bola2',	'Ndufa5',	'Ndufc1',	'Mrps21',	'Snrpe',	'Uqcr11',	'Hist1h2ae',	'Sp3os',	'Ubl5',	'Rpl39',	'Tuba1c',	'Rpl37',	'Prelp',	'Mrpl33',	'Snrpd2',	'Naa38',	'Polr2k',	'Uqcr10',	'Cox17',	'Smim22',	'Rplp2',	'Prc1',	'Rpl37a',	'Rpl31',	'Sgo1',	'Rpl21',	'Fmc1',	'Ifitm3',	'Ppih',	'Nusap1',	'Trappc10',	'Zfas1',	'Rpl36',	'Ncapg',	'Polr2i',	'Tmsb15b2',	'Cox20',	'Ndufv3',	'Gm11808',	'Atp5e',	'Arpp19',	'Adra2a',	'Fxyd3',	'Mt1',	'Esco2',	'2410006H16Rik',	'Rpl35a',	'Ndufa1',	'Elob',	'Hist1h1b',	'Ndufs6',	'2310009B15Rik',	'Cit',	'Rpl41',	'Hspe1-rs1',	'S100a11',	'Ndufa7',	'Ndufb4',	'S100g',	'Srek1ip1',	'Sem1',	'Hmgcs2',	'Mif',	'Aurka',	'Diaph3',	'Cd44',	'Knl1',	'Snhg3',	'Top2a',	'Hells',	'Rrp8',	'Plk1',	'Ddx18',	'Cep164',	'Kif20b',	'Cdca2',	'Pmepa1',	'Jaml',	'Incenp',	'Smim26',	'Trip13',	'Mis18bp1',	'Rbmx2',	'Lsm6',	'Spc25',	'Rrp1b',	'Rps17',	'Lsm5',	'S100a13',	'Rpgrip1',	'Cdca8',	'Gas5',	'Ddx24',	'Hist1h1e',	'Cenpf',	'Hmmr',	'Atp5j2',	'Acot1',	'Kif22',	'Neil3',	'Rpl34',	'Mis12',	'Cox7a1',	'Tuba1b',	'Cbx6',	'Tubb2b',	'Ccdc34',	'Rpl36a',	'Rbmxl1',	'Ndc80',	'2310009A05Rik',	'Cep295',	'Bub1',	'Spink4',	'Hirip3',	'H2afj',	'Dlgap5',	'Atp5l',	'Smim11',	'Cox6c',	'Dnajb1',	'Kif11',	'Ubap2',	'Tpx2',	'Cebpz',	'Ctc1',	'Crip1',	'Sycn',	'Telo2',	'Rbis',	'Rassf4',	'Shcbp1',	'Cox7c',	'Gm43813',	'Rgcc',	'Notch1',	'Serf1',	'Ndufa13',	'Aurkb',	'Pbk',	'Sec61b',	'Lrrc31',	'Chordc1',	'Rps27l',	'1810037I17Rik',	'Rpp21',	'2200002D01Rik',	'E030030I06Rik',	'Fbxo5',	'Rad51ap1',	'Uqcrq',	'Polr2f',	'Rpl35',	'Ttr',	'Nup210',	'Krtcap2',	'Mcph1',	'Rps27',	'Ndufaf8',	'Lig1',	'Cd3eap',	'Ppwd1',	'1700097N02Rik',	'2810408I11Rik',	'Ckap2',	'2810004N23Rik',	'Rps2',	'Prpf4',	'Snhg18',	'Kif4',	'Sgo2a',	'Higd1a',	'Ndufs5',	'Xpo5',	'Soat1',	'Rps15',	'Anln',	'Rps26',	'Grwd1',	'Lrrc14',	'Cox7a2',	'Ncapd2',	'Tcof1',	'Rad54l',	'Ier3ip1',	'Cenpe',	'Apobec3',	'Lyrm2',	'Bub1b',	'Ccdc85c',	'Cep192',	'Prpf6',	'Cox7b',	'Rps25',	'Pinx1',	'Tex10',	'Gpatch4',	'Gnl3')          
AA_induced_stem2 <- c('Plekhm3',	'H2-Eb1',	'Cd74',	'H2-Ab1',	'Nucb2',	'Map7d1',	'Sik1',	'Upb1',	'Arid3a',	'Cyb561a3',	'Ankrd44',	'Ctsl',	'Tbc1d8',	'Slc29a3',	'Orai3',	'Dusp1',	'Ifnar2',	'Tcf4',	'Irf7',	'Unc93b1',	'Dusp6',	'Ptpn6',	'Tubgcp5',	'S100a6',	'Ptpre',	'Rell1',	'Tom1',	'Crlf2',	'Snx29',	'Gm37529',	'Gsn',	'Snx18',	'Gmfg',	'Mob3a',	'Mfsd12',	'Unc119',	'Uvrag',	'Cybc1',	'Cmtm7',	'Cnp',	'Npc1',	'Map4k2',	'Stat2',	'H2-Q4',	'Gm43329',	'Tmem123',	'Cdkn2d',	'Nek7',	'Slc15a4',	'Gns',	'Anxa5',	'Bcl11a',	'Ppp1r15a',	'Sema4b',	'Itpr1',	'Aldh3b1',	'Igtp',	'Gng10',	'Dock8',	'Fam174a',	'Pmepa1',	'Rogdi',	'Bbc3',	'Rab33b',	'Ftx',	'Pde7a',	'Lyn',	'Dclre1c',	'Atp13a2',	'Ccl5',	'Cbx4',	'Tgfbr1',	'Gpcpd1',	'Pqlc3',	'Rpgrip1',	'Arid3b',	'Ablim1',	'Gm31718',	'Dipk1a',	'Lgmn',	'Grina',	'Rhobtb2',	'Snx30',	'Nudt18',	'Trim12a',	'Slc49a4',	'Myo9b',	'Sgk3',	'Tor3a',	'Bicd2',	'Rwdd2a',	'Dtx2',	'Ctss',	'Slc25a4',	'5031439G07Rik',	'Grn',	'Tex2',	'Vamp1',	'A430035B10Rik',	'Inpp5k',	'Ccnd3',	'Irf8',	'Stat5a',	'Rab2b',	'Tap1',	'Pacc1',	'Ifnar1',	'Ier5',	'Pml',	'Acox3',	'Arl10',	'Il4ra',	'Tap2',	'Apobec3',	'Litaf',	'Malt1',	'Il17ra',	'Bmyc',	'Arhgap17',	'Zfp319',	'Bnip2',	'Pkd1',	'Tpst2',	'H2-T22',	'Calcoco1',	'Pip5k1c',	'Gmip',	'Crtc3',	'Ptpn1',	'L3mbtl3',	'Tcirg1',	'Prcp',	'Stat1',	'Cyth1',	'Nabp1',	'Nek6',	'Bin1',	'Arhgap27',	'Wdr81',	'Ptprj')

All_Genes <- ARDKO_atac_obj@assays$RNA@data@Dimnames[[1]]
AA_induced_invivo <- intersect(All_Genes, AA_induced_invivo)
AA_induced_stem2 <- intersect(All_Genes, AA_induced_stem2)
Fetal <- intersect(All_Genes, Fetal)
Granuloma <- intersect(All_Genes, Granuloma)
Radiation <- intersect(All_Genes, Radiation)
Regeneration <- intersect(All_Genes, Regeneration)

mean.exp <- zscore(colMeans(x = ARDKO_atac_obj@assays$RNA@data[AA_induced_stem2, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = ARDKO_atac_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  ARDKO_atac_obj@meta.data$Stem2_Markers <- mean.exp
}

mean.exp <- zscore(colMeans(x = ARDKO_atac_obj@assays$RNA@data[AA_induced_invivo, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = ARDKO_atac_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  ARDKO_atac_obj@meta.data$AA_induced_invivo <- mean.exp
}
score.list <- c('AA_induced_invivo', 'Stem2_Markers')
DefaultAssay(ARDKO_atac_obj) <- 'RNA'
for(i in score.list) {
  selected_cells <- names(ARDKO_atac_obj$Cell_Type)
  vln_data <- FetchData(ARDKO_atac_obj,
                        vars = c(i,"Treatment", "Cell_Type"),
                        cells = selected_cells,
                        layer = "data")
  vln_data <- melt(vln_data)
  
  All <- VlnPlot(ARDKO_atac_obj, split.by = "Treatment",group.by = 'Cell_Type', features = i, pt.size = 0, assay = "RNA",  cols = c('#d95f02' ,'#01949A'), log = FALSE, split.plot = TRUE) + 
    theme(legend.position = 'none') + 
    geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd= .2) + xlab('') + 
    theme(text = element_text(size=9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size =9)) +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Treatment), method = "wilcox.test",  label = "p.signif", size = 4, hide.ns = TRUE) #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  All$layers[[1]]$aes_params$size = .15
  All
  ggsave(file = paste0('ARDKO_ATAC_percelltype_', i, '.pdf'), plot=All, width=3.5, height=3.5, units="in")
}

#housekeeping
#saveRDS(AYKO_obj, file = "/Users/vyom/CSHL Dropbox Team Dropbox/Vyom Shah/Seurat_Objects/jun20_AA_YAP_ko_new_n2_scRNA.rds")
#AYKO_obj <- readRDS('/Users/vyom/CSHL Dropbox Team Dropbox/Vyom Shah/Seurat_Objects/jun20_AA_YAP_ko_new_n2_scRNA.rds', refhook = NULL)





