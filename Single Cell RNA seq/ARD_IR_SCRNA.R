#ARD + IR scRNA analysis
#Spring 2024
#Vyom Shah
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
IR1 <- Read10X("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/ARD_IR/Beyaz_OE04_SB0027_Control_IR/outs/filtered_feature_bc_matrix")
IR2 <- Read10X("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/ARD_IR/Beyaz_OE04_SB0028_Control_IR/outs/filtered_feature_bc_matrix")

ARD1 <- Read10X( "/Users/vyom/cloud_computing/data_norepl/sequencing_backup/ARD_IR/Beyaz_OE04_SB0042_ARD_IR/outs/filtered_feature_bc_matrix")
ARD2 <- Read10X( "/Users/vyom/cloud_computing/data_norepl/sequencing_backup/ARD_IR/Beyaz_OE04_SB0043_ARD_IR/outs/filtered_feature_bc_matrix")


IR1 <- CreateSeuratObject(IR1, project = "IR1")
IR2 <- CreateSeuratObject(IR2, project = "IR2")

ARD1 <- CreateSeuratObject(ARD1, project = "ARD1")
ARD2 <- CreateSeuratObject(ARD2, project = "ARD2")

d <- merge(IR1, y = c(IR2, ARD1, ARD2), add.cell.ids = c("IR1", "IR2","ARD1", "ARD2"), project = "ARD_intest")

#QC and filtering
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "mt-")
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0 ,ncol = 3)
d <- subset(d, subset = nCount_RNA > 100 & nCount_RNA < 7500 & nFeature_RNA > 100 & nFeature_RNA < 3000 & percent.mt < 10)
d
plot1 <- FeatureScatter(d, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(d, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Normalization
Data.list <- SplitObject(d, split.by = "ident")
Data.list <- Data.list[c("IR1", "IR2","ARD1", "ARD2")]
for (i in 1:length(Data.list)) {
  
  Data.list[[i]] <- SCTransform(Data.list[[i]], verbose = FALSE)
}

#select highly variable genes 
Data.features <- SelectIntegrationFeatures(object.list = Data.list, nfeatures = 3000)
Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                verbose = FALSE)
Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                       anchor.features = Data.features, verbose = FALSE)
ARD_obj <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                          verbose = TRUE)

# Visulization and Clustering
ARD_obj <- RunPCA(ARD_obj)
VizDimLoadings(ARD_obj, dims = 1:2, reduction = "pca")

DimPlot(ARD_obj, reduction = "pca")
#determine number of PCs to use for clustering
ElbowPlot(ARD_obj, ndims = 50, reduction = "pca")

ARD_obj <- RunUMAP(ARD_obj, dims = 1:30, n.epochs = 500)

ARD_obj <- FindNeighbors(ARD_obj, dims = 1:30, n.trees = 200)

ARD_obj <- FindClusters(ARD_obj, resolution = 2)
DimPlot(ARD_obj, reduction = "umap", label = TRUE)
DimPlot(ARD_obj, reduction = "umap", label = TRUE, group.by = 'Treatment')

new.cluster.ids <- c("Immune", "Immune", "Immune", "Immune", "Immune", "Immune", "Immune", "Immune", "Immune", "Immune", "Immune", "Immune", "Immune", "Immune", "Ent.", "Ent.", "Immune", "Immune", "Immune", "Immune", "Ent.", "Ent.", "Immune", "Ent.", "Immune", "Immune", "Immune", "Immune", "Immune", "Immune", "RevSC", "Ent.", "ArdSC", "Immune", "Entend.", "Sec.")

ARD_obj[["CellType"]] <- Idents(ARD_obj)
names(new.cluster.ids) <- levels(ARD_obj)
ARD_obj <- RenameIdents(ARD_obj, new.cluster.ids)
ARD_obj[["CellType"]] <- Idents(ARD_obj)
DimPlot(ARD_obj, reduction = "umap", group.by= 'CellType')

#continue downstream analysis on epithelial lineages only
Idents(ARD_obj) <- ARD_obj$CellType
IR_imm_obj <- subset(ARD_obj, idents = 'Immune')
ARDIR_obj <- subset(ARD_obj, idents = c('Ent.', 'RevSC', 'ArdSC', 'Entend.', 'Sec.'))


#MAGIC imputation:
#normalize and scale all
use_python("/Users/vyom/miniconda3/bin/python")
DefaultAssay(ARDIR_obj) <- "RNA"
ARDIR_obj <- NormalizeData(ARDIR_obj, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(ARDIR_obj)
ARDIR_obj <- ScaleData(ARDIR_obj, features = all.genes)
ARDIR_obj[["RNA"]] <- as(object = ARDIR_obj[["RNA"]], Class = "Assay")
ARDIR_obj <- magic(ARDIR_obj)


#Differential Expression
ARDIR_obj <- PrepSCTFindMarkers(ARDIR_obj)
Idents(ARDIR_obj) <- ARDIR_obj$Cell_Type
DE_all <- FindMarkers(ARDIR_obj, ident.1 = "Stem 2", test.use = "MAST", logfc.threshold = 1, min.pct = .25, assay = 'SCT')
DE_all$genes <- rownames(DE_all)

write.csv(DE_all, paste0('All_ARD_IR_scrna_DE.csv')) 

if (max(-log10(DE_all$p_val_adj)) < 321){
  ylim_num <- max(-log(DE_all$p_val_adj))} else {
    ylim_num <- 320}
EnhancedVolcano(DE_all, lab = rownames(DE_all), x = 'avg_log2FC', y = 'p_val_adj', title = 'ARD_IR vs IR',
                pCutoff = .05, FCcutoff = 0.25, xlim = c(min(DE_all$avg_log2FC)-.05,max(DE_all$avg_log2FC)+.05),ylim = c(0,ylim_num+5),
                subtitle = 'ALL CELLS', gridlines.minor = FALSE, gridlines.major = FALSE)
ggsave(file = paste0('All_DE_peak_volcano_SCRNA.pdf'), width=6, height=6, units="in")


Cell_Types <- c(levels(as.factor(ARDIR_obj$Treatment)))
Idents(ARDIR_obj) <- ARDIR_obj$Treatment
for(i in Cell_Types){
  subset_cell <- subset(ARDIR_obj,  idents = i)
  subset_cell <- PrepSCTFindMarkers(subset_cell, verbose = TRUE)
  Idents(subset_cell) <- subset_cell$Cell_Type
  DE_subset <- FindMarkers(subset_cell, ident.1 = "ArdSC", ident.2 = "RevSC",  test.use = "MAST",recorrect_umi=FALSE, logfc.threshold = .1, min.pct = .01, assay = 'SCT')
  write.csv(DE_subset, paste0(i,'_ARD_IR_scrna_DE.csv')) 
  if (max(-log10(DE_subset$p_val_adj)) < 321){
    ylim_num <- max(-log10(DE_subset$p_val_adj))} else {
      ylim_num <- 320}
  EnhancedVolcano(DE_subset, lab = rownames(DE_subset), x = 'avg_log2FC', y = 'p_val_adj', title = 'ARD_IR vs IR',
                  pCutoff = .05, FCcutoff = 0.1, xlim = c(min(DE_subset$avg_log2FC)-.05,max(DE_subset$avg_log2FC)+.05),ylim = c(0,ylim_num+5),
                  subtitle = i, gridlines.minor = FALSE, gridlines.major = FALSE)
  ggsave(file = paste0(i,'_DE_peak_volcano_SCRNA.pdf'), width=6, height=6, units="in")
}

#score:
AA_induced_organoid_scrna <- c('Fabp1',	'S100a6',	'Zg16',	'Dmbt1',	'Lgals3',	'Agr2',	'Reg1',	'Akr1b8',	'Anxa2',	'Cbr3',	'Ldha',	'Hmgcs2',	'2210407C18Rik',	'Gsta1',	'Krt19',	'Pkm',	'Acadl',	'AA467197',	'Apoa4',	'Pgam1',	'Gm20594',	'Sis',	'Apoa1',	'Reg4',	'Krt8',	'Slc11a2',	'Gna11',	'Anxa13',	'Gsta4',	'Guca2b',	'Prdx6',	'Serpinb6a',	'Nme2',	'Tubb4b',	'Krt18',	'Lypd8',	'Rbp2',	'Mttp',	'Nfkbia',	'Cldn3',	'Rpl15',	'Dstn',	'Tuba1b',	'Ckmt1',	'Gpx2',	'Krt7',	'Prss32',	'Ube2m',	'Alpi',	'Gsta3',	'Wfdc2',	'Slc16a3',	'Anxa3',	'Capg',	'Aldob',	'Cbr1',	'Gstm1',	'Arpc4',	'Ddah1',	'Egln3',	'Tomm22',	'Tfrc',	'Mdh2',	'Actg1',	'Eif5a',	'Akr1c19',	'Pgk1',	'Ptgr1',	'Cmas',	'Aldh1a1',	'Ezr',	'Cct7',	'Aprt',	'Tm4sf4',	'2610528A11Rik',	'Eif6',	'Ankrd37',	'Tuba4a',	'Ddx3x',	'Acaa2',	'Galk1',	'Tm4sf5',	'Lars2',	'Nfe2l2',	'Oit1',	'Aldoa',	'Ddx39',	'Plpp2',	'Arhgdig',	'Vdac2',	'Ctsl',	'Gclm',	'Tpi1',	'Prdx1',	'Gsr',	'Gpd1',	'Tspan8',	'Hspd1',	'Lgals4',	'Phb',	'Ctsz',	'Actb',	'Me1',	'St3gal4',	'Tmem54',	'Gm10116',	'Acp5',	'Cldn15',	'Ero1l',	'Akr1c13',	'Hras',	'Sri',	'Tmem237',	'Eif4a1',	'Id1',	'Areg',	'Arf1',	'Pdcd6',	'Tomm20',	'Pgp',	'Rpl7l1',	'Pfkp',	'Fabp2',	'Uqcrfs1',	'Sfn',	'Muc2',	'Plxnb2',	'Gmds',	'Snrpb',	'1110008F13Rik',	'Pfkl',	'Lman2',	'Taldo1',	'Suclg1',	'Cyc1',	'Ppp1ca',	'Cda',	'Tspan1',	'Slc25a3',	'Cybrd1',	'Ech1',	'Ppa1',	'Uba52',	'Eif4g1',	'Srsf7',	'Esd',	'Rpl6',	'Cyba',	'Lrrc59',	'Ier3',	'Ckb',	'Uqcrc1',	'Ran',	'Polr2e',	'Cd9',	'Gipc1',	'Gstp1',	'Tomm40',	'Cldn2',	'Aldh1b1',	'Ctsd',	'Stoml2',	'Rac1',	'Syngr2',	'Eci2',	'Slc25a4',	'Psmd13',	'Lamp1',	'Plcb3',	'Gsto1',	'Ak2',	'Ccnd1',	'Tmem45b',	'Lmna',	'Rab25',	'Acadvl',	'Plac8',	'Maoa',	'Oxct1',	'Sumo1',	'Adh6a',	'Map2k2',	'Ugt2b34',	'Gpi1',	'Mapk13',	'Tmprss2',	'Kif5b',	'Txnrd1',	'Ctsh',	'Taf10',	'Grpel1',	'Cyp4f14',	'Fbp2',	'Actr3',	'Clic1',	'Twf1',	'Tmem14c',	'Hsp90aa1',	'Atp1a1',	'Junb',	'Wdr18',	'Elf3',	'Akr1b3',	'Cyp2d26',	'Ahsa1',	'Psma1',	'Ctnnb1',	'Hspa4',	'Txn2',	'Psmd3',	'Mrps10',	'Bag1',	'Mcm3',	'Cct2',	'Sdc1',	'Prmt1',	'Ddb1',	'Tsta3',	'Ccnd2',	'Eps8l3',	'Bsg' )      
AA_induced_invivo <- c('Hspa1b',	'Gm2000',	'Defa24',	'Banp',	'Lockd',	'Dpm3',	'Rps27rt',	'Snrpg',	'mt-Nd3',	'Cenpw',	'Slirp',	'Pet100',	'Tomm7',	'Tmem256',	'Cops9',	'Atp5k',	'Snhg6',	'Smim4',	'Mrpl52',	'Ndufb1-ps',	'Epb41l4aos',	'Rpl38',	'Pin4',	'Snhg8',	'Snrpf',	'Ndufa2',	'Psmg4',	'Atp5mpl',	'Romo1',	'Rps29',	'Rps21',	'Rps28',	'Anapc13',	'Atp5md',	'Lsm7',	'Tmem258',	'Uqcc2',	'Plp2',	'S100a6',	'Phgdh',	'Ndufa3',	'Eif3j1',	'Nop10',	'Sec61g',	'Aspm',	'Atox1',	'Polr2l',	'Ndufb2',	'Ccdc167',	'Fkbp5',	'Bola2',	'Ndufa5',	'Ndufc1',	'Mrps21',	'Snrpe',	'Uqcr11',	'Hist1h2ae',	'Sp3os',	'Ubl5',	'Rpl39',	'Tuba1c',	'Rpl37',	'Prelp',	'Mrpl33',	'Snrpd2',	'Naa38',	'Polr2k',	'Uqcr10',	'Cox17',	'Smim22',	'Rplp2',	'Prc1',	'Rpl37a',	'Rpl31',	'Sgo1',	'Rpl21',	'Fmc1',	'Ifitm3',	'Ppih',	'Nusap1',	'Trappc10',	'Zfas1',	'Rpl36',	'Ncapg',	'Polr2i',	'Tmsb15b2',	'Cox20',	'Ndufv3',	'Gm11808',	'Atp5e',	'Arpp19',	'Adra2a',	'Fxyd3',	'Mt1',	'Esco2',	'2410006H16Rik',	'Rpl35a',	'Ndufa1',	'Elob',	'Hist1h1b',	'Ndufs6',	'2310009B15Rik',	'Cit',	'Rpl41',	'Hspe1-rs1',	'S100a11',	'Ndufa7',	'Ndufb4',	'S100g',	'Srek1ip1',	'Sem1',	'Hmgcs2',	'Mif',	'Aurka',	'Diaph3',	'Cd44',	'Knl1',	'Snhg3',	'Top2a',	'Hells',	'Rrp8',	'Plk1',	'Ddx18',	'Cep164',	'Kif20b',	'Cdca2',	'Pmepa1',	'Jaml',	'Incenp',	'Smim26',	'Trip13',	'Mis18bp1',	'Rbmx2',	'Lsm6',	'Spc25',	'Rrp1b',	'Rps17',	'Lsm5',	'S100a13',	'Rpgrip1',	'Cdca8',	'Gas5',	'Ddx24',	'Hist1h1e',	'Cenpf',	'Hmmr',	'Atp5j2',	'Acot1',	'Kif22',	'Neil3',	'Rpl34',	'Mis12',	'Cox7a1',	'Tuba1b',	'Cbx6',	'Tubb2b',	'Ccdc34',	'Rpl36a',	'Rbmxl1',	'Ndc80',	'2310009A05Rik',	'Cep295',	'Bub1',	'Spink4',	'Hirip3',	'H2afj',	'Dlgap5',	'Atp5l',	'Smim11',	'Cox6c',	'Dnajb1',	'Kif11',	'Ubap2',	'Tpx2',	'Cebpz',	'Ctc1',	'Crip1',	'Sycn',	'Telo2',	'Rbis',	'Rassf4',	'Shcbp1',	'Cox7c',	'Gm43813',	'Rgcc',	'Notch1',	'Serf1',	'Ndufa13',	'Aurkb',	'Pbk',	'Sec61b',	'Lrrc31',	'Chordc1',	'Rps27l',	'1810037I17Rik',	'Rpp21',	'2200002D01Rik',	'E030030I06Rik',	'Fbxo5',	'Rad51ap1',	'Uqcrq',	'Polr2f',	'Rpl35',	'Ttr',	'Nup210',	'Krtcap2',	'Mcph1',	'Rps27',	'Ndufaf8',	'Lig1',	'Cd3eap',	'Ppwd1',	'1700097N02Rik',	'2810408I11Rik',	'Ckap2',	'2810004N23Rik',	'Rps2',	'Prpf4',	'Snhg18',	'Kif4',	'Sgo2a',	'Higd1a',	'Ndufs5',	'Xpo5',	'Soat1',	'Rps15',	'Anln',	'Rps26',	'Grwd1',	'Lrrc14',	'Cox7a2',	'Ncapd2',	'Tcof1',	'Rad54l',	'Ier3ip1',	'Cenpe',	'Apobec3',	'Lyrm2',	'Bub1b',	'Ccdc85c',	'Cep192',	'Prpf6',	'Cox7b',	'Rps25',	'Pinx1',	'Tex10',	'Gpatch4',	'Gnl3')          
AA_induced_stem2_marker <- c('Plekhm3',	'H2-Eb1',	'Cd74',	'H2-Ab1',	'Nucb2',	'Map7d1',	'Sik1',	'Upb1',	'Arid3a',	'Cyb561a3',	'Ankrd44',	'Ctsl',	'Tbc1d8',	'Slc29a3',	'Orai3',	'Dusp1',	'Ifnar2',	'Tcf4',	'Irf7',	'Unc93b1',	'Dusp6',	'Ptpn6',	'Tubgcp5',	'S100a6',	'Ptpre',	'Rell1',	'Tom1',	'Crlf2',	'Snx29',	'Gm37529',	'Gsn',	'Snx18',	'Gmfg',	'Mob3a',	'Mfsd12',	'Unc119',	'Uvrag',	'Cybc1',	'Cmtm7',	'Cnp',	'Npc1',	'Map4k2',	'Stat2',	'H2-Q4',	'Gm43329',	'Tmem123',	'Cdkn2d',	'Nek7',	'Slc15a4',	'Gns',	'Anxa5',	'Bcl11a',	'Ppp1r15a',	'Sema4b',	'Itpr1',	'Aldh3b1',	'Igtp',	'Gng10',	'Dock8',	'Fam174a',	'Pmepa1',	'Rogdi',	'Bbc3',	'Rab33b',	'Ftx',	'Pde7a',	'Lyn',	'Dclre1c',	'Atp13a2',	'Ccl5',	'Cbx4',	'Tgfbr1',	'Gpcpd1',	'Pqlc3',	'Rpgrip1',	'Arid3b',	'Ablim1',	'Gm31718',	'Dipk1a',	'Lgmn',	'Grina',	'Rhobtb2',	'Snx30',	'Nudt18',	'Trim12a',	'Slc49a4',	'Myo9b',	'Sgk3',	'Tor3a',	'Bicd2',	'Rwdd2a',	'Dtx2',	'Ctss',	'Slc25a4',	'5031439G07Rik',	'Grn',	'Tex2',	'Vamp1',	'A430035B10Rik',	'Inpp5k',	'Ccnd3',	'Irf8',	'Stat5a',	'Rab2b',	'Tap1',	'Pacc1',	'Ifnar1',	'Ier5',	'Pml',	'Acox3',	'Arl10',	'Il4ra',	'Tap2',	'Apobec3',	'Litaf',	'Malt1',	'Il17ra',	'Bmyc',	'Arhgap17',	'Zfp319',	'Bnip2',	'Pkd1',	'Tpst2',	'H2-T22',	'Calcoco1',	'Pip5k1c',	'Gmip',	'Crtc3',	'Ptpn1',	'L3mbtl3',	'Tcirg1',	'Prcp',	'Stat1',	'Cyth1',	'Nabp1',	'Nek6',	'Bin1',	'Arhgap27',	'Wdr81',	'Ptprj')

All_Genes <- rownames(ARDIR_obj@assays$MAGIC_RNA)
AA_induced_invivo <- intersect(All_Genes, AA_induced_invivo)
AA_induced_stem2_marker <- intersect(All_Genes, AA_induced_stem2_marker)
AA_induced_organoid_scrna <- intersect(All_Genes, AA_induced_organoid_scrna)

mean.exp <- zscore(colMeans(x = ARDIR_obj@assays$MAGIC_RNA@data[AA_induced_organoid_scrna, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = ARDIR_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  ARDIR_obj@meta.data$AA_organoid <- mean.exp
}
mean.exp <- zscore(colMeans(x = ARDIR_obj@assays$MAGIC_RNA@data[AA_induced_invivo, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = ARDIR_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  ARDIR_obj@meta.data$AA_induced <- mean.exp
}
mean.exp <- zscore(colMeans(x = ARDIR_obj@assays$MAGIC_RNA@data[AA_induced_stem2_marker, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = ARDIR_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  ARDIR_obj@meta.data$AA_stem_marker <- mean.exp
}

#pseudotime analysis
library(SeuratWrappers)
library(monocle3)
#install.packages("Signac")
library(Signac)
library(org.Mm.eg.db)
gene_symbol <- as.list(org.Mm.egSYMBOL)
Idents(ARDIR_obj) <- ARDIR_obj$Cell_Type

cds <- as.cell_data_set(ARDIR_obj)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(ARDIR_obj[["SCT"]])
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE,   close_loop = FALSE)
cds <- order_cells(cds, reduction_method = "UMAP")


plot_cells(cds = cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, cell_size = 2)
pseudo_umap <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = FALSE, label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, trajectory_graph_color = "#A8A8A8", graph_label_size=1.5,  cell_size = 2) 
pseudo_umap + scale_fill_gradientn(colors = annColors$Pseudotime)  + geom_polygon(data = hulls2, aes(x = x, y = y, group=CellType), fill=NA, color="black", size=0.3, alpha = 0.3) + geom_shadowtext(data = clusterMedian, aes(x = UMAP_1,y= UMAP_2, label = CellType, group = as.factor(CellType)), size = 3.5, bg.colour="black") + labs(fill = "Pseudotime")
ggsave(file="ARDIR_obj_Pseudotime_Umap_no_trajectory.svg", width=10, height=10)
plot_cells(cds,
           color_cells_by = "Cell_Type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

ARDIR_obj <- AddMetaData(
  object = ARDIR_obj,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)
ARDIR_obj$Pseudotime
Idents(ARDIR_obj) <- ARDIR_obj$Cell_Type
FeaturePlot(ARDIR_obj, c("Pseudotime"), pt.size = .0001, label = TRUE, repel = TRUE) +scale_fill_gradientn(colors = rev(brewer.pal(11,'Spectral'))) + scale_color_gradientn(colors =  rev(brewer.pal(11,'Spectral')))
ggsave(file = 'ARDIR_obj_pseudotime_umap.pdf', width=4, height=4, units="in")


DimPlot(ARDIR_obj)
#housekeeping
#saveRDS(ARD_obj, file = "/Users/vyom/data/Seurat_objects/ARD_IR_scRNA_all.rds")
#saveRDS(ARDIR_obj, file = "/Users/vyom/data/Seurat_objects/ARD_IR_scRNA_epith.rds")
#ARDIR_obj <- readRDS('/Users/vyom/CSHL Dropbox Team Dropbox/Vyom Shah/Seurat_Objects/ARD_IR_scRNA_epith.rds', refhook = NULL)

levels(factor(ARDIR_obj@meta.data$Cell_Type))
ARDIR_obj$Cell_Type <- factor(ARDIR_obj$Cell_Type, levels = c("ArdSC","RevSC","Ent.","Entend.", "Sec."))
Idents(ARDIR_obj) <- ARDIR_obj$Cell_Type
ARDIR_obj[["Cell_Type"]] <- Idents(ARDIR_obj)
new.cluster.ids <- c("Stem 2", "Clu/Atoh1+", "Ent.", "End.", "Sec. Prg.")
names(new.cluster.ids) <- levels(ARDIR_obj)
ARDIR_obj <- RenameIdents(ARDIR_obj, new.cluster.ids)
ARDIR_obj[["Cell_Type"]] <- Idents(ARDIR_obj)
DimPlot(ARDIR_obj, group.by = 'Cell_Type')
levels(factor(ARDIR_obj@meta.data$Treatment))

ARDIR_obj$Pseudotime

