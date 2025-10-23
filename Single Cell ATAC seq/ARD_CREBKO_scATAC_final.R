#Vyom Shah - Spring 2024
#ARD + ARDKO scATAC
#scATAC analysis
library(plyr)
library(Seurat)
library(Signac)
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
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(GenomicRanges)
library(future)
library(chromVARmotifs)

# load the RNA and ATAC data
ARD1_counts <- Read10X_h5("/Users/vyom/data/OE_5_ARD_creb_ATAC/Beyaz_OE06_7964_Creb1WT/outs/filtered_peak_bc_matrix.h5")
ARD1_fragpath <- "/Users/vyom/data/OE_5_ARD_creb_ATAC/Beyaz_OE06_7964_Creb1WT/outs/fragments.tsv.gz"

KO_ARD1_counts <- Read10X_h5("/Users/vyom/data/OE_5_ARD_creb_ATAC/Beyaz_OE06_7962_Creb1KO/outs/filtered_peak_bc_matrix.h5")
KO_ARD1_fragpath <- "/Users/vyom/data/OE_5_ARD_creb_ATAC/Beyaz_OE06_7962_Creb1KO/outs/fragments.tsv.gz"

KO_ARD2_counts <- Read10X_h5("/Users/vyom/data/OE_5_ARD_creb_ATAC/Beyaz_OE06_7963_Creb1KO/outs/filtered_peak_bc_matrix.h5")
KO_ARD2_fragpath <- "/Users/vyom/data/OE_5_ARD_creb_ATAC/Beyaz_OE06_7963_Creb1KO/outs/fragments.tsv.gz"

ARD_counts <- Read10X_h5("/Users/vyom/data/ARD_CREBKO_scATAC/Beyaz_OE02_7699_ARD_WT_ATAC/outs/filtered_peak_bc_matrix.h5")
ARD_fragpath <- "/Users/vyom/data/ARD_CREBKO_scATAC/Beyaz_OE02_7699_ARD_WT_ATAC/outs/fragments.tsv.gz"

KO_ARD_counts <- Read10X_h5("/Users/vyom/data/ARD_CREBKO_scATAC/Beyaz_OE02_7700_ARD_Creb1KO_ATAC/outs/filtered_peak_bc_matrix.h5")
KO_ARD_fragpath <- "/Users/vyom/data/ARD_CREBKO_scATAC/Beyaz_OE02_7700_ARD_Creb1KO_ATAC/outs/fragments.tsv.gz"

# read in peak sets
peaks.ARD1 <- read.table(
  file = "/Users/vyom/data/OE_5_ARD_creb_ATAC/Beyaz_OE06_7964_Creb1WT/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.KO_ARD1 <- read.table(
  file = "/Users/vyom/data/OE_5_ARD_creb_ATAC/Beyaz_OE06_7962_Creb1KO/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.KO_ARD2 <- read.table(
  file = "/Users/vyom/data/OE_5_ARD_creb_ATAC/Beyaz_OE06_7963_Creb1KO/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.ARD <- read.table(
  file = "/Users/vyom/data/ARD_CREBKO_scATAC/Beyaz_OE02_7699_ARD_WT_ATAC/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.KO_ARD <- read.table(
  file = "/Users/vyom/data/ARD_CREBKO_scATAC/Beyaz_OE02_7700_ARD_Creb1KO_ATAC/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)


# convert to genomic ranges
gr.ARD1 <- makeGRangesFromDataFrame(peaks.ARD1)
gr.ARD <- makeGRangesFromDataFrame(peaks.ARD)
gr.KO_ARD1 <- makeGRangesFromDataFrame(peaks.KO_ARD1)
gr.KO_ARD2 <- makeGRangesFromDataFrame(peaks.KO_ARD2)

gr.KO_ARD <- makeGRangesFromDataFrame(peaks.KO_ARD)

combined.peaks <- reduce(x = c(gr.ARD1, gr.KO_ARD2, gr.KO_ARD1, gr.ARD, gr.KO_ARD))


peakwidths <- width(combined.peaks)
plot(density(peakwidths))
combined.peaks <- combined.peaks[peakwidths  < 1650 & peakwidths > 100]
combined.peaks

#create fragment objects
#load metadata
md.ARD1 <- read.table(
  file = "/Users/vyom/data/OE_5_ARD_creb_ATAC/Beyaz_OE06_7964_Creb1WT/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.KO_ARD1 <- read.table(
  file = "/Users/vyom/data/OE_5_ARD_creb_ATAC/Beyaz_OE06_7962_Creb1KO/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.KO_ARD2 <- read.table(
  file = "/Users/vyom/data/OE_5_ARD_creb_ATAC/Beyaz_OE06_7963_Creb1KO/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.ARD <- read.table(
  file = "/Users/vyom/data/ARD_CREBKO_scATAC/Beyaz_OE02_7699_ARD_WT_ATAC/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.KO_ARD <- read.table(
  file = "/Users/vyom/data/ARD_CREBKO_scATAC/Beyaz_OE02_7700_ARD_Creb1KO_ATAC/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

boxplot(md.ARD1$passed_filters)
boxplot(md.ARD$passed_filters)
boxplot(md.KO_ARD1$passed_filters)
boxplot(md.KO_ARD2$passed_filters)
boxplot(md.KO_ARD$passed_filters)

md.ARD1 <- md.ARD1[md.ARD1$passed_filters > 1, ]
md.ARD <- md.ARD[md.ARD$passed_filters > 1, ]
md.KO_ARD1 <- md.KO_ARD1[md.KO_ARD1$passed_filters > 1, ]
md.KO_ARD2 <- md.KO_ARD2[md.KO_ARD2$passed_filters > 1, ]
md.KO_ARD <- md.KO_ARD[md.KO_ARD$passed_filters > 1, ]

frags.ARD1 <- CreateFragmentObject(
  path = ARD1_fragpath,
  cells = rownames(md.ARD1)
)
frags.KO_ARD1 <- CreateFragmentObject(
  path = KO_ARD1_fragpath,
  cells = rownames(md.KO_ARD1)
)
frags.KO_ARD2 <- CreateFragmentObject(
  path = KO_ARD2_fragpath,
  cells = rownames(md.KO_ARD2)
)

frags.ARD <- CreateFragmentObject(
  path = ARD_fragpath,
  cells = rownames(md.ARD)
)

frags.KO_ARD <- CreateFragmentObject(
  path = KO_ARD_fragpath,
  cells = rownames(md.KO_ARD)
)


ARD1.counts <- FeatureMatrix(
  fragments = frags.ARD1,
  features = combined.peaks,
  cells = rownames(md.ARD1)
)
KO_ARD1.counts <- FeatureMatrix(
  fragments = frags.KO_ARD1,
  features = combined.peaks,
  cells = rownames(md.KO_ARD1)
)
KO_ARD2.counts <- FeatureMatrix(
  fragments = frags.KO_ARD2,
  features = combined.peaks,
  cells = rownames(md.KO_ARD2)
)

ARD.counts <- FeatureMatrix(
  fragments = frags.ARD,
  features = combined.peaks,
  cells = rownames(md.ARD)
)
KO_ARD.counts <- FeatureMatrix(
  fragments = frags.KO_ARD,
  features = combined.peaks,
  cells = rownames(md.KO_ARD)
)

ARD1_assay <- CreateChromatinAssay(ARD1.counts, fragments = frags.ARD1)
ARD1_obj <- CreateSeuratObject(ARD1_assay, assay = "ATAC", meta.data=md.ARD1)

KO_ARD1_assay <- CreateChromatinAssay(KO_ARD1.counts, fragments = frags.KO_ARD1)
KO_ARD1_obj <- CreateSeuratObject(KO_ARD1_assay, assay = "ATAC", meta.data=md.KO_ARD1)

KO_ARD2_assay <- CreateChromatinAssay(KO_ARD2.counts, fragments = frags.KO_ARD2)
KO_ARD2_obj <- CreateSeuratObject(KO_ARD2_assay, assay = "ATAC", meta.data=md.KO_ARD2)

ARD_assay <- CreateChromatinAssay(ARD.counts, fragments = frags.ARD)
ARD_atac_obj <- CreateSeuratObject(ARD_assay, assay = "ATAC", meta.data=md.ARD)

KO_ARD_assay <- CreateChromatinAssay(KO_ARD.counts, fragments = frags.KO_ARD)
KO_ARD_atac_obj <- CreateSeuratObject(KO_ARD_assay, assay = "ATAC", meta.data=md.KO_ARD)

ARD1_obj$dataset <- 'ARD1'
ARD_atac_obj$dataset <-  'ARD'
KO_ARD1_obj$dataset <- 'KO_ARD1'
KO_ARD_atac_obj$dataset <-  'KO_ARD'
KO_ARD2_obj$dataset <- 'KO_ARD2'

ARDKO_atac_obj <- merge(
  x = ARD1_obj,
  y = list(ARD_atac_obj, KO_ARD1_obj, KO_ARD_atac_obj, KO_ARD2_obj),
  add.cell.ids = c("ARD1", "ARD",  "KO_ARD1", "KO_ARD", 'KO_ARD2'))


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
seqlevelsStyle(annotations) <- "UCSC"

  # add the gene information to the object
  Annotation(ARDKO_atac_obj) <- annotations
  rownames(ARDKO_atac_obj)
  # compute nucleosome signal score per cell
  ARDKO_atac_obj <- NucleosomeSignal(object = ARDKO_atac_obj)
  
  # compute TSS enrichment score per cell
  ARDKO_atac_obj <- TSSEnrichment(object = ARDKO_atac_obj, fast = FALSE)
  
  # add blacklist ratio and fraction of reads in peaks
  ARDKO_atac_obj$pct_reads_in_peaks <- ARDKO_atac_obj$peak_region_fragments / ARDKO_atac_obj$passed_filters * 100
  ARDKO_atac_obj$blacklist_ratio <- ARDKO_atac_obj$blacklist_region_fragments / ARDKO_atac_obj$peak_region_fragments
  
  ARDKO_atac_obj$high.tss <- ifelse(ARDKO_atac_obj$TSS.enrichment > 3, 'High', 'Low')
  TSSPlot(ARDKO_atac_obj, group.by = 'dataset') + NoLegend()
  ARDKO_atac_obj$blacklist_region_fragments
  VlnPlot(
    object = ARDKO_atac_obj,
    features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal', 'pct_reads_in_peaks','passed_filters'),
    pt.size = 0.0,
    ncol = 5,
    group.by = 'dataset'
  )
  ARDKO_atac_obj@assays$ATAC
  ARDKO_atac_obj <- subset(
    x = ARDKO_atac_obj,
    subset = 
      nCount_ATAC > 100 &
      nCount_ATAC < 20000 &
      nucleosome_signal < 5 &
      nucleosome_signal > 0 &
      TSS.enrichment > 0 &
      TSS.enrichment < 10 &
      pct_reads_in_peaks > 25 &
      passed_filters > 1000)

  ARDKO_atac_obj <- FindTopFeatures(ARDKO_atac_obj, min.cutoff = 10)
  ARDKO_atac_obj <- RunTFIDF(ARDKO_atac_obj)
  ARDKO_atac_obj <- RunSVD(ARDKO_atac_obj)

DepthCor(ARDKO_atac_obj, n = 50, reduction = 'lsi')
ARDKO_atac_obj <- RunUMAP(ARDKO_atac_obj, reduction = "lsi", dims = 2:50)
p1 <- DimPlot(ARDKO_atac_obj, group.by = "dataset", label = TRUE)
p1

#clustering and annotation
ARDKO_atac_obj <- RunUMAP(ARDKO_atac_obj, reduction = "lsi", dims = 2:50)
ARDKO_atac_obj <- FindNeighbors(object = ARDKO_atac_obj, reduction = 'lsi', dims = 2:50)
ARDKO_atac_obj <- FindClusters(object = ARDKO_atac_obj, verbose = FALSE, algorithm = 3, resolution = 1)

DotPlot_Sig <- unique(c('Epcam','Ptprc','Cd3g','Cd3e','Cd8a','Cd8b1','Cd4','Trac','Tcrg-C1','Lag3','Pdcd1','Havcr2','Tox','Tcf7','Gzmb','Tbx21','Ifng','Gata3','Il5','Il13','Rorc','Il17a','Il17f','Foxp3','Il10','Il2rb','Il2ra','Klrd1','Cd19','Jchain','Ighm','Ighg1','S100a8','S100a9','Ly6a','Cd14','Cd68','Cd74','Ciita','Nrc1','Klre1','Itgam','Itgax','H2-Eb1','H2-Ab1','Arg1','Mrc1','Tgfbi','Ccr2','Vegfa','Prdx1','Clec4d','Ccl5','Cd83','Ccr7','Fcn1','Msrb1','Ly6g','Col3a1','Sparc', 'Cd11b', 'Csf1r', 'Cxcl9','Cxcl10','Ccl5','Cxcl13','Il10','Ccl1','Ccl17','Ifng'))     
DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3","S100a6","Ly6a","Anxa3", "Areg","Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a",'Tmigd1', 'Fabp6', 'Slc51b', 'Slc51a', "Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Lyz1","Defa17","Defa24","Ang4","Pou2f3","Avil","Tuba1a","Adh1",'Epcam','Ptprc','Cd8a','Cd4','Cd74','Itgax') 

DotPlot(ARDKO_atac_obj, features = DotPlot_Sig, assay = 'RNA') + labs(y= "Cell Treatment", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 5)) +
  theme(text = element_text(size=5), axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())

cluster_ids <- c("Ent.", "Stem 1", "TA", "EP", "Ent.", "T", "EP", "TA", "T", "T", "T", "T", "Gob.", "Ent.", "Gob.", "Stem 1", "Entend.", "Myeloid", "T", "Entend.", "Tuft", "Entend.", "Stem 2")

ARDKO_atac_obj[["Cell_Type"]] <- Idents(ARDKO_atac_obj)
names(cluster_ids) <- levels(ARDKO_atac_obj)
ARDKO_atac_obj <- RenameIdents(ARDKO_atac_obj, cluster_ids)
ARDKO_atac_obj[["Cell_Type"]] <- Idents(ARDKO_atac_obj)
DimPlot(ARDKO_atac_obj, reduction = "umap", group.by= 'Cell_Type')

DimPlot(ARDKO_atac_obj, reduction = "umap", group.by= 'Treatment')
unique(ARDKO_atac_obj$Cell_Type)
unique(ARDKO_atac_obj$Treatment)

my_levels <- c("Stem 1", "Stem 2", "TA", "EP", "Ent.",'SP', "End.", "Gob.","Tuft",'T','DP T', 'Myeloid')
my_levels1 <- c( "ARD", "KO_ARD")
ARDKO_atac_obj$Cell_Type <- factor(x = ARDKO_atac_obj$Cell_Type, levels = my_levels)
ARDKO_atac_obj$Treatment <- factor(x = ARDKO_atac_obj$Treatment, levels = my_levels1)

#motif analysis
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2024)
#library(JASPAR2020)
library(TFBSTools)
library(JASPAR2024) 
library(TFBSTools)
jaspar <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
pfm <- TFBSTools::getMatrixSet(sq24, list(tax_group = "vertebrates", collection = "CORE"))

JASPAR2024_obj <- JASPAR2024()
class(JASPAR2024_obj) <- 'JASPAR2020'
pfm1 <- getMatrixSet(
  x = JASPAR2024_obj,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
library(BSgenome.Mmusculus.UCSC.mm10)

main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
keep.peaks <- as.logical(seqnames(granges(ARDKO_atac_obj)) %in% main.chroms)
ARDKO_atac_obj <- ARDKO_atac_obj[keep.peaks, ]

ARDKO_atac_obj <- AddMotifs(
  object = ARDKO_atac_obj,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

Idents(ARDKO_atac_obj) <- ARDKO_atac_obj$Treatment
DE_subset <- FindMarkers(ARDKO_atac_obj, ident.1 = "KO_ARD", ident.2 = "ARD",  test.use = 'MAST', logfc.threshold = .1, min.pct = .05, latent.vars = 'nCount_ATAC', assay = 'ATAC')
DE_subset$peaks <- rownames(DE_subset)
Open_Regions <- rownames(DE_subset)
Clostest_genes <- ClosestFeature(ARDKO_atac_obj, regions = Open_Regions)
DE_subset$gene <-  Clostest_genes$gene_name
write.csv(DE_subset, paste0('All_ARDvsKO_ARD_ATAC_DE_no_imm_rev.csv')) 

celltypes <- c(levels(as.factor(ARDKO_atac_obj$Cell_Type)))

Idents(ARDKO_atac_obj) <- ARDKO_atac_obj$Cell_Type
for(i in celltypes){
  subset_cell <- subset(ARDKO_atac_obj,  idents = i)
  Idents(subset_cell) <- subset_cell$Treatment
  DE_subset <- FindMarkers(subset_cell, ident.1 = "ARD", ident.2 = "KO_ARD",test.use = 'MAST',latent.vars = 'nCount_ATAC', logfc.threshold = .1, min.pct = .05, assay = 'ATAC')
  DE_subset$peaks <- rownames(DE_subset)
  Open_Regions <- rownames(DE_subset)
  Clostest_genes <- ClosestFeature(subset_cell, regions = Open_Regions)
  DE_subset$gene <-  Clostest_genes$gene_name
  write.csv(DE_subset, paste0(i,'_ARD_KO_ATAC_DE.csv')) 
}

DefaultAssay(ARDKO_atac_obj) <- 'ATAC'
FeaturePlot(ARDKO_atac_obj, features = 'S100a6')

ararev.bed <- granges(reversal_obj)
df <- data.frame(seqnames=seqnames(ararev.bed),
                 starts=start(ararev.bed)-1,
                 ends=end(ararev.bed),
                 names=c(rep(".", length(ararev.bed))),
                 scores=c(rep(".", length(ararev.bed))),
                 strands=strand(ararev.bed))
write.table(df, file="ARD_wt_atac_combined.bed", quote=F, sep="\t", row.names=F, col.names=F)


#homer code for command line: annotatePeaks.pl /Users/vyom/data/HOMER/ARD_KOactual_atac_combined.bed mm10 -gtf  /Users/vyom/data/reversal_scATAC/genes.gtf -annStats combined_stats.txt > ARD_KO_Homer_combined.txt
Homer_peaks <- read.delim("~/data/HOMER/ARD_WT_Homer_combined.txt")
Homer_peaks <- read.delim("~/data/HOMER/ARD_KO_Homer_combined.txt")
Homer_peaks$peak <- paste0(Homer_peaks$Chr,'-',Homer_peaks$Start,'-', Homer_peaks$End)
Homer_peaks$Annotation_simple <- gsub("([A-Za-z]+).*", "\\1", Homer_peaks$Annotation)

celltypes <- levels(as.factor(ARDKO_atac_obj$Cell_Type))
for(i in celltypes) 
  {
  
  DE_subset_ARD_KO <- read.csv(paste0('~/All_ARDvsKO_ARD_ATAC_DE_no_imm_rev.csv'))
  DE_subset_ARD_KO <- DE_subset_ARD_KO[DE_subset_ARD_KO$p_val_adj < 0.05 ,]
  #DE_subset_ARD_KO <- read.csv(paste0(i,'_ARD_KO_ATAC_DE.csv'))
  #DE_subset_ARD_KO <- DE_subset1
  #DE_subset_ARD_KO <- read.csv(paste0('ARD_Vs_control_DE.csv'))
  #DE_subset_ARD_KO <- read.csv(paste0('reversal_Vs_control_DE.csv'))

  Up_matched_ARD_KO <- intersect(DE_subset_ARD_KO[DE_subset_ARD_KO$avg_log2FC > 0 & DE_subset_ARD_KO$p_val_adj < 0.05 ,]$X, Homer_peaks$peak)
  Down_matched_ARD_KO <- intersect(DE_subset_ARD_KO[DE_subset_ARD_KO$avg_log2FC < 0 & DE_subset_ARD_KO$p_val_adj < 0.05 ,]$X, Homer_peaks$peak)
  
  cat_anno <- c('exon', 'Intergenic', 'intron', 'promoter', 'TTS')
  
  down_ARD_KO <- Homer_peaks[Homer_peaks$peak %in% Down_matched_ARD_KO,]$Annotation_simple %>% factor(levels = cat_anno) %>% tabulate()
  up_ARD_KO <- Homer_peaks[Homer_peaks$peak %in% Up_matched_ARD_KO,]$Annotation_simple %>% factor(levels = cat_anno) %>% tabulate()
  if (length(up_ARD_KO) < 5) {
    up_ARD_KO <- c(up_ARD_KO,0)
  }
  annotations <- Homer_peaks[Homer_peaks$peak %in% Down_matched_ARD_KO,]$Annotation_simple %>% as.factor() %>% levels()
  Annotation_homer <- data.frame(down_ARD_KO, up_ARD_KO)
  Annotation_homer_prop <- Annotation_homer %>% mutate(across(where(is.numeric), ~ ./sum(.)))
  Annotation_homer_prop <- cbind(annotations,Annotation_homer_prop)
  
  
  Annotation_homer_prop <- melt(Annotation_homer_prop)
  Annotation_homer_prop$annotations <- factor(Annotation_homer_prop$annotations,levels = annotations)
  ggplot(Annotation_homer_prop, aes(x=variable, y= value, fill = annotations))+theme_vyom  + geom_col(colour = "black", position = "fill") + scale_y_continuous(labels = scales::percent) + scale_fill_brewer(palette = "Pastel2")
  ggsave(file = paste0(i,'_scATAC_ARD_KO_per_accessible_annotations_proportion.pdf'), width=3, height=3, units="in")
  
  Annotation_homer_sum <- cbind(annotations,Annotation_homer)
  Annotation_homer_sum <- melt(Annotation_homer_sum)
  Annotation_homer_sum$annotations <- factor(Annotation_homer_sum$annotations,levels = annotations)
  ggplot(Annotation_homer_sum, aes(x=variable, y= value, fill = annotations))+theme_vyom   + geom_col(colour = "black")  + scale_fill_brewer(palette = "Pastel2")
  ggsave(file = paste0(i,'sARD_KO_scATAC_per_accessible_annotations_sum_rev.pdf'), width=3, height=3, units="in")
  
  Up_ARD_KO_PEAKS_promoter <- Homer_peaks[Homer_peaks$peak %in% Up_matched_ARD_KO & Homer_peaks$Annotation_simple == 'promoter',]$peak
  Down_ARD_KO_PEAKS_promoter <- Homer_peaks[Homer_peaks$peak %in% Down_matched_ARD_KO & Homer_peaks$Annotation_simple == 'promoter',]$peak
  
  motifs_ARD_KO_up <- FindMotifs(
    object = ARDKO_atac_obj,
    features = Up_ARD_KO_PEAKS_promoter)
  motifs_ARD_KO_up <- motifs_ARD_KO_up[motifs_ARD_KO_up$pvalue < .05,]
  
  write.csv(motifs_ARD_KO_up, paste0(i,'_promoter_ARD_KO_up_motifs.csv'))
  
  motifs_ARD_KO_up <- FindMotifs(
    object = ARDKO_atac_obj,
    features = Down_ARD_KO_PEAKS_promoter)
  motifs_ARD_KO_up <- motifs_ARD_KO_up[motifs_ARD_KO_up$pvalue < .05,]
  
  write.csv(motifs_ARD_KO_up, paste0(i,'_promoter_ARD_KO_down_motifs.csv'))
  
  #Up_ARD_KO_genes_promoter <- DE_subset_ARD_KO[DE_subset_ARD_KO$X %in% Homer_peaks[Homer_peaks$peak %in% Homer_peaks$Annotation_simple == 'promoter',]$peak,]
  #Up_ARD_KO_genes_promoter$peak <- Up_ARD_KO_genes_promoter$X
  #Up_ARD_KO_genes_promoter1 <- merge(Up_ARD_KO_genes_promoter, Homer_peaks[Homer_peaks$peak %in% Up_ARD_KO_genes_promoter$peak,] %>% dplyr::select(Gene.Name,peak), by= 'peak', all.x=TRUE)
  
  #write.csv(Up_ARD_KO_genes_promoter, paste0(i,'_promoter_ARD_KO_up_PEAKS.csv'))
  
  matched_ard <- intersect(DE_subset_ARD_KO[DE_subset_ARD_KO$p_val_adj < 0.05,]$X, Homer_peaks$peak)
  DE_subset_ARD_KO1 <- DE_subset_ARD_KO[DE_subset_ARD_KO$X %in% Homer_peaks[Homer_peaks$peak %in% matched_ard,]$peak,]
  homer_subset <- Homer_peaks[Homer_peaks$peak %in% DE_subset_ARD_KO1$X,] %>% dplyr::select(peak, Annotation_simple, Gene.Name )
  DE_subset_ARD_KO1$peak <- DE_subset_ARD_KO1$X
  DE_subset_ARD_KO1 <- merge(DE_subset_ARD_KO1, homer_subset, by="peak")
  write.csv(DE_subset_ARD_KO1, paste0(i, '_ARD_KO_promoters_DE_peak.csv'))
  
  if (max(-log10(DE_subset_ARD_KO$p_val_adj)) < 321){
    ylim_num <- max(-log10(DE_subset_ARD_KO$p_val_adj))} else {
      ylim_num <- 320}
  EnhancedVolcano(DE_subset_ARD_KO, lab = DE_subset_ARD_KO$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'ARD_KO vs Control',
                  pCutoff = .05, FCcutoff = 0.25, xlim = c(min(DE_subset_ARD_KO$avg_log2FC)-.05,max(DE_subset_ARD_KO$avg_log2FC)+.05),ylim = c(0,ylim_num+5),
                  subtitle = 'ALL CELLS', gridlines.minor = FALSE, gridlines.major = FALSE)
  ggsave(file = paste0(i,'_DE_peak_volcano_promoters.pdf'), width=6, height=6, units="in")
  
}
gene.list <- c('Fgfbp1','S100a6', 'Lgr5','chr10-115580000-115590000', 'Ascl2','Anxa3','Anxa2',  'Clu', 'Yes1', 'Lef1', 'Ly6a','Creb1','Olfm4','Sox9','Mki67')

DefaultAssay(ARDKO_atac_obj) <- 'ATAC'
for(i in unique(gene.list)){
  plot <- CoveragePlot(
    object = ARDKO_atac_obj,
    region = i,
    group.by = 'Treatment_celltype',
    expression.assay = "ATAC",
    extend.upstream = 2000,
    extend.downstream = 2000,
    annotation = TRUE,
    peaks = FALSE,
    tile = FALSE,
    links = FALSE
  ) + scale_fill_manual(values = rep(c('#d95f02' ,'#01949A'), 12)) 
  plot <- plot & scale_fill_manual(values = rep(c('#d95f02' ,'#01949A'), 12)) 
  ggsave(file = paste0(i,'_scatac_coverage_track_ARD_KO.pdf'),plot = plot, width=5, height=7, units="in")
}

#housekeeping
# saveRDS(ARDKO_atac_obj, file = "./data/Seurat_Objects/ARD_KO_ATAC_fin_5_23_Sobj.rds")
# ARDKO_atac_obj <- readRDS("/Users/vyom/CSHL Dropbox Team Dropbox/Vyom Shah/Seurat_Objects/ARD_KO_ATAC_fin_5_23_Sobj.rds", refhook = NULL)

