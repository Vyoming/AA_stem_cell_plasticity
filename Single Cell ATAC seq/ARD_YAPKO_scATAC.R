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
ARD1_counts <- Read10X_h5("/Users/vyom/data/yapko_ovy_new/Beyaz_OE07_ATAC_6812_YAP_WT_ARD/outs/filtered_peak_bc_matrix.h5")
ARD1_fragpath <- "/Users/vyom/data/yapko_ovy_new/Beyaz_OE07_ATAC_6812_YAP_WT_ARD/outs/fragments.tsv.gz"

KO_ARD1_counts <- Read10X_h5("/Users/vyom/data/yapko_ovy_new/Beyaz_OE07_ATAC_6813_YAP_KO_ARD/outs/filtered_peak_bc_matrix.h5")
KO_ARD1_fragpath <- "/Users/vyom/data/yapko_ovy_new/Beyaz_OE07_ATAC_6813_YAP_KO_ARD/outs/fragments.tsv.gz"

KO_ARD2_counts <- Read10X_h5("/Users/vyom/data/yapko_ovy_new/Beyaz_OE07_ATAC_6406_YAP_KO_ARD/outs/filtered_peak_bc_matrix.h5")
KO_ARD2_fragpath <- "/Users/vyom/data/yapko_ovy_new/Beyaz_OE07_ATAC_6406_YAP_KO_ARD/outs/fragments.tsv.gz"

ARD_counts <- Read10X_h5("/Users/vyom/data/ARD_YAPKO_scATAC/Beyaz_OE05_6800_ARD_WT/outs/filtered_peak_bc_matrix.h5")
ARD_fragpath <- "/Users/vyom/data/ARD_YAPKO_scATAC/Beyaz_OE05_6800_ARD_WT/outs/fragments.tsv.gz"

KO_ARD_counts <- Read10X_h5("/Users/vyom/data/ARD_YAPKO_scATAC/Beyaz_OE05_6801_ARD_YAP_KO/outs/filtered_peak_bc_matrix.h5")
KO_ARD_fragpath <- "/Users/vyom/data/ARD_YAPKO_scATAC/Beyaz_OE05_6801_ARD_YAP_KO/outs/fragments.tsv.gz"

# read in peak sets
peaks.ARD1 <- read.table(
  file = "/Users/vyom/data/yapko_ovy_new/Beyaz_OE07_ATAC_6812_YAP_WT_ARD/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.KO_ARD1 <- read.table(
  file = "/Users/vyom/data/yapko_ovy_new/Beyaz_OE07_ATAC_6813_YAP_KO_ARD/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.KO_ARD2 <- read.table(
  file = "/Users/vyom/data/yapko_ovy_new/Beyaz_OE07_ATAC_6406_YAP_KO_ARD/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.ARD <- read.table(
  file = "/Users/vyom/data/ARD_YAPKO_scATAC/Beyaz_OE05_6800_ARD_WT/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.KO_ARD <- read.table(
  file = "/Users/vyom/data/ARD_YAPKO_scATAC/Beyaz_OE05_6801_ARD_YAP_KO/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)


# convert to genomic ranges
gr.ARD1 <- makeGRangesFromDataFrame(peaks.ARD1)
gr.ARD <- makeGRangesFromDataFrame(peaks.ARD)
gr.KO_ARD1 <- makeGRangesFromDataFrame(peaks.KO_ARD1)
gr.KO_ARD2 <- makeGRangesFromDataFrame(peaks.KO_ARD2)
gr.KO_ARD <- makeGRangesFromDataFrame(peaks.KO_ARD)

combined.peaks <- GenomicRanges::reduce(x = c(gr.ARD1, gr.KO_ARD2, gr.KO_ARD1, gr.ARD, gr.KO_ARD))


peakwidths <- width(combined.peaks)
plot(density(peakwidths))
combined.peaks <- combined.peaks[peakwidths  < 1650 & peakwidths > 100]
combined.peaks

#create fragment objects
#load metadata
md.ARD1 <- read.table(
  file = "/Users/vyom/data/yapko_ovy_new/Beyaz_OE07_ATAC_6812_YAP_WT_ARD/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.KO_ARD1 <- read.table(
  file = "/Users/vyom/data/yapko_ovy_new/Beyaz_OE07_ATAC_6813_YAP_KO_ARD/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.KO_ARD2 <- read.table(
  file = "/Users/vyom/data/yapko_ovy_new/Beyaz_OE07_ATAC_6406_YAP_KO_ARD/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.ARD <- read.table(
  file = "/Users/vyom/data/ARD_YAPKO_scATAC/Beyaz_OE05_6800_ARD_WT/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.KO_ARD <- read.table(
  file = "/Users/vyom/data/ARD_YAPKO_scATAC/Beyaz_OE05_6801_ARD_YAP_KO/outs/singlecell.csv",
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

md.ARD1 <- md.ARD1[md.ARD1$passed_filters > 100, ]
md.ARD <- md.ARD[md.ARD$passed_filters > 100, ]
md.KO_ARD1 <- md.KO_ARD1[md.KO_ARD1$passed_filters > 100, ]
md.KO_ARD2 <- md.KO_ARD2[md.KO_ARD2$passed_filters > 100, ]
md.KO_ARD <- md.KO_ARD[md.KO_ARD$passed_filters > 100, ]

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

ARDKO_atac_obj <- subset(
  x = ARDKO_atac_obj,
  subset = 
    nCount_ATAC > 100 &
    nCount_ATAC < 10000 &
    nucleosome_signal < 10 &
    nucleosome_signal > 0 &
    TSS.enrichment > 0 &
    TSS.enrichment < 10 &
    pct_reads_in_peaks > 25 &
    passed_filters > 1000)

ARDKO_atac_obj <- FindTopFeatures(ARDKO_atac_obj, min.cutoff = 10)
ARDKO_atac_obj <- RunTFIDF(ARDKO_atac_obj)
ARDKO_atac_obj <- RunSVD(ARDKO_atac_obj)

DefaultAssay(ARDKO_atac_obj)
DepthCor(ARDKO_atac_obj, n = 50, reduction = 'lsi')
ARDKO_atac_obj <- RunUMAP(ARDKO_atac_obj, reduction = "lsi", dims = 2:50)
p1 <- DimPlot(ARDKO_atac_obj, group.by = "dataset", label = TRUE)
p1
table(ARDKO_atac_obj$dataset)
#umap ge
#ARDKO_atac_obj <- ARDko_atac_obj
table(ARDKO_atac_obj$dataset)

library(harmony)
ARDKO_atac_obj
ARDKO_atac_obj$Treatment <- ARDKO_atac_obj$dataset
ARDKO_atac_obj <- RunHarmony(
  object = ARDKO_atac_obj,
  group.by.vars = 'dataset',
  reduction.use = 'lsi',
  assay.use = 'ATAC',
  epsilon.harmony = 1e-10,
  epsilon.cluster = 1e-10,
  plot_convergence = TRUE,
  project.dim = FALSE
)
ARDKO_atac_obj@assays$ATAC
table(ARDKO_atac_obj$dataset)
ARDKO_atac_obj <- RunUMAP(ARDKO_atac_obj, dims = 1:50, reduction = 'harmony')

DimPlot(ARDKO_atac_obj, group.by = "dataset", label = TRUE)

ARDKO_atac_obj <- FindNeighbors(object = ARDKO_atac_obj, reduction = 'harmony', dims = 2:50)
ARDKO_atac_obj <- FindClusters(object = ARDKO_atac_obj, verbose = FALSE, algorithm = 2, resolution =1)
DimPlot(object = ARDKO_atac_obj, label = TRUE, group.by = 'ATAC_snn_res.1') 
FeaturePlot(ARDKO_atac_obj, features = 'passed_filters')
VlnPlot(ARDKO_atac_obj, features = 'passed_filters', group.by = 'ATAC_snn_res.1', pt.size = 0)

#transfer reference to query: annotate
ARDKO_atac_obj <- readRDS('/Users/vyom/data/Seurat_Objects/scRNA.rds', refhook = NULL)
DefaultAssay(ARDKO_atac_obj) <- 'integrated'

#create gene activity matrix
gene.activities <- GeneActivity(ARDKO_atac_obj)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
ARDKO_atac_obj[['RNA']] <- CreateAssayObject(counts = gene.activities)
ARDKO_atac_obj <- NormalizeData(
  object = ARDKO_atac_obj,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(ARDKO_atac_obj$nCount_RNA)
)

cluster_ids <- c("TA", "Stem 1", "Ent.", "EP", "T", "T", "TA", "Gob.", "EP", "Stem 1", "TA", "End.", "TA", "Stem 1", "Tuft", "End.", "Stem 1", "Pan.", "Myeloid","Stem 2",  "Stem 1",  "End.",'TA', "EP", "TA", "TA", "TA", "Ent.")

Idents(ARDKO_atac_obj) <- ARDKO_atac_obj$seurat_clusters
ARDKO_atac_obj[["Cell_Type"]] <- Idents(ARDKO_atac_obj)
names(cluster_ids) <- levels(ARDKO_atac_obj)
ARDKO_atac_obj <- RenameIdents(ARDKO_atac_obj, cluster_ids)
ARDKO_atac_obj[["Cell_Type"]] <- Idents(ARDKO_atac_obj)
DimPlot(ARDKO_atac_obj, reduction = "umap", group.by= 'Cell_Type')

DimPlot(ARDKO_atac_obj, reduction = "umap", group.by= 'Treatment')
unique(ARDKO_atac_obj$Cell_Type)
unique(ARDKO_atac_obj$Treatment)

my_levels <- c("Stem 1", "Stem 2", "TA", "EP", "Ent.", "End.", "Gob.",'Pan.',"Tuft",'T', 'Myeloid')
my_levels1 <- c( "ARD", "KO_ARD")
ARDKO_atac_obj$Cell_Type <- factor(x = ARDKO_atac_obj$Cell_Type, levels = my_levels)
ARDKO_atac_obj$Treatment <- factor(x = ARDKO_atac_obj$Treatment, levels = my_levels1)
DimPlot(ARDKO_atac_obj, reduction = "umap", group.by= 'Cell_Type')
DimPlot(ARDKO_atac_obj, reduction = "umap", group.by= 'Treatment')

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
DE_subset <- FindMarkers(ARDKO_atac_obj, ident.1 = "ARD", ident.2 = "KO_ARD",  test.use = 'MAST', logfc.threshold = .1, min.pct = .05, latent.vars = 'nCount_ATAC', assay = 'ATAC')
DE_subset$peaks <- rownames(DE_subset)
DE_subset1 <- DE_subset
DE_subset <- DE_subset[!grepl(paste(c('GL456233.1','JH584304.1','GL456216.1', 'JH584295.1'), collapse="|"),DE_subset$peaks),]
Open_Regions <- rownames(DE_subset)
Clostest_genes <- ClosestFeature(ARDKO_atac_obj, regions = Open_Regions)
DE_subset$gene <-  Clostest_genes$gene_name
write.csv(DE_subset, paste0('All_ARDvsKO_ARD_ATAC_DE.csv')) 

celltypes <- c(levels(as.factor(ARDKO_atac_obj$Cell_Type)))
celltypes <- c('Stem 1', 'Stem 2','TA','EP','T','Myeloid','Enterocyte', 'Tuft')
celltypes <- c('Stem 1', 'Stem 2')

Idents(ARDKO_atac_obj) <- ARDKO_atac_obj$Cell_Type
i = 'Stem 2'
i = 'TA'
subset_cell$nCount_ATAC
for(i in celltypes){
  subset_cell <- subset(ARDKO_atac_obj,  idents = i)
  Idents(subset_cell) <- subset_cell$Treatment
  DE_subset <- FindMarkers(subset_cell, ident.1 = "ARD", ident.2 = "KO_ARD",test.use = 'MAST',latent.vars = 'nCount_ATAC', logfc.threshold = .1, min.pct = .05, assay = 'ATAC')
  DE_subset$peaks <- rownames(DE_subset)
  DE_subset <- DE_subset[!grepl(paste(c('GL456233.1','JH584304.1','GL456216.1', 'JH584295.1'), collapse="|"),DE_subset$peaks),]
  Open_Regions <- rownames(DE_subset)
  Clostest_genes <- ClosestFeature(subset_cell, regions = Open_Regions)
  DE_subset$gene <-  Clostest_genes$gene_name
  write.csv(DE_subset, paste0(i,'_ARD_KO_ATAC_DE.csv')) 
}

DefaultAssay(ARDKO_atac_obj) <- 'ATAC'
FeaturePlot(ARDKO_atac_obj, features = 'S100a6')

ararev.bed <- granges(ARDKO_atac_obj)
df <- data.frame(seqnames=seqnames(ararev.bed),
                 starts=start(ararev.bed)-1,
                 ends=end(ararev.bed),
                 names=c(rep(".", length(ararev.bed))),
                 scores=c(rep(".", length(ararev.bed))),
                 strands=strand(ararev.bed))
write.table(df, file="ARD_KO_atac_combined.bed", quote=F, sep="\t", row.names=F, col.names=F)


#homer code for command line: annotatePeaks.pl /Users/vyom/data/HOMER/ARD_KO_atac_combined.bed mm10 -gtf  /Users/vyom/data/reversal_scATAC/genes.gtf -annStats combined_stats.txt > ARD_KO_Homer_combined.txt
Homer_peaks <- read.delim("~/data/HOMER/ARD_WT_Homer_combined.txt")
Homer_peaks$peak <- paste0(Homer_peaks$Chr,'-',Homer_peaks$Start,'-', Homer_peaks$End)
Homer_peaks$Annotation_simple <- gsub("([A-Za-z]+).*", "\\1", Homer_peaks$Annotation)
i = 'All'
celltypes <- levels(as.factor(ARDKO_atac_obj$Cell_Type))
for(i in celltypes) 
{
  
  DE_subset_ARD_KO <- read.csv(paste0('~/Stem 2_ARD_KO_ATAC_DE.csv'))
  DE_subset_ARD_KO <- DE_subset_ARD_KO[DE_subset_ARD_KO$p_val < 0.05 ,]
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
  ggsave(file = paste0(i,'sARD_KO_scATAC_per_accessible_annotations_sum_new.pdf'), width=3, height=3, units="in")
  
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
  
  if (max(-log(DE_subset_ARD_KO$p_val_adj)) < 321){
    ylim_num <- max(-log(DE_subset_ARD_KO$p_val_adj))} else {
      ylim_num <- 320}
  EnhancedVolcano(DE_subset_ARD_KO, lab = DE_subset_ARD_KO$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'ARD_KO vs Control',
                  pCutoff = .05, FCcutoff = 0.25, xlim = c(min(DE_subset_ARD_KO$avg_log2FC)-.05,max(DE_subset_ARD_KO$avg_log2FC)+.05),ylim = c(0,ylim_num+5),
                  subtitle = 'ALL CELLS', gridlines.minor = FALSE, gridlines.major = FALSE)
  ggsave(file = paste0(i,'_DE_peak_volcano_promoters.pdf'), width=6, height=6, units="in")
  
}

ARDKO_atac_obj$Treatment_celltype <- paste0(ARDKO_atac_obj$Cell_Type, "_", ARDKO_atac_obj$Treatment)
levels(as.factor(ARDKO_atac_obj$Treatment_celltype))
levels_idents <- c("Stem 1_ARD1", "Stem 1_ARD", "Stem 1_KO", "Stem 1_KO_ARD", 
                   "Stem 2_ARD1", "Stem 2_ARD", "Stem 2_KO", "Stem 2_KO_ARD", 
                   "TA_ARD1", "TA_ARD", "TA_KO", "TA_KO_ARD", 
                   "EP_ARD1",  "EP_ARD", "EP_KO", "EP_KO_ARD", 
                   "Ent._ARD1", "Ent._ARD", "Ent._KO", "Ent._KO_ARD", 
                   "Entend._ARD1", "Entend._ARD", "Entend._KO", "Entend._KO_ARD", 
                   "Gob._ARD1", "Gob._ARD", "Gob._KO",  "Gob._KO_ARD", 
                   "Tuft_ARD1", "Tuft_ARD", "Tuft_KO", "Tuft_KO_ARD", 
                   "T_ARD1",  "T_ARD", "T_KO", "T_KO_ARD",   
                   "Myeloid_ARD1", "Myeloid_ARD", "Myeloid_KO", "Myeloid_KO_ARD")

ARDKO_atac_obj$Treatment_celltype <- factor(x = ARDKO_atac_obj$Treatment_celltype, levels = levels_idents)

CoverageBrowser(ARDKO_atac_obj, region = 'S100a6')
CoveragePlot(
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
)

# gather the footprinting information for sets of motifs
ARDKO_atac_obj <- Footprint(
  object = ARDKO_atac_obj,
  motif.name = c('PPARA::RXRA','Pparg::Rxra',"PPARG", "PPARD", "CREB1"),
  genome = BSgenome.Mmusculus.UCSC.mm10,
  in.peaks = TRUE
)

#housekeeping
# saveRDS(ARDKO_atac_obj, file = "/Users/vyom/CSHL Dropbox Team Dropbox/Vyom Shah/Seurat_Objects/ARD_YAPKO_ATAC_YAP_fin_6_17_Sobj.rds")
# ARDKO_atac_obj <- readRDS("/Users/vyom/CSHL Dropbox Team Dropbox/Vyom Shah/Seurat_Objects/ARD_YAPKO_ATAC_YAP_fin_6_17_Sobj.rds", refhook = NULL)


