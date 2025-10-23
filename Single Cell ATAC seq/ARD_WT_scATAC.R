#Vyom Shah - Summer 2023
#AA Reversal project
#scATAC/scRNA analysis

library(Rmagic)
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
Control_counts <- Read10X_h5("/Users/vyom/data/reversal_scATAC/count/Beyaz_SB44_Control/outs/filtered_peak_bc_matrix.h5")
Control_fragpath <- "/Users/vyom/data/reversal_scATAC/count/Beyaz_SB44_Control/outs/fragments.tsv.gz"

ARD_counts <- Read10X_h5("/Users/vyom/data/reversal_scATAC/count/Beyaz_SB44_ARD/outs/filtered_peak_bc_matrix.h5")
ARD_fragpath <- "/Users/vyom/data/reversal_scATAC/count/Beyaz_SB44_ARD/outs/fragments.tsv.gz"

# read in peak sets
peaks.control <- read.table(
  file = "/Users/vyom/data/reversal_scATAC/count/Beyaz_SB44_Control/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.ard <- read.table(
  file = "/Users/vyom/data/reversal_scATAC/count/Beyaz_SB44_ARD/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
gr.control <- makeGRangesFromDataFrame(peaks.control)
gr.ard <- makeGRangesFromDataFrame(peaks.ard)

combined.peaks <- reduce(x = c(gr.control, gr.ard))

peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

#create fragment objects
#load metadata

md.Control <- read.table(
  file = "/Users/vyom/data/reversal_scATAC/count/Beyaz_SB44_Control/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.Ard <- read.table(
  file = "/Users/vyom/data/reversal_scATAC/count/Beyaz_SB44_ARD/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]


boxplot(md.Ard$passed_filters)
md.Control <- md.Control[md.Control$passed_filters > 50, ]
md.Ard <- md.Ard[md.Ard$passed_filters > 35, ]

frags.Control <- CreateFragmentObject(
  path = Control_fragpath,
  cells = rownames(md.Control)
)

frags.Ard <- CreateFragmentObject(
  path = ARD_fragpath,
  cells = rownames(md.Ard)
)

Control.counts <- FeatureMatrix(
  fragments = frags.Control,
  features = combined.peaks,
  cells = rownames(md.Control)
)

Ard.counts <- FeatureMatrix(
  fragments = frags.Ard,
  features = combined.peaks,
  cells = rownames(md.Ard)
)

Control_assay <- CreateChromatinAssay(Control.counts, fragments = frags.Control)
Control_obj <- CreateSeuratObject(Control_assay, assay = "ATAC", meta.data=md.Control)

ARD_assay <- CreateChromatinAssay(Ard.counts, fragments = frags.Ard)
ARD_obj <- CreateSeuratObject(ARD_assay, assay = "ATAC", meta.data=md.Ard)

Control_obj$dataset <- 'Control'
ARD_obj$dataset <-  'ARD'

ARD_obj <- merge(
  x = Control_obj,
  y = list(ARD_obj, REV_obj),
  add.cell.ids = c("Control", "ARD")
)
ARD_obj[["ATAC"]]

CoveragePlot(
  object = ARD_obj,
  group.by = 'dataset',
  region = ""
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
# add the gene information to the object
Annotation(ARD_obj) <- annotations


# compute nucleosome signal score per cell
ARD_obj <- NucleosomeSignal(object = ARD_obj)

# compute TSS enrichment score per cell
ARD_obj <- TSSEnrichment(object = ARD_obj, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
ARD_obj$pct_reads_in_peaks <- ARD_obj$peak_region_fragments / ARD_obj$passed_filters * 100
ARD_obj$blacklist_ratio <- ARD_obj$blacklist_region_fragments / ARD_obj$peak_region_fragments

ARD_obj$high.tss <- ifelse(ARD_obj$TSS.enrichment > 3, 'High', 'Low')
TSSPlot(ARD_obj, group.by = 'dataset') + NoLegend()

VlnPlot(
  object = ARD_obj,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.0,
  ncol = 4,
  group.by = 'dataset'
)
ARD_obj <- subset(
  x = ARD_obj,
  subset = nCount_ATAC > 750 &
    nCount_ATAC < 15000 &
    pct_reads_in_peaks > 35 &
    nucleosome_signal < 3 &
    TSS.enrichment > 1 &
    TSS.enrichment < 10
)

ARD_obj <- RunTFIDF(ARD_obj)
ARD_obj <- FindTopFeatures(ARD_obj, min.cutoff = 5)
ARD_obj <- RunSVD(ARD_obj)
DepthCor(ARD_obj)
ARD_obj@assays$ATAC
library(harmony)

ARD_obj$Treatment <- ARD_obj$dataset
ARD_obj <- RunHarmony(
  object = ARD_obj,
  group.by.vars = 'Treatment',
  reduction = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE,
  epsilon.harmony = 1e-10,
  epsilon.cluster = 1e-10,
  plot_convergence = TRUE
  
)

ARD_obj <- RunUMAP(ARD_obj, dims = 1:50,n.neighbors = 50,min.dist = .3,n.components = 2,reduction = 'harmony')

ARD_obj <- FindNeighbors(object = ARD_obj, reduction = 'lsi', dims = 1:25)
ARD_obj <- FindClusters(object = ARD_obj, verbose = FALSE, algorithm = 2, resolution = .5)
DimPlot(object = ARD_obj, label = TRUE) 


Prop_table <- prop.table(x = table(ARD_obj$Treatment, ARD_obj$seurat_clusters))

DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3","S100a6","Ly6a","Anxa3", "Areg","Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a",'Tmigd1', 'Fabp6', 'Slc51b', 'Slc51a', "Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Lyz1","Defa17","Defa24","Ang4","Pou2f3","Avil","Tuba1a","Adh1",'Ptprc','Cd8a','Cd4','Cd3e','Cd3d') 
DotPlot_Sig <- unique(c('Ptprc',"S100a6","Ly6a","Anxa3", "Areg",'Cd3g','Cd3e','Cd8a','Cd4','Trac','Tcrg-C1','Lag3','Pdcd1','Havcr2','Tox','Tcf7','Gzmb','Tbx21','Ifng','Gata3','Il4','Rorc','Il17a','Il17f','Foxp3','Il10','Il2rb','Il2ra','Klrd1','Cd19','Ighm','Ighg1','Cd74','Ciita','Nrc1','Klre1','Itgam','Itgax','H2-Eb1','H2-Ab1','Arg1','Mrc1','Tgfbi','Ccr2','Vegfa','Prdx1','Clec4d','Ccl5','Cd83','Ccr7','Fcn1','Msrb1','Ly6g','Col3a1','Sparc'))

DotPlot(ARD_obj, features = DotPlot_Sig, assay = 'RNA') + labs(y= "Cell Type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1)) +
  theme(text = element_text(size=5), axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())

ARD_obj[["Cell_Type"]] <- Idents(ARD_obj)
new.cluster.ids <- c( "Stem 1", "Stem 2", "Transit Amplifying", "Enterocyte Progenitor", "Enterocyte", "Enteroendocrine", "Goblet", "Tuft", "Macrophage", "T cell", "T cell", "NK cell")
names(new.cluster.ids) <- levels(ARD_obj)
ARD_obj <- RenameIdents(ARD_obj, new.cluster.ids)
ARD_obj[["Cell_Type"]] <- Idents(ARD_obj)
DimPlot(ARD_obj, group.by = 'Cell_Type')
rev_filter_obj$Cell_Type <- factor(x = rev_filter_obj$Cell_Type, levels = new.cluster.ids)
rev_filter_obj$Treatmnet_celltype <- paste0(rev_filter_obj$Cell_Type, "_", rev_filter_obj$Treatment)
levels(as.factor(rev_filter_obj$Treatmnet_celltype))
levels_idents <- c( "Stem1_Control","Stem1_ARD",  "Stem1_Reversal",
                    "Stem2_Control", "Stem2_ARD",  "Stem2_Reversal",
                    "TA_Control", "TA_ARD",  "TA_Reversal",
                    "EnterocyteProgenitor_Control", "EnterocyteProgenitor_ARD",  "EnterocyteProgenitor_Reversal",
                    "Enterocyte_Control", "Enterocyte_ARD",  "Enterocyte_Reversal",
                    "Goblet_Control", "Goblet_ARD",  "Goblet_Reversal",
                    "Tuft_Control","Tuft_ARD",  "Tuft_Reversal",
                    "Enteroendocrine_Control", "Enteroendocrine_ARD",  "Enteroendocrine_Reversal")
rev_filter_obj$Treatmnet_celltype <- factor(x = rev_filter_obj$Treatmnet_celltype, levels = levels_idents)

#DE per cell type
for(i in celltypes){
  subset_cell <- subset(ARD_obj,  idents = i)
  Idents(subset_cell) <- subset_cell$Treatment
  for(j in Treatment_groups) {
  DE_subset <- FindMarkers(subset_cell, ident.1 = j, ident.2 = "Control",  test.use = 'MAST', logfc.threshold = .1, min.pct = .01, latent.vars = 'nCount_ATAC', assay = 'ATAC')
  DE_subset$peaks <- rownames(DE_subset)
  Open_Regions <- rownames(DE_subset)
  Clostest_genes <- ClosestFeature(subset_cell, regions = Open_Regions)
  DE_subset$gene <-  Clostest_genes$gene_name
  write.csv(DE_subset, paste0(i,'_',j,'_ATAC_DE.csv')) 
  }}

#MOTIF and Annotations per cell type
#create gene activity matrix
gene.activities <- GeneActivity(ARD_obj)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
ARD_obj[['RNA']] <- CreateAssayObject(counts = gene.activities)
ARD_obj <- NormalizeData(
  object = ARD_obj,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(ARD_obj$nCount_RNA)
)

ararev.bed <- granges(ARD_obj)
df <- data.frame(seqnames=seqnames(ararev.bed),
                 starts=start(ararev.bed)-1,
                 ends=end(ararev.bed),
                 names=c(rep(".", length(ararev.bed))),
                 scores=c(rep(".", length(ararev.bed))),
                 strands=strand(ararev.bed))
write.table(df, file="atac_combined.bed", quote=F, sep="\t", row.names=F, col.names=F)
#homer code for command line: annotatePeaks.pl /Users/vyom/data/HOMER/rev_atac_combined.bed mm10 -gtf  /Users/vyom/data/reversal_scATAC/genes.gtf -annStats combined_stats.txt > Homer_combined.txt


Homer_peaks <- read.delim("~/data/HOMER/Homer_combined.txt")
Homer_peaks$peak <- paste0(Homer_peaks$Chr,'-',Homer_peaks$Start,'-', Homer_peaks$End)
Homer_peaks$Annotation_simple <- gsub("([A-Za-z]+).*", "\\1", Homer_peaks$Annotation)
i = 'All'
for(i in celltypes){
  
  DE_subset_P300 <- read.csv(paste0(i,'_ARD_ATAC_DE.csv'))

  intersect(DE_subset_P300$X, Homer_peaks$peak)
  
  Up_matched_P300 <- intersect(DE_subset_P300[DE_subset_P300$avg_log2FC > 0 & DE_subset_P300$p_val < 0.05 ,]$X, Homer_peaks$peak)
  Down_matched_P300 <- intersect(DE_subset_P300[DE_subset_P300$avg_log2FC < 0 & DE_subset_P300$p_val < 0.05 ,]$X, Homer_peaks$peak)
  
  cat_anno <- c('exon', 'Intergenic', 'intron', 'promoter', 'TTS')
  unique(Homer_peaks[Homer_peaks$peak %in% Down_matched_P300,]$Annotation_simple)
  unique(Homer_peaks[Homer_peaks$peak %in% Up_matched_P300,]$Annotation_simple)
  down_P300 <- Homer_peaks[Homer_peaks$peak %in% Down_matched_P300,]$Annotation_simple %>% factor(levels = cat_anno) %>% tabulate()
  up_P300 <- Homer_peaks[Homer_peaks$peak %in% Up_matched_P300,]$Annotation_simple %>% factor(levels = cat_anno) %>% tabulate()
  if (length(up_P300) < 5) {
    up_P300 <- c(up_P300,0)
  }
  if (length(down_P300) < 5) {
    down_P300 <- c(down_P300,0)
  }
  annotations <- Homer_peaks[Homer_peaks$peak %in% Up_matched_P300,]$Annotation_simple %>% as.factor() %>% levels()
  Annotation_homer <- data.frame(down_P300, up_P300)
  Annotation_homer_prop <- Annotation_homer %>% mutate(across(where(is.numeric), ~ ./sum(.)))
  Annotation_homer_prop <- cbind(annotations,Annotation_homer_prop)
  
  Annotation_homer_prop <- melt(Annotation_homer_prop)
  Annotation_homer_prop$annotations <- factor(Annotation_homer_prop$annotations,levels = annotations)
  ggplot(Annotation_homer_prop, aes(x=variable, y= value, fill = annotations))+theme_vyom  + geom_col(colour = "black", position = "fill") + scale_y_continuous(labels = scales::percent) + scale_fill_brewer(palette = "Pastel2")
  ggsave(file = paste0(i,'_scATAC_ARD_per_accessible_annotations_proportion.pdf'), width=3, height=3, units="in")
  
  Annotation_homer_sum <- cbind(annotations,Annotation_homer)
  Annotation_homer_sum <- melt(Annotation_homer_sum)
  Annotation_homer_sum$annotations <- factor(Annotation_homer_sum$annotations,levels = annotations)
  ggplot(Annotation_homer_sum, aes(x=variable, y= value, fill = annotations))+theme_vyom   + geom_col(colour = "black")  + scale_fill_brewer(palette = "Pastel2")
  ggsave(file = paste0(i,'_ARD_scATAC_per_accessible_annotations_sum.pdf'), width=3, height=3, units="in")
  
  Up_P300_PEAKS_promoter <- Homer_peaks[Homer_peaks$peak %in% Up_matched_P300 & Homer_peaks$Annotation_simple == 'promoter',]$peak
  Down_P300_PEAKS_promoter <- Homer_peaks[Homer_peaks$peak %in% Down_matched_P300 & Homer_peaks$Annotation_simple == 'promoter',]$peak
  
  motifs_P300_up <- FindMotifs(
    object = ARD_obj,
    features = Up_P300_PEAKS_promoter)
  motifs_P300_up <- motifs_P300_up[motifs_P300_up$pvalue < .05,]
  
  write.csv(motifs_P300_up, paste0(i,'_promoter_ARD_up_motifs.csv'))
  
  motifs_P300_up <- FindMotifs(
    object = ARD_obj,
    features = Down_P300_PEAKS_promoter)
  motifs_P300_up <- motifs_P300_up[motifs_P300_up$pvalue < .05,]
  
  write.csv(motifs_P300_up, paste0(i,'_promoter_ARD_down_motifs.csv'))
  
  Up_P300_genes_promoter <- DE_subset_P300[DE_subset_P300$X %in% Homer_peaks[Homer_peaks$peak %in% Homer_peaks$Annotation_simple == 'promoter',]$peak,]
  Up_P300_genes_promoter$peak <- Up_P300_genes_promoter$X
  Up_P300_genes_promoter1 <- merge(Up_P300_genes_promoter, Homer_peaks[Homer_peaks$peak %in% Up_P300_genes_promoter$peak,] %>% dplyr::select(Gene.Name,peak), by= 'peak', all.x=TRUE)
  
  write.csv(Up_P300_genes_promoter, paste0(i,'_promoter_ARD_up_PEAKS.csv'))
  
  matched_ard <- intersect(DE_subset_P300[DE_subset_P300$p_val_adj < 0.05,]$X, Homer_peaks$peak)
  DE_subset_P3001 <- DE_subset_P300[DE_subset_P300$X %in% Homer_peaks[Homer_peaks$peak %in% matched_ard,]$peak,]
  homer_subset <- Homer_peaks[Homer_peaks$peak %in% DE_subset_P3001$X,] %>% dplyr::select(peak, Annotation_simple, Gene.Name )
  DE_subset_P3001$peak <- DE_subset_P3001$X
  DE_subset_P3001 <- merge(DE_subset_P3001, homer_subset, by="peak")
  write.csv(DE_subset_P3001, paste0(i, '_ARD_promoters_DE_peak.csv'))
  
  if (max(-log(DE_subset_P300$p_val_adj)) < 321){
    ylim_num <- max(-log(DE_subset_P300$p_val_adj))} else {
      ylim_num <- 320}
  EnhancedVolcano(DE_subset_P300, lab = DE_subset_P300$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'REV vs WT',
                  pCutoff = .05, FCcutoff = 0.1, xlim = c(min(DE_subset_P300$avg_log2FC)-.05,max(DE_subset_P300$avg_log2FC)+.05),ylim = c(0,ylim_num+5),
                  subtitle = i, gridlines.minor = FALSE, gridlines.major = FALSE)
  ggsave(file = paste0(i,'_ARD_DE_peak_volcano_promoters.pdf'), width=6, height=6, units="in")
  
}




#motif analysis
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)
#library(JASPAR2022)
library(TFBSTools)
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
library(BSgenome.Mmusculus.UCSC.mm10)

main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
keep.peaks <- as.logical(seqnames(granges(ARD_obj)) %in% main.chroms)
ARD_obj <- ARD_obj[keep.peaks, ]

ARD_obj <- AddMotifs(
  object = ARD_obj,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)
Idents(ARD_obj) <- ARD_obj$Treatment
DE_output <- FindMarkers(ARD_obj, ident.1 = "Control", test.use = 'LR', logfc.threshold = .1, min.pct = .01, latent.vars = 'nCount_ATAC', assay = 'ATAC')

enriched.motifs <- FindMotifs(
  object = ARD_obj,
  features = rownames(DE_output[DE_output$p_val_adj < 0.05 & DE_output$avg_log2FC > 0, ])
)

#dive into stem2 specific
Idents(ARD_obj) <- ARD_obj$Cell_Type
subset_obj <- subset(ARD_obj, idents = c('Stem2'))

Idents(subset_obj) <- subset_obj$Treatment
DE_ARD_peaklist_stem2 <- FindMarkers(subset_obj, ident.1 = "ARD", ident.2 = "Control",  test.use = 'LR', logfc.threshold = .1, min.pct = .01, latent.vars = 'nCount_ATAC', assay = 'ATAC')

#DE_REV_peaklist1 <- DE_REV_peaklist
#DE_ARD_peaklist1 <- DE_ARD_peaklist
rownames(ARD_obj)
Atac_upset_Stem2 <- list()
Atac_upset_Stem2$Ard_up <- rownames(DE_ARD_peaklist_stem2[DE_ARD_peaklist_stem2$p_val < 0.05 & DE_ARD_peaklist_stem2$avg_log2FC > 0, ])
#Atac_upset$Ard_ns <- rownames(ARD_obj)[!(rownames(ARD_obj) %in% rownames(DE_ARD_peaklist[DE_ARD_peaklist$p_val_adj < 0.01, ]))]
Atac_upset_Stem2$Ard_down <- rownames(DE_ARD_peaklist_stem2[DE_ARD_peaklist_stem2$p_val < 0.05 & DE_ARD_peaklist_stem2$avg_log2FC < 0, ])

#housekeeping
#ararev_obj <- readRDS('./data/Seurat_Objects/reversal_scrna.Rds', refhook = NULL) # this is santhi's analysis scRNA
# saveRDS(ARD_obj, file = "./data/Seurat_Objects/Ard_Rev_ATAC_Sobj.rds")
# ARD_obj <- readRDS("/Users/vyom/CSHL Dropbox Team Dropbox/Vyom Shah/Seurat_Objects/Ard_Rev_ATAC_Sobj.rds", refhook = NULL)


