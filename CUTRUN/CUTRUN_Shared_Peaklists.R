#find shared peaks
library(DiffBind)
library(tidyverse)
library(edgeR)
library(DESeq2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPpeakAnno)
library(org.Mm.eg.db)
library(GenomicRanges)
library(ChIPseeker)
library(clusterProfiler)
library(ReactomePA)
library(biomaRt)
library(ChIPpeakAnno)
library(readxl)
#PAN-TEAD Target Definitions
{
#load pan-TEAD data
m4_ard_pantead_gr <- toGRanges('/Users/vyom/cloud_computing/home/vshah/analysis/YAP_CUTRUN/results/03_peak_calling/05_consensus_peaks/M4-ARD-panTEAD.consensus.peaks.awk.bed', format="narrowPeak", header=FALSE)
m4_control_pantead_gr <- toGRanges('/Users/vyom/cloud_computing/home/vshah/analysis/YAP_CUTRUN/results/03_peak_calling/05_consensus_peaks/M4-Control-panTEAD.consensus.peaks.awk.bed', format="narrowPeak", header=FALSE)

m3_ard_pantead_gr <- toGRanges('/Users/vyom/cloud_computing/home/vshah/analysis/YAP_CUTRUN/results/03_peak_calling/05_consensus_peaks/M3-ARD-panTEAD.consensus.peaks.awk.bed', format="narrowPeak", header=FALSE)
m3_control_pantead_gr <- toGRanges('/Users/vyom/cloud_computing/home/vshah/analysis/YAP_CUTRUN/results/03_peak_calling/05_consensus_peaks/M3-Control-panTEAD.consensus.peaks.awk.bed', format="narrowPeak", header=FALSE)

m4_ard_pantead_gr$pValue <- NULL
m4_ard_pantead_gr$qValue <- NULL

m4_control_pantead_gr$pValue <- NULL
m4_control_pantead_gr$qValue <- NULL

m3_ard_pantead_gr$pValue <- NULL
m3_ard_pantead_gr$qValue <- NULL

m3_control_pantead_gr$pValue <- NULL
m3_control_pantead_gr$qValue <- NULL

data(TSS.mouse.GRCm38)
m3_ard_pantead_anno <- annotatePeakInBatch(m3_ard_pantead_gr, AnnotationData=TSS.mouse.GRCm38)

m3_ard_pantead_anno <- addGeneIDs(annotatedPeak=m3_ard_pantead_anno,
                        orgAnn="org.Mm.eg.db",
                        IDs2Add="symbol")
m3_ard_pantead_anno <- data.frame(m3_ard_pantead_anno)

m4_ard_pantead_anno <- annotatePeakInBatch(m4_ard_pantead_gr, AnnotationData=TSS.mouse.GRCm38)

m4_ard_pantead_anno <- addGeneIDs(annotatedPeak=m4_ard_pantead_anno,
                                  orgAnn="org.Mm.eg.db",
                                  IDs2Add="symbol")
m4_ard_pantead_anno <- data.frame(m4_ard_pantead_anno)

m4_control_pantead_anno <- annotatePeakInBatch(m4_control_pantead_gr, AnnotationData=TSS.mouse.GRCm38)

m4_control_pantead_anno <- addGeneIDs(annotatedPeak=m4_control_pantead_anno,
                                  orgAnn="org.Mm.eg.db",
                                  IDs2Add="symbol")
m4_control_pantead_anno <- data.frame(m4_control_pantead_anno)

m3_control_pantead_anno <- annotatePeakInBatch(m3_control_pantead_gr, AnnotationData=TSS.mouse.GRCm38)

m3_control_pantead_anno <- addGeneIDs(annotatedPeak=m3_control_pantead_anno,
                                  orgAnn="org.Mm.eg.db",
                                  IDs2Add="symbol")
m3_control_pantead_anno <- data.frame(m3_control_pantead_anno)

m3_control_pantead_anno <- m3_control_pantead_anno[complete.cases(m3_control_pantead_anno),]
m4_control_pantead_anno <- m4_control_pantead_anno[complete.cases(m4_control_pantead_anno),]
m4_ard_pantead_anno <- m4_ard_pantead_anno[complete.cases(m4_ard_pantead_anno),]
m3_ard_pantead_anno <- m3_ard_pantead_anno[complete.cases(m3_ard_pantead_anno),]
unique(m3_ard_pantead_anno$symbol)
#library(BRGenomics)
#pan_tead_peaks <- mergeGRangesData(m4_control_pantead_gr, m4_ard_pantead_gr, m3_ard_pantead_gr, m3_control_pantead_gr , field = "peak", ncores = 1)
m3_ard_pantead_anno[m3_ard_pantead_anno$insideFeature %in% c("upstream",'overlapStart') & m3_ard_pantead_anno$shortestDistance < 500, ]$symbol
#load CREB data
Pan_tead_broad_peak_genes <- unique(c(m3_ard_pantead_anno$symbol, m4_ard_pantead_anno$symbol,  m4_control_pantead_anno$symbol, m3_control_pantead_anno$symbol))
Pan_tead_broad_peak_genes <- unique(c(m3_ard_pantead_anno[m3_ard_pantead_anno$insideFeature == "upstream", ]$symbol,
                                      m4_ard_pantead_anno[m4_ard_pantead_anno$insideFeature == "upstream", ]$symbol,
                                      m4_control_pantead_anno[m4_control_pantead_anno$insideFeature == "upstream", ]$symbol,
                                      m3_control_pantead_anno[m3_control_pantead_anno$insideFeature == "upstream", ]$symbol))

Pan_tead_broad_peak_genes_promoter <- unique(c(m3_ard_pantead_anno[m3_ard_pantead_anno$insideFeature %in% c("upstream",'overlapStart', 'inside') & m3_ard_pantead_anno$shortestDistance < 1000 & m3_ard_pantead_anno$signalValue > 20, ]$symbol,
                                      m4_ard_pantead_anno[m4_ard_pantead_anno$insideFeature %in% c("upstream",'overlapStart', 'inside') & m4_ard_pantead_anno$shortestDistance < 1000 & m4_ard_pantead_anno$signalValue > 20, ]$symbol,
                                      m4_control_pantead_anno[m4_control_pantead_anno$insideFeature %in% c("upstream",'overlapStart', 'inside') & m4_control_pantead_anno$shortestDistance < 1000 & m4_control_pantead_anno$signalValue > 20, ]$symbol,
                                      m3_control_pantead_anno[m3_control_pantead_anno$insideFeature %in% c("upstream",'overlapStart', 'inside') & m3_control_pantead_anno$shortestDistance < 1000 & m3_control_pantead_anno$signalValue > 20, ]$symbol))

Pan_tead_broad_peak_genes_enhancer <- unique(c(m3_ard_pantead_anno[!(m3_ard_pantead_anno$insideFeature %in% c("upstream",'overlapStart', 'inside') & m3_ard_pantead_anno$shortestDistance < 1000 )& m3_ard_pantead_anno$signalValue > 20, ]$symbol,
                                      m4_ard_pantead_anno[!(m4_ard_pantead_anno$insideFeature %in% c("upstream",'overlapStart', 'inside') & m4_ard_pantead_anno$shortestDistance < 1000 )& m4_ard_pantead_anno$signalValue > 20, ]$symbol,
                                      m4_control_pantead_anno[!(m4_control_pantead_anno$insideFeature %in% c("upstream",'overlapStart', 'inside') & m4_control_pantead_anno$shortestDistance < 1000 )& m4_control_pantead_anno$signalValue > 20, ]$symbol,
                                      m3_control_pantead_anno[!(m3_control_pantead_anno$insideFeature %in% c("upstream",'overlapStart', 'inside') & m3_control_pantead_anno$shortestDistance < 1000 )& m3_control_pantead_anno$signalValue > 20, ]$symbol))


YAP_Targets <- read_excel("Users/vyom/analysis/Finalized_Signatures.xlsx", sheet = "Yap_Targets")
Yap_refine <- intersect(YAP_Targets$Gene_Name, Pan_tead_broad_peak_genes)
write.csv(Pan_tead_broad_peak_genes_promoter,'/Users/vyom/YAP_Peaks_Promoters.csv' )
write.csv(Pan_tead_broad_peak_genes_enhancer,'/Users/vyom/YAP_Peaks_Enhancers.csv' )
}

#Creb Target Definitions
{
  #load pan-TEAD data
  m1_ard_CREB_gr <- toGRanges('/Users/vyom/cloud_computing/home/vshah/analysis/NAT_ARD_CUTRUN/results/03_peak_calling/05_consensus_peaks/M1-ARD-CREB.consensus.peaks.awk.bed', format="narrowPeak", header=FALSE)
  m1_control_CREB_gr <- toGRanges('/Users/vyom/cloud_computing/home/vshah/analysis/NAT_ARD_CUTRUN/results/03_peak_calling/05_consensus_peaks/M1-Control-CREB.consensus.peaks.awk.bed', format="narrowPeak", header=FALSE)
  
  m2_ard_CREB_gr <- toGRanges('/Users/vyom/cloud_computing/home/vshah/analysis/NAT_ARD_CUTRUN/results/03_peak_calling/05_consensus_peaks/M2-ARD-CREB.consensus.peaks.awk.bed', format="narrowPeak", header=FALSE)
  m2_control_CREB_gr <- toGRanges('/Users/vyom/cloud_computing/home/vshah/analysis/NAT_ARD_CUTRUN/results/03_peak_calling/05_consensus_peaks/M2-Control-CREB.consensus.peaks.awk.bed', format="narrowPeak", header=FALSE)
  
  m1_ard_CREB_gr$pValue <- NULL
  m1_ard_CREB_gr$qValue <- NULL
  
  m1_control_CREB_gr$pValue <- NULL
  m1_control_CREB_gr$qValue <- NULL
  
  m2_ard_CREB_gr$pValue <- NULL
  m2_ard_CREB_gr$qValue <- NULL
  
  m2_control_CREB_gr$pValue <- NULL
  m2_control_CREB_gr$qValue <- NULL
  
  data(TSS.mouse.GRCm38)
  m2_ard_CREB_anno <- annotatePeakInBatch(m2_ard_CREB_gr, AnnotationData=TSS.mouse.GRCm38)
  
  m2_ard_CREB_anno <- addGeneIDs(annotatedPeak=m2_ard_CREB_anno,
                                    orgAnn="org.Mm.eg.db",
                                    IDs2Add="symbol")
  m2_ard_CREB_anno <- data.frame(m2_ard_CREB_anno)
  
  m1_ard_CREB_anno <- annotatePeakInBatch(m1_ard_CREB_gr, AnnotationData=TSS.mouse.GRCm38)
  
  m1_ard_CREB_anno <- addGeneIDs(annotatedPeak=m1_ard_CREB_anno,
                                    orgAnn="org.Mm.eg.db",
                                    IDs2Add="symbol")
  m1_ard_CREB_anno <- data.frame(m1_ard_CREB_anno)
  
  m1_control_CREB_anno <- annotatePeakInBatch(m1_control_CREB_gr, AnnotationData=TSS.mouse.GRCm38)
  
  m1_control_CREB_anno <- addGeneIDs(annotatedPeak=m1_control_CREB_anno,
                                        orgAnn="org.Mm.eg.db",
                                        IDs2Add="symbol")
  m1_control_CREB_anno <- data.frame(m1_control_CREB_anno)
  
  m2_control_CREB_anno <- annotatePeakInBatch(m2_control_CREB_gr, AnnotationData=TSS.mouse.GRCm38)
  
  m2_control_CREB_anno <- addGeneIDs(annotatedPeak=m2_control_CREB_anno,
                                        orgAnn="org.Mm.eg.db",
                                        IDs2Add="symbol")
  m2_control_CREB_anno <- data.frame(m2_control_CREB_anno)
  
  m2_control_CREB_anno <- m2_control_CREB_anno[complete.cases(m2_control_CREB_anno),]
  m1_control_CREB_anno <- m1_control_CREB_anno[complete.cases(m1_control_CREB_anno),]
  m1_ard_CREB_anno <- m1_ard_CREB_anno[complete.cases(m1_ard_CREB_anno),]
  m2_ard_CREB_anno <- m2_ard_CREB_anno[complete.cases(m2_ard_CREB_anno),]
  unique(m2_ard_CREB_anno$symbol)
  #library(BRGenomics)
  #Creb_peaks <- mergeGRangesData(m1_control_CREB_gr, m1_ard_CREB_gr, m2_ard_CREB_gr, m2_control_CREB_gr , field = "peak", ncores = 1)
  m2_ard_CREB_anno[m2_ard_CREB_anno$insideFeature %in% c("upstream",'overlapStart') & m2_ard_CREB_anno$shortestDistance < 500, ]$symbol
  #load CREB data
  Creb_broad_peak_genes <- unique(c(m2_ard_CREB_anno$symbol, m1_ard_CREB_anno$symbol,  m1_control_CREB_anno$symbol, m2_control_CREB_anno$symbol))
  Creb_broad_peak_genes <- unique(c(m2_ard_CREB_anno[m2_ard_CREB_anno$insideFeature == "upstream", ]$symbol,
                                        m1_ard_CREB_anno[m1_ard_CREB_anno$insideFeature == "upstream", ]$symbol,
                                        m1_control_CREB_anno[m1_control_CREB_anno$insideFeature == "upstream", ]$symbol,
                                        m2_control_CREB_anno[m2_control_CREB_anno$insideFeature == "upstream", ]$symbol))
  
  Creb_broad_peak_genes_promoter <- unique(c(m2_ard_CREB_anno[m2_ard_CREB_anno$insideFeature %in% c("upstream",'overlapStart', 'inside') & m2_ard_CREB_anno$shortestDistance < 1000 & m2_ard_CREB_anno$signalValue > 20, ]$symbol,
                                                 m1_ard_CREB_anno[m1_ard_CREB_anno$insideFeature %in% c("upstream",'overlapStart', 'inside') & m1_ard_CREB_anno$shortestDistance < 1000 & m1_ard_CREB_anno$signalValue > 20, ]$symbol,
                                                 m1_control_CREB_anno[m1_control_CREB_anno$insideFeature %in% c("upstream",'overlapStart', 'inside') & m1_control_CREB_anno$shortestDistance < 1000 & m1_control_CREB_anno$signalValue > 20, ]$symbol,
                                                 m2_control_CREB_anno[m2_control_CREB_anno$insideFeature %in% c("upstream",'overlapStart', 'inside') & m2_control_CREB_anno$shortestDistance < 1000 & m2_control_CREB_anno$signalValue > 20, ]$symbol))
  
  Creb_broad_peak_genes_enhancer <- unique(c(m2_ard_CREB_anno[!(m2_ard_CREB_anno$insideFeature %in% c("upstream",'overlapStart', 'inside') & m2_ard_CREB_anno$shortestDistance < 1000 )& m2_ard_CREB_anno$signalValue > 20, ]$symbol,
                                                 m1_ard_CREB_anno[!(m1_ard_CREB_anno$insideFeature %in% c("upstream",'overlapStart','inside') & m1_ard_CREB_anno$shortestDistance < 1000 )& m1_ard_CREB_anno$signalValue > 20, ]$symbol,
                                                 m1_control_CREB_anno[!(m1_control_CREB_anno$insideFeature %in% c("upstream",'overlapStart', 'inside') & m1_control_CREB_anno$shortestDistance < 1000 )& m1_control_CREB_anno$signalValue > 20, ]$symbol,
                                                 m2_control_CREB_anno[!(m2_control_CREB_anno$insideFeature %in% c("upstream",'overlapStart','inside') & m2_control_CREB_anno$shortestDistance < 1000 )& m2_control_CREB_anno$signalValue > 20, ]$symbol))
  
  
  write.csv(Creb_broad_peak_genes_promoter,'/Users/vyom/CREB_Peaks_Promoters.csv' )
  write.csv(Creb_broad_peak_genes_enhancer,'/Users/vyom/CREB_Peaks_Enhancers.csv' )
}

#Enhancer Target Definition:
{
  #load pan-TEAD data
  m4_ard_H3K4me1_gr <- toGRanges('/Users/vyom/cloud_computing/home/vshah/analysis/YAP_CUTRUN/results/03_peak_calling/05_consensus_peaks/M4-ARD-H3K4me1.consensus.peaks.awk.bed', format="narrowPeak", header=FALSE)
  m4_control_H3K4me1_gr <- toGRanges('/Users/vyom/cloud_computing/home/vshah/analysis/YAP_CUTRUN/results/03_peak_calling/05_consensus_peaks/M4-Control-H3K4me1.consensus.peaks.awk.bed', format="narrowPeak", header=FALSE)
  
  m3_ard_H3K4me1_gr <- toGRanges('/Users/vyom/cloud_computing/home/vshah/analysis/YAP_CUTRUN/results/03_peak_calling/05_consensus_peaks/M3-ARD-H3K4me1.consensus.peaks.awk.bed', format="narrowPeak", header=FALSE)
  m3_control_H3K4me1_gr <- toGRanges('/Users/vyom/cloud_computing/home/vshah/analysis/YAP_CUTRUN/results/03_peak_calling/05_consensus_peaks/M3-Control-H3K4me1.consensus.peaks.awk.bed', format="narrowPeak", header=FALSE)
  
  m4_ard_H3K4me1_gr$pValue <- NULL
  m4_ard_H3K4me1_gr$qValue <- NULL
  
  m4_control_H3K4me1_gr$pValue <- NULL
  m4_control_H3K4me1_gr$qValue <- NULL
  
  m3_ard_H3K4me1_gr$pValue <- NULL
  m3_ard_H3K4me1_gr$qValue <- NULL
  
  m3_control_H3K4me1_gr$pValue <- NULL
  m3_control_H3K4me1_gr$qValue <- NULL
  
  data(TSS.mouse.GRCm38)
  m3_ard_H3K4me1_anno <- annotatePeakInBatch(m3_ard_H3K4me1_gr, AnnotationData=TSS.mouse.GRCm38)
  
  m3_ard_H3K4me1_anno <- addGeneIDs(annotatedPeak=m3_ard_H3K4me1_anno,
                                    orgAnn="org.Mm.eg.db",
                                    IDs2Add="symbol")
  m3_ard_H3K4me1_anno <- data.frame(m3_ard_H3K4me1_anno)
  
  m4_ard_H3K4me1_anno <- annotatePeakInBatch(m4_ard_H3K4me1_gr, AnnotationData=TSS.mouse.GRCm38)
  
  m4_ard_H3K4me1_anno <- addGeneIDs(annotatedPeak=m4_ard_H3K4me1_anno,
                                    orgAnn="org.Mm.eg.db",
                                    IDs2Add="symbol")
  m4_ard_H3K4me1_anno <- data.frame(m4_ard_H3K4me1_anno)
  
  m4_control_H3K4me1_anno <- annotatePeakInBatch(m4_control_H3K4me1_gr, AnnotationData=TSS.mouse.GRCm38)
  
  m4_control_H3K4me1_anno <- addGeneIDs(annotatedPeak=m4_control_H3K4me1_anno,
                                        orgAnn="org.Mm.eg.db",
                                        IDs2Add="symbol")
  m4_control_H3K4me1_anno <- data.frame(m4_control_H3K4me1_anno)
  
  m3_control_H3K4me1_anno <- annotatePeakInBatch(m3_control_H3K4me1_gr, AnnotationData=TSS.mouse.GRCm38)
  
  m3_control_H3K4me1_anno <- addGeneIDs(annotatedPeak=m3_control_H3K4me1_anno,
                                        orgAnn="org.Mm.eg.db",
                                        IDs2Add="symbol")
  m3_control_H3K4me1_anno <- data.frame(m3_control_H3K4me1_anno)
  
  m3_control_H3K4me1_anno <- m3_control_H3K4me1_anno[complete.cases(m3_control_H3K4me1_anno),]
  m4_control_H3K4me1_anno <- m4_control_H3K4me1_anno[complete.cases(m4_control_H3K4me1_anno),]
  m4_ard_H3K4me1_anno <- m4_ard_H3K4me1_anno[complete.cases(m4_ard_H3K4me1_anno),]
  m3_ard_H3K4me1_anno <- m3_ard_H3K4me1_anno[complete.cases(m3_ard_H3K4me1_anno),]
  unique(m3_ard_H3K4me1_anno$symbol)
  #library(BRGenomics)
  #H3K4me1_peaks <- mergeGRangesData(m4_control_H3K4me1_gr, m4_ard_H3K4me1_gr, m3_ard_H3K4me1_gr, m3_control_H3K4me1_gr , field = "peak", ncores = 1)
  m3_ard_H3K4me1_anno[m3_ard_H3K4me1_anno$insideFeature %in% c("upstream",'overlapStart') & m3_ard_H3K4me1_anno$shortestDistance < 500, ]$symbol
  #load CREB data
  #H3K4me1_broad_peak_genes <- unique(c(m3_ard_H3K4me1_anno$symbol, m4_ard_H3K4me1_anno$symbol,  m4_control_H3K4me1_anno$symbol, m3_control_H3K4me1_anno$symbol))
  H3K4me1_broad_peak_genes <- unique(c(m3_ard_H3K4me1_anno[m3_ard_H3K4me1_anno$signalValue > 20, ]$symbol,
                                        m4_ard_H3K4me1_anno[m4_ard_H3K4me1_anno$signalValue > 20, ]$symbol,
                                        m4_control_H3K4me1_anno[m4_control_H3K4me1_anno$signalValue > 20, ]$symbol,
                                        m3_control_H3K4me1_anno[m3_control_H3K4me1_anno$signalValue > 20, ]$symbol))
  
  H3K4me1_broad_peak_genes_promoter <- unique(c(m3_ard_H3K4me1_anno[m3_ard_H3K4me1_anno$insideFeature %in% c("upstream",'overlapStart') & m3_ard_H3K4me1_anno$shortestDistance < 1000 & m3_ard_H3K4me1_anno$signalValue > 20, ]$symbol,
                                                 m4_ard_H3K4me1_anno[m4_ard_H3K4me1_anno$insideFeature %in% c("upstream",'overlapStart') & m4_ard_H3K4me1_anno$shortestDistance < 1000 & m4_ard_H3K4me1_anno$signalValue > 20, ]$symbol,
                                                 m4_control_H3K4me1_anno[m4_control_H3K4me1_anno$insideFeature %in% c("upstream",'overlapStart') & m4_control_H3K4me1_anno$shortestDistance < 1000 & m4_control_H3K4me1_anno$signalValue > 20, ]$symbol,
                                                 m3_control_H3K4me1_anno[m3_control_H3K4me1_anno$insideFeature %in% c("upstream",'overlapStart') & m3_control_H3K4me1_anno$shortestDistance < 1000 & m3_control_H3K4me1_anno$signalValue > 20, ]$symbol))
  
  H3K4me1_broad_peak_genes_enhancer <- unique(c(m3_ard_H3K4me1_anno[!(m3_ard_H3K4me1_anno$insideFeature %in% c("upstream",'overlapStart') & m3_ard_H3K4me1_anno$shortestDistance < 1000 & m3_ard_H3K4me1_anno$signalValue > 20), ]$symbol,
                                                 m4_ard_H3K4me1_anno[!(m4_ard_H3K4me1_anno$insideFeature %in% c("upstream",'overlapStart') & m4_ard_H3K4me1_anno$shortestDistance < 1000 & m4_ard_H3K4me1_anno$signalValue > 20), ]$symbol,
                                                 m4_control_H3K4me1_anno[!(m4_control_H3K4me1_anno$insideFeature %in% c("upstream",'overlapStart') & m4_control_H3K4me1_anno$shortestDistance < 1000 & m4_control_H3K4me1_anno$signalValue > 20), ]$symbol,
                                                 m3_control_H3K4me1_anno[!(m3_control_H3K4me1_anno$insideFeature %in% c("upstream",'overlapStart') & m3_control_H3K4me1_anno$shortestDistance < 1000 & m3_control_H3K4me1_anno$signalValue > 20), ]$symbol))
  
  write.csv(H3K4me1_broad_peak_genes_enhancer,'/Users/vyom/H3K4ME3_Peaks_Enhancers.csv' )
}

#promoter target definition
{
  #load pan-TEAD data
  m1_ard_H3K4me3_gr <- toGRanges('/Users/vyom/cloud_computing/home/vshah/analysis/NAT_ARD_CUTRUN/results/03_peak_calling/05_consensus_peaks/M1-ARD-H3K4me3.consensus.peaks.awk.bed', format="narrowPeak", header=FALSE)
  m1_control_H3K4me3_gr <- toGRanges('/Users/vyom/cloud_computing/home/vshah/analysis/NAT_ARD_CUTRUN/results/03_peak_calling/05_consensus_peaks/M1-Control-H3K4me3.consensus.peaks.awk.bed', format="narrowPeak", header=FALSE)
  
  m2_ard_H3K4me3_gr <- toGRanges('/Users/vyom/cloud_computing/home/vshah/analysis/NAT_ARD_CUTRUN/results/03_peak_calling/05_consensus_peaks/M2-ARD-H3K4me3.consensus.peaks.awk.bed', format="narrowPeak", header=FALSE)
  m2_control_H3K4me3_gr <- toGRanges('/Users/vyom/cloud_computing/home/vshah/analysis/NAT_ARD_CUTRUN/results/03_peak_calling/05_consensus_peaks/M2-Control-H3K4me3.consensus.peaks.awk.bed', format="narrowPeak", header=FALSE)
  
  
  m1_ard_H3K4me3_gr$pValue <- NULL
  m1_ard_H3K4me3_gr$qValue <- NULL
  
  m1_control_H3K4me3_gr$pValue <- NULL
  m1_control_H3K4me3_gr$qValue <- NULL
  
  m2_ard_H3K4me3_gr$pValue <- NULL
  m2_ard_H3K4me3_gr$qValue <- NULL
  
  m2_control_H3K4me3_gr$pValue <- NULL
  m2_control_H3K4me3_gr$qValue <- NULL
  
  data(TSS.mouse.GRCm38)
  m2_ard_H3K4me3_anno <- annotatePeakInBatch(m2_ard_H3K4me3_gr, AnnotationData=TSS.mouse.GRCm38)
  
  m2_ard_H3K4me3_anno <- addGeneIDs(annotatedPeak=m2_ard_H3K4me3_anno,
                                    orgAnn="org.Mm.eg.db",
                                    IDs2Add="symbol")
  m2_ard_H3K4me3_anno <- data.frame(m2_ard_H3K4me3_anno)
  
  m1_ard_H3K4me3_anno <- annotatePeakInBatch(m1_ard_H3K4me3_gr, AnnotationData=TSS.mouse.GRCm38)
  
  m1_ard_H3K4me3_anno <- addGeneIDs(annotatedPeak=m1_ard_H3K4me3_anno,
                                    orgAnn="org.Mm.eg.db",
                                    IDs2Add="symbol")
  m1_ard_H3K4me3_anno <- data.frame(m1_ard_H3K4me3_anno)
  
  m1_control_H3K4me3_anno <- annotatePeakInBatch(m1_control_H3K4me3_gr, AnnotationData=TSS.mouse.GRCm38)
  
  m1_control_H3K4me3_anno <- addGeneIDs(annotatedPeak=m1_control_H3K4me3_anno,
                                        orgAnn="org.Mm.eg.db",
                                        IDs2Add="symbol")
  m1_control_H3K4me3_anno <- data.frame(m1_control_H3K4me3_anno)
  
  m2_control_H3K4me3_anno <- annotatePeakInBatch(m2_control_H3K4me3_gr, AnnotationData=TSS.mouse.GRCm38)
  
  m2_control_H3K4me3_anno <- addGeneIDs(annotatedPeak=m2_control_H3K4me3_anno,
                                        orgAnn="org.Mm.eg.db",
                                        IDs2Add="symbol")
  m2_control_H3K4me3_anno <- data.frame(m2_control_H3K4me3_anno)
  
  m2_control_H3K4me3_anno <- m2_control_H3K4me3_anno[complete.cases(m2_control_H3K4me3_anno),]
  m1_control_H3K4me3_anno <- m1_control_H3K4me3_anno[complete.cases(m1_control_H3K4me3_anno),]
  m1_ard_H3K4me3_anno <- m1_ard_H3K4me3_anno[complete.cases(m1_ard_H3K4me3_anno),]
  m2_ard_H3K4me3_anno <- m2_ard_H3K4me3_anno[complete.cases(m2_ard_H3K4me3_anno),]
  unique(m2_ard_H3K4me3_anno$symbol)
  #library(BRGenomics)
  #H3K4me1_peaks <- mergeGRangesData(m1_control_H3K4me3_gr, m1_ard_H3K4me3_gr, m2_ard_H3K4me3_gr, m2_control_H3K4me3_gr , field = "peak", ncores = 1)
  m2_ard_H3K4me3_anno[m2_ard_H3K4me3_anno$insideFeature %in% c("upstream",'overlapStart') & m2_ard_H3K4me3_anno$shortestDistance < 500, ]$symbol
  #load H3K4me3 data
  #H3K4me1_broad_peak_genes <- unique(c(m2_ard_H3K4me3_anno$symbol, m1_ard_H3K4me3_anno$symbol,  m1_control_H3K4me3_anno$symbol, m2_control_H3K4me3_anno$symbol))
  H3K4me3_broad_peak_genes_promoter <- unique(c(m2_ard_H3K4me3_anno[m2_ard_H3K4me3_anno$signalValue > 20, ]$symbol,
                                       m1_ard_H3K4me3_anno[m1_ard_H3K4me3_anno$signalValue > 20, ]$symbol,
                                       m1_control_H3K4me3_anno[m1_control_H3K4me3_anno$signalValue > 20, ]$symbol,
                                       m2_control_H3K4me3_anno[m2_control_H3K4me3_anno$signalValue > 20, ]$symbol))

  write.csv(H3K4me3_broad_peak_genes_promoter,'/Users/vyom/H3K4ME3_Peaks_promoters.csv' )
}
m3_ard_H3K4me1_anno$peak

write.csv(Yap_refine,'Users/vyom/refined_YAP_targets.csv' )
read.csv('Users/vyom/refined_YAP_targets.csv')
intersect_genes_refined <- intersect(CREB_Targets_Vyom$`Creb Targets`, Yap_refine)

intersect_genes <- intersect(CREB_Targets_Vyom$`Creb Targets`, Pan_tead_broad_peak_genes)
write.csv(intersect_genes, 'Users/vyom/Shared_targets_YAP_CREB.csv')
read.csv('Users/vyom/Shared_targets_YAP_CREB.csv')
AA_induced_organoid_scrna <- c('Fabp1',	'S100a6',	'Zg16',	'Dmbt1',	'Lgals3',	'Agr2',	'Reg1',	'Akr1b8',	'Anxa2',	'Cbr3',	'Ldha',	'Hmgcs2',	'2210407C18Rik',	'Gsta1',	'Krt19',	'Pkm',	'Acadl',	'AA467197',	'Apoa4',	'Pgam1',	'Gm20594',	'Sis',	'Apoa1',	'Reg4',	'Krt8',	'Slc11a2',	'Gna11',	'Anxa13',	'Gsta4',	'Guca2b',	'Prdx6',	'Serpinb6a',	'Nme2',	'Tubb4b',	'Krt18',	'Lypd8',	'Rbp2',	'Mttp',	'Nfkbia',	'Cldn3',	'Rpl15',	'Dstn',	'Tuba1b',	'Ckmt1',	'Gpx2',	'Krt7',	'Prss32',	'Ube2m',	'Alpi',	'Gsta3',	'Wfdc2',	'Slc16a3',	'Anxa3',	'Capg',	'Aldob',	'Cbr1',	'Gstm1',	'Arpc4',	'Ddah1',	'Egln3',	'Tomm22',	'Tfrc',	'Mdh2',	'Actg1',	'Eif5a',	'Akr1c19',	'Pgk1',	'Ptgr1',	'Cmas',	'Aldh1a1',	'Ezr',	'Cct7',	'Aprt',	'Tm4sf4',	'2610528A11Rik',	'Eif6',	'Ankrd37',	'Tuba4a',	'Ddx3x',	'Acaa2',	'Galk1',	'Tm4sf5',	'Lars2',	'Nfe2l2',	'Oit1',	'Aldoa',	'Ddx39',	'Plpp2',	'Arhgdig',	'Vdac2',	'Ctsl',	'Gclm',	'Tpi1',	'Prdx1',	'Gsr',	'Gpd1',	'Tspan8',	'Hspd1',	'Lgals4',	'Phb',	'Ctsz',	'Actb',	'Me1',	'St3gal4',	'Tmem54',	'Gm10116',	'Acp5',	'Cldn15',	'Ero1l',	'Akr1c13',	'Hras',	'Sri',	'Tmem237',	'Eif4a1',	'Id1',	'Areg',	'Arf1',	'Pdcd6',	'Tomm20',	'Pgp',	'Rpl7l1',	'Pfkp',	'Fabp2',	'Uqcrfs1',	'Sfn',	'Muc2',	'Plxnb2',	'Gmds',	'Snrpb',	'1110008F13Rik',	'Pfkl',	'Lman2',	'Taldo1',	'Suclg1',	'Cyc1',	'Ppp1ca',	'Cda',	'Tspan1',	'Slc25a3',	'Cybrd1',	'Ech1',	'Ppa1',	'Uba52',	'Eif4g1',	'Srsf7',	'Esd',	'Rpl6',	'Cyba',	'Lrrc59',	'Ier3',	'Ckb',	'Uqcrc1',	'Ran',	'Polr2e',	'Cd9',	'Gipc1',	'Gstp1',	'Tomm40',	'Cldn2',	'Aldh1b1',	'Ctsd',	'Stoml2',	'Rac1',	'Syngr2',	'Eci2',	'Slc25a4',	'Psmd13',	'Lamp1',	'Plcb3',	'Gsto1',	'Ak2',	'Ccnd1',	'Tmem45b',	'Lmna',	'Rab25',	'Acadvl',	'Plac8',	'Maoa',	'Oxct1',	'Sumo1',	'Adh6a',	'Map2k2',	'Ugt2b34',	'Gpi1',	'Mapk13',	'Tmprss2',	'Kif5b',	'Txnrd1',	'Ctsh',	'Taf10',	'Grpel1',	'Cyp4f14',	'Fbp2',	'Actr3',	'Clic1',	'Twf1',	'Tmem14c',	'Hsp90aa1',	'Atp1a1',	'Junb',	'Wdr18',	'Elf3',	'Akr1b3',	'Cyp2d26',	'Ahsa1',	'Psma1',	'Ctnnb1',	'Hspa4',	'Txn2',	'Psmd3',	'Mrps10',	'Bag1',	'Mcm3',	'Cct2',	'Sdc1',	'Prmt1',	'Ddb1',	'Tsta3',	'Ccnd2',	'Eps8l3',	'Bsg' )      
AA_induced_invivo <- c('Hspa1b',	'Gm2000',	'Defa24',	'Banp',	'Lockd',	'Dpm3',	'Rps27rt',	'Snrpg',	'mt-Nd3',	'Cenpw',	'Slirp',	'Pet100',	'Tomm7',	'Tmem256',	'Cops9',	'Atp5k',	'Snhg6',	'Smim4',	'Mrpl52',	'Ndufb1-ps',	'Epb41l4aos',	'Rpl38',	'Pin4',	'Snhg8',	'Snrpf',	'Ndufa2',	'Psmg4',	'Atp5mpl',	'Romo1',	'Rps29',	'Rps21',	'Rps28',	'Anapc13',	'Atp5md',	'Lsm7',	'Tmem258',	'Uqcc2',	'Plp2',	'S100a6',	'Phgdh',	'Ndufa3',	'Eif3j1',	'Nop10',	'Sec61g',	'Aspm',	'Atox1',	'Polr2l',	'Ndufb2',	'Ccdc167',	'Fkbp5',	'Bola2',	'Ndufa5',	'Ndufc1',	'Mrps21',	'Snrpe',	'Uqcr11',	'Hist1h2ae',	'Sp3os',	'Ubl5',	'Rpl39',	'Tuba1c',	'Rpl37',	'Prelp',	'Mrpl33',	'Snrpd2',	'Naa38',	'Polr2k',	'Uqcr10',	'Cox17',	'Smim22',	'Rplp2',	'Prc1',	'Rpl37a',	'Rpl31',	'Sgo1',	'Rpl21',	'Fmc1',	'Ifitm3',	'Ppih',	'Nusap1',	'Trappc10',	'Zfas1',	'Rpl36',	'Ncapg',	'Polr2i',	'Tmsb15b2',	'Cox20',	'Ndufv3',	'Gm11808',	'Atp5e',	'Arpp19',	'Adra2a',	'Fxyd3',	'Mt1',	'Esco2',	'2410006H16Rik',	'Rpl35a',	'Ndufa1',	'Elob',	'Hist1h1b',	'Ndufs6',	'2310009B15Rik',	'Cit',	'Rpl41',	'Hspe1-rs1',	'S100a11',	'Ndufa7',	'Ndufb4',	'S100g',	'Srek1ip1',	'Sem1',	'Hmgcs2',	'Mif',	'Aurka',	'Diaph3',	'Cd44',	'Knl1',	'Snhg3',	'Top2a',	'Hells',	'Rrp8',	'Plk1',	'Ddx18',	'Cep164',	'Kif20b',	'Cdca2',	'Pmepa1',	'Jaml',	'Incenp',	'Smim26',	'Trip13',	'Mis18bp1',	'Rbmx2',	'Lsm6',	'Spc25',	'Rrp1b',	'Rps17',	'Lsm5',	'S100a13',	'Rpgrip1',	'Cdca8',	'Gas5',	'Ddx24',	'Hist1h1e',	'Cenpf',	'Hmmr',	'Atp5j2',	'Acot1',	'Kif22',	'Neil3',	'Rpl34',	'Mis12',	'Cox7a1',	'Tuba1b',	'Cbx6',	'Tubb2b',	'Ccdc34',	'Rpl36a',	'Rbmxl1',	'Ndc80',	'2310009A05Rik',	'Cep295',	'Bub1',	'Spink4',	'Hirip3',	'H2afj',	'Dlgap5',	'Atp5l',	'Smim11',	'Cox6c',	'Dnajb1',	'Kif11',	'Ubap2',	'Tpx2',	'Cebpz',	'Ctc1',	'Crip1',	'Sycn',	'Telo2',	'Rbis',	'Rassf4',	'Shcbp1',	'Cox7c',	'Gm43813',	'Rgcc',	'Notch1',	'Serf1',	'Ndufa13',	'Aurkb',	'Pbk',	'Sec61b',	'Lrrc31',	'Chordc1',	'Rps27l',	'1810037I17Rik',	'Rpp21',	'2200002D01Rik',	'E030030I06Rik',	'Fbxo5',	'Rad51ap1',	'Uqcrq',	'Polr2f',	'Rpl35',	'Ttr',	'Nup210',	'Krtcap2',	'Mcph1',	'Rps27',	'Ndufaf8',	'Lig1',	'Cd3eap',	'Ppwd1',	'1700097N02Rik',	'2810408I11Rik',	'Ckap2',	'2810004N23Rik',	'Rps2',	'Prpf4',	'Snhg18',	'Kif4',	'Sgo2a',	'Higd1a',	'Ndufs5',	'Xpo5',	'Soat1',	'Rps15',	'Anln',	'Rps26',	'Grwd1',	'Lrrc14',	'Cox7a2',	'Ncapd2',	'Tcof1',	'Rad54l',	'Ier3ip1',	'Cenpe',	'Apobec3',	'Lyrm2',	'Bub1b',	'Ccdc85c',	'Cep192',	'Prpf6',	'Cox7b',	'Rps25',	'Pinx1',	'Tex10',	'Gpatch4',	'Gnl3')          
AA_induced_stem2 <- c('S100a6','Ly6a','Malat1','Tmsb4x','Pigr','Ftl1','Gnas','Ubc')          
AA_induced_stem2 <- c('Plekhm3',	'H2-Eb1',	'Cd74',	'H2-Ab1',	'Nucb2',	'Map7d1',	'Sik1',	'Upb1',	'Arid3a',	'Cyb561a3',	'Ankrd44',	'Ctsl',	'Tbc1d8',	'Slc29a3',	'Orai3',	'Dusp1',	'Ifnar2',	'Tcf4',	'Irf7',	'Unc93b1',	'Dusp6',	'Ptpn6',	'Tubgcp5',	'S100a6',	'Ptpre',	'Rell1',	'Tom1',	'Crlf2',	'Snx29',	'Gm37529',	'Gsn',	'Snx18',	'Gmfg',	'Mob3a',	'Mfsd12',	'Unc119',	'Uvrag',	'Cybc1',	'Cmtm7',	'Cnp',	'Npc1',	'Map4k2',	'Stat2',	'H2-Q4',	'Gm43329',	'Tmem123',	'Cdkn2d',	'Nek7',	'Slc15a4',	'Gns',	'Anxa5',	'Bcl11a',	'Ppp1r15a',	'Sema4b',	'Itpr1',	'Aldh3b1',	'Igtp',	'Gng10',	'Dock8',	'Fam174a',	'Pmepa1',	'Rogdi',	'Bbc3',	'Rab33b',	'Ftx',	'Pde7a',	'Lyn',	'Dclre1c',	'Atp13a2',	'Ccl5',	'Cbx4',	'Tgfbr1',	'Gpcpd1',	'Pqlc3',	'Rpgrip1',	'Arid3b',	'Ablim1',	'Gm31718',	'Dipk1a',	'Lgmn',	'Grina',	'Rhobtb2',	'Snx30',	'Nudt18',	'Trim12a',	'Slc49a4',	'Myo9b',	'Sgk3',	'Tor3a',	'Bicd2',	'Rwdd2a',	'Dtx2',	'Ctss',	'Slc25a4',	'5031439G07Rik',	'Grn',	'Tex2',	'Vamp1',	'A430035B10Rik',	'Inpp5k',	'Ccnd3',	'Irf8',	'Stat5a',	'Rab2b',	'Tap1',	'Pacc1',	'Ifnar1',	'Ier5',	'Pml',	'Acox3',	'Arl10',	'Il4ra',	'Tap2',	'Apobec3',	'Litaf',	'Malt1',	'Il17ra',	'Bmyc',	'Arhgap17',	'Zfp319',	'Bnip2',	'Pkd1',	'Tpst2',	'H2-T22',	'Calcoco1',	'Pip5k1c',	'Gmip',	'Crtc3',	'Ptpn1',	'L3mbtl3',	'Tcirg1',	'Prcp',	'Stat1',	'Cyth1',	'Nabp1',	'Nek6',	'Bin1',	'Arhgap27',	'Wdr81',	'Ptprj')

intersect(AA_induced_stem2, CREB_Targets_Vyom$`Creb Targets`)

#make upset plot for all groups
ARD_upset <- list()
#ARD_upset$AA_organoid <- AA_induced_organoid_scrna
ARD_upset$AA_induced <- AA_induced_invivo
#ARD_upset$AA_induced_stem2 <- AA_induced_stem2
Creb_enhancer <- intersect(Creb_broad_peak_genes_enhancer, H3K4me1_broad_peak_genes )
Pan_Tead_enhancer <- intersect(H3K4me1_broad_peak_genes, Pan_tead_broad_peak_genes_enhancer )

Creb_promoter<- intersect(Creb_broad_peak_genes_promoter, H3K4me3_broad_peak_genes_promoter )
Pan_Tead_promter <- intersect(Pan_tead_broad_peak_genes_promoter, H3K4me3_broad_peak_genes_promoter )
ARD_upset$Creb_promoter <- Creb_promoter
ARD_upset$Creb_enhancer <- Creb_enhancer
ARD_upset$Pan_Tead_promoter <- Pan_Tead_promter
ARD_upset$Pan_Tead_enhancer <- Pan_Tead_enhancer

write.csv(Creb_promoter, 'Creb_promoter.csv')
write.csv(Creb_enhancer, 'Creb_enhancer.csv')
write.csv(Pan_Tead_promter, 'Pan_Tead_promter.csv')
write.csv(Pan_Tead_enhancer, 'Pan_Tead_enhancer.csv')

library(UpSetR)
#upset(fromList(ARD_upset),nsets = 7,  order.by = "freq", keep.order = T,  sets = rev(c("AA_organoid", "AA_induced", "AA_induced_stem2", "Creb_promoter", "Creb_enhancer", "Pan_Tead_promoter", "Pan_Tead_enhancer")))
upset(fromList(ARD_upset),nsets = 5, nintersects = 20,  order.by = "freq", keep.order = T,number.angles = 0, mb.ratio = c(0.7, 0.3),  sets = rev(c("AA_induced", "Creb_promoter", "Creb_enhancer", "Pan_Tead_promoter", "Pan_Tead_enhancer")))

Reduce(intersect, list(AA_induced_invivo, Creb_promoter, Creb_enhancer, Pan_Tead_promter, Pan_Tead_enhancer))

Reduce(intersect, list(Creb_promoter, Creb_enhancer, Pan_Tead_promter, Pan_Tead_enhancer))

#make a venn diagram comparing the things
write.csv(Creb_promoter, 'Creb_promoter.csv')
write.csv(Creb_enhancer, 'Creb_enhancer.csv')
write.csv(Pan_Tead_promter, 'Pan_Tead_promter.csv')
write.csv(Pan_Tead_enhancer, 'Pan_Tead_enhancer.csv')




CREB_Targets <- CREB_Targets[CREB_Targets$insideFeature %in% c("upstream",'overlapStart') & CREB_Targets$shortestDistance < 500, ]$symbol
YAP_targets <- YAP_targets$x
YAP_targets <-  YAP_targets %>% na.omit()
CREB_Targets <-  CREB_Targets %>% na.omit()

pastel_colors <- c( "#B0E2FF", "#CAB8FF", "#FDFD96" ,  "#A6FBB2", "#FF9999")
devtools::install_github("nicolash2/ggvenn")

library(ggvenn)
ggplot() + geom_venn(ARD_upset, fill_color = pastel_colors)+ theme_void()
ggsave(file = paste0('ARD_YAP_CREB_Shared_Targets.pdf'), width=5, height=4, units="in")
names(ARD_upset)
library("ggVennDiagram")
# Default plot
ggVennDiagram(ARD_upset) + scale_fill_gradientn(colors = pastel_colors)
intersect_genes <- Reduce(intersect, list(AA_induced_invivo, Creb_promoter, Creb_enhancer, Pan_Tead_promter, Pan_Tead_enhancer))
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023", "MSigDB_Hallmark_2020")
enriched <- enrichr(intersect_genes, dbs)
Enriched1 <- rbind(enriched[[1]],enriched[[2]],enriched[[3]],enriched[[4]] )
write.csv(Enriched1, paste0('Shared_Yap_creb_targets_ENRICHR.csv'))
Enriched1
Enriched_filter <- Enriched1[Enriched1$Adjusted.P.value < .05,]
plotEnrich(Enriched_filter, showTerms = 100, numChar = 100, y = "Count", orderBy = "Adjusted.P.value")

#get appropriate peaklists
Creb_promoter <- read.csv("~/analysis/Creb_promoter.csv")$x
Creb_enhancer <- read.csv("~/analysis/Creb_enhancer.csv")$x
Pan_Tead_promter <- read.csv("~/analysis/Pan_Tead_promter.csv")$x
Pan_Tead_enhancer <- read.csv("~/analysis/Pan_Tead_enhancer.csv")$x
AA_induced_invivo <- c('Hspa1b',	'Gm2000',	'Defa24',	'Banp',	'Lockd',	'Dpm3',	'Rps27rt',	'Snrpg',	'mt-Nd3',	'Cenpw',	'Slirp',	'Pet100',	'Tomm7',	'Tmem256',	'Cops9',	'Atp5k',	'Snhg6',	'Smim4',	'Mrpl52',	'Ndufb1-ps',	'Epb41l4aos',	'Rpl38',	'Pin4',	'Snhg8',	'Snrpf',	'Ndufa2',	'Psmg4',	'Atp5mpl',	'Romo1',	'Rps29',	'Rps21',	'Rps28',	'Anapc13',	'Atp5md',	'Lsm7',	'Tmem258',	'Uqcc2',	'Plp2',	'S100a6',	'Phgdh',	'Ndufa3',	'Eif3j1',	'Nop10',	'Sec61g',	'Aspm',	'Atox1',	'Polr2l',	'Ndufb2',	'Ccdc167',	'Fkbp5',	'Bola2',	'Ndufa5',	'Ndufc1',	'Mrps21',	'Snrpe',	'Uqcr11',	'Hist1h2ae',	'Sp3os',	'Ubl5',	'Rpl39',	'Tuba1c',	'Rpl37',	'Prelp',	'Mrpl33',	'Snrpd2',	'Naa38',	'Polr2k',	'Uqcr10',	'Cox17',	'Smim22',	'Rplp2',	'Prc1',	'Rpl37a',	'Rpl31',	'Sgo1',	'Rpl21',	'Fmc1',	'Ifitm3',	'Ppih',	'Nusap1',	'Trappc10',	'Zfas1',	'Rpl36',	'Ncapg',	'Polr2i',	'Tmsb15b2',	'Cox20',	'Ndufv3',	'Gm11808',	'Atp5e',	'Arpp19',	'Adra2a',	'Fxyd3',	'Mt1',	'Esco2',	'2410006H16Rik',	'Rpl35a',	'Ndufa1',	'Elob',	'Hist1h1b',	'Ndufs6',	'2310009B15Rik',	'Cit',	'Rpl41',	'Hspe1-rs1',	'S100a11',	'Ndufa7',	'Ndufb4',	'S100g',	'Srek1ip1',	'Sem1',	'Hmgcs2',	'Mif',	'Aurka',	'Diaph3',	'Cd44',	'Knl1',	'Snhg3',	'Top2a',	'Hells',	'Rrp8',	'Plk1',	'Ddx18',	'Cep164',	'Kif20b',	'Cdca2',	'Pmepa1',	'Jaml',	'Incenp',	'Smim26',	'Trip13',	'Mis18bp1',	'Rbmx2',	'Lsm6',	'Spc25',	'Rrp1b',	'Rps17',	'Lsm5',	'S100a13',	'Rpgrip1',	'Cdca8',	'Gas5',	'Ddx24',	'Hist1h1e',	'Cenpf',	'Hmmr',	'Atp5j2',	'Acot1',	'Kif22',	'Neil3',	'Rpl34',	'Mis12',	'Cox7a1',	'Tuba1b',	'Cbx6',	'Tubb2b',	'Ccdc34',	'Rpl36a',	'Rbmxl1',	'Ndc80',	'2310009A05Rik',	'Cep295',	'Bub1',	'Spink4',	'Hirip3',	'H2afj',	'Dlgap5',	'Atp5l',	'Smim11',	'Cox6c',	'Dnajb1',	'Kif11',	'Ubap2',	'Tpx2',	'Cebpz',	'Ctc1',	'Crip1',	'Sycn',	'Telo2',	'Rbis',	'Rassf4',	'Shcbp1',	'Cox7c',	'Gm43813',	'Rgcc',	'Notch1',	'Serf1',	'Ndufa13',	'Aurkb',	'Pbk',	'Sec61b',	'Lrrc31',	'Chordc1',	'Rps27l',	'1810037I17Rik',	'Rpp21',	'2200002D01Rik',	'E030030I06Rik',	'Fbxo5',	'Rad51ap1',	'Uqcrq',	'Polr2f',	'Rpl35',	'Ttr',	'Nup210',	'Krtcap2',	'Mcph1',	'Rps27',	'Ndufaf8',	'Lig1',	'Cd3eap',	'Ppwd1',	'1700097N02Rik',	'2810408I11Rik',	'Ckap2',	'2810004N23Rik',	'Rps2',	'Prpf4',	'Snhg18',	'Kif4',	'Sgo2a',	'Higd1a',	'Ndufs5',	'Xpo5',	'Soat1',	'Rps15',	'Anln',	'Rps26',	'Grwd1',	'Lrrc14',	'Cox7a2',	'Ncapd2',	'Tcof1',	'Rad54l',	'Ier3ip1',	'Cenpe',	'Apobec3',	'Lyrm2',	'Bub1b',	'Ccdc85c',	'Cep192',	'Prpf6',	'Cox7b',	'Rps25',	'Pinx1',	'Tex10',	'Gpatch4',	'Gnl3')          

AA_Target_bound <- Reduce(intersect, list(AA_induced_invivo, Creb_promoter, Creb_enhancer, Pan_Tead_promter, Pan_Tead_enhancer))

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
data(TSS.mouse.GRCm38)
gr <- transcripts(txdb)
gr <- annotatePeakInBatch(gr, AnnotationData=TSS.mouse.GRCm38)

gr1 <- addGeneIDs(annotatedPeak=gr,
                  orgAnn="org.Mm.eg.db",
                  IDs2Add="symbol")
gr2 <- gr1[gr1$symbol %in% AA_Target_bound]

df <- data.frame(seqnames=seqnames(gr2),
                 starts=start(gr2)-1,
                 ends=end(gr2),
                 names=c(rep(".", length(gr2))),
                 scores=c(rep(".", length(gr2))),
                 strands=strand(gr2))

write.table(df, file="~/analysis/AA_Target_Bount_induced.bed", quote=F, sep="\t", row.names=F, col.names=F)






