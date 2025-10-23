library(DiffBind)
library(profileplyr)
theme_linedraw2 = theme_linedraw() + theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), strip.text.x=element_text(colour="black"), strip.text.y=element_text(colour="black"))
theme_vyom = theme_linedraw2 + theme(legend.position="right", legend.title=element_text(size=15), legend.text=element_text(size=14), axis.text.x = element_text(size=12, angle=-90, hjust=0, vjust=0.5), axis.text.y=element_text(size=12), axis.title=element_text(size=15), axis.title.y=element_text(vjust=1), plot.title = element_text(size=18, vjust=1.5), strip.background = element_rect(fill="#EEEEEE"), strip.text = element_text(size = 11), panel.grid.major = element_line(colour = "grey98"), panel.grid.minor = element_blank())

#creb
proplyrObject <- import_deepToolsMat("CREB_proper_ARD_induced.bed_computematrix.mat.gz") 
#proplyrObject <- import_deepToolsMat("ARDKO_induced_creb_computematrix.mat.gz") 

AR2 <- proplyrObject@assays@data$`M1-ARD-CREB_R2`
AR3 <- proplyrObject@assays@data$`M1-ARD-CREB_R3`
AR1 <- proplyrObject@assays@data$`M2-ARD-CREB_R1`

ARD <- (AR1 + AR2 + AR3) / 3
#ARD <- (log10(AR1) + log10(AR2) + log10(AR3)) / 3
ARD <- log2(ARD)
ARD[!is.finite(ARD)] <- 0

C1 <- proplyrObject@assays@data$`M1-Control-CREB_R1`
C2 <- proplyrObject@assays@data$`M2-Control-CREB_R2`
C3 <- proplyrObject@assays@data$`M2-Control-CREB_R3`

WT <- (C1 + C2 + C3) / 3
#WT <- (log10(C1) + log10(C2) + log10(C3)) / 3
WT <- log2(WT)
WT[!is.finite(WT)] <- 0

vector1 <- c(as.vector(AR1),as.vector(AR2), as.vector(AR3))
vector2 <- c(as.vector(C1),as.vector(C2), as.vector(C3))

# Perform the Wilcoxon rank-sum test
ks_Res <- ks.test(log2(vector1), log2(vector1), alternative = c("two.sided"))
ks_Res
result <- wilcox.test(vector1, vector2, paired = FALSE, B=2000)
result

proplyrObject@assays@data$ARD <- ARD
proplyrObject@assays@data$WT <- WT

proplyrObject@assays@data$`M1-ARD-CREB_R2` = NULL
proplyrObject@assays@data$`M1-ARD-CREB_R3`= NULL
proplyrObject@assays@data$`M2-ARD-CREB_R1`= NULL
proplyrObject@assays@data$`M1-Control-CREB_R1`= NULL
proplyrObject@assays@data$`M2-Control-CREB_R2`= NULL
proplyrObject@assays@data$`M2-Control-CREB_R3`= NULL



sampleData <- proplyrObject@sampleData[1:2,]

rownames(sampleData) <- c('ARD','WT')
sampleData$sample_labels <-c('ARD','WT')
proplyrObject@sampleData <- sampleData

output_path <- file.path("ARD_merged_cutrun_creb.gz")
export_deepToolsMat(proplyrObject, con = output_path, overwrite=TRUE)

dba.plotProfile(proplyrObject, raster_quality = 1.5, normalize=TRUE)

#yap

proplyrObject <- import_deepToolsMat("YAP_ARD_induced_computematrix.mat.gz") 
#proplyrObject <- import_deepToolsMat("ARDKO_induced_yap_computematrix.mat.gz") 


proplyrObject
names(proplyrObject@assays@data)
AR2 <- proplyrObject@assays@data$`M3-ARD-panTEAD_R1`
AR3 <- proplyrObject@assays@data$`M3-ARD-panTEAD_R2`
AR1 <- proplyrObject@assays@data$`M4-ARD-panTEAD_R1`
AR4 <- proplyrObject@assays@data$`M4-ARD-panTEAD_R2`

ARD <- (AR1 + AR2 + AR3 + AR4) / 4


C1 <- proplyrObject@assays@data$`M3-Control-panTEAD_R1`
C2 <- proplyrObject@assays@data$`M3-Control-panTEAD_R2`
C3 <- proplyrObject@assays@data$`M4-Control-panTEAD_R1`
C4 <- proplyrObject@assays@data$`M4-Control-panTEAD_R2`

WT <- (C1 + C2 + C3 + C4) / 4

proplyrObject@assays@data$ARD <- ARD
proplyrObject@assays@data$WT <- WT

proplyrObject@assays@data$`M3-Control-panTEAD_R1` = NULL
proplyrObject@assays@data$`M3-Control-panTEAD_R2`= NULL
proplyrObject@assays@data$`M4-Control-panTEAD_R1`= NULL
proplyrObject@assays@data$`M4-Control-panTEAD_R2`= NULL
proplyrObject@assays@data$`M3-ARD-panTEAD_R1`= NULL
proplyrObject@assays@data$`M3-ARD-panTEAD_R2`= NULL
proplyrObject@assays@data$`M4-ARD-panTEAD_R1`= NULL
proplyrObject@assays@data$`M4-ARD-panTEAD_R2`= NULL

sampleData <- proplyrObject@sampleData[1:2,]

rownames(sampleData) <- c('ARD','WT')
sampleData$sample_labels <-c('ARD','WT')
proplyrObject@sampleData <- sampleData

dba.plotProfile(proplyrObject, raster_quality = 1.5)

output_path <- file.path("ARD_merged_cutrun_YAP.gz")
export_deepToolsMat(proplyrObject, con = output_path, overwrite=TRUE)
#do signficiance test:

# Sample matrices
matrix1 <- ARD
matrix2 <- WT

# Check dimensions to ensure they are the same
if(!all(dim(matrix1) == dim(matrix2))) {
  stop("Matrices must have the same dimensions.")
}

# Flatten the matrices into vectors
vector1 <- as.vector(matrix1)
vector2 <- as.vector(matrix2)

# Perform the Wilcoxon rank-sum test
ks_Res <- ks.test(log2(vector1), log2(vector1), alternative = c("two.sided"))
ks_Res
result <- wilcox.test(vector1, vector2, paired = FALSE, B=2000)

# Display the result
print(result)

CREBres <- result
YAPres <- result

#make boxplot for max values: CREB
ARD <- (AR1 + AR2 + AR3) / 3
WT <- (C1 + C2 + C3) / 3

c(max(AR1), max(AR2), max(AR3))
c(max(C1), max(C2), max(C3))

#kmeans clustering:


set.seed(0)
kmeans <- clusterRanges(proplyrObject, 
                        fun = rowMeans, 
                        kmeans_k = 2, 
                        silent = FALSE)
hclust <- clusterRanges(proplyrObject, 
                        fun = rowMeans, 
                        cutree_rows = 2, 
                        silent = FALSE)

generateEnrichedHeatmap(hclust)
library(magrittr)
library(SummarizedExperiment)
proplyrObject %>%
  clusterRanges(fun = rowMeans, cutree_rows = 2, silent = TRUE) %>% 
  profileplyr::summarize(fun = rowMeans, output = "long") %>% 
  ggplot(aes(x = Sample, y = Signal)) + 
  geom_boxplot() + 
  facet_grid(~cluster) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
median(AR1)
proplyrObject_long <- profileplyr::summarize(proplyrObject, 
                                fun = rowMeans, 
                                output = "long") 
# 1 for creb 5 for yap
n <- 50
rows_use <- c()
groups <- unique(proplyrObject_long$Sample)

for(i in groups){
  prof_subset <- proplyrObject_long[proplyrObject_long$Sample == i,]
  rows_use <- c(rows_use,rownames(prof_subset[prof_subset$Signal > quantile(prof_subset$Signal,prob=1-n/100),]))                    
  
}
box_signal <- proplyrObject_long
box_signal <- proplyrObject_long[rows_use, ]

#creb Specific analysis 
box_signal_rename <- box_signal %>% mutate(Sample = recode(Sample, `M1-ARD-CREB_R2` = 'ARD 1', `M1-ARD-CREB_R3`= 'ARD 2', `M1-Control-CREB_R1`= 'WT 1', `M2-ARD-CREB_R1`= 'ARD 3', `M2-Control-CREB_R2`= 'WT 2', `M2-Control-CREB_R3`= 'WT 3'))                
box_signal_rename$Sample <- factor(box_signal_rename$Sample, levels = c('WT 1', 'WT 2', 'WT 3', 'ARD 1', 'ARD 2', 'ARD 3'))

ggplot(box_signal_rename, aes(x = Sample, y = Signal)) + geom_boxplot(notch = FALSE, outlier.size = .25) + theme_vyom
ggsave(file = paste0('ARD_induced_Creb_cutrun_boxplot_by_sample.pdf'),  width=2.25, height=3, units="in")



min(box_signal_merge$Signal)
box_signal_merge <- box_signal %>% mutate(Sample = recode(Sample, `M1-ARD-CREB_R2` = 'ARD', `M1-ARD-CREB_R3`= 'ARD', `M1-Control-CREB_R1`= 'WT', `M2-ARD-CREB_R1`= 'ARD', `M2-Control-CREB_R2`= 'WT', `M2-Control-CREB_R3`= 'WT'))                
box_signal_merge$Sample <- factor(box_signal_merge$Sample, levels = c('WT', 'ARD'))

ggplot(box_signal_merge, aes(x = Sample, y = Signal)) + geom_boxplot(notch = TRUE, outlier.shape = NA) + theme_vyom +
  ggpubr::geom_pwc(data= box_signal_merge, aes(x = Sample, y = Signal),na.rm = TRUE,  method = "wilcox.test", hjust = .5,  tip.length = 0.01, vjust = .5, y.position = 50, label.size = 7, label = "p.signif", bracket.nudge.y = .03, step.increase = 0.05)
ggsave(file = paste0('ARD_induced_Creb_cutrun_boxplot_combined_break.pdf'),  width=1.5, height=3, units="in")

box_signal_merge$Signal1 <- log2(box_signal_merge$Signal)
ggplot(box_signal_merge, aes(x = Sample, y = Signal1)) + geom_violin(notch = TRUE, outlier.size = .25) + theme_vyom

#YAP Specific analysis 
box_signal_rename <- box_signal %>% mutate(Sample = recode(Sample, `M3-ARD-panTEAD_R1`= 'ARD 1', `M3-ARD-panTEAD_R2`= 'ARD 2', `M3-Control-panTEAD_R1`= 'WT 1',  `M3-Control-panTEAD_R2` = 'WT 2',  `M4-ARD-panTEAD_R1`= 'ARD 3', `M4-ARD-panTEAD_R2`= 'ARD 4',  `M4-Control-panTEAD_R1`= 'WT 3', `M4-Control-panTEAD_R2` = 'WT 4'))                
box_signal_rename$Sample <- factor(box_signal_rename$Sample, levels = c('WT 1', 'WT 2', 'WT 3', 'WT 4', 'ARD 1', 'ARD 2', 'ARD 3', 'ARD 4'))

ggplot(box_signal_rename, aes(x = Sample, y = Signal)) + geom_boxplot(notch = FALSE, outlier.size = .25) + theme_vyom
ggsave(file = paste0('ARD_induced_YAP_bound_new_boxplot_by_sample.pdf'),  width=2.45, height=3, units="in")
ggplot(box_signal_rename, aes(x = Sample, y = Signal)) + geom_violin(notch = FALSE, outlier.size = .25) + theme_vyom
ggsave(file = paste0('ARD_induced_YAP_bound_new_violin_by_sample.pdf'),  width=2.45, height=3, units="in")


box_signal_merge <- box_signal %>% mutate(Sample = recode(Sample, `M3-ARD-panTEAD_R1`= 'ARD', `M3-ARD-panTEAD_R2`= 'ARD', `M3-Control-panTEAD_R1`= 'WT',  `M3-Control-panTEAD_R2` = 'WT',  `M4-ARD-panTEAD_R1`= 'ARD', `M4-ARD-panTEAD_R2`= 'ARD',  `M4-Control-panTEAD_R1`= 'WT', `M4-Control-panTEAD_R2` = 'WT'))                
box_signal_merge$Sample <- factor(box_signal_merge$Sample, levels = c('WT', 'ARD'))

ggplot(box_signal_merge, aes(x = Sample, y = Signal)) + geom_boxplot(notch = TRUE, outlier.shape = NA) + theme_vyom + ylim(0, 80)
  ggpubr::geom_pwc(data= box_signal_merge, aes(x = Sample, y = Signal),na.rm = TRUE,  method = "t_test", hjust = .5,  tip.length = 0.01, vjust = .5, y.position = 190, label.size = 7, label = "p.signif", bracket.nudge.y = .03, step.increase = 0.05)
ggsave(file = paste0('ARD_induced_YAP_new_bound_boxplot_combined_break.pdf'),  width=1.5, height=3, units="in")




#annotation:
anno_great <- annotateRanges_great(proplyrObject, species = "mm10")
peak <- rowRanges(anno_great)
peak <- annotatePeakInBatch(peak, AnnotationData=TSS.mouse.GRCm38)
peak <- addGeneIDs(annotatedPeak=peak,
                       orgAnn="org.Mm.eg.db",
                       IDs2Add="symbol")
rowRanges(anno_great) <- peak 
anno_great@assays@data$ARD

#creb
{
CREB_ARD_Mat <- ARD-WT
CREB_ARD_Mat <- as.data.frame(CREB_ARD_Mat)
CREB_ARD_Mat$enrich <- rowMeans(CREB_ARD_Mat)
CREB_ARD_Mat$gene <- peak$symbol

CREB_ARD_summ <- CREB_ARD_Mat %>% 
  dplyr::select(enrich, gene) %>% 
  group_by(gene) %>% 
  dplyr::summarize(enrich = mean(enrich))

CREB_ARD_summ
CREB_ARD_summ <- CREB_ARD_summ[rev(order(CREB_ARD_summ$enrich)),]
}

#yap
{
  YAP_ARD_Mat <- ARD-WT
  YAP_ARD_Mat <- as.data.frame(YAP_ARD_Mat)
  YAP_ARD_Mat$enrich <- rowMeans(YAP_ARD_Mat)
  YAP_ARD_Mat$gene <- peak$symbol
  
  YAP_ARD_summ <- YAP_ARD_Mat %>% 
    dplyr::select(enrich, gene) %>% 
    group_by(gene) %>% 
    dplyr::summarize(enrich = mean(enrich))
  
  YAP_ARD_summ
  YAP_ARD_summ <- YAP_ARD_summ[rev(order(YAP_ARD_summ$enrich)),]
}
CREB_ARD_summ
YAP_ARD_summ
intersect(CREB_ARD_summ[CREB_ARD_summ$enrich >0,]$gene, YAP_ARD_summ[YAP_ARD_summ$enrich >0,]$gene)

intersect(CREB_ARD_summ[CREB_ARD_summ$enrich <0,]$gene, YAP_ARD_summ[YAP_ARD_summ$enrich <0,]$gene)

#plot concordance
# concordance plot new 
cond1 <- CREB_ARD_summ
cond2 <- YAP_ARD_summ

cond1_unique_genes <- c(cond1$gene[!(cond1$gene %in% cond2$gene)])
cond2_unique_genes <- c(cond2$gene[!(cond2$gene %in% cond1$gene)])

add_cond1 <- data.frame(gene = cond2_unique_genes, enrich = rep(0, length(cond2_unique_genes)))
add_cond2 <- data.frame(gene = cond1_unique_genes, enrich = rep(0, length(cond1_unique_genes)))

cond1 <- bind_rows(cond1, add_cond1)
cond2 <- bind_rows(cond2, add_cond2)
DE_merge <- merge(cond1,cond2, by=c("gene"))

p_cal <- cor.test(DE_merge$enrich.x, DE_merge$enrich.y)
p_cal$p.value
r_squared <- cor(DE_merge$enrich.x, DE_merge$enrich.y)

ggplot(DE_merge, aes(x=enrich.x, y=enrich.y)) + geom_vline(xintercept = 0, linetype="dashed") + geom_hline(yintercept = 0, linetype="dashed") +
  ggrepel::geom_text_repel(data=DE_merge,
                           aes(x=enrich.x,
                               y=enrich.y,
                               label=gene),
                           size = 2, color='Black',
                           point.padding=0.1,
                           segment.size = 0.2,
                           min.segment.length=0.2) +
  geom_point(size = 1, color=c('Black')) + theme_vyom + labs(y = 'YAP enrichment difference (ARD - WT)', x = 'CREB enrichment difference (ARD - WT)')
ggsave(file = paste0('CREB_YAP_Co_occupancy_ARD_induced.pdf'),  width=5, height=5, units="in")


ARD_upset <- list()
ARD_upset$CREB_Up <- CREB_ARD_summ[CREB_ARD_summ$enrich >0,]$gene
ARD_upset$CREB_Down <- CREB_ARD_summ[CREB_ARD_summ$enrich <0,]$gene
ARD_upset$YAP_Up <- YAP_ARD_summ[YAP_ARD_summ$enrich >0,]$gene
ARD_upset$YAP_Down <- YAP_ARD_summ[YAP_ARD_summ$enrich <0,]$gene



library(UpSetR)
#upset(fromList(ARD_upset),nsets = 7,  order.by = "freq", keep.order = T,  sets = rev(c("AA_organoid", "AA_induced", "AA_induced_stem2", "Creb_promoter", "Creb_enhancer", "Pan_Tead_promoter", "Pan_Tead_enhancer")))
upset(fromList(ARD_upset),nsets = 5, nintersects = 20,  order.by = "freq", keep.order = T,number.angles = 0, mb.ratio = c(0.7, 0.3),  sets = rev(c("CREB_Up", "CREB_Down", "YAP_Up", "YAP_Down")))

Reduce(intersect, list(AA_induced_invivo, Creb_promoter, Creb_enhancer, Pan_Tead_promter, Pan_Tead_enhancer))




