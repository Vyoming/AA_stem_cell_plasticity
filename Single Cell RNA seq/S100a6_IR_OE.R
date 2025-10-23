#LGR5 Sort DESEQ2 analysis
library("tidyverse")
library("AnnotationHub")
library("DESeq2")
library('apeglm')
library('reshape')
library('EnhancedVolcano')
library('fgsea')
library('biomaRt')

theme_linedraw2 = theme_linedraw() + theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), strip.text.x=element_text(colour="black"), strip.text.y=element_text(colour="black"))
theme_vyom = theme_linedraw2 + theme(legend.position="right", legend.title=element_text(size=15), legend.text=element_text(size=14), axis.text.x = element_text(size=12, angle=-90, hjust=0, vjust=0.5), axis.text.y=element_text(size=12), axis.title=element_text(size=15), axis.title.y=element_text(vjust=1), plot.title = element_text(size=18, vjust=1.5), strip.background = element_rect(fill="#EEEEEE"), strip.text = element_text(size = 11), panel.grid.major = element_line(colour = "grey98"), panel.grid.minor = element_blank())

load("~/data/S100a6_OE_IR_bulk/deseq2.dds.RData")
dds
rld <- rlog(dds, blind=TRUE)
pcaData <- plotPCA(rld, intgroup=c("sample"), returnData = TRUE) # vsd and plotPCA are part of DESeq2 package, nothing with my example below. 
percentVar <- round(100 * attr(pcaData, "percentVar")) 

ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(sample))) + 
  geom_point(size =3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_vyom
ggsave(file="ppar_agonist_Bulk_PCA_per_sample.pdf",  width=6, height=3, units="in")

dds@colData$Treatment <- dds@colData$sample
dds@colData$Treatment <- factor(c('WT',	'OE',	'WT',	'WT',	'OE',	'OE',	'OE',	'WT',	'WT',	'OE'))


colnames(dds)
colData(dds)

levels(dds@colData$Treatment)
my_levels <- c( "WT","OE")


dds@colData$Treatment <- factor(dds@colData$Treatment,levels = my_levels)


design(dds) <- formula(~Treatment)
dds <- DESeq(dds)
resultsNames(dds)
result_levels <- c("Treatment_OE_vs_WT")

for (i in result_levels){
  results(dds)
  
  rownames(dds) <- gsub("\\.\\d*", "", rownames(dds))
  hub <- AnnotationHub()
  hubid <- "AH7799"
  anno <- hub[[hubid]]
  genemap <- tibble(gene_id=anno$gene_id,
                    symbol=anno$gene_name) %>%
    distinct()
  
  featureData <- tibble(gene_id=rownames(dds)) %>%
    left_join(genemap, by="gene_id") %>%
    mutate(symbol=case_when(is.na(symbol) ~ gene_id,
                            TRUE ~ symbol)) %>%
    dplyr::select(symbol)
  mcols(dds) <- DataFrame(mcols(dds), featureData)
  
  dds <- nbinomWaldTest(dds)
  res <- results(dds, name=i)
  res <- lfcShrink(dds, coef=i, res=res, type='apeglm')
  
  res$gene_name <- mcols(dds)$symbol
  res <- res[order(res$log2FoldChange),]
  res_format <- res %>%
    as.data.frame() %>%
    rownames_to_column(var="gene_id") %>%
    as_tibble() %>%
    rename_at(vars(-gene_id, -gene_name), ~ paste0(., ""))
  res_format <- res_format[complete.cases(res_format),]
  res_format$log2FoldChange
  #res_format <- res_format[order(-res_format$log2FoldChange),]
  write.csv(res_format, paste0(i,'_deseq2.csv'))
  
  
  #Deseq_results <- mcols(dds) 
  #Normalized_mean_counts <- assays(dds)[['counts']]
  
  #Normalized_mean_counts <- assays(dds)[['rlog']]
  #Normalized_mean_counts <- data.frame(Normalized_mean_counts)
  #Normalized_mean_counts$Gene_Name <- mcols(dds)$symbol
  #Normalized_mean_counts$Zscore <- mcols(dds)[[i]]
  #Normalized_mean_counts[complete.cases(Normalized_mean_counts),]
  #Normalized_mean_counts <- Normalized_mean_counts[order(-Normalized_mean_counts$Zscore),]
  
  #write.csv(Normalized_mean_counts, paste0(i,'_deseq2_counts.csv'))
}

#rld <- rlog(dds, blind=FALSE)
rld <- rlog(dds, blind=TRUE)
pcaData <- plotPCA(rld, intgroup=c("Treatment"), returnData = TRUE) # vsd and plotPCA are part of DESeq2 package, nothing with my example below. 
percentVar <- round(100 * attr(pcaData, "percentVar")) 

ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(Treatment))) + 
  geom_point(size =3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_vyom
ggsave(file="S100a6_OE_IR_bulk.pdf",  width=5, height=3, units="in")
#Normalized_mean_counts <- assays(rld)[['rlog']]
EnhancedVolcano(res_format, lab = res_format$gene_name, x = 'log2FoldChange', y = 'padj', title = 'Naive_Experimental_vs_Naive_Control', pCutoff = .01, FCcutoff = .5, xlim = c(-5,11), ylim = c(0,18) )
tail(res_format)
res_format <- read.csv('Treatment_PPAR3_vs_Control_deseq2.csv')

#Gene set enrichment analysis
res_format$gene <- convert_mouse_to_human_symbols(res_format$gene_name)
res_format$gene <- res_format$gene_name
res_filter <- res_format[complete.cases(res_format),]
res_filter <- res_filter[res_filter$padj<.05,]
res2 <- res_filter %>% 
  dplyr::select(gene, log2FoldChange) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene) %>% 
  summarize(avg_log2FC=mean(log2FoldChange))
res2

ranks <- deframe(res2)
head(ranks, 20)
#c2.all.v2022.1.Hs.symbols.gmt - most comprehensive all
#c5.all.v2022.1.Hs.symbols.gmt - second most comprehensive
#h.all.v2022.1.Hs.symbols.gmt - overall summary
CREB_Targets <- read_csv("CSHL Dropbox Team Dropbox/Vyom Shah/Vyom Shah/CREB1_CBP_P300/CREB_Targets_Vyom.csv")$`Creb Targets`
CREB_Targets
#WT IR : Creb Targets | NES: 3.29 | Padj: 5.20e-22
#KO IR: Creb Targets | NES: 2.78 | Padj: 9.80e-18
All <- plotEnrichment(CREB_Targets,
                      ranks) + labs(title="Creb Targets | NES: 2.78 | Padj: 9.80e-18") + theme_vyom +
  theme(plot.title = element_text(size=6),axis.title.x = element_text(size = 6), axis.title.y = element_text(size = 6), axis.text.x = element_text(size=6, angle = 45, hjust = 1, vjust= 1), axis.text.y = element_text(size=6))
All
ggsave(file = 'KO_IR_vs_WT_CREB_TARGET_enrichment.pdf', plot=All, width=3, height=3, units="in")


pathways.hallmark <- list(CREB_Targets)
names(pathways.hallmark) <- c('CREB_Targets')

pathways.hallmark <- gmtPathways("Documents/c5.all.v2022.1.Hs.symbols.gmt")
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -log2err) %>% 
  arrange(padj) %>% 
  DT::datatable()
fgseaResTidier <- fgseaResTidy[which(fgseaResTidy$pval < .05 ),]
fgseaResTidier %>% 
  dplyr::select(-leadingEdge, -ES, -log2err) %>% 
  arrange(padj) %>% 
  DT::datatable()
ggplot(fgseaResTidier, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
gsea_output <- as.data.frame(fgseaResTidier)
gsea_output$leadingEdge <- NULL
write.csv(gsea_output,'PPAR3_v_contr_GSEA.csv')
colnames(gsea_output)
#pathways <- c('GOCC_KERATIN_FILAMENT', 'GOBP_HUMORAL_IMMUNE_RESPONSE', 'GOBP_REGULATION_OF_VASCULATURE_DEVELOPMENT', 'GOBP_RESPONSE_TO_LIPID', 'GOBP_ANATOMICAL_STRUCTURE_FORMATION_INVOLVED_IN_MORPHOGENESIS', 'GOBP_DEFENSE_RESPONSE', 'GOBP_CELL_ADHESION', 'GOBP_TISSUE_DEVELOPMENT', 'GOBP_REGULATION_OF_CELL_DEATH', 'GOCC_CATION_CHANNEL_COMPLEX')      
pathways <- c('GOCC_EXTERNAL_ENCAPSULATING_STRUCTURE', 'GOBP_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND', 'GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT', 'GOBP_HUMORAL_IMMUNE_RESPONSE', 'GOBP_REGULATION_OF_VASCULATURE_DEVELOPMENT', 'GOMF_STRUCTURAL_MOLECULE_ACTIVITY', 'GOBP_DEFENSE_RESPONSE', 'GOBP_CELL_ADHESION', 'GOCC_KERATIN_FILAMENT', 'GOBP_ANATOMICAL_STRUCTURE_FORMATION_INVOLVED_IN_MORPHOGENESIS')      

pathwayColors <- c( "#FDC926FF", "#FA9E3BFF", "#ED7953FF", "#D8576BFF", "#BD3786FF", "#9C179EFF", "#7301A8FF", "#47039FFF")
gsea_output1 <- gsea_output[order(-gsea_output$NES),]
gsea_output1 <- gsea_output1[ gsea_output1$pathway %in% pathways,]
gsea_output1$pathway <- factor(gsea_output1$pathway, levels = gsea_output1$pathway)
All <- ggplot(data = gsea_output1,aes(x=pathway,y=NES)) +  theme_vyom +
  theme(axis.text.y = element_text(size = 3) ,axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black", size = .75, linetype = 'solid')) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=NES, color = padj))+
  geom_point(color = 'black' ,size = 1.2 ) +
  geom_point(aes(color = padj) ,size = 1 ) + geom_vline(xintercept = 0) +
  scale_color_gradientn(colors=pathwayColors, breaks = c(1,.05, .01, 0),labels =  c('1','.05', '.01', '0'),limits = c(0,1), trans = scales::boxcox_trans(.35)) +
  coord_flip() + 
  scale_y_continuous(limits = c(0,2.2),expand = expansion(mult = c(0, 0)), breaks = scales::breaks_extended(n = 3)) +
  labs(y= "Normalized Enrichment Score", x="Pathway") 
All + geom_vline(xintercept = 0)
ggsave(file = paste0('PPAR_agonist_GSEA.pdf'), plot=All, width=5, height=3, units="in")



#volcano plots

#Gene set enrichment analysis
comparisons_files <- c("Treatment_OE_vs_WT_deseq2.csv")
#volcano
remove_specific_strings <- function(input_string, strings_to_remove) {
  pattern <- paste(strings_to_remove, collapse = "|")
  cleaned_string <- gsub(pattern, "", input_string)
  return(cleaned_string)
}
res_format$padj
for (i in comparisons_files) {
  res_format <- read.csv(paste0('~/',i))
  res_format <- res_format[complete.cases(res_format),]
  if (max(-log10(res_format$padj)) < 321){
    ylim_num <- max(-log10(res_format$padj))} else {
      ylim_num <- 320}
  strings_to_remove <- c("_deseq2.csv")
  k <- remove_specific_strings(i, strings_to_remove)
  
  EnhancedVolcano(res_format, lab = res_format$gene_name, x = 'log2FoldChange', y = 'padj', title = k,
                  pCutoff = .05, FCcutoff = 0.1, xlim = c(min(res_format$log2FoldChange)-.05,max(res_format$log2FoldChange)+.05),ylim = c(0,ylim_num+5),
                  gridlines.minor = FALSE, gridlines.major = FALSE)
  ggsave(file = paste0(k,'_DE_volcano.pdf'), width=6, height=6, units="in")
}
max(-log10(res_format$padj), na.rm = TRUE)





#GEO Term enrichment
library(enrichR)
dbs <- listEnrichrDbs()
comparisons_files <- c('Treatment_OE_vs_WT_deseq2.csv')

for (i in comparisons_files) {
  res_format <- read.csv(paste0('~/',i))
  res_format <- res_format[complete.cases(res_format),]
  strings_to_remove <- c("_deseq2.csv")
  k <- remove_specific_strings(i, strings_to_remove)
  dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023", "MSigDB_Hallmark_2020")
  enriched <- enrichr(unique(res_format[(res_format$padj < .05 & res_format$log2FoldChange > 0.01 ),]$gene_name), dbs)
  Enriched1 <- rbind(enriched[[1]],enriched[[2]],enriched[[3]],enriched[[4]] )
  Enriched_filter <- Enriched1[Enriched1$Adjusted.P.value < .05,]
  plotEnrich(Enriched_filter, showTerms = 100, numChar = 100, y = "Count", orderBy = "Adjusted.P.value")
  write.csv(Enriched_filter, paste0(k,'_enrichR_UP_gsea.csv'))
}


#lets make some dotplots:
comparisons_files <- c('Treatment_OE_vs_WT_deseq2.csv')
res_format <- read.csv(paste0('~/',i))
res_format <- res_format[complete.cases(res_format),]


Heatmap_data <- data.frame(res_format)
OE_gene <- c('Clu', 'Ptgs2', 'Anxa2', 'Msln', 'Myof', 'Sox9')
Heatmap_data <- Heatmap_data[Heatmap_data$gene_name %in% OE_gene,]
Heatmap_data <- Heatmap_data[,c('log2FoldChange', 'gene_name', 'padj')]

pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
Heatmap_data$dataset <- 'OE'
Heatmap_data$Log_P <- -log10(Heatmap_data$padj)

ggplot(Heatmap_data, aes(x = dataset, y = gene_name, size = Log_P ,color= log2FoldChange)) + geom_point() + theme_vyom + 
  scale_color_gradientn(name = "Log2 FC", colors = pathwayColorsDiff, limits = c(-2,2)) + labs(size="-log(P)") +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

ggsave(file="S100a6_OE_IR_gene.pdf",  width=2.15, height=2.5, units="in")
 
library(reshape)
Heatmap_data <-melt(Heatmap_data)
pathwayColors =rev(viridis::magma(10))
All_genes <- c('H2-q1','H2-T3','H2-DMb2','H2-T10','H2-Ea-ps','H2-Q10','H2-D1','Klra1','Klra4','Klrk1', 'Ctsh','Ctse','Lyz2','Hist2h3c2','Hist1h2ag','Hist4h4','Hist1h2ak','Hist1h2ap','Hist2h2aa2')
Heatmap_data$gene_name <- factor(Heatmap_data$gene_name,levels = All_genes)
ggplot(Heatmap_data, aes(gene_name, variable, fill= value)) + geom_tile() + scale_fill_gradientn(name = "Log2FC", colors = pathwayColors) + theme(panel.background = element_blank())

gsave(file="Bulk_logfc_scale.pdf",  width=15, height=2, units="in")
IR_KO_vs_IR_WT_enrichR_gsea.csv
Enriched_filter <- read.csv(paste0('~/Treatment_OE_vs_WT_enrichR_UP_gsea.csv'))
pathways <-  c('Cell-Cell Junction (GO:0005911)', 'p53 Pathway', 'Androgen Response', 'Reactive Oxygen Species Pathway', 'PI3K/AKT/mTOR  Signaling', 'Myc Targets V1', 'KRAS Signaling Up')

pathwayColors <- rev(c( "#FDC926FF", "#FA9E3BFF", "#ED7953FF", "#D8576BFF", "#BD3786FF", "#9C179EFF", "#7301A8FF", "#47039FFF"))
#pathwayColors <- c("#440154", "#482475", "#414487", "#355f8d", "#2a788e", "#21918c", "#22a884", "#44bf70", "#7ad151", "#bddf26", "#fde725")
#pathwayColors <- c("#ff0000", "#ff9700", "#d1ff00", "#3aff00", "#00ff5c", "#00fff3", "#0074ff", "#2200ff", "#b900ff", "#47039FFF")
Enriched_filter <- Enriched_filter[order(-Enriched_filter$Combined.Score),]
Enriched_filter <- Enriched_filter[ Enriched_filter$Term %in% pathways,]
Enriched_filter$Term <- factor(Enriched_filter$Term, levels = Enriched_filter$Term)
Enriched_filter$p_adj <- -log(Enriched_filter$Adjusted.P.value)
Enriched_filter$Overlap
#Enriched_filter$num_genes <- c(8,80,22,9,13,13,60,24,16,45)
#Enriched_filter$num_genes <- c(32,32,9,20,5)

Enriched_filter <- Enriched_filter[order(Enriched_filter$Odds.Ratio),]
Enriched_filter$Term <- factor(Enriched_filter$Term, levels = Enriched_filter$Term)
ggplot(Enriched_filter, aes(x = Term, y = Odds.Ratio, size = p_adj ,color= num_genes)) + geom_point() + theme_vyom + 
  scale_color_gradientn(name = "# of Genes", colors = pathwayColors, limits = c(-0,100), trans = scales::boxcox_trans(.6), breaks = c(0, 10, 25, 50, 100),labels =  c('0', '10', '25', '50', '100')) + coord_flip() + labs(size="-log(p)") +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =  element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.y=element_blank()) + ylab('Odds Ratio')       
ggsave(file="S100a6OE_IR_GSEA_dotplot.pdf",  width=6, height=3, units="in")


