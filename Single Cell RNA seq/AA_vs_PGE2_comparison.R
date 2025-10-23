#PGE2 vs AA comparison, In-Vitro and In-Vivo

ARD_DE <- read.csv('ARAsco_DE.csv')
PGE_DE <- read.csv('PGE2_DE.csv')
AA_DE <- read.csv('AA_DE.csv')


ARD_DE <- ARD_DE[ARD_DE$p_val_adj <.05,]
PGE_DE <- PGE_DE[PGE_DE$p_val_adj <.05,]
AA_DE <- AA_DE[AA_DE$p_val_adj <.05,]

AA_PGE2_up <- intersect(AA_DE[AA_DE$avg_logFC>.1,]$X, PGE_DE[PGE_DE$avg_logFC>.1,]$X)


ARD_upset <- list()
ARD_upset$AA_invitro <- AA_DE[AA_DE$avg_logFC>.1,]$X
ARD_upset$AA_invivo <- ARD_DE[ARD_DE$avg_logFC>.1,]$X
ARD_upset$PGE2_invitro <- PGE_DE[PGE_DE$avg_logFC>.1,]$X

pastel_colors <- c( "#B0E2FF", "#CAB8FF", "#FDFD96" ,  "#A6FBB2", "#FF9999")

library(ggvenn)
ggplot() + geom_venn(ARD_upset, fill_color = pastel_colors)+ theme_void()
ggsave(file = paste0('ARD_PGE2_Shared_Targets.pdf'), width=5, height=4, units="in")
ggsave(file = paste0('ARD_invivo_vitro_Shared_Targets.pdf'), width=5, height=4, units="in")

names(ARD_upset)
library("ggVennDiagram")
# Default plot
ggVennDiagram(ARD_upset) + scale_fill_gradientn(colors = pastel_colors)
PGE_unique <- setdiff(ARD_upset$PGE2_invitro, ARD_upset$AA_invitro)


Reduce(intersect, list(ARD_upset$PGE2_invitro, ARD_upset$AA_invitro, ARD_upset$AA_invivo))

ARD_invivo_unique <- Reduce(setdiff, list(ARD_upset$AA_invivo, ARD_upset$PGE2_invitro, ARD_upset$AA_invitro))

dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023", "MSigDB_Hallmark_2020")
enriched <- enrichr(ARD_invivo_unique, dbs)
Enriched1 <- rbind(enriched[[1]],enriched[[2]],enriched[[3]],enriched[[4]] )
write.csv(Enriched1, paste0('AA_invivo_up_ENRICHR.csv'))
Enriched1
Enriched_filter <- Enriched1[Enriched1$Adjusted.P.value < .05,]
plotEnrich(Enriched_filter, showTerms = 100, numChar = 100, y = "Count", orderBy = "Adjusted.P.value")
ARD_enrichr <- Enriched_filter

setdiff(PGE2_enrichr$Term,AA_enrichr$Term)
setdiff(AA_enrichr$Term, PGE2_enrichr$Term)
pathways <- c("Interferon Gamma Response", "TNF-alpha Signaling via NF-kB", "Respiratory Chain Complex III (GO:0045275)", "Hydrogen Peroxide Catabolic Process (GO:0042744)", "Superoxide Metabolic Process (GO:0006801)")             
pathways <- c("Positive Regulation Of Mitotic Cell Cycle Phase Transition (GO:1901992)", "Protein Localization To Nucleus (GO:0034504)", "Pyruvate Metabolic Process (GO:0006090)", "Epiboly Involved In Wound Healing (GO:0090505)", "Long-Chain fatty-acyl-CoA Metabolic Process (GO:0035336)", "Fatty Acid Catabolic Process (GO:0009062)", "Fatty Acid Beta-Oxidation (GO:0006635)", "Intestinal Absorption (GO:0050892)", "Clathrin-Dependent Endocytosis (GO:0072583)")                
pathways <- c('Aerobic Respiration (GO:0009060)', 'Adipogenesis', 'DNA Repair', 'G2-M Checkpoint', 'Myc Targets V1')

Enriched_filter <- PGE2_enrichr
Enriched_filter <- AA_enrichr
Enriched_filter <- ARD_enrichr
pathwayColors <- rev(c( "#FDC926FF", "#FA9E3BFF", "#ED7953FF", "#D8576BFF", "#BD3786FF", "#9C179EFF", "#7301A8FF", "#47039FFF"))
Enriched_filter <- Enriched_filter[order(-Enriched_filter$Combined.Score),]
Enriched_filter <- Enriched_filter[ Enriched_filter$Term %in% pathways,]
Enriched_filter$Term <- factor(Enriched_filter$Term, levels = Enriched_filter$Term)
Enriched_filter$p_adj <- -log(Enriched_filter$Adjusted.P.value)
Enriched_filter$Overlap
#ched_filter$num_genes <- c( 9, 7, 11, 20, 40, 22, 17, 20, 10)
Enriched_filter$num_genes <- c( 20, 8, 8, 7, 7)
Enriched_filter <- Enriched_filter[order(Enriched_filter$Odds.Ratio),]
Enriched_filter$Term <- factor(Enriched_filter$Term, levels = Enriched_filter$Term)
ggplot(Enriched_filter, aes(x = Term, y = Odds.Ratio, size = p_adj ,color= num_genes)) + geom_point() + theme_vyom + 
  scale_color_gradientn(name = "# of Genes", colors = pathwayColors, limits = c(-0,20)) + coord_flip() + labs(size="-log(p)") +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =  element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.y=element_blank()) + ylab('Odds Ratio')       
ggsave(file="AA_invivo_unique_GSEA_dotplot.pdf",  width=6, height=3, units="in")







