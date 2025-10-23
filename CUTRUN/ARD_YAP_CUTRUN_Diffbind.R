#Vyom Shah 4/2024 - ARD Project: YAP/TEAD Cut&Run
#CUT&RUN downstream analysis
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

samples <- read.csv('~/analysis/YAP_ARD_CUTRUN_Samplesheet.csv')
samples

res = dba(sampleSheet=samples, config=data.frame(RunParallel=TRUE), scoreCol=5)

# counting reads mapping to intervals (peaks)
res.cnt = dba.count(res)
res.cnt1 <- res.cnt
# applying  TMM normalisation
res.norm=dba.normalize(res.cnt, normalize=DBA_NORM_TMM)

# inspecting the object: notice the FRiP values!
res.norm

#update samplesheet
res.cnt$samples <- samples
res.cnt
res.cnt1
# setting the contrast
res.cnt2 = dba.contrast(res.norm, design = '~Condition')

# inspecting the object: how many contrasts were set in the previous step
res.cnt2$DESeq2

# plotting the correlation of libraries based on normalised counts of reads in peaks
pdf("correlation_libraries_normalised.pdf")
plot(res.cnt)
dev.off()

# PCA scores plot: data overview
pdf("PCA_normalised_libraries.pdf")
dba.plotPCA(res.cnt,DBA_CONDITION,label=DBA_CONDITION)
dev.off()

# we will skip generating greylists (regions of high signal in input samples) because of time - it is recommended to perform this step in your own analyses though!
res.cnt2$config$doGreylist=FALSE

# performing analysis of differential binding

res.norm=dba.normalize(res.cnt, normalize=DBA_NORM_TMM)
res.cnt2 = dba.contrast(res.norm, categories=DBA_CONDITION, minMembers=2)

res.cnt3 = dba.analyze(res.cnt2)
dds <- dba.analyze(res.cnt3, bRetrieveAnalysis=DBA_DESEQ2)

res.cnt3 = dba.analyze(res.cnt)
dds1 <- dba.analyze(res.cnt3, bRetrieveAnalysis=DBA_DESEQ2)

# inspecting the object: which condition are most alike, which are most different, is this expected?
dba.show(res.cnt3, bContrasts = T)

# correlation heatmap  using only significantly differentially bound sites
# choose the contrast of interest
pdf("correlation_YAP1_ARD.pdf")
plot(res.cnt3, contrast=1)
dev.off()

# boxplots to view how read distributions differ between classes of binding sites
# are reads distributed evenly between those that increase binding affinity?
pdf("Boxplot_YAP1_ARD.pdf")
dba.plotBox(res.cnt3, contrast=1)
dev.off()

# plotting overlaps of sites bound by REST in different cell types
pdf("binding_site_overlap_ARD.pdf")
dba.plotVenn(res.cnt3, contrast=1)
dev.off()

pdf("YAP1_ARD_PCA.pdf")
dba.plotPCA(res.cnt3,contrast=1,label=DBA_CONDITION)
dev.off()

pdf("YAP1_ARD_plotprofile.pdf")
dba.plotProfile(res.cnt3,contrast=1)
dev.off()

#saveRDS(res.cnt3, file = "./data/YAPKO_Processing.rds")
#res.cnt3 <- readRDS('~/data/YAPKO_Processing.rds', refhook = NULL)
res.cnt3$samples
res.cnt3
# extracting differentially occupied sites in a GRanges object
res.db1 = dba.report(res.cnt3)
head(res.db1)
gr <- res.db1
data(TSS.mouse.GRCm38)
macs.anno <- annotatePeakInBatch(gr, AnnotationData=TSS.mouse.GRCm38)

macs.anno <- addGeneIDs(annotatedPeak=macs.anno,
                        orgAnn="org.Mm.eg.db",
                        IDs2Add="symbol")
Anno_output <- data.frame(macs.anno)
write.csv(Anno_output, 'YAP1_ARDvsWT_DiffPeak.csv')
Anno_output$Fold
#Conduct Gene Set enrichment analysis for this dataset
library(nichenetr)
library(fgsea)

Anno_output$gene <- convert_mouse_to_human_symbols(Anno_output$symbol)
res_filter <- Anno_output[complete.cases(Anno_output),]
res_filter <- res_filter[res_filter$p.value <.05,]
res2 <- res_filter %>% 
  dplyr::select(gene, Fold) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene) %>% 
  summarize(Fold=mean(Fold))
res2

ranks <- deframe(res2)
head(ranks, 20)
#c2.all.v2022.1.Hs.symbols.gmt - most comprehensive all
#c5.go.v2023.1.Hs.symbols.gmt
#c2.cp.v2023.1.Hs.symbols.gmt 
#c5.all.v2022.1.Hs.symbols.gmt - second most comprehensive
#h.all.v2022.1.Hs.symbols.gmt - overall summary

pathways.hallmark <- gmtPathways("Documents/c2.cp.v2023.1.Hs.symbols.gmt")
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
fgseaResTidier <- fgseaResTidy[which(fgseaResTidy$pval < .01 ),]
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
write.csv(gsea_output, paste0('YAPKO_CUTRUN_gsea_pathways.csv'))

pathways <- c('WP_NOTCH_SIGNALING', 'WNT_SIGNALING', 'WP_PROSTAGLANDIN_SIGNALING', 'WP_FATTY_ACID_BETAOXIDATION', 'BIOCARTA_MONOCYTE_PATHWAY', 'BIOCARTA_LYMPHOCYTE_PATHWAY', 'REACTOME_PKA_MEDIATED_PHOSPHORYLATION_OF_YAP', 'REACTOME_GLYCOGEN_METABOLISM', 'REACTOME_GPER1_SIGNALING', 'PID_PS1_PATHWAY')      

pathwayColors <- c( "#FDC926FF", "#FA9E3BFF", "#ED7953FF", "#D8576BFF", "#BD3786FF", "#9C179EFF", "#7301A8FF", "#47039FFF")
gsea_output1 <- gsea_output[order(-gsea_output$NES),]
gsea_output1 <- gsea_output1[ gsea_output1$pathway %in% pathways,]
gsea_output1$pathway <- factor(gsea_output1$pathway, levels = gsea_output1$pathway)
All <- ggplot(data = gsea_output1,aes(x=pathway,y=NES)) +  theme_vyom +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size = 6) ,axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black", size = .75, linetype = 'solid')) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=NES, color = pval))+
  geom_point(color = 'black' ,size = 1.2 ) +
  geom_point(aes(color = pval) ,size = 1 ) + geom_vline(xintercept = 0) +
  scale_color_gradientn(colors=rev(pathwayColors), breaks = c(1,.05, .01, 0),labels =  c('1','.05', '.01', '0'),limits = c(0,1), trans = scales::boxcox_trans(.35)) +
  coord_flip() + 
  scale_y_continuous(limits = c(-2.5,0),expand = expansion(mult = c(0, 0)), breaks = scales::breaks_extended(n = 3)) +
  labs(y= "Normalized Enrichment Score", x="Pathway") 
All <- All + geom_hline(yintercept = 0, linewidth= .25, linetype = 'dashed' )
ggsave(file = paste0('YAPko_GSEA.pdf'), plot=All, width=4.75, height=3, units="in")



#use MonaLisa for motif enrichment
library(monaLisa)

#Use CHIPSeeker for further downstream analysis

dba.show(res.cnt3)
dba.show(res.cnt3)
peaks.consensus <- dba.peakset(res.cnt3, bRetrieve = T)

peaks.KO_YAP1_R1 <- peaks.consensus[res.cnt3$called[,1]==1]
peaks.KO_YAP1_R2 <- peaks.consensus[res.cnt3$called[,2]==1]

peaks.WT_YAP1_R1 <- peaks.consensus[res.cnt3$called[,5]==1]
peaks.WT_YAP1_R2 <- peaks.consensus[res.cnt3$called[,6]==1]


# adding an unified affinity scores column (re-formatting data)
peaks.KO_YAP1_R1$Score <- peaks.KO_YAP1_R1$KO_YAP1_R1
peaks.KO_YAP1_R2$Score <- peaks.KO_YAP1_R2$KO_YAP1_R2
peaks.WT_YAP1_R1$Score <- peaks.WT_YAP1_R1$WT_YAP1_R1
peaks.WT_YAP1_R2$Score <- peaks.WT_YAP1_R2$WT_YAP1_R2

# plotting coverage for replicate 1, using affinity scores as a weight for peaks height
covplot(peaks.KO_YAP1_R1, weightCol = "Score")
covplot(peaks.WT_YAP1_R1, weightCol = "Score")

# extracting peaks for each cell line present across replicates
peaks.KO_YAP1 <- peaks.consensus[res.cnt3$called[,1]==1 & res.cnt3$called[,2]==1]
peaks.WT_YAP1 <- peaks.consensus[res.cnt3$called[,5]==1 & res.cnt3$called[,6]==1]

res.cnt3$samples
# getting TSS regions
Tuft_genes <- c('Alox5ap',	'Hck',	'Lrmp',	'Avil',	'Trpm5',	'Spib',	'Rgs13',	'Ltc4s',	'Pygl',	'Sh2d7',	'Dclk1',	'Alox5',	'Pik3r5',	'Fyb',	'Vav1',	'Matk',	'Tspan6',	'Strip2',	'Pou2f3',
                '1810046K07Rik',	'Ptpn6',	'Bmx',	'Tuba1a',	'Espn',	'Plcb2',	'Ffar3',	'Ccdc109b',	'Plcg2',	'Ly6g6f',	'Hpgds',	'Pea15a',	'Ly6g6d',	'Pik3cg',	'Inpp5d',	'Ccdc28b',	'Snrnp25',	'Kctd12',	'Siglec5',	'Skap2',	'Ccdc129',	'Nebl',	'Gprc5c',	'Rgs22',	'Gfi1b',	'Hmx3',	'Cbr3',	'Pfkfb3',	'Prss53',	'Itpr2',	'Limd2',	'Cd300lf',	'Chn2',	'Smpx',	'Ptgs1',	'A4galt',	'Rac2',	'Csk',	'Slco4a1',	'Ptpn18',	'Chat',	'Hebp1',	'Ppp1r14c',	'Dgki',	'Inpp5j',	'Tppp3',	'Gng13',	'Ildr1',	'Cwh43',	'Il17rb',	'Ncf2',	'Fut2',	'Coprs',	'Ddah1',	'Tmem116',	'Sucnr1',	'Tmem176a',	'Ccrl1',	'1110007C09Rik',	'Adcy5',	'Fnbp1',	'Plk2',	'Hmx2',	'Tmem141',	'Krt23',	'Gprc5a',	'Rgs2',	'Camk2b',	'Fes',	'Bpgm',	'Acacb',	'Il13ra1',	'Zfp428',	'Ppp1r3b',	'Ccnj',	'Bcl2l14',	'Tmem229a',	'Ethe1',	'Runx1',	'Gga2',	'Apobec1',	'Serpini1',	'St6galnac6',	'Fbxl21',	'9030624J02Rik',	'Inpp5b',	'Samd14',	'Pgm2l1',	'Pla2g4a',	'Ptprc',	'Aldh2',	'Ifi27l1',	'Pnpla3',	'Jarid2',	'Rgs19',	'Reep5',	'Tiparp',	'Gnai2',	'Fam49a',	'Cacna2d2',	'Ypel2',	'Cd24a',	'Acot7',	'Svil',	'Abhd16a',	'Fam101a',	'Trim40',	'Trak1',	'Sec14l1',	'4930539E08Rik',	'Smtn',	'Galk1',	'Tbc1d1',	'Tmem176b',	'Fcna',	'Abhd2',	'Hsbp1l1',	'Slc4a8',	'Myo1b',	'Tmem38b',	'Hk1',	'Neurl1a',
                'Dmxl2',	'Bub3',	'Ptprj',	'Trib2',	'Stard5',	'Ubtd1',	'Slc41a3',	'Plekhg5',	'Rbm38',	'Fam57a',	'Eef2k',	'Cables2',	'Fbxo25',	'Ap1s2',	'1300002K09Rik',	'Ero1lb',	'Clmn',	'Fam49b',	'Cpvl',	'Prr15',	'Lpcat4',	'Tmem74b',	'Mn1',	'Eppk1',	'Samd9l',	'Tmem245',	'Glyctk',	'Aldh3a2',	'Ppp3ca',	'Cpne3',	'Slc4a7',	'Nfatc1',	'Kit',	'Fam117b',	'Nradd',	'Tmem121',	'Cpm',	'Asah1',	'Slc9a9',	'Ubl7',	'Abca3',	'Pde6d',	'Bmp2',	'Kdm4a',	'Camkk2',	'Arhgap8',	'Agt',	'Ptpra',	'Adh1',	'Dusp14',	'Clic4',	'Gimap1',	'Cpne5',	'Ceacam2',	'Zfp710',	'Gcnt1',	'B4galt5',	'Suco',	'Pim3',	'Ogdhl',	'Oas1g',	'Dcp1b',	'Myzap',	'Cdkn1a',	'Cd37',	'Brms1',	'Lrrc42',	'Pld2',	'Tmem9',	'Cpeb4',	'Ssx2ip',	'Ddah2',	'Tmem65',	'5430417L22Rik',	'2210016L21Rik',	'Msi2',	'B4galt4',	'Rabgap1l',	'Pik3r3',	'Nt5c3',	'Palld',	'AA467197',	'Pip5k1b',	'Krt18',	'Map1a',	'Lmf1',	'Arhgef28',	'Nsfl1c',	'Txndc16',	'Pstpip2',	'Ttll11',	'Exph5',	'2700086A05Rik',	'Gadd45a',	'Plekhs1',	'Fam188a',	'Jmy',	'Atat1',	'Arhgef2',	'Lmbr1',	'Rhoc',	'Card10',	'Kcnj16',	'Arhgap4',	'Acsl4',	'Rhog',	'Fam221a',	'Dynlt1b',	'C2',	'Zbtb41',	'Socs1',	'Atp6ap1',	'Fam171a1',	'Wnk2',	'Kcnd3',	'Slc27a1',	'Atxn1',	'Rabgap1',	'Myrfl',	'Crot',	'Tm4sf4',	'Ube2j1',	'Sort1',	'Lima1',	'Mov10',	'Lca5',	'Gimap9',	'Mlip',	'1110008P14Rik',	'Ckap4',	'Tor4a',	'Rmdn1',	'Oas2',	'Dsp',	'Sox9',	'Osbpl3',	'Kif21b',	'Tbcb',	'Arap2',	'Casp3',	'Enc1',	'Il25',	'Lman2l',	'Zmiz1',	'Nav2',	'Atp2a3',	'Gimap8',	'Folr1',	'Fn1',	'Hspa4l',	'Sufu',	'Atp8a1',	'Vps53',	'Rgs14',	'Gm17660',	'Pdcl',	'Shkbp1',	'Oas1a',	'Pkp1',	'Ccdc23',	'Il4ra',	'1700112E06Rik',	'Dvl1',	'Zfhx3',	'Adam22',	'Gramd1c',	'Tmem45b',	'Unc5b',	'Mical3',	'Kctd13',	'Ak7',	'Tcta',	'Nek7',	'D730039F16Rik',	'Plekho2',	'Myo6',	'Chdh',	'Opn3',	'Tle3',	'Ttll10',	'Strada',	'Ypel3',	'Cmip',	'Cachd1',	'Pigc',	'Atp6v1d',	'Rdx',	'S100a11',	'Spa17',	'Gimap5',	'Cystm1',	'Zdhhc17',	'Lect2',	'Vdac3',	'Hspb11',	'Gm4952',	'Slc16a2',	'Abhd5',	'Rhbdf1',	'Cblb',	'Nfe2l3',	'Pla2g16',	'Sept8',	'Gpcpd1',	'Psd3',	'Anxa11',	'Slc25a12',	'Ehf',	'Akr1b10',	'Dapp1',	'Vmn2r26',	'Esyt1',	'Ppt1',	'Cd47',	'Chi3l1',	'Mical1',	'Gna14',	'Pacs2',	'Lyn',	'Rmnd5a',	'Ankrd12',	'BC022687',	'Rit1',	'Camta2',	'Mocs2',	'Usp49',	'Nrbp2',	'Ifnar2',	'Epha4',	'Arl5a',	'Rgl2',	'St18',	'BC016579',	'Tead1',	'Enpp4',	'Tmem158',	'Tnfaip3',	'Gys1',	'Hivep2',	'Cap1',	'Slc4a2',	'Map4k4',	'Desi1',	'H2-D1',	'Man2a1',	'Cyp17a1',	'Cyhr1',	'Morf4l1',	'Mllt4',	'Phf17',	'Stox2',	'Hist3h2a',	'Hdac6',	'Prox1',	'Dtnb',	'Lrch4',	'Spire2',	'Klf6',	'Rab5b',	'Anxa4',	'Rab4b',	'Iqsec1',	'Pdpk1',	'Stk40',	'Gde1',	'Mtmr11',	'Cib2',	'March2',	'Capg',	'Narf',	'Mgst3',	'Angel1',	'Bicd1',	'Ifitm1',	'Stx3',	'S100a1',	'Omd',	'0610040J01Rik',	'Arpc5',	'Homer3',	'Cdc42se1',	'Abcc3',	'Hsf2',	'Pnpla6',	'Ccdc68',	'Fryl',	'Lmtk2',	'Tas1r3',	'4931406H21Rik',	'Uspl1',	'Ajuba',	'Kalrn',	'Basp1',	'Pip5kl1',	'Slc26a2',	'Atp2b2',	'Smug1',	'Myadm',	'D330041H03Rik',	'Wdfy2',	'Trim38',	'Arf3',	'Scand1',	'Dpysl2',	'Ndufaf3',	'Sik1',	'Wdr7',	'Sfxn3',	'Kcnq4',	'Mll1',	'Hsbp1',	'Calml4',	'Atf7ip',	'Gpr137b-ps',	'Hap1',	'Kctd15',	'Prcp',	'9430023L20Rik',	'Gmip',	'Cmtm3',	'Madd',	'Krt222',	'Nsf',	'Klhl28',	'Pparg',	'Eml3',	'Phlda1',	'P2rx1',	'Pde9a',	'Otud7b',	'Tfpi2',	'Rilpl2',	'Klf3',	'Gyg',	'4930455F23Rik',	'Armcx1',	'Lzts2',	'Plek',	'Vamp8',	'Stat2',	'Znf512b',	'Ptplad1',	'1110058L19Rik',	'Tmem160',	'Tmem51',	'Cdhr5',	'Stk38',	'Atp13a2',	'Nptn',	'Sirt5',	'Gabarapl2',	'Nudt14',	'2010111I01Rik',	'Alkbh7',	'Slc18a3',	'4930427A07Rik',	'Ttll7',	'Acss2',	'Siae')  
tuft_genes <- c('Lrmp',	'Alox5ap',	'Rgs13',	'Sh2d6',	'Ltc4s',	'Avil',	'Hck',	'Dclk1',	'Snrnp25',	'Cd24a',	'Trpm5',	'Kctd12',	'Aldh2',	'Il13ra1',	'Gng13',	'Tmem176a',	'Skap2',	'Ptpn6',	'Ly6g6f',	'Fyb',	'Adh1',	'Tmem176b',	'Hpgds',	'Reep5',	'Ptpn18',	'Spib',	'Bpgm',	'Galk1',	'Matk',	'Tuba1a',	'1810046K07Rik',	'Hmx2',	'Ccdc28b',	'Ethe1',	'Limd2',	'Sh2d7',	'Ccdc109b',	'Tspan6',	'Smpx',	'Vav1',	'Ly6g6d',	'Pik3r5',	'Nebl',	'Plcg2',	'Rbm38',	'Vdac3',	'Krt18',	'Asah1',	'Cd47',	'Krt23',	'Bcl2l14',	'Lima1',	'Pygl',	'Itpr2',	'Inpp5j',	'Pea15a',	'Rac2',	'Pou2f3',	'Atp2a3',	'Bmx',	'Acot7',	'Gnai2',	'Alox5',	'Ppp3ca',	'Ptgs1',	'Calm2',	'Zfp428',	'Tmem141',	'Myo1b',	'Siglecf',	'Pla2g4a',	'Inpp5b',	'Fam221a',	'Bub3',	'Arpc5',	'Pla2g16',	'1110007C09Rik',	'Gimap1',	'Coprs',	'Lect2',	'Nrgn',	'Agt',	'Ffar3',	'Tmem45b',	'Ccdc23',	'Rgs2',	'Mlip',	'Csk',	'2210016L21Rik',	'St6galnac2',	'Ildr1',	'Gprc5c',	'Mocs2',	'Nrep',	'Pik3cg',	'Malat1',	'Sec14l1',	'Ndufaf3',	'Inpp5d',	'Pim3',	'Tmem9',	'Gga2',	'Nt5c3')
tuft_genes <- c("Trpm5" ,"Gfi1b" ,"Il25" ,"Alox5ap" ,"Lrmp" ,"Rgs13" ,"Ltc4s" ,"Adh1", 'Dclk1', 'Dusp4','Pten','Ltc4s','Stub1','Vil1')
# Example R vector of MHC related mouse genes
mhc1_genes <- c("H2-K1", "H2-K2", "H2-D1", "H2-D2", "H2-L", "H2-Q1", "H2-Q2",'H2-Q3','H2-Q4','H2-Q6', 'H2-M3', "H2-T3", 'H2-DMa',
                "H2-Q10", "H2-M2", "H2-M3", "H2-Oa", "H2-Ob", 'H2-T25', "H2-T23",'H2-T27', 'H2-DMb1', 'H2-DMb2', 'H2-Ob',
                "H2-T24", 'H2-T10', 'H2-Q10')
mhc2_genes <- c('H2-Ea', "H2-Eb1", "H2-Eb2", "H2-Ab1", "H2-Ab2","H2-Aa", 'H2-D1','Ciita','Cd74')
Stem_genes <- c('Lgr5',	'Gkn3',	'Ascl2',	'Olfm4',	'Rgmb',	'Igfbp4',	'2210407C18Rik',	'Jun',	'Pdgfa',	'Soat1',	'Tnfrsf19',	'Cyp2e1',	'Fstl1',	'H2-Eb1',	'Ifitm3',	'Prelp',	'Scn2b',	'A930009A15Rik',	'H2-Ab1',	'Slc1a2',	'Cd74',	'Sp5',	'Noxa1',	'Rgcc',	'Sorbs2',	'Sectm1b',	'H2-Aa',	'Cdo1',	'Slc14a1',	'Clca2',	'Tifa',	'Pls3',	'Hmgcs2',	'Arid5b',	'Agr3',	'Slc12a2',	'Rassf5',	'Rnf43',	'Nrn1',	'Lamb3',	'Cd44',	'Axin2',	'Slc27a2',	'Afap1l1',	'Ccdc3',	'Lrig1',	'Noxo1',	'Cdk6',	'Amica1',	'Tgif1',	'Tns3',	'Nr2e3',	'Efna4',	'Rnf32',	'Prss23',	'2010009K17Rik',	'Smoc2',	'Mecom',	'Esrrg',	'Aqp1',	'Znrf3',	'Grb7',	'Phgdh',	'2410004N09Rik',	'Clca4',	'Aqp4',	'Lcp1',	'E030011O05Rik',	'Snhg1',	'BC064078',	'Car12',	'Zbtb38',	'Cdca7',	'Fam13a',	'Shisa2',	'Dtx4',	'Slc19a2',	'Fam115c',	'Mir703',	'Cd14',	'Mettl20',	'Myo9a',	'App',	'Clic6',	'Wee1',	'2410006H16Rik',	'Lancl1',	'1500012F01Rik',	'Casp12',	'Sh3rf1',	'Lrp4',	'Arhgef26',	'Etv6',	'1700024F13Rik',	'Cttnbp2',	'Slc16a13',	'Htr4',	'Pdxk',	'Immp2l',	'Rps15a-ps6',	'Rps15a-ps4',	'Nap1l1',	'Sdc4',	'Epn3',	'Sipa1l1',	'Wfdc15b',	'Zfp341',	'Ngef',	'Nrg4',	'Csad',	'Rpl34-ps1',	'Rin2',	'Cd81',	'Irf2bp2',	'Sesn3',	'Phlpp1',	'Yap1',	'Mfge8',	'Zfp825',	'Itga1',	'Pcdh8',	'Vdr',	'Kcnq1',	'Slc28a2',	'Zfp36l1',	'Urod',	'Rgs12',	'Nfib',	'Sdsl',	'Nfia')
Stem2_genes <- c('S100a6', 'Clu', 'Ly6a', 'Anxa3','Anxa2','Msi1')


promoter <- getPromoters(TxDb=txdb, upstream=5000, downstream=5000)
macs.anno <- annotatePeakInBatch(promoter, AnnotationData=TSS.mouse.GRCm38)

macs.anno <- addGeneIDs(annotatedPeak=macs.anno,
                        orgAnn="org.Mm.eg.db",
                        IDs2Add="symbol")
macs.anno <- macs.anno[macs.anno$symbol %in% mhc2_genes]

promoter <- promoter[promoter@ranges %in% macs.anno@ranges]

# calculating tagMatrix
newStyle <- mapSeqlevels(seqlevels(peaks.KO_YAP1_R1), "UCSC")
peaks.KO_YAP1 <- renameSeqlevels(peaks.KO_YAP1, newStyle)
newStyle <- mapSeqlevels(seqlevels(peaks.WT_YAP1), "UCSC")
peaks.WT_YAP1 <- renameSeqlevels(peaks.WT_YAP1, newStyle)

tagMatrix.1 <- getTagMatrix(peaks.KO_YAP1, windows=promoter)
tagMatrix.2 <- getTagMatrix(peaks.WT_YAP1, windows=promoter)

tagMatrixList <- list(KO=tagMatrix.1, WT=tagMatrix.2)

#tagHeatmap(tagMatrixList)
plotAvgProf(tagMatrixList, xlim=c(-5000, 5000))

#write promoter granges as bed
AA_induced_organoid_scrna <- c('Fabp1',	'S100a6',	'Zg16',	'Dmbt1',	'Lgals3',	'Agr2',	'Reg1',	'Akr1b8',	'Anxa2',	'Cbr3',	'Ldha',	'Hmgcs2',	'2210407C18Rik',	'Gsta1',	'Krt19',	'Pkm',	'Acadl',	'AA467197',	'Apoa4',	'Pgam1',	'Gm20594',	'Sis',	'Apoa1',	'Reg4',	'Krt8',	'Slc11a2',	'Gna11',	'Anxa13',	'Gsta4',	'Guca2b',	'Prdx6',	'Serpinb6a',	'Nme2',	'Tubb4b',	'Krt18',	'Lypd8',	'Rbp2',	'Mttp',	'Nfkbia',	'Cldn3',	'Rpl15',	'Dstn',	'Tuba1b',	'Ckmt1',	'Gpx2',	'Krt7',	'Prss32',	'Ube2m',	'Alpi',	'Gsta3',	'Wfdc2',	'Slc16a3',	'Anxa3',	'Capg',	'Aldob',	'Cbr1',	'Gstm1',	'Arpc4',	'Ddah1',	'Egln3',	'Tomm22',	'Tfrc',	'Mdh2',	'Actg1',	'Eif5a',	'Akr1c19',	'Pgk1',	'Ptgr1',	'Cmas',	'Aldh1a1',	'Ezr',	'Cct7',	'Aprt',	'Tm4sf4',	'2610528A11Rik',	'Eif6',	'Ankrd37',	'Tuba4a',	'Ddx3x',	'Acaa2',	'Galk1',	'Tm4sf5',	'Lars2',	'Nfe2l2',	'Oit1',	'Aldoa',	'Ddx39',	'Plpp2',	'Arhgdig',	'Vdac2',	'Ctsl',	'Gclm',	'Tpi1',	'Prdx1',	'Gsr',	'Gpd1',	'Tspan8',	'Hspd1',	'Lgals4',	'Phb',	'Ctsz',	'Actb',	'Me1',	'St3gal4',	'Tmem54',	'Gm10116',	'Acp5',	'Cldn15',	'Ero1l',	'Akr1c13',	'Hras',	'Sri',	'Tmem237',	'Eif4a1',	'Id1',	'Areg',	'Arf1',	'Pdcd6',	'Tomm20',	'Pgp',	'Rpl7l1',	'Pfkp',	'Fabp2',	'Uqcrfs1',	'Sfn',	'Muc2',	'Plxnb2',	'Gmds',	'Snrpb',	'1110008F13Rik',	'Pfkl',	'Lman2',	'Taldo1',	'Suclg1',	'Cyc1',	'Ppp1ca',	'Cda',	'Tspan1',	'Slc25a3',	'Cybrd1',	'Ech1',	'Ppa1',	'Uba52',	'Eif4g1',	'Srsf7',	'Esd',	'Rpl6',	'Cyba',	'Lrrc59',	'Ier3',	'Ckb',	'Uqcrc1',	'Ran',	'Polr2e',	'Cd9',	'Gipc1',	'Gstp1',	'Tomm40',	'Cldn2',	'Aldh1b1',	'Ctsd',	'Stoml2',	'Rac1',	'Syngr2',	'Eci2',	'Slc25a4',	'Psmd13',	'Lamp1',	'Plcb3',	'Gsto1',	'Ak2',	'Ccnd1',	'Tmem45b',	'Lmna',	'Rab25',	'Acadvl',	'Plac8',	'Maoa',	'Oxct1',	'Sumo1',	'Adh6a',	'Map2k2',	'Ugt2b34',	'Gpi1',	'Mapk13',	'Tmprss2',	'Kif5b',	'Txnrd1',	'Ctsh',	'Taf10',	'Grpel1',	'Cyp4f14',	'Fbp2',	'Actr3',	'Clic1',	'Twf1',	'Tmem14c',	'Hsp90aa1',	'Atp1a1',	'Junb',	'Wdr18',	'Elf3',	'Akr1b3',	'Cyp2d26',	'Ahsa1',	'Psma1',	'Ctnnb1',	'Hspa4',	'Txn2',	'Psmd3',	'Mrps10',	'Bag1',	'Mcm3',	'Cct2',	'Sdc1',	'Prmt1',	'Ddb1',	'Tsta3',	'Ccnd2',	'Eps8l3',	'Bsg' )      
AA_induced_invivo <- c('Hspa1b',	'Gm2000',	'Defa24',	'Banp',	'Lockd',	'Dpm3',	'Rps27rt',	'Snrpg',	'mt-Nd3',	'Cenpw',	'Slirp',	'Pet100',	'Tomm7',	'Tmem256',	'Cops9',	'Atp5k',	'Snhg6',	'Smim4',	'Mrpl52',	'Ndufb1-ps',	'Epb41l4aos',	'Rpl38',	'Pin4',	'Snhg8',	'Snrpf',	'Ndufa2',	'Psmg4',	'Atp5mpl',	'Romo1',	'Rps29',	'Rps21',	'Rps28',	'Anapc13',	'Atp5md',	'Lsm7',	'Tmem258',	'Uqcc2',	'Plp2',	'S100a6',	'Phgdh',	'Ndufa3',	'Eif3j1',	'Nop10',	'Sec61g',	'Aspm',	'Atox1',	'Polr2l',	'Ndufb2',	'Ccdc167',	'Fkbp5',	'Bola2',	'Ndufa5',	'Ndufc1',	'Mrps21',	'Snrpe',	'Uqcr11',	'Hist1h2ae',	'Sp3os',	'Ubl5',	'Rpl39',	'Tuba1c',	'Rpl37',	'Prelp',	'Mrpl33',	'Snrpd2',	'Naa38',	'Polr2k',	'Uqcr10',	'Cox17',	'Smim22',	'Rplp2',	'Prc1',	'Rpl37a',	'Rpl31',	'Sgo1',	'Rpl21',	'Fmc1',	'Ifitm3',	'Ppih',	'Nusap1',	'Trappc10',	'Zfas1',	'Rpl36',	'Ncapg',	'Polr2i',	'Tmsb15b2',	'Cox20',	'Ndufv3',	'Gm11808',	'Atp5e',	'Arpp19',	'Adra2a',	'Fxyd3',	'Mt1',	'Esco2',	'2410006H16Rik',	'Rpl35a',	'Ndufa1',	'Elob',	'Hist1h1b',	'Ndufs6',	'2310009B15Rik',	'Cit',	'Rpl41',	'Hspe1-rs1',	'S100a11',	'Ndufa7',	'Ndufb4',	'S100g',	'Srek1ip1',	'Sem1',	'Hmgcs2',	'Mif',	'Aurka',	'Diaph3',	'Cd44',	'Knl1',	'Snhg3',	'Top2a',	'Hells',	'Rrp8',	'Plk1',	'Ddx18',	'Cep164',	'Kif20b',	'Cdca2',	'Pmepa1',	'Jaml',	'Incenp',	'Smim26',	'Trip13',	'Mis18bp1',	'Rbmx2',	'Lsm6',	'Spc25',	'Rrp1b',	'Rps17',	'Lsm5',	'S100a13',	'Rpgrip1',	'Cdca8',	'Gas5',	'Ddx24',	'Hist1h1e',	'Cenpf',	'Hmmr',	'Atp5j2',	'Acot1',	'Kif22',	'Neil3',	'Rpl34',	'Mis12',	'Cox7a1',	'Tuba1b',	'Cbx6',	'Tubb2b',	'Ccdc34',	'Rpl36a',	'Rbmxl1',	'Ndc80',	'2310009A05Rik',	'Cep295',	'Bub1',	'Spink4',	'Hirip3',	'H2afj',	'Dlgap5',	'Atp5l',	'Smim11',	'Cox6c',	'Dnajb1',	'Kif11',	'Ubap2',	'Tpx2',	'Cebpz',	'Ctc1',	'Crip1',	'Sycn',	'Telo2',	'Rbis',	'Rassf4',	'Shcbp1',	'Cox7c',	'Gm43813',	'Rgcc',	'Notch1',	'Serf1',	'Ndufa13',	'Aurkb',	'Pbk',	'Sec61b',	'Lrrc31',	'Chordc1',	'Rps27l',	'1810037I17Rik',	'Rpp21',	'2200002D01Rik',	'E030030I06Rik',	'Fbxo5',	'Rad51ap1',	'Uqcrq',	'Polr2f',	'Rpl35',	'Ttr',	'Nup210',	'Krtcap2',	'Mcph1',	'Rps27',	'Ndufaf8',	'Lig1',	'Cd3eap',	'Ppwd1',	'1700097N02Rik',	'2810408I11Rik',	'Ckap2',	'2810004N23Rik',	'Rps2',	'Prpf4',	'Snhg18',	'Kif4',	'Sgo2a',	'Higd1a',	'Ndufs5',	'Xpo5',	'Soat1',	'Rps15',	'Anln',	'Rps26',	'Grwd1',	'Lrrc14',	'Cox7a2',	'Ncapd2',	'Tcof1',	'Rad54l',	'Ier3ip1',	'Cenpe',	'Apobec3',	'Lyrm2',	'Bub1b',	'Ccdc85c',	'Cep192',	'Prpf6',	'Cox7b',	'Rps25',	'Pinx1',	'Tex10',	'Gpatch4',	'Gnl3')          
AA_induced_stem2 <- c('S100a6','Ly6a','Malat1','Tmsb4x','Pigr','Ftl1','Gnas','Ubc')          
AA_induced_stem2 <- c('Plekhm3',	'H2-Eb1',	'Cd74',	'H2-Ab1',	'Nucb2',	'Map7d1',	'Sik1',	'Upb1',	'Arid3a',	'Cyb561a3',	'Ankrd44',	'Ctsl',	'Tbc1d8',	'Slc29a3',	'Orai3',	'Dusp1',	'Ifnar2',	'Tcf4',	'Irf7',	'Unc93b1',	'Dusp6',	'Ptpn6',	'Tubgcp5',	'S100a6',	'Ptpre',	'Rell1',	'Tom1',	'Crlf2',	'Snx29',	'Gm37529',	'Gsn',	'Snx18',	'Gmfg',	'Mob3a',	'Mfsd12',	'Unc119',	'Uvrag',	'Cybc1',	'Cmtm7',	'Cnp',	'Npc1',	'Map4k2',	'Stat2',	'H2-Q4',	'Gm43329',	'Tmem123',	'Cdkn2d',	'Nek7',	'Slc15a4',	'Gns',	'Anxa5',	'Bcl11a',	'Ppp1r15a',	'Sema4b',	'Itpr1',	'Aldh3b1',	'Igtp',	'Gng10',	'Dock8',	'Fam174a',	'Pmepa1',	'Rogdi',	'Bbc3',	'Rab33b',	'Ftx',	'Pde7a',	'Lyn',	'Dclre1c',	'Atp13a2',	'Ccl5',	'Cbx4',	'Tgfbr1',	'Gpcpd1',	'Pqlc3',	'Rpgrip1',	'Arid3b',	'Ablim1',	'Gm31718',	'Dipk1a',	'Lgmn',	'Grina',	'Rhobtb2',	'Snx30',	'Nudt18',	'Trim12a',	'Slc49a4',	'Myo9b',	'Sgk3',	'Tor3a',	'Bicd2',	'Rwdd2a',	'Dtx2',	'Ctss',	'Slc25a4',	'5031439G07Rik',	'Grn',	'Tex2',	'Vamp1',	'A430035B10Rik',	'Inpp5k',	'Ccnd3',	'Irf8',	'Stat5a',	'Rab2b',	'Tap1',	'Pacc1',	'Ifnar1',	'Ier5',	'Pml',	'Acox3',	'Arl10',	'Il4ra',	'Tap2',	'Apobec3',	'Litaf',	'Malt1',	'Il17ra',	'Bmyc',	'Arhgap17',	'Zfp319',	'Bnip2',	'Pkd1',	'Tpst2',	'H2-T22',	'Calcoco1',	'Pip5k1c',	'Gmip',	'Crtc3',	'Ptpn1',	'L3mbtl3',	'Tcirg1',	'Prcp',	'Stat1',	'Cyth1',	'Nabp1',	'Nek6',	'Bin1',	'Arhgap27',	'Wdr81',	'Ptprj')

Tuft_genes <- c('Alox5ap',	'Hck',	'Lrmp',	'Avil',	'Trpm5',	'Spib',	'Rgs13',	'Ltc4s',	'Pygl',	'Sh2d7',	'Dclk1',	'Alox5',	'Pik3r5',	'Fyb',	'Vav1',	'Matk',	'Tspan6',	'Strip2',	'Pou2f3',
                '1810046K07Rik',	'Ptpn6',	'Bmx',	'Tuba1a',	'Espn',	'Plcb2',	'Ffar3',	'Ccdc109b',	'Plcg2',	'Ly6g6f',	'Hpgds',	'Pea15a',	'Ly6g6d',	'Pik3cg',	'Inpp5d',	'Ccdc28b',	'Snrnp25',	'Kctd12',	'Siglec5',	'Skap2',	'Ccdc129',	'Nebl',	'Gprc5c',	'Rgs22',	'Gfi1b',	'Hmx3',	'Cbr3',	'Pfkfb3',	'Prss53',	'Itpr2',	'Limd2',	'Cd300lf',	'Chn2',	'Smpx',	'Ptgs1',	'A4galt',	'Rac2',	'Csk',	'Slco4a1',	'Ptpn18',	'Chat',	'Hebp1',	'Ppp1r14c',	'Dgki',	'Inpp5j',	'Tppp3',	'Gng13',	'Ildr1',	'Cwh43',	'Il17rb',	'Ncf2',	'Fut2',	'Coprs',	'Ddah1',	'Tmem116',	'Sucnr1',	'Tmem176a',	'Ccrl1',	'1110007C09Rik',	'Adcy5',	'Fnbp1',	'Plk2',	'Hmx2',	'Tmem141',	'Krt23',	'Gprc5a',	'Rgs2',	'Camk2b',	'Fes',	'Bpgm',	'Acacb',	'Il13ra1',	'Zfp428',	'Ppp1r3b',	'Ccnj',	'Bcl2l14',	'Tmem229a',	'Ethe1',	'Runx1',	'Gga2',	'Apobec1',	'Serpini1',	'St6galnac6',	'Fbxl21',	'9030624J02Rik',	'Inpp5b',	'Samd14',	'Pgm2l1',	'Pla2g4a',	'Ptprc',	'Aldh2',	'Ifi27l1',	'Pnpla3',	'Jarid2',	'Rgs19',	'Reep5',	'Tiparp',	'Gnai2',	'Fam49a',	'Cacna2d2',	'Ypel2',	'Cd24a',	'Acot7',	'Svil',	'Abhd16a',	'Fam101a',	'Trim40',	'Trak1',	'Sec14l1',	'4930539E08Rik',	'Smtn',	'Galk1',	'Tbc1d1',	'Tmem176b',	'Fcna',	'Abhd2',	'Hsbp1l1',	'Slc4a8',	'Myo1b',	'Tmem38b',	'Hk1',	'Neurl1a',
                'Dmxl2',	'Bub3',	'Ptprj',	'Trib2',	'Stard5',	'Ubtd1',	'Slc41a3',	'Plekhg5',	'Rbm38',	'Fam57a',	'Eef2k',	'Cables2',	'Fbxo25',	'Ap1s2',	'1300002K09Rik',	'Ero1lb',	'Clmn',	'Fam49b',	'Cpvl',	'Prr15',	'Lpcat4',	'Tmem74b',	'Mn1',	'Eppk1',	'Samd9l',	'Tmem245',	'Glyctk',	'Aldh3a2',	'Ppp3ca',	'Cpne3',	'Slc4a7',	'Nfatc1',	'Kit',	'Fam117b',	'Nradd',	'Tmem121',	'Cpm',	'Asah1',	'Slc9a9',	'Ubl7',	'Abca3',	'Pde6d',	'Bmp2',	'Kdm4a',	'Camkk2',	'Arhgap8',	'Agt',	'Ptpra',	'Adh1',	'Dusp14',	'Clic4',	'Gimap1',	'Cpne5',	'Ceacam2',	'Zfp710',	'Gcnt1',	'B4galt5',	'Suco',	'Pim3',	'Ogdhl',	'Oas1g',	'Dcp1b',	'Myzap',	'Cdkn1a',	'Cd37',	'Brms1',	'Lrrc42',	'Pld2',	'Tmem9',	'Cpeb4',	'Ssx2ip',	'Ddah2',	'Tmem65',	'5430417L22Rik',	'2210016L21Rik',	'Msi2',	'B4galt4',	'Rabgap1l',	'Pik3r3',	'Nt5c3',	'Palld',	'AA467197',	'Pip5k1b',	'Krt18',	'Map1a',	'Lmf1',	'Arhgef28',	'Nsfl1c',	'Txndc16',	'Pstpip2',	'Ttll11',	'Exph5',	'2700086A05Rik',	'Gadd45a',	'Plekhs1',	'Fam188a',	'Jmy',	'Atat1',	'Arhgef2',	'Lmbr1',	'Rhoc',	'Card10',	'Kcnj16',	'Arhgap4',	'Acsl4',	'Rhog',	'Fam221a',	'Dynlt1b',	'C2',	'Zbtb41',	'Socs1',	'Atp6ap1',	'Fam171a1',	'Wnk2',	'Kcnd3',	'Slc27a1',	'Atxn1',	'Rabgap1',	'Myrfl',	'Crot',	'Tm4sf4',	'Ube2j1',	'Sort1',	'Lima1',	'Mov10',	'Lca5',	'Gimap9',	'Mlip',	'1110008P14Rik',	'Ckap4',	'Tor4a',	'Rmdn1',	'Oas2',	'Dsp',	'Sox9',	'Osbpl3',	'Kif21b',	'Tbcb',	'Arap2',	'Casp3',	'Enc1',	'Il25',	'Lman2l',	'Zmiz1',	'Nav2',	'Atp2a3',	'Gimap8',	'Folr1',	'Fn1',	'Hspa4l',	'Sufu',	'Atp8a1',	'Vps53',	'Rgs14',	'Gm17660',	'Pdcl',	'Shkbp1',	'Oas1a',	'Pkp1',	'Ccdc23',	'Il4ra',	'1700112E06Rik',	'Dvl1',	'Zfhx3',	'Adam22',	'Gramd1c',	'Tmem45b',	'Unc5b',	'Mical3',	'Kctd13',	'Ak7',	'Tcta',	'Nek7',	'D730039F16Rik',	'Plekho2',	'Myo6',	'Chdh',	'Opn3',	'Tle3',	'Ttll10',	'Strada',	'Ypel3',	'Cmip',	'Cachd1',	'Pigc',	'Atp6v1d',	'Rdx',	'S100a11',	'Spa17',	'Gimap5',	'Cystm1',	'Zdhhc17',	'Lect2',	'Vdac3',	'Hspb11',	'Gm4952',	'Slc16a2',	'Abhd5',	'Rhbdf1',	'Cblb',	'Nfe2l3',	'Pla2g16',	'Sept8',	'Gpcpd1',	'Psd3',	'Anxa11',	'Slc25a12',	'Ehf',	'Akr1b10',	'Dapp1',	'Vmn2r26',	'Esyt1',	'Ppt1',	'Cd47',	'Chi3l1',	'Mical1',	'Gna14',	'Pacs2',	'Lyn',	'Rmnd5a',	'Ankrd12',	'BC022687',	'Rit1',	'Camta2',	'Mocs2',	'Usp49',	'Nrbp2',	'Ifnar2',	'Epha4',	'Arl5a',	'Rgl2',	'St18',	'BC016579',	'Tead1',	'Enpp4',	'Tmem158',	'Tnfaip3',	'Gys1',	'Hivep2',	'Cap1',	'Slc4a2',	'Map4k4',	'Desi1',	'H2-D1',	'Man2a1',	'Cyp17a1',	'Cyhr1',	'Morf4l1',	'Mllt4',	'Phf17',	'Stox2',	'Hist3h2a',	'Hdac6',	'Prox1',	'Dtnb',	'Lrch4',	'Spire2',	'Klf6',	'Rab5b',	'Anxa4',	'Rab4b',	'Iqsec1',	'Pdpk1',	'Stk40',	'Gde1',	'Mtmr11',	'Cib2',	'March2',	'Capg',	'Narf',	'Mgst3',	'Angel1',	'Bicd1',	'Ifitm1',	'Stx3',	'S100a1',	'Omd',	'0610040J01Rik',	'Arpc5',	'Homer3',	'Cdc42se1',	'Abcc3',	'Hsf2',	'Pnpla6',	'Ccdc68',	'Fryl',	'Lmtk2',	'Tas1r3',	'4931406H21Rik',	'Uspl1',	'Ajuba',	'Kalrn',	'Basp1',	'Pip5kl1',	'Slc26a2',	'Atp2b2',	'Smug1',	'Myadm',	'D330041H03Rik',	'Wdfy2',	'Trim38',	'Arf3',	'Scand1',	'Dpysl2',	'Ndufaf3',	'Sik1',	'Wdr7',	'Sfxn3',	'Kcnq4',	'Mll1',	'Hsbp1',	'Calml4',	'Atf7ip',	'Gpr137b-ps',	'Hap1',	'Kctd15',	'Prcp',	'9430023L20Rik',	'Gmip',	'Cmtm3',	'Madd',	'Krt222',	'Nsf',	'Klhl28',	'Pparg',	'Eml3',	'Phlda1',	'P2rx1',	'Pde9a',	'Otud7b',	'Tfpi2',	'Rilpl2',	'Klf3',	'Gyg',	'4930455F23Rik',	'Armcx1',	'Lzts2',	'Plek',	'Vamp8',	'Stat2',	'Znf512b',	'Ptplad1',	'1110058L19Rik',	'Tmem160',	'Tmem51',	'Cdhr5',	'Stk38',	'Atp13a2',	'Nptn',	'Sirt5',	'Gabarapl2',	'Nudt14',	'2010111I01Rik',	'Alkbh7',	'Slc18a3',	'4930427A07Rik',	'Ttll7',	'Acss2',	'Siae')  
tuft_genes <- c('Lrmp',	'Alox5ap',	'Rgs13',	'Sh2d6',	'Ltc4s',	'Avil',	'Hck',	'Dclk1',	'Snrnp25',	'Cd24a',	'Trpm5',	'Kctd12',	'Aldh2',	'Il13ra1',	'Gng13',	'Tmem176a',	'Skap2',	'Ptpn6',	'Ly6g6f',	'Fyb',	'Adh1',	'Tmem176b',	'Hpgds',	'Reep5',	'Ptpn18',	'Spib',	'Bpgm',	'Galk1',	'Matk',	'Tuba1a',	'1810046K07Rik',	'Hmx2',	'Ccdc28b',	'Ethe1',	'Limd2',	'Sh2d7',	'Ccdc109b',	'Tspan6',	'Smpx',	'Vav1',	'Ly6g6d',	'Pik3r5',	'Nebl',	'Plcg2',	'Rbm38',	'Vdac3',	'Krt18',	'Asah1',	'Cd47',	'Krt23',	'Bcl2l14',	'Lima1',	'Pygl',	'Itpr2',	'Inpp5j',	'Pea15a',	'Rac2',	'Pou2f3',	'Atp2a3',	'Bmx',	'Acot7',	'Gnai2',	'Alox5',	'Ppp3ca',	'Ptgs1',	'Calm2',	'Zfp428',	'Tmem141',	'Myo1b',	'Siglecf',	'Pla2g4a',	'Inpp5b',	'Fam221a',	'Bub3',	'Arpc5',	'Pla2g16',	'1110007C09Rik',	'Gimap1',	'Coprs',	'Lect2',	'Nrgn',	'Agt',	'Ffar3',	'Tmem45b',	'Ccdc23',	'Rgs2',	'Mlip',	'Csk',	'2210016L21Rik',	'St6galnac2',	'Ildr1',	'Gprc5c',	'Mocs2',	'Nrep',	'Pik3cg',	'Malat1',	'Sec14l1',	'Ndufaf3',	'Inpp5d',	'Pim3',	'Tmem9',	'Gga2',	'Nt5c3')
tuft_genes <- c("Trpm5" ,"Gfi1b" ,"Il25" ,"Alox5ap" ,"Lrmp" ,"Rgs13" ,"Ltc4s" ,"Adh1", 'Dclk1', 'Dusp4','Pten','Ltc4s','Stub1','Vil1')
mhc1_genes <- c("H2-K1", "H2-K2", "H2-D1", "H2-D2", "H2-L", "H2-Q1", "H2-Q2",'H2-Q3','H2-Q4','H2-Q6', 'H2-M3', "H2-T3", 'H2-DMa',
                "H2-Q10", "H2-M2", "H2-M3", "H2-Oa", "H2-Ob", 'H2-T25', "H2-T23",'H2-T27', 'H2-DMb1', 'H2-DMb2', 'H2-Ob',
                "H2-T24", 'H2-T10', 'H2-Q10')
mhc2_genes <- c('H2-Ea', "H2-Eb1", "H2-Eb2", "H2-Ab1", "H2-Ab2","H2-Aa", 'H2-D1','Ciita','Cd74')
Stem_genes <- c('Lgr5',	'Gkn3',	'Ascl2',	'Olfm4',	'Rgmb',	'Igfbp4',	'2210407C18Rik',	'Jun',	'Pdgfa',	'Soat1',	'Tnfrsf19',	'Cyp2e1',	'Fstl1',	'H2-Eb1',	'Ifitm3',	'Prelp',	'Scn2b',	'A930009A15Rik',	'H2-Ab1',	'Slc1a2',	'Cd74',	'Sp5',	'Noxa1',	'Rgcc',	'Sorbs2',	'Sectm1b',	'H2-Aa',	'Cdo1',	'Slc14a1',	'Clca2',	'Tifa',	'Pls3',	'Hmgcs2',	'Arid5b',	'Agr3',	'Slc12a2',	'Rassf5',	'Rnf43',	'Nrn1',	'Lamb3',	'Cd44',	'Axin2',	'Slc27a2',	'Afap1l1',	'Ccdc3',	'Lrig1',	'Noxo1',	'Cdk6',	'Amica1',	'Tgif1',	'Tns3',	'Nr2e3',	'Efna4',	'Rnf32',	'Prss23',	'2010009K17Rik',	'Smoc2',	'Mecom',	'Esrrg',	'Aqp1',	'Znrf3',	'Grb7',	'Phgdh',	'2410004N09Rik',	'Clca4',	'Aqp4',	'Lcp1',	'E030011O05Rik',	'Snhg1',	'BC064078',	'Car12',	'Zbtb38',	'Cdca7',	'Fam13a',	'Shisa2',	'Dtx4',	'Slc19a2',	'Fam115c',	'Mir703',	'Cd14',	'Mettl20',	'Myo9a',	'App',	'Clic6',	'Wee1',	'2410006H16Rik',	'Lancl1',	'1500012F01Rik',	'Casp12',	'Sh3rf1',	'Lrp4',	'Arhgef26',	'Etv6',	'1700024F13Rik',	'Cttnbp2',	'Slc16a13',	'Htr4',	'Pdxk',	'Immp2l',	'Rps15a-ps6',	'Rps15a-ps4',	'Nap1l1',	'Sdc4',	'Epn3',	'Sipa1l1',	'Wfdc15b',	'Zfp341',	'Ngef',	'Nrg4',	'Csad',	'Rpl34-ps1',	'Rin2',	'Cd81',	'Irf2bp2',	'Sesn3',	'Phlpp1',	'Yap1',	'Mfge8',	'Zfp825',	'Itga1',	'Pcdh8',	'Vdr',	'Kcnq1',	'Slc28a2',	'Zfp36l1',	'Urod',	'Rgs12',	'Nfib',	'Sdsl',	'Nfia')
Stem2_genes <- c('S100a6', 'Clu', 'Ly6a', 'Anxa3','Anxa2','Msi1')


gr <- transcripts(txdb)
gr <- annotatePeakInBatch(gr, AnnotationData=TSS.mouse.GRCm38)

gr1 <- addGeneIDs(annotatedPeak=gr,
                  orgAnn="org.Mm.eg.db",
                  IDs2Add="symbol")
gr2 <- gr1[gr1$symbol %in% AA_induced_stem2]

df <- data.frame(seqnames=seqnames(gr2),
                 starts=start(gr2)-1,
                 ends=end(gr2),
                 names=c(rep(".", length(gr2))),
                 scores=c(rep(".", length(gr2))),
                 strands=strand(gr2))

write.table(df, file="~/analysis/stem2_markers_induced.bed", quote=F, sep="\t", row.names=F, col.names=F)



