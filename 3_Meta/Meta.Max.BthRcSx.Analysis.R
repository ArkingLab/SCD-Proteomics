####Description####
## Analyze Meta Max BthRcSx dataset for SCA Proteomics
## Min = Primary model
## RskFct = Full model

####Load packages####
library(magrittr)
#library(metafor)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(cowplot)
library(colorspace)
library(ggforestplot)
library(VennDiagram)
library(tidyr)

####Load data####
aric.chs.meta.results <- readRDS("aric.prot.real.archive/vthuyduo/scd.assoc/5_meta/with.10032022.CHS/BthRcSx.Meta.Results.RDS")
load("aric.prot.real.archive/vthuyduo/scd.assoc/5_meta/with.10032022.CHS/PtnRefList.rds")

####Volcano plots####
make.simp.col.volc2 <- function(data, annotate, pcolumn, betacolumn, seqidcolumn, targetcolumn, entrezsymbcolumn, wordlength, sigcutoff, overlaps, xstart, xend){
  #Added options to choose SeqId, target, entrez columns
  #Seq Ids must be in the XXXX-XX format

  colors <- brewer.pal(9, "Set1")
  colors <- colors[-6]
  colors2 <- brewer.pal(12, "Set3")

  bonf <- 0.05/4955
  bonf_log <- -log10(bonf)

  ## Label assignment
  data$pcol <- NA
  data$pcol <- ifelse(data[,pcolumn] < bonf, "Significant", data$pcol)
  #data$pcol <- ifelse(data[,5] < bonf & data[,11] < bonf, "Significant in Risk Factor model", data$pcol)
  data$pcol <- ifelse(data[,pcolumn] >= bonf, "Non-significant", data$pcol)

  data$label <- data[, targetcolumn]
  data$label <- ifelse(data[, targetcolumn] == "Deprecated", data[, entrezsymbcolumn], data$label)
  data$label <- ifelse(data$pcol != "Non-significant", data[, targetcolumn], NA)
  data$label <- ifelse(nchar(data$label)>wordlength | as.numeric(-log10(data[,pcolumn]))<=sigcutoff, NA, data$label)
  data$label <- ifelse(data[, seqidcolumn] == "7655-11", "N-terminal pro-BNP", data$label)

  ## Count
  data <- data %>% add_count(pcol) %>% mutate(znn = paste0(pcol, ' (', n, ')'))

  data$HR <- exp(data[,betacolumn])
  data$LogP <- -log10(data[,pcolumn])

  ## Plot
  plot <- ggplot(data, aes(x=HR, y=LogP, col=pcol, label=label)) +
    geom_point(size = 2, alpha = 0.8) + geom_segment(x=(xstart-0.1), xend=(xend+0.1), y=bonf_log, yend = bonf_log,linetype = "dashed", color = "black", size = 0.2) +
    geom_text_repel(size = 2.5, max.overlaps = overlaps, col = "black", min.segment.length = 0, hjust = 0.5, nudge_y=0.3, nudge_x = 0.01, segment.size = 0.1) + ylim(0,40) + ylab("-log10(P)") + xlab("Hazard Ratio") +
    annotate(geom = "label", label = annotate, x = -Inf, y = Inf, hjust = 0, vjust = 1) +
    scale_x_continuous(breaks = seq(xstart,xend,0.1)) + coord_cartesian(xlim = c(xstart,xend)) + scale_color_manual(values = c("royalblue", colors2[9]), breaks = c("Significant","Non-significant")) +
    theme_classic(base_size = 12) + theme(legend.title = element_blank(), legend.position = "None")
  return(plot)
}

pdf(file = "aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/AllCohorts.Max.BthRcSx.AllMods.AllPtns_Volcano.pdf", width = 15, height = 15)
plot_grid(
  subset(aric.chs.meta.results, model.cohort == "M2.meta") %>%
  make.simp.col.volc2(., annotate = "Meta Max BthRcSx Min", pcolumn = "P", betacolumn = "Beta", seqidcolumn = "SeqId.use", targetcolumn = "ARIC.target", 
                      entrezsymbcolumn = "ARIC.entrezgenesymbol", wordlength = 15, sigcutoff = 7, overlaps = 10, xstart = 0.4, xend = 2.6),
  subset(aric.chs.meta.results, model.cohort == "M1.meta") %>%
    make.simp.col.volc2(., annotate = "Meta Max BthRcSx Dz", pcolumn = "P", betacolumn = "Beta", seqidcolumn = "SeqId.use", targetcolumn = "ARIC.target",
                        entrezsymbcolumn = "ARIC.entrezgenesymbol", wordlength = 15, sigcutoff = 7, overlaps = 10, xstart = 0.4, xend = 2.6),
  subset(aric.chs.meta.results, model.cohort == "M3.meta") %>%
    make.simp.col.volc2(., annotate = "Meta Max BthRcSx RskFct", pcolumn = "P", betacolumn = "Beta", seqidcolumn = "SeqId.use", targetcolumn = "ARIC.target",
                        entrezsymbcolumn = "ARIC.entrezgenesymbol", wordlength = 15, sigcutoff = 5, overlaps = 20, xstart = 0.4, xend = 2.6),

  subset(aric.chs.meta.results, model.cohort == "M2.aric") %>%
    make.simp.col.volc2(., annotate = "ARIC Max BthRcSx Min", pcolumn = "P", betacolumn = "Beta", seqidcolumn = "SeqId.use", targetcolumn = "ARIC.target",
                        entrezsymbcolumn = "ARIC.entrezgenesymbol", wordlength = 15, sigcutoff = 4, overlaps = 10, xstart = 0.4, xend = 2.6),
  subset(aric.chs.meta.results, model.cohort == "M1.aric") %>%
    make.simp.col.volc2(., annotate = "ARIC Max BthRcSx Dz", pcolumn = "P", betacolumn = "Beta", seqidcolumn = "SeqId.use", targetcolumn = "ARIC.target",
                        entrezsymbcolumn = "ARIC.entrezgenesymbol", wordlength = 15, sigcutoff = 4, overlaps = 10, xstart = 0.4, xend = 2.6),
  subset(aric.chs.meta.results, model.cohort == "M3.aric") %>%
    make.simp.col.volc2(., annotate = "ARIC Max BthRcSx RskFct", pcolumn = "P", betacolumn = "Beta", seqidcolumn = "SeqId.use", targetcolumn = "ARIC.target",
                        entrezsymbcolumn = "ARIC.entrezgenesymbol", wordlength = 15, sigcutoff = 4, overlaps = 20, xstart = 0.4, xend = 2.6),

  subset(aric.chs.meta.results, model.cohort == "M2.chs") %>%
    make.simp.col.volc2(., annotate = "CHS Max BthRcSx Min", pcolumn = "P", betacolumn = "Beta", seqidcolumn = "SeqId.use", targetcolumn = "ARIC.target",
                        entrezsymbcolumn = "ARIC.entrezgenesymbol", wordlength = 15, sigcutoff = 4, overlaps = 10, xstart = 0.4, xend = 2.6),
  subset(aric.chs.meta.results, model.cohort == "M1.chs") %>%
    make.simp.col.volc2(., annotate = "CHS Max BthRcSx Dz", pcolumn = "P", betacolumn = "Beta", seqidcolumn = "SeqId.use", targetcolumn = "ARIC.target",
                        entrezsymbcolumn = "ARIC.entrezgenesymbol", wordlength = 15, sigcutoff = 4, overlaps = 10, xstart = 0.4, xend = 2.6),
  subset(aric.chs.meta.results, model.cohort == "M3.chs") %>%
    make.simp.col.volc2(., annotate = "CHS Max BthRcSx RskFct", pcolumn = "P", betacolumn = "Beta", seqidcolumn = "SeqId.use", targetcolumn = "ARIC.target",
                        entrezsymbcolumn = "ARIC.entrezgenesymbol", wordlength = 15, sigcutoff = 4, overlaps = 20, xstart = 0.4, xend = 2.6),
  nrow = 3, ncol = 3)
dev.off()

####Defining protein subsets####

#Proteins that meet Bonf cutoff in the min model
meta.norm.min.bonf.seqid <- subset(aric.chs.meta.results, model == "M2" & PASSED == "YES" & cohort == "meta")$SeqId.use

#Min proteins that meet Bonf cutoff for rskfct model as well
meta.norm.min.rskfct.bonf.seqid <- subset(aric.chs.meta.results, SeqId.use %in% meta.norm.min.bonf.seqid) %>% 
  subset(., model == "M3" & PASSED == "YES" & cohort == "meta") %>% .[, "SeqId.use"]

#Proteins that meet Bonf cutoff across all models
meta.norm.all.mods.bonf.seqid <- c(subset(aric.chs.meta.results, model == "M2" & PASSED == "YES" & cohort == "meta")$SeqId.use,
                                   subset(aric.chs.meta.results, model == "M1" & PASSED == "YES" & cohort == "meta")$SeqId.use,
                                   subset(aric.chs.meta.results, model == "M3" & PASSED == "YES" & cohort == "meta")$SeqId.use) %>% unique()

####Heatmap####

## Add asterisk for proteins that pass Bonf
aric.chs.meta.results$PASSED_PLOT <- ifelse(aric.chs.meta.results$PASSED == "YES", "*", "")

## ORDER MADE - add column to order proteins by Meta min beta and P value [column name: order_meta.M2.p.beta]
aric.chs.meta.results <- as.data.frame(rbind(subset(aric.chs.meta.results, Beta > 0 & model == "M2" & cohort == "meta") %>% .[order(.$P),] %>% 
                                               mutate(order_meta.M2.p.beta = 1:nrow(.)) %>% .[, c("SeqId.use", "order_meta.M2.p.beta")],
                                             subset(aric.chs.meta.results, Beta < 0 & model == "M2" & cohort == "meta") %>% .[order(.$P),] %>% 
                                               mutate(order_meta.M2.p.beta = 2662:4955) %>% .[, c("SeqId.use", "order_meta.M2.p.beta")])) %>% 
  merge(aric.chs.meta.results, ., by = "SeqId.use", sort = F, all.x = TRUE)

## Add * to GOLM1, SVEP1
aric.chs.meta.results$target_use <- aric.chs.meta.results$ARIC.target
aric.chs.meta.results$target_use[grep("11178-21", aric.chs.meta.results$SeqId.use)] <- paste0(aric.chs.meta.results$target_use[grep("11178-21", aric.chs.meta.results$SeqId.use)], "*")
aric.chs.meta.results$target_use[grep("8983-7", aric.chs.meta.results$SeqId.use)] <- paste0(aric.chs.meta.results$target_use[grep("8983-7", aric.chs.meta.results$SeqId.use)], "*")

## Fix deprecated naming
aric.chs.meta.results$target_use <- ifelse(aric.chs.meta.results$target_use == "Deprecated", aric.chs.meta.results$ARIC.entrezgenesymbol, aric.chs.meta.results$target_use)

## Edit model column
aric.chs.meta.results$model_use <- aric.chs.meta.results$model
aric.chs.meta.results$model_use <- ifelse(aric.chs.meta.results$model == "M2", "Minimum", aric.chs.meta.results$model_use)
aric.chs.meta.results$model_use <- ifelse(aric.chs.meta.results$model == "M1", "Disease", aric.chs.meta.results$model_use)
aric.chs.meta.results$model_use <- ifelse(aric.chs.meta.results$model == "M3", "Risk Factor", aric.chs.meta.results$model_use)
aric.chs.meta.results$model_use <- factor(aric.chs.meta.results$model_use, ordered = T, levels = c("Minimum", "Disease", "Risk Factor"))

## Min Bonf proteins only
ggplot(subset(aric.chs.meta.results, SeqId.use %in% meta.norm.min.bonf.seqid & cohort == "meta"), aes(y = model_use, x = reorder(target_use, order_meta.M2.p.beta), fill=HR)) +
  geom_tile() + geom_text(aes(label = PASSED_PLOT), color = "white", size = 2.5, vjust="center", hjust = "middle") +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 1, l3 = 0, p3 = .4, p4 = .6, p2 = 1, rev = TRUE) + coord_fixed(ratio = 4) +
  scale_y_discrete(limits = rev) +  theme_classic() + labs(fill = "Hazard\nRatio") +
  ylab("Model") + xlab("Protein") + theme(axis.text.y = element_text(size = 6.8), axis.text.x = element_text(size = 6.8, angle = 60, hjust=1), plot.caption = element_text(hjust = 1, size = 7)) +
  labs(caption = "Meta Max BthRcSx MinBonfPtns in order of Min Beta P")
ggsave("aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/Meta.Max.BthRcSx.AllMods.MinBonfPtns_HeatMap.Horiz.png", height = 5, width = 14, units = "in")

## Bonf proteins across all models
ggplot(subset(aric.chs.meta.results, SeqId.use %in% meta.norm.all.mods.bonf.seqid & cohort == "meta"), aes(y = model_use, x = reorder(target_use, order_meta.M2.p.beta), fill=HR)) +
  geom_tile() + geom_text(aes(label = PASSED_PLOT), color = "white", size = 2.5, vjust="center", hjust = "middle") +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 1, l3 = 0, p3 = .4, p4 = .6, p2 = 1, rev = TRUE) + coord_fixed(ratio = 4) +
  scale_y_discrete(limits = rev) +  theme_classic() + labs(fill = "Hazard\nRatio") +
  ylab("Model") + xlab("Protein") + theme(axis.text.y = element_text(size = 6.8), axis.text.x = element_text(size = 6.8, angle = 60, hjust=1), plot.caption = element_text(hjust = 1, size = 7)) +
  labs(caption = "Meta Max BthRcSx AllModsBonfPtns in order of Min Beta P")
ggsave("aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/Meta.Max.BthRcSx.AllMods.AllModsBonfPtns_HeatMap.Horiz.png", height = 5, width = 14, units = "in")

####Forest plots####
quince.forest <- function(data, column){
  forestplot(df = data, name = target_use, estimate = Beta, se = SE,
             colour = column, pvalue = P, xlab = "Hazard Ratio", logodds = TRUE, psignif = (0.05/4955)) + theme_bw() +
    theme(text = element_text(size = 10), legend.direction = "horizontal", legend.margin=margin(c(-10,1,1,1)),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank())
}

## ORDER MADE - add column to order proteins by Meta Min P value [column name: order_meta.M2.p]
aric.chs.meta.results <- subset(aric.chs.meta.results, model == "M2" & cohort == "meta") %>% 
  .[order(.$P),] %>% mutate(order_meta.M2.p = 1:nrow(.)) %>% .[, c("SeqId.use", "order_meta.M2.p")] %>% 
  merge(aric.chs.meta.results, ., by = "SeqId.use", all.x = TRUE, sort = F)

aric.chs.meta.results$cohort_use <- aric.chs.meta.results$cohort
aric.chs.meta.results$cohort_use <- toupper(aric.chs.meta.results$cohort_use)
aric.chs.meta.results$cohort_use <- factor(aric.chs.meta.results$cohort_use, ordered = T, levels = c("META", "CHS", "ARIC"))
aric.chs.meta.results$model.cohort_use <- paste0(aric.chs.meta.results$model_use, " ", aric.chs.meta.results$cohort_use)
aric.chs.meta.results$model.cohort_use <- factor(aric.chs.meta.results$model.cohort_use, ordered = T, 
                                                 levels = c("Minimum ARIC", "Minimum CHS", "Minimum META", "Disease ARIC", "Disease CHS", "Disease META",  "Risk Factor ARIC", "Risk Factor CHS", "Risk Factor META"))

## ORDER MADE - add column to order Meta Min Bonf proteins by P value (only proteins that pass Bonf for Min model will have a number) [column name: order_meta.M2.p.bonf]
aric.chs.meta.results <- subset(aric.chs.meta.results, model != "M3") %>% subset(., model != "M1") %>% subset(., cohort == "meta") %>% subset(., SeqId.use %in% meta.norm.min.bonf.seqid) %>% .[order(.$P),] %>%
  mutate(order_meta.M2.p.bonf = 1:nrow(.)) %>% .[, c("SeqId.use", "order_meta.M2.p.bonf")] %>% merge(aric.chs.meta.results, ., by = "SeqId.use", sort = F, all.x = TRUE)

## ORDER MADE - add column to randomly plot Meta Min Bonf proteins [column name: order_rand]
set.seed(10)
aric.chs.meta.results <- subset(aric.chs.meta.results, model != "M3") %>% subset(., model != "M1") %>% subset(., cohort == "meta") %>% 
  mutate(order_rand = sample(1:nrow(.))) %>% .[, c("SeqId.use", "order_rand")] %>% 
  merge(aric.chs.meta.results, ., by = "SeqId.use", sort = F, all.x = TRUE)

## Plot all models in one .pdf  
min.forest.df <- subset(aric.chs.meta.results, model == "M2" & SeqId.use %in% meta.norm.min.bonf.seqid) %>% .[order(.$order_meta.M2.p.bonf),]
dz.forest.df <- subset(aric.chs.meta.results, model == "M1" & SeqId.use %in% meta.norm.min.bonf.seqid) %>% .[order(.$order_meta.M2.p.bonf),]
rskfct.forest.df <- subset(aric.chs.meta.results, model == "M3" & SeqId.use %in% meta.norm.min.bonf.seqid) %>% .[order(.$order_meta.M2.p.bonf),]
meta.forest.df <- subset(aric.chs.meta.results, cohort == "meta" & SeqId.use %in% meta.norm.min.bonf.seqid) %>% .[order(.$order_meta.M2.p.bonf),]
  
pdf(file = "aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/AllCohorts.Max.BthRcSx.SepByMods.MinBonfPtns_Forest.pdf", height = 13, width = 16)
plot_grid(
  subset(min.forest.df, order_meta.M2.p.bonf <=15) %>% quince.forest(., .$cohort_use) + theme(legend.position = "none"),
  subset(min.forest.df, order_meta.M2.p.bonf >15 & order_meta.M2.p.bonf<=30) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(min.forest.df, order_meta.M2.p.bonf >30 & order_meta.M2.p.bonf<=45) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(min.forest.df, order_meta.M2.p.bonf >45 & order_meta.M2.p.bonf<=60) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(min.forest.df, order_meta.M2.p.bonf >60 & order_meta.M2.p.bonf<=75) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(min.forest.df, order_meta.M2.p.bonf >75 & order_meta.M2.p.bonf<=90) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(min.forest.df, order_meta.M2.p.bonf >91 & order_meta.M2.p.bonf<=105) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(min.forest.df, order_meta.M2.p.bonf >106 & order_meta.M2.p.bonf<=120) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(min.forest.df, order_meta.M2.p.bonf >120 & order_meta.M2.p.bonf<=129) %>% quince.forest(., .$cohort_use)  +
    theme(legend.position = "bottom", plot.caption = element_text(hjust = 1, size = 9)) + labs(caption = "Max BthRcSx MinBonfPtns Min HR in order of Min P"), nrow = 3, ncol = 3)

plot_grid(
  subset(dz.forest.df, order_meta.M2.p.bonf <=15) %>% quince.forest(., .$cohort_use) + theme(legend.position = "none"),
  subset(dz.forest.df, order_meta.M2.p.bonf >15 & order_meta.M2.p.bonf<=30) %>% quince.forest(., .$cohort_use) + theme(legend.position = "none"),
  subset(dz.forest.df, order_meta.M2.p.bonf >30 & order_meta.M2.p.bonf<=45) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(dz.forest.df, order_meta.M2.p.bonf >45 & order_meta.M2.p.bonf<=60) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(dz.forest.df, order_meta.M2.p.bonf >60 & order_meta.M2.p.bonf<=75) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(dz.forest.df, order_meta.M2.p.bonf >75 & order_meta.M2.p.bonf<=90) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(dz.forest.df, order_meta.M2.p.bonf >91 & order_meta.M2.p.bonf<=105) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(dz.forest.df, order_meta.M2.p.bonf >106 & order_meta.M2.p.bonf<=120) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(dz.forest.df, order_meta.M2.p.bonf >120 & order_meta.M2.p.bonf<=129) %>% quince.forest(., .$cohort_use) +
    theme(legend.position = "bottom", plot.caption = element_text(hjust = 1, size = 9)) + labs(caption = "Max BthRcSx MinBonfPtns Dz HR in order of Min P"), nrow = 3, ncol = 3)

plot_grid(
  subset(rskfct.forest.df, order_meta.M2.p.bonf <=15) %>% quince.forest(., .$cohort_use) + theme(legend.position = "none"),
  subset(rskfct.forest.df, order_meta.M2.p.bonf >15 & order_meta.M2.p.bonf<=30) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(rskfct.forest.df, order_meta.M2.p.bonf >30 & order_meta.M2.p.bonf<=45) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(rskfct.forest.df, order_meta.M2.p.bonf >45 & order_meta.M2.p.bonf<=60) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(rskfct.forest.df, order_meta.M2.p.bonf >60 & order_meta.M2.p.bonf<=75) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(rskfct.forest.df, order_meta.M2.p.bonf >75 & order_meta.M2.p.bonf<=90) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(rskfct.forest.df, order_meta.M2.p.bonf >91 & order_meta.M2.p.bonf<=105) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(rskfct.forest.df, order_meta.M2.p.bonf >106 & order_meta.M2.p.bonf<=120) %>% quince.forest(., .$cohort_use)  + theme(legend.position = "none"),
  subset(rskfct.forest.df, order_meta.M2.p.bonf >120 & order_meta.M2.p.bonf<=129) %>% quince.forest(., .$cohort_use) +
    theme(legend.position = "bottom", plot.caption = element_text(hjust = 1, size = 9)) + labs(caption = "Max BthRcSx MinBonfPtns RskFct HR in order of Min P"), nrow = 3, ncol = 3)
dev.off()

pdf(file = "aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/Meta.Max.BthRcSx.SepByMods.MinBonfPtns_Forest.pdf", height = 13, width = 16)
plot_grid(
  subset(meta.forest.df, order_meta.M2.p.bonf <=15) %>% quince.forest(., .$model_use) + theme(legend.position = "none"),
  subset(meta.forest.df, order_meta.M2.p.bonf >15 & order_meta.M2.p.bonf<=30) %>% quince.forest(., .$model_use)  + theme(legend.position = "none"),
  subset(meta.forest.df, order_meta.M2.p.bonf >30 & order_meta.M2.p.bonf<=45) %>% quince.forest(., .$model_use)  + theme(legend.position = "none"),
  subset(meta.forest.df, order_meta.M2.p.bonf >45 & order_meta.M2.p.bonf<=60) %>% quince.forest(., .$model_use)  + theme(legend.position = "none"),
  subset(meta.forest.df, order_meta.M2.p.bonf >60 & order_meta.M2.p.bonf<=75) %>% quince.forest(., .$model_use)  + theme(legend.position = "none"),
  subset(meta.forest.df, order_meta.M2.p.bonf >75 & order_meta.M2.p.bonf<=90) %>% quince.forest(., .$model_use)  + theme(legend.position = "none"),
  subset(meta.forest.df, order_meta.M2.p.bonf >91 & order_meta.M2.p.bonf<=105) %>% quince.forest(., .$model_use)  + theme(legend.position = "none"),
  subset(meta.forest.df, order_meta.M2.p.bonf >106 & order_meta.M2.p.bonf<=120) %>% quince.forest(., .$model_use)  + theme(legend.position = "none"),
  subset(meta.forest.df, order_meta.M2.p.bonf >120 & order_meta.M2.p.bonf<=129) %>% quince.forest(., .$model_use) +
    theme(legend.position = "bottom", plot.caption = element_text(hjust = 1, size = 9)) + labs(caption = "Meta Max BthRcSx MinBonfPtns HR in order of Min P"), nrow = 3, ncol = 3)
dev.off()

####Heatmap - by cohort & model####
ggplot(subset(aric.chs.meta.results, SeqId.use %in% meta.norm.min.bonf.seqid), aes(y = model.cohort_use, x = reorder(target_use, order_meta.M2.p.beta), fill=HR)) +
  geom_tile() + geom_text(aes(label = PASSED_PLOT), color = "white", size = 2.5, vjust="center", hjust = "middle") +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 1, l3 = 0, p3 = .4, p4 = .6, p2 = 1, rev = TRUE) + coord_fixed(ratio = 4) +
  scale_y_discrete(limits = rev) +  theme_classic() + labs(fill = "Hazard\nRatio") +
  ylab("Model") + xlab("Protein") + theme(axis.text.y = element_text(size = 6.8), axis.text.x = element_text(size = 6.8, angle = 60, hjust=1), plot.caption = element_text(hjust = 1, size = 9)) +
  labs(caption = "Max BthRcSx MinBonfPtns in order of Min Beta P")
ggsave("aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/AllCohorts.Max.BthRcSx.AllMods.MinBonfPtns_HeatMap.Horiz.png", height = 5, width = 14, units = "in")

ggplot(subset(aric.chs.meta.results, SeqId.use %in% meta.norm.all.mods.bonf.seqid), aes(y = model.cohort_use, x = reorder(target_use, order_meta.M2.p.beta), fill=HR)) +
  geom_tile() + geom_text(aes(label = PASSED_PLOT), color = "white", size = 2.5, vjust="center", hjust = "middle") +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 1, l3 = 0, p3 = .4, p4 = .6, p2 = 1, rev = TRUE) + coord_fixed(ratio = 4) +
  scale_y_discrete(limits = rev) +  theme_classic() + labs(fill = "Hazard\nRatio") +
  ylab("Model") + xlab("Protein") + theme(axis.text.y = element_text(size = 6.8), axis.text.x = element_text(size = 6.8, angle = 60, hjust=1), plot.caption = element_text(hjust = 1, size = 9)) +
  labs(caption = "Max BthRcSx AllModsBonfPtns in order of Min Beta P")
ggsave("aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/AllCohorts.Max.BthRcSx.AllMods.AllBonfPtns_HeatMap.Horiz.png", height = 5, width = 14, units = "in")

####Dropouts & Addins from/to ARIC min model####
aric.norm.min.bonf.seqid <- subset(aric.chs.meta.results, model == "M2" & PASSED == "YES" & cohort == "aric")$SeqId.use
table(aric.norm.min.bonf.seqid %in% meta.norm.min.bonf.seqid)
#114/131 are in common between ARIC and META min
#17 are not

## Dropouts between Min ARIC and Min Meta, random plot order
meta.min.dropouts <- aric.norm.min.bonf.seqid[which(!(aric.norm.min.bonf.seqid %in% meta.norm.min.bonf.seqid))]
dropouts.df <- subset(aric.chs.meta.results, SeqId.use %in% meta.min.dropouts)

## Addins between Min ARIC and Min Meta, 
meta.min.addins <- meta.norm.min.bonf.seqid[which(!(meta.norm.min.bonf.seqid %in% aric.norm.min.bonf.seqid))]
addins.df <- subset(aric.chs.meta.results, SeqId.use %in% meta.min.addins)

####Venn diagram for meta models####
#meta.norm.min.bonf.seqid
meta.norm.dz.bonf.seqid <- subset(aric.chs.meta.results, model == "M1" & PASSED == "YES" & cohort == "meta")$SeqId.use
meta.norm.rskfct.bonf.seqid <- subset(aric.chs.meta.results, model == "M3" & PASSED == "YES" & cohort == "meta")$SeqId.use

venn.diagram(
  x = list(meta.norm.min.bonf.seqid, meta.norm.dz.bonf.seqid, meta.norm.rskfct.bonf.seqid),
  category.names = c("Min", "Dz", "RskFct"),
  filename = 'aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/Meta.Max.BthRcSx.AllMods.BonfPtns_VennDiag.png',
  output=TRUE,
  disable.logging = TRUE,
  main = "Meta Max BthRc",
  main.cex = 0.4,

  # Output features
  imagetype="png" ,
  height = 480 ,
  width = 480 ,
  resolution = 300,
  compression = "lzw",

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = brewer.pal(3, "Pastel2"),

  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.065, 0.025, 0.045),
  cat.fontfamily = "sans",
  rotation = 1
)

####Proteins that pass meta min and meta rskfct models####
meta.norm.min.rskfct.bonf.seqid <- meta.norm.min.bonf.seqid[which(meta.norm.min.bonf.seqid %in% meta.norm.rskfct.bonf.seqid)]

## Forest plot
forestplot(
  df = subset(aric.chs.meta.results, model != "M1" & SeqId.use %in% meta.norm.min.rskfct.bonf.seqid) %>% .[order(.$order_meta.M2.p.beta),],
  name = target_use,
  estimate = Beta,
  se = SE,
  pvalue = P,
  psignif = (0.05/4955),
  colour = cohort_use,
  xlab = "Hazard Ratio",
  logodds = TRUE
) +
  ggforce::facet_col(
    facets = ~model_use,
    scales = "free_y",
    space = "free"
  ) + scale_x_continuous(breaks = seq(0.5,2.0,0.1)) + theme_classic() +
  theme(text = element_text(size = 10), legend.position = "bottom", legend.direction = "horizontal", legend.margin=margin(c(-10,1,1,1)), legend.title = element_blank(),
        plot.caption = element_text(hjust = 1, size = 9)) + labs(caption = "Meta Max BthRcSx Min & RskFct Bonf Ptns in order of Min Beta P")
ggsave("aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/Meta.Max.BthRcSx.MinRskFctBonfPtns_Forest.png", width = 7, height = 14, units = "in")

## Table
rotunda <- function(num){
  numa <- signif(num, digits = 3) %>% sprintf("%.2f",.)
  return(numa)
}
rndP <- function(num){
  numa <- signif(num, digits = 2)
  return(numa)
}

write.table(
  cbind(
    subset(aric.chs.meta.results, model == "M2" & SeqId.use %in% meta.norm.min.rskfct.bonf.seqid & cohort == "meta") %>% .[order(.$order_meta.M2.p.bonf),] %>% 
      .[, c("target_use", "HR", "HR.CI.upper", "HR.CI.lower", "P", "model.cohort")] %>% 
      mutate(HR.CI = paste0(rotunda(.$HR), " [", rotunda(.$HR.CI.lower), "-", rotunda(.$HR.CI.upper), "]"), P = rndP(.$P)) %>% .[, c("target_use", "HR.CI", "P")],
    subset(aric.chs.meta.results, model == "M3" & SeqId.use %in% meta.norm.min.rskfct.bonf.seqid & cohort == "meta") %>% .[order(.$order_meta.M2.p.bonf),] %>% 
      .[, c("target_use", "HR", "HR.CI.upper", "HR.CI.lower", "P", "model.cohort")] %>%
      mutate(HR.CI = paste0(rotunda(.$HR), " [", rotunda(.$HR.CI.lower), "-", rotunda(.$HR.CI.upper), "]"), P = rndP(.$P)) %>% .[, c("HR.CI", "P")]) %>% 
    as.data.frame() %>% setNames(c("target", "Meta.HR.CI.Min", "Meta.P.Min", "Meta.HR.CI.RskFct", "Meta.P.RskFct")),
  file = "aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/Meta.Max.BthRcSx.MinRskFctBonfPtns_Table.txt",
  row.names = F, quote = F, sep = "\t")

####Dropouts & Addins from/to ARIC Min & RskFct Hits####
#ARIC Min & RskFct hits - 32 hits found to remain after adjusting for HD and risk factors

## Dropouts
aric.norm.min.rskfct.bonf.seqid <- subset(aric.chs.meta.results, cohort == "aric" & SeqId.use %in% aric.norm.min.bonf.seqid & model == "M3" & P<(0.05/4955)) %>% .$SeqId.use
aric.hits.dropouts.df <- subset(aric.chs.meta.results, SeqId.use %in% aric.norm.min.rskfct.bonf.seqid[which(!(aric.norm.min.rskfct.bonf.seqid %in% meta.norm.min.rskfct.bonf.seqid))])

## Addins
aric.hits.addins.df <- subset(aric.chs.meta.results, SeqId.use %in% meta.norm.min.rskfct.bonf.seqid[which(!(meta.norm.min.rskfct.bonf.seqid %in% aric.norm.min.rskfct.bonf.seqid))])

####Absolute Beta Difference####
calc.min.abs.beta.diff <- function(df){
  ## Calculate difference
  df$M2.M1.BetaDiff <- abs(df[,"M2.Beta"] - df[,"M1.Beta"])
  df$M2.M3.BetaDiff <- abs(df[,"M2.Beta"] - df[,"M3.Beta"])
  df$M1.M3.BetaDiff <- abs(df[,"M1.Beta"] - df[,"M3.Beta"])
  
  ## Make columns to order by differences
  df <- df[order(df$M2.M1.BetaDiff),] %>% mutate(M2.M1.BetaDiff.Order = nrow(df):1)
  df <- df[order(df$M2.M3.BetaDiff),] %>% mutate(M2.M3.BetaDiff.Order = nrow(df):1)
  df <- df[order(df$M1.M3.BetaDiff),] %>% mutate(M1.M3.BetaDiff.Order = nrow(df):1)
  
  df <- df %>% gather(comparison, Abs.Beta.Diff,c(M2.M1.BetaDiff, M2.M3.BetaDiff, M1.M3.BetaDiff))
  
  ## Rename difference category
  df$comparison_use <- ifelse(df$comparison == "M2.M1.BetaDiff", "Min-Dz", df$comparison)
  df$comparison_use <- ifelse(df$comparison == "M2.M3.BetaDiff", "Min-RskFct", df$comparison_use)
  df$comparison_use <- ifelse(df$comparison == "M1.M3.BetaDiff", "Dz-RskFct", df$comparison_use)
  df$comparison_use <- factor(df$comparison_use, ordered = T, levels = c("Min-Dz", "Min-RskFct", "Dz-RskFct"))
  
  return(df)
}

meta.results.wide <- cbind(
  subset(aric.chs.meta.results, model.cohort == "M2.meta") %>% .[, c("SeqId.use", "target_use", "Beta", "SE", "P")],
  subset(aric.chs.meta.results, model.cohort == "M1.meta") %>% .[, c("Beta", "SE", "P")],
  subset(aric.chs.meta.results, model.cohort == "M3.meta") %>% .[, c("Beta", "SE", "P")]) %>% as.data.frame() %>% 
  setNames(c("SeqId.use", "target_use", "M2.Beta", "M2.SE", "M2.P", "M1.Beta", "M1.SE", "M1.P", "M3.Beta", "M3.SE", "M3.P"))

meta.abs.beta.diff <- calc.min.abs.beta.diff(meta.results.wide)

## Heatmap
ggplot(subset(meta.abs.beta.diff, comparison != "M1.M3.BetaDiff" & SeqId.use %in% meta.norm.min.bonf.seqid),
       aes(x = comparison_use, y = reorder(target_use, M2.M3.BetaDiff.Order), fill=Abs.Beta.Diff)) +
  geom_tile() + scale_fill_continuous_divergingx(palette = 'RdBu') + scale_y_discrete(limits = rev) + theme_classic() +
  labs(fill = "Beta\nDiff", caption = "Meta Max BthRcSx Abs Beta Diff") +
  xlab("Model Comp") + ylab("Protein") + theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7), plot.caption = element_text(hjust = 1, size = 9))
ggsave("aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/Meta.Max.BthRcSx.MinBonfPtns.AbsBetaDiff_Heatmap.png", width = 5, height = 11, units = "in")

## Forest plot for top 10
meta.abs.beta.diff.top10.seqid <- subset(meta.abs.beta.diff, SeqId.use %in% meta.norm.min.bonf.seqid) %>% .[order(.$M2.M3.BetaDiff.Order), c("SeqId.use", "M2.M3.BetaDiff.Order")] %>% 
  unique() %>% .[1:10,] %>% as.data.frame() %>% setNames(c("SeqId.use", "M2.M3.BetaDiff.Order"))

forestplot(
  df = merge(aric.chs.meta.results, meta.abs.beta.diff.top10.seqid, by = "SeqId.use") %>% .[order(.$M2.M3.BetaDiff.Order),] %>% subset(., cohort == "meta"),
  name = target_use,
  estimate = Beta,
  se = SE,
  pvalue = P,
  psignif = (0.05/4955),
  colour = model_use,
  xlab = "Hazard Ratio",
  logodds = TRUE
) + scale_x_continuous(breaks = seq(0.5,2.0,0.1)) + theme_classic() +
  theme(text = element_text(size = 10), legend.position = "bottom", legend.direction = "horizontal", legend.margin=margin(c(-10,1,1,1)), legend.title = element_blank(), plot.caption = element_text(hjust = 1, size = 9)) +
  labs(caption = "Meta Max BthRcSx Abs Beta Diff Top 10")
ggsave("aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/Meta.Max.BthRcSx.MinBonfPtns.AbsBetaDiff.Top10_Forest.png", width = 7, height = 7, units = "in")

####Max Risk Beta Scatter Plots#####
plot_grid(
  ## Min vs. Dz
  ggplot(meta.results.wide, aes(x = M2.Beta, y = M1.Beta)) + geom_point(col = "royalblue") +
    geom_abline(intercept = 0, linetype = "dashed", col = "red") + xlim(-0.5,0.6) + ylim(-0.5,0.6) +
    annotate(geom = "label", label = paste0("r = ", round(cor(meta.results.wide$M2.Beta, meta.results.wide$M1.Beta), digits = 2)), x = -0.4, y = 0.6) +
    theme_bw() + xlab("Minimum Beta") + ylab("Disease Beta"),

  ## Min vs. RskFct
  ggplot(meta.results.wide, aes(x = M2.Beta, y = M3.Beta)) + geom_point(col = "royalblue") +
    geom_abline(intercept = 0, linetype = "dashed", col = "red") + xlim(-0.5,0.6) + ylim(-0.5,0.6) +
    annotate(geom = "label", label = paste0("r = ", round(cor(meta.results.wide$M2.Beta, meta.results.wide$M3.Beta), digits = 2)), x = -0.4, y = 0.6) +
    theme_bw() + xlab("Minimum Beta") + ylab("Risk Factor Beta"),

  ## RskFct vs. Dz
  ggplot(meta.results.wide, aes(x = M1.Beta, y = M3.Beta)) + geom_point(col = "royalblue") +
    geom_abline(intercept = 0, linetype = "dashed", col = "red") + xlim(-0.5,0.6) + ylim(-0.5,0.6) +
    annotate(geom = "label", label = paste0("r = ", round(cor(meta.results.wide$M1.Beta, meta.results.wide$M3.Beta), digits = 2)), x = -0.4, y = 0.6) +
    theme_bw() + xlab("Disease Beta") + ylab("Risk Factor Beta") + theme(plot.caption = element_text(hjust = 1, size = 7)) +
    labs(caption = "Meta Max BthRcSx"),
  nrow = 1, ncol = 3)
ggsave("aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/Meta.Max.BthRcSx.AllPtns.BetaComp_Scatter.png", width = 12, height = 5, units = "in")

####Save spot####
save.image("aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/Meta.Max.BthRcSx.Analysis.RData")
saveRDS(aric.chs.meta.results, "aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/Meta.Max.BthRcSx.Processed.RDS")
saveRDS(meta.results.wide, "aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/MetaOnly.Max.BthRcSx.Wide.RDS")
saveRDS(meta.abs.beta.diff, "aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/Meta.Max.BthRcSx.AbsBetaDiff.RDS")
