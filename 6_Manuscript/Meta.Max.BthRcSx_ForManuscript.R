####Description####
## Creates figures and data of meta-analysis for research letter
## Panels C-E (ver19) are generated
## Note that for tables, model names were kept as is written in the script but manually changed in Excel
## Dz+RskFct = Full model
## Min = Primary model

####Load packages####
library(magrittr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(forestplot)
library(cowplot)
library(writexl)
library(ggstats)
# library(ggforestplot)
# library(colorspace)
# library(ggforestplot)
# library(VennDiagram)
# library(tidyr)
# library(gridExtra)
# library(grid)

####Load data####
alldata <- readRDS("aric.prot.real.archive/scd.assoc/9_other/combine.data/All.SCAxProteomics.Results_10072022.RDS")
heatmap.color.limits <- readRDS("aric.prot.real.archive/scd.assoc/9_other/combine.data/Heatmap.Color.Limits_10072022.RDS")
aric.chs.meta.results <- subset(alldata, strat == "all")
meta.results <- subset(aric.chs.meta.results, cohort == "meta")
meta.results.wide <- readRDS("aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/MetaOnly.Max.BthRcSx.Wide.RDS")
meta.results.wide$target_use <- NULL
#Second column
target.ref <- readRDS("aric.prot.real.archive/scd.assoc/9_other/combine.data/Unique.SeqId.Target.RDS")
meta.results.wide <- merge(target.ref, meta.results.wide, by = "SeqId.use")

####Get unique protein number####
meta.m2.bonf <- subset(meta.results, P <(0.05/4955) & model == "M2")
length(unique(meta.m2.bonf$ARIC.uniprot_id))
#125, doubles: Q8NBJ4 (GOLM1), Q4LDE5 (SVEP1),  P16860 (N-terminal pro-BNP, BNP), P51858 (HDGF)
length(unique(meta.m2.bonf$ARIC.target))
#126 if use target, identical targets: 2 SVEP1, 2 HDGF, 2 GOLM1; N-terminal pro-BNP share same UniProt but are different proteins after processing
#FINAL DECISION: stick to 126 because BNP and NPPB are biologically different even if they come from the same gene

####[VER 3] Volcano plot####
#Modified from VER 2 code, takes out line in make.volc.data where HR is calculated since HR is already calculated in alldata

make.volc.data <- function(data, pcolumn, betacolumn, seqidcolumn, targetcolumn, wordlength, sigcutoff){
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
  #data$label <- ifelse(data[, targetcolumn] == "Deprecated", data[, entrezsymbcolumn], data$label)
  data$label <- ifelse(data$pcol != "Non-significant", data[, targetcolumn], NA)
  data$label <- ifelse(nchar(data$label)>wordlength | as.numeric(-log10(data[,pcolumn]))<=sigcutoff, NA, data$label)
  data$label <- ifelse(data[, seqidcolumn] == "7655-11", "N-terminal pro-BNP", data$label)
  
  ## Count
  #data <- data %>% add_count(pcol) %>% mutate(znn = paste0(pcol, ' (', n, ')')) 
  
  #data$HR <- exp(data[,betacolumn])
  data$LogP <- -log10(data[,pcolumn])
  return(data)
}

volc.man.data <- subset(aric.chs.meta.results, model.cohort == "M2.meta") %>% 
  make.volc.data(., pcolumn = "P", betacolumn = "Beta", seqidcolumn = "SeqId.use", 
                 targetcolumn = "target_use", wordlength = 15, sigcutoff = 6.5)

set.seed(10)

volc.man.plot <-
  ggplot(volc.man.data, aes(x=HR, y=LogP, col=pcol, label=label)) +
  geom_point(size = 1, alpha = 0.8) + 
  geom_segment(x=(0.6-0.1), xend=(1.8+0.2), y=-log10(0.05/4955), yend = -log10(0.05/4955), 
               linetype = "dashed", color = "black", linewidth = 0.2) +
  geom_text_repel(size = 1.85, max.overlaps = 35, col = "black", 
                  min.segment.length = 0, 
                  hjust = 0.5, nudge_y=0.3, nudge_x = 0.015, 
                  segment.size = 0.05, seed = 10, direction = "both") + 
  ylim(0,40) + ylab("-log10(P)") + xlab("Hazard Ratio") +
  scale_x_continuous(breaks = seq(0.5,2.0,0.25)) + coord_cartesian(xlim = c(0.5,2.0)) + 
  scale_color_manual(values = c("royalblue", "#D9D9D9"), breaks = c("Significant","Non-significant")) + 
  theme_classic(base_size = 9) + 
  theme(legend.title = element_blank(), 
        plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), units="line"),
        legend.position = c(0.18,0.97), 
        legend.key.size = unit(0.35, 'lines'), legend.margin=margin(c(1,1,1,0)),
        axis.title=element_text(size=7),
        axis.text=element_text(size=6), legend.text=element_text(size=6.5))

tiff("aric.prot.real.archive/scd.assoc/6_manuscript/ver17/Meta.Max.AllMods.AllPtns_VolcanoSingleCol.tiff", units = "in", res = 1600, width = 2.8, height = 2.5)
volc.man.plot
dev.off()
#Panel C in ver19

####[VER 3] Forest plot by cohort, with table####
## Updated to work with forestplot_3.1.3 and dplyr_1.1.3

tree <- rbind(subset(meta.m2.bonf, cohort == "meta" & Beta>0) %>% 
          .[order(.$P),] %>% .[1:5,] %>% mutate(order = 1:nrow(.)) %>% .[1:5,c("SeqId.use", "order")],
        subset(meta.m2.bonf, cohort == "meta" & Beta<0) %>% 
          .[order(.$P),] %>% .[1:5,] %>% mutate(order = 6:10) %>% .[1:5,c("SeqId.use", "order")]) %>% as.data.frame() %>%
  merge(aric.chs.meta.results, ., by = "SeqId.use")
tree <- tree[order(tree$order),] %>% subset(., model == "M2")

tree$cohort_prop <- ifelse(tree$cohort == "meta", "Meta-analysis", tree$cohort)
tree$cohort_prop <- ifelse(tree$cohort == "aric", "ARIC", tree$cohort_prop)
tree$cohort_prop <- ifelse(tree$cohort == "chs", "CHS", tree$cohort_prop)
#tree$cohort_prop <- factor(tree$cohort_prop, ordered = T, levels = c("ARIC", "CHS", "Meta-analysis"))
tree$cohort_prop <- factor(tree$cohort_prop, ordered = T, levels = c("Meta-analysis", "CHS", "ARIC"))
tree$P_plot <- tree$P
tree$P_plot <- ifelse(tree$cohort != "meta", NA, tree$P_plot)
tree$HR.CI <- paste0(sprintf("%0.2f", round(tree$HR, digits = 3)), " [", 
                     sprintf("%0.2f", round(tree$HR.CI.lower, digits = 3)),  "-", 
                     sprintf("%0.2f", round(tree$HR.CI.upper, digits = 3)), "]")

tree.plot <- tree[, c("target_use", "HR.CI", "P_plot", "HR", "HR.CI.lower", "HR.CI.upper", "cohort_prop", "cohort")]
protein.order <- tree.plot[, c("target_use")] %>% as.data.frame() %>% unique()
names(protein.order)[1] <- "target_use"
protein.order$prot.order <- seq(2,11, by = 1)
tree.plot <- merge(tree.plot, protein.order, by = "target_use")

forest.table <- subset(tree.plot, cohort == "meta") %>% .[, c("target_use", "HR.CI", "P_plot")]
names(forest.table) <- c("Protein", "HR [95% CI]", "P")
forest.table$P <- signif(forest.table$P, digits = 2)

temp <- as.data.frame(matrix(nrow = 1, ncol = 3))
temp[1,] <- c("Protein", "HR [95% CI]", "P")
names(temp) <- c("Protein", "HR [95% CI]", "P")
forest.table <- rbind(temp, forest.table) %>% as.data.frame()
forest.table$P[which(forest.table$Protein == "CILP2")] <- "6.0e-21"

temp <- as.data.frame(matrix(nrow = 3, ncol = 9))
names(temp) <- c("target_use", "HR.CI", "P_plot", "HR", "HR.CI.lower", "HR.CI.upper", "cohort_prop", "cohort", "prot.order")
temp$cohort <- c("meta", "aric", "chs")
temp$cohort_prop <- c("Meta-analysis", "ARIC", "CHS")
tree.plot <- rbind(temp, tree.plot)
tree.plot$prot.order[1:3]<- "1"

#tree.plot$order <- 1:nrow(tree.plot)
names(forest.table)[2] <- "HR.CI.Print"
tree.plot <- merge(tree.plot, forest.table, by.x = "target_use", by.y = "Protein", all.x = T)
tree.plot$prot.order <- as.numeric(tree.plot$prot.order)
tree.plot <- tree.plot[order(tree.plot$prot.order, tree.plot$cohort_prop),]

## Plot
# forest.table2 <- forest.table
# forest.table2[1,] <- c("Protein", "", "P")
forestplot <-
  subset(tree.plot, prot.order != 1) %>% group_by(cohort_prop) %>% 
  forestplot(mean = HR, lower = HR.CI.lower, upper = HR.CI.upper, 
             labeltext = c(target_use, HR.CI.Print, P), graph.pos = 2,
             align = c("l","l","l"), grid = F, zero = 1, cex = 0.03, 
             lineheight = "auto", xticks = seq(0.5,2.0,by=0.25), colgap=unit(3, "mm"), xticks.digits = 2,
             hrzl_lines=list("2" = gpar(lwd=1, col="#444444", lty = 2),
                             "3" = gpar(lwd=1, lineend="butt", col="#99999922"),
                             "4" = gpar(lwd=1, lineend="butt", col="#99999922"),
                             "5" = gpar(lwd=1, lineend="butt", col="#99999922"),
                             "6" = gpar(lwd=1, lineend="butt", col="#99999922"),
                             "7" = gpar(lwd=1, lineend="butt", col="#99999922"),
                             "8" = gpar(lwd=1, lineend="butt", col="#99999922"),
                             "9" = gpar(lwd=1, lineend="butt", col="#99999922"),
                             "10" = gpar(lwd=1, lineend="butt", col="#99999922"),
                             "11" = gpar(lwd=1, lineend="butt", col="#99999922")),
             col = fpColors(box = c("#F89F5B", "#7F00FF", "black"), 
                            lines = c("#F89F5B", "#7C5295", "black")), 
             xlab = "Hazard Ratio (HR)", boxsize = 0.12,
             clip = c(0.5, 2), fn.ci_norm = c(fpDrawCircleCI, fpDrawCircleCI, fpDrawDiamondCI), 
             legend_args = fpLegend(pos = list(x = 0.5, y = -0.14, "align" = "horizontal")),
             txt_gp = fpTxtGp(ticks = gpar(cex = 0.38),
                              xlab = gpar(cex = 0.40),
                              label = gpar(cex = 0.44))) |>
  fp_add_header(target_use = ("Protein"),
                HR.CI.Print = ("HR [95% CI]"),
                P = ("P"))

tiff("aric.prot.real.archive/scd.assoc/6_manuscript/ver17/Meta.Max.VolcanoSingle.Forestplot.Combined_withTable.tiff", units = "in", 
     res = 1600, width =  2.8, height = 2.8)
forestplot
dev.off()
#Panel D (ver19)

####Heat map####
meta.results$model.abrv_use <- as.character(meta.results$model.abrv_use)
meta.results$model.abrv_use <- factor(meta.results$model.abrv_use, ordered = T, levels = c("Dz+RskFct", "Dz", "Min"))

heatmap.horiz <- plot_grid(
  subset(meta.results, SeqId.use %in% meta.m2.bonf$SeqId.use & cohort == "meta" & Beta>0) %>%
    ggplot(., aes(y = model.abrv_use, x = reorder(target_use, order_meta.M2.p.beta), fill=HR)) +
      geom_tile() + geom_text(aes(label = PASSED_PLOT), color = "white", size = 2.5, vjust="center", hjust = "middle") +
    scale_fill_gradient2(low="white", mid="#DA012D", high="#960018", 
                         midpoint = 1.5, breaks = seq(1,2,0.2), limits = c(1,heatmap.color.limits[2])) +
      theme_classic() + labs(fill = "Hazard Ratio") + 
      ylab("Model") + xlab("Protein") + 
      theme(legend.position = "top", legend.direction = "horizontal", legend.justification = "right", 
            legend.key.height = unit(1.5, 'mm'), legend.title = element_text(size = 4.5), 
            legend.text = element_text(size = 4.3), legend.box.margin=margin(1,1,-10,1),
            axis.text.y = element_text(size = 6.5), axis.text.x = element_text(size = 5, angle = 60, hjust=1), 
            axis.title = element_text(size = 6.8)) + guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
    scale_y_discrete(labels = c("Dz+RskFct" = "Full", "Dz" = "Disease", "Min" = "Primary")),
  
  subset(meta.results, SeqId.use %in% meta.m2.bonf$SeqId.use & cohort == "meta" & Beta<0) %>%
    ggplot(., aes(y = model.abrv_use, x = reorder(target_use, order_meta.M2.p.beta), fill=HR)) +
    geom_tile() + geom_text(aes(label = PASSED_PLOT), color = "white", size = 2.5, vjust="center", hjust = "middle") +
    scale_fill_gradient2(low="#1034A6", mid="#1E90FF", high="white", 
                         midpoint = 0.75, breaks = seq(0.5,1,0.1), limits = c(heatmap.color.limits[1],1))+
    theme_classic() + labs(fill = "Hazard Ratio") + 
    ylab("Model") + xlab("Protein") + 
    theme(legend.position = "top", legend.direction = "horizontal", legend.justification = "right", 
          legend.key.height = unit(1.5, 'mm'), legend.title = element_text(size = 4.5), 
          legend.text = element_text(size = 4.3), legend.box.margin=margin(1,1,-10,1),
          axis.text.y = element_text(size = 6.5), axis.text.x = element_text(size = 5, angle = 60, hjust=1), 
          axis.title = element_text(size = 6.8)) + guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
    scale_y_discrete(labels = c("Dz+RskFct" = "Full", "Dz" = "Disease", "Min" = "Primary")),
  
  nrow = 2, ncol = 1, labels = "AUTO", label_size = 9, rel_heights = c(0.9,0.85))

tiff("aric.prot.real.archive/scd.assoc/6_manuscript/ver17/Meta.Max.AllMods.MinBonfPtns_Heatmap.Horiz.tiff", units = "in", res = 1600, width = 9, height = 4.6)
heatmap.horiz
dev.off()

####Full results, supplemental table####
ci.upper <- function(data, beta, SE){
  beta.data <- data[, beta]
  se.data <- data[, SE] 
  ciu <- exp(beta.data + (1.96*se.data))
  return(ciu)
}
ci.lower <- function(data, beta, SE){
  beta.data <- data[, beta]
  se.data <- data[, SE] 
  ciu <- exp(beta.data - (1.96*se.data))
  return(ciu)
}

## Merge to get aptamer information
full.results.supp <- merge(meta.results.wide, subset(meta.results, model == "M2") %>% 
                             .[, c("SeqId.use", "ARIC.targetfullname", "ARIC.uniprot_id", "ARIC.entrezgenesymbol", "ARIC.type", "ARIC.multi.apt")], by = "SeqId.use")
full.results.supp <- full.results.supp[order(full.results.supp$M2.P),]
full.results.supp <- full.results.supp[, c("target_use", "ARIC.targetfullname", "ARIC.uniprot_id", "ARIC.entrezgenesymbol", "SeqId.use", 
                                           "ARIC.type", "ARIC.multi.apt", "M2.Beta", "M2.SE", "M2.P", "M1.Beta", "M1.SE", "M1.P", "M3.Beta", "M3.SE", "M3.P")]

## Change column noting multi apt to YES/NO
full.results.supp$ARIC.multi.apt <- ifelse(full.results.supp$ARIC.multi.apt == 1, "YES", "NO")

## Add HR and CI
full.results.supp$M2.HR <- exp(full.results.supp$M2.Beta)
full.results.supp$M1.HR <- exp(full.results.supp$M1.Beta)
full.results.supp$M3.HR <- exp(full.results.supp$M3.Beta)

full.results.supp$M2.CI.lower <- ci.lower(full.results.supp, "M2.Beta", "M2.SE")
full.results.supp$M1.CI.lower <- ci.lower(full.results.supp, "M1.Beta", "M1.SE")
full.results.supp$M3.CI.lower <- ci.lower(full.results.supp, "M3.Beta", "M3.SE")

full.results.supp$M2.CI.upper <- ci.upper(full.results.supp, "M2.Beta", "M2.SE")
full.results.supp$M1.CI.upper <- ci.upper(full.results.supp, "M1.Beta", "M1.SE")
full.results.supp$M3.CI.upper <- ci.upper(full.results.supp, "M3.Beta", "M3.SE")

## Format numbers
rndRup <- function(data){
  HR <- grep("[.]HR", names(data))
  CI <- grep("[.]CI[.]", names(data))
  betas <- grep("[.]Beta", names(data))
  SE <- grep("[.]SE", names(data))
  nonP <- c(HR, CI, betas, SE) %>% unique()
  P <- grep("[.]P", names(data))
  
  data[,nonP] <- lapply(data[,nonP],sprintf,fmt="%.2f")
  data[,P] <- lapply(data[,P],format,scientific = TRUE,digits=2)
  return(data)
}

full.results.supp <- rndRup(full.results.supp)
full.results.supp$M2.HR.CI <- paste0(full.results.supp$M2.HR, " [", full.results.supp$M2.CI.lower,"-",full.results.supp$M2.CI.upper, "]")
full.results.supp$M1.HR.CI <- paste0(full.results.supp$M1.HR, " [", full.results.supp$M1.CI.lower,"-",full.results.supp$M1.CI.upper, "]")
full.results.supp$M3.HR.CI <- paste0(full.results.supp$M3.HR, " [", full.results.supp$M3.CI.lower,"-",full.results.supp$M3.CI.upper, "]")

## Select columns and rename
full.results.supp.final <- full.results.supp[, c("SeqId.use","target_use",
                                                 "M2.HR.CI", "M2.P",
                                                 "M1.HR.CI", "M1.P",
                                                 "M3.HR.CI", "M3.P")]

names(full.results.supp.final) <- c("SeqId", "Target", "Min HR [95% CI]", "Min P", "Dz HR [95% CI]", "Dz P", "Dz+RskFct HR [95% CI]", "Dz+RskFct P")

## Write out file, note that proteins are in order of Meta Max Min P value
#write_xlsx(full.results.supp.final[, c(1,5,8:19)], "Meta.BthRcSx.AllMods_FullResults.xlsx")
write_xlsx(full.results.supp.final, "aric.prot.real.archive/scd.assoc/6_manuscript/ver17/Meta.BthRcSx.AllMods_FullResults.xlsx")
write_xlsx(subset(full.results.supp.final, SeqId %in% meta.m2.bonf$SeqId.use), "aric.prot.real.archive/scd.assoc/6_manuscript/ver17/Meta.BthRcSx.AllMods_129FullResults.xlsx")

## Write out dictionary
dictionary <- full.results.supp[, c("SeqId.use", "target_use", "ARIC.uniprot_id", "ARIC.entrezgenesymbol", "ARIC.type", "ARIC.multi.apt")]
dictionary$ARIC.type <- ifelse(dictionary$ARIC.type == "Protein", "Ptn", 'Dep')
dictionary$ARIC.multi.apt <- ifelse(dictionary$ARIC.multi.apt == "NO", "N", 'Y')
names(dictionary) <- c("SeqId", "Target", "UniProt", "Gene", "Type", "Multi")
# write_xlsx(dictionary, "Protein.Dictionary_ForManuscript.xlsx")
           
####Table of min and rskfct sig ptns####
rotunda <- function(num){
  numa <- signif(num, digits = 3) %>% sprintf("%.2f",.)
  return(numa)
}
rndP <- function(num){
  numa <- signif(num, digits = 2)
  numa <- format(numa, scientific = TRUE)
  return(numa)
}

meta.max.min.rskfct.bonf.seqid <- subset(meta.results.wide, M2.P<(0.05/4955) & M3.P<(0.05/4955))$SeqId.use

meta.max.min.rskfct.bonf.table <- 
  cbind(
  subset(aric.chs.meta.results, model == "M2" & SeqId.use %in% meta.max.min.rskfct.bonf.seqid & cohort == "meta") %>% 
    .[order(.$order_meta.M2.p.beta),] %>% .[, c("target_use", "HR", "HR.CI.upper", "HR.CI.lower", "P", "model.cohort")] %>% 
    mutate(HR.CI = paste0(rotunda(.$HR), " [", rotunda(.$HR.CI.lower), "-", rotunda(.$HR.CI.upper), "]"), P = rndP(.$P)) %>% .[, c("target_use", "HR.CI", "P")],
  subset(aric.chs.meta.results, model == "M3" & SeqId.use %in% meta.max.min.rskfct.bonf.seqid & cohort == "meta") %>% 
    .[order(.$order_meta.M2.p.beta),] %>% .[, c("target_use", "HR", "HR.CI.upper", "HR.CI.lower", "P", "model.cohort")] %>%
    mutate(HR.CI = paste0(rotunda(.$HR), " [", rotunda(.$HR.CI.lower), "-", rotunda(.$HR.CI.upper), "]"), P = rndP(.$P)) %>% .[, c("HR.CI", "P")]) %>% 
  as.data.frame() %>% 
  setNames(c("Target", "Min HR [95% CI]", "Min P", "Dz+RskFct HR [95% CI]", "Dz+RskFct P"))

## Write out file, note that proteins are in order of Meta Max Min P Beta value
write_xlsx(meta.max.min.rskfct.bonf.table, "aric.prot.real.archive/scd.assoc/6_manuscript/ver17/Meta.Max.BthRcSx.MinRskFctBonfPtns_Table.xlsx")

####Plot for min and rskfct sig ptns####
meta.max.min.rskfct.bonf.table.plot <- 
  rbind(
    subset(aric.chs.meta.results, model == "M2" & 
             SeqId.use %in% meta.max.min.rskfct.bonf.seqid & cohort == "meta") %>% 
      .[order(.$order_meta.M2.p.beta),] %>% 
      .[, c("target_use", "Beta", "SE", "HR", "HR.CI.upper", "HR.CI.lower", "P", "model.cohort")],
    subset(aric.chs.meta.results, model == "M3" & 
             SeqId.use %in% meta.max.min.rskfct.bonf.seqid & cohort == "meta") %>% 
      .[order(.$order_meta.M2.p.beta),] %>% 
      .[, c("target_use", "Beta", "SE", "HR", "HR.CI.upper", "HR.CI.lower", "P", "model.cohort")]) %>%
  as.data.frame()

meta.max.min.rskfct.bonf.table.plot$model_use <- ifelse(meta.max.min.rskfct.bonf.table.plot$model.cohort == "M2.meta", "Primary", "Full")
meta.max.min.rskfct.bonf.table.plot$model_use <- factor(meta.max.min.rskfct.bonf.table.plot$model_use, ordered = T, levels = c("Full", "Primary"))
cutoff <- (0.05/4955)

meta.max.min.rskfct.bonf.table.plot$p_label <- NA
meta.max.min.rskfct.bonf.table.plot$p_label <- ifelse(meta.max.min.rskfct.bonf.table.plot$P<cutoff, "*", meta.max.min.rskfct.bonf.table.plot$p_label)
meta.max.min.rskfct.bonf.table.plot$p_label <- ifelse(meta.max.min.rskfct.bonf.table.plot$P<1e-10, "**", meta.max.min.rskfct.bonf.table.plot$p_label)
meta.max.min.rskfct.bonf.table.plot$p_label <- ifelse(meta.max.min.rskfct.bonf.table.plot$P<1e-15, "***", meta.max.min.rskfct.bonf.table.plot$p_label)
# meta.max.min.rskfct.bonf.table.plot$p_label <- ifelse(meta.max.min.rskfct.bonf.table.plot$P<1e-20, "+", meta.max.min.rskfct.bonf.table.plot$p_label)

meta.max.min.rskfct.bonf.table.plot$x_label <- meta.max.min.rskfct.bonf.table.plot$HR.CI.upper+0.04
meta.max.min.rskfct.bonf.table.plot$x_label <- ifelse(meta.max.min.rskfct.bonf.table.plot$p_label == "*", meta.max.min.rskfct.bonf.table.plot$HR.CI.upper+0.02, meta.max.min.rskfct.bonf.table.plot$x_label)
# meta.max.min.rskfct.bonf.table.plot$x2_label <- meta.max.min.rskfct.bonf.table.plot$HR.CI.lower-0.05
# meta.max.min.rskfct.bonf.table.plot$p_label2 <- ifelse(meta.max.min.rskfct.bonf.table.plot$model.cohort == "M2.meta", NA, meta.max.min.rskfct.bonf.table.plot$p_label)

protein.order2 <- meta.max.min.rskfct.bonf.table.plot[, c("target_use")] %>% as.data.frame() %>% setNames("target_use") %>% mutate(order = 1:nrow(.)) %>% subset(., order <=24)
meta.max.min.rskfct.bonf.table.plot <- merge(meta.max.min.rskfct.bonf.table.plot, protein.order2, by = "target_use")
# meta.max.min.rskfct.bonf.table.plot <- meta.max.min.rskfct.bonf.table.plot[order(as.numeric(meta.max.min.rskfct.bonf.table.plot$order)),]

tiff("aric.prot.real.archive/scd.assoc/6_manuscript/ver17/Meta.Max.BthRcSx.MinRskFctBonfPtns_Forestplot.tiff", units = "in", res = 1600, width = 3.2, height = 3.8)
ggplot(meta.max.min.rskfct.bonf.table.plot, 
       aes(y=reorder(target_use, -order),
           x=HR, xmin=HR.CI.lower, xmax=HR.CI.upper, col=model_use, fill=model_use)) + 
  geom_linerange(linewidth=0.8, position = position_dodge(width = 0.8)) +
  geom_point(size=1.5, shape=21, colour="white", stroke = 0.5, position=position_dodge(width = 0.8)) +
  geom_vline(xintercept=1, lty=2) +
  geom_text(aes(label = p_label, x = x_label), 
            size = 2.8, colour = "black", vjust = 0.8, position = position_dodge(width = 1)) +
  scale_color_manual(values = c("#508CC8", "goldenrod"), breaks = c("Primary", "Full")) + 
  scale_fill_manual(values = c("#508CC8", "goldenrod"), breaks = c("Primary", "Full")) +
  scale_x_continuous(breaks = seq(0.5,2.0,0.25)) + coord_cartesian(xlim = c(0.5,2.0)) + 
  theme_bw(base_size = 12) + xlab("Hazard Ratio") +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 5.5),
        axis.title.x = element_text(size = 7),
        legend.position = "top", legend.direction = "horizontal",
        plot.margin=unit(c(-0.4, 0.25, 0.1, 0.1), units="line"),
        legend.box.spacing = unit(0, "pt"), 
        legend.margin=margin(0,0,0,0),
        legend.justification = "right",
        legend.title=element_blank(),
        legend.text=element_text(size=5)) + 
  geom_stripped_rows(odd = "white", even = "lightgray", colour = "gray", alpha = 0.1) 
dev.off()
#Panel E (ver19)

session.info <- sessionInfo()
writeLines(capture.output(sessionInfo()), "aric.prot.real.archive/scd.assoc/6_manuscript/ver17/sessionInfo.txt")

####Save spot####
save.image("aric.prot.real.archive/scd.assoc/6_manuscript/ver17/BthRcSx.Meta.Results/BthRcSx.ForManuscript.RData")
