####Description####
## Creates figures and data of MR for research letter
## Panel F (ver19) is generated

####Load packages###
library(magrittr)
library(ggplot2)
#Version 3.4.0

library(TwoSampleMR)
library(tidyr)
library(writexl)

####Load function####
convert_mr_results <- function(data){
  data_wide <- data[, -c(1:3)]
  ivw <- subset(data_wide, method == "Inverse variance weighted") %>% .[, c("exposure", "b", "se", "pval")]
  ivw$adj_pval <- p.adjust(ivw$pval, method = "fdr")
  names(ivw)[2:ncol(ivw)] <- paste0("ivw.", names(ivw)[2:ncol(ivw)])
  
  simp.med <- subset(data_wide, method == "Simple median") %>% .[, c("exposure", "b", "se", "pval")]
  names(simp.med)[2:ncol(simp.med)] <- paste0("simp.med.", names(simp.med)[2:ncol(simp.med)])
  
  weight.med <- subset(data_wide, method == "Weighted median") %>% .[, c("exposure", "b", "se", "pval")]
  names(weight.med)[2:ncol(weight.med)] <- paste0("weight.med.", names(weight.med)[2:ncol(weight.med)])
  
  egger <- subset(data_wide, method == "MR Egger") %>% .[, c("exposure", "b", "se", "pval")]
  names(egger)[2:ncol(egger)] <- paste0("egger.", names(egger)[2:ncol(egger)])
  
  merged <- merge(data_wide[, c("exposure", "nsnp")] %>% unique(), ivw, by = "exposure", all.x = T)
  merged <- merge(merged, egger, by = "exposure", all.x = T, all.y = T)
  merged <- merge(merged, simp.med, by = "exposure", all.x = T, all.y = T)
  merged <- merge(merged, weight.med, by = "exposure", all.x = T, all.y = T)
  merged <- merge(harmonized.df[, c("exposure", "SeqId.use")] %>% unique(), merged, by = "exposure")
  merged$SeqId.use <- gsub("_", "-", merged$SeqId.use)
  
  key <- readxl::read_xlsx("aric.prot.real.archive/scd.assoc/6_manuscript/ver6/BthRcSx.Meta.Results/Meta.BthRcSx.AllMods_FullResults.xlsx") %>% as.data.frame() %>% .[, c("SeqId", "Target")]
  merged <- merge(key, merged, by.y = "SeqId.use", by.x = "SeqId")
  merged <- merged[order(merged$ivw.adj_pval),]
  
  merged$ivw.CI.lower <- ci.lower.noexp(merged, "ivw.b", "ivw.se")
  merged$ivw.CI.upper <- ci.upper.noexp(merged, "ivw.b", "ivw.se")
  
  merged$egger.CI.lower <- ci.lower.noexp(merged, "egger.b", "egger.se")
  merged$egger.CI.upper <- ci.upper.noexp(merged, "egger.b", "egger.se")
  
  merged$simp.med.CI.lower <- ci.lower.noexp(merged, "simp.med.b", "simp.med.se")
  merged$simp.med.CI.upper <- ci.upper.noexp(merged, "simp.med.b", "simp.med.se")
  
  merged$weight.med.CI.lower <- ci.lower.noexp(merged, "weight.med.b", "weight.med.se")
  merged$weight.med.CI.upper <- ci.upper.noexp(merged, "weight.med.b", "weight.med.se")
  
  names(merged) <- gsub(".b", ".Beta", names(merged))
  names(merged) <- gsub(".se", ".SE", names(merged))
  names(merged) <- gsub("pval", "P", names(merged))
  names(merged) <- gsub("_P", ".P", names(merged))
  
  merged.formatted <- rndRup(merged)
  merged.formatted.print <- merged.formatted[, c("Target", "SeqId", "nsnp", "ivw.Beta", "ivw.CI.lower", "ivw.CI.upper", "ivw.P", "ivw.adj.P",
                                                 "egger.Beta", "egger.CI.lower", "egger.CI.upper", "egger.P",
                                                 "simp.med.Beta", "simp.med.CI.lower", "simp.med.CI.upper", "simp.med.P",
                                                 "weight.med.Beta", "weight.med.CI.lower", "weight.med.CI.upper", "weight.med.P")]
  merged.formatted.print <- merged.formatted.print[merged.formatted.print$ivw.Beta!= "NA",]
  names(merged.formatted.print) <- c("Target", "SeqId", "nSNP", "IVW Beta", "IVW Lower 95% CI", "IVW Upper 95% CI", "IVW P", "IVW Adj. P",
                                     "MR-Egger Beta", "MR-Egger Lower 95% CI", "MR-Egger Upper 95% CI", "MR-Egger P",
                                     "Simple Median Beta", "Simple Median Lower 95% CI", "Simple Median Upper 95% CI", "Simple Median P",
                                     "Median Weighted Beta", "Median Weighted Lower 95% CI", "Median Weighted Upper 95% CI", "Median Weighted P")
  return(merged.formatted.print)
}
ci.upper.noexp <- function(data, beta, SE){
  beta.data <- data[, beta]
  se.data <- data[, SE] 
  ciu <- beta.data + (1.96*se.data)
  return(ciu)
}
ci.lower.noexp <- function(data, beta, SE){
  beta.data <- data[, beta]
  se.data <- data[, SE] 
  ciu <- beta.data - (1.96*se.data)
  return(ciu)
}
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
convert_wald_results <- function(data){
  data_wide <- data[, -c(1:3)]
  wald <- subset(data_wide, method == "Wald ratio") %>% .[, c("exposure", "b", "se", "pval")]
  names(wald)[2:ncol(wald)] <- paste0("wald.", names(wald)[2:ncol(wald)])
  wald.print <- merge(harmonized.df[, c("exposure", "SeqId.use")] %>% unique(), wald, by = "exposure")
  wald.print$SeqId <- gsub("_", "-", wald.print$SeqId.use)
  key <- readxl::read_xlsx("aric.prot.real.archive/scd.assoc/6_manuscript/ver6/BthRcSx.Meta.Results/Meta.BthRcSx.AllMods_FullResults.xlsx") %>% as.data.frame() %>% .[, c("SeqId", "Target")]
  wald.print <- merge(key, wald.print, by="SeqId")
  wald.print <- wald.print[, -c(3:4)]
  names(wald.print) <- c("SeqId", "Target", "Wald Beta", "Wald SE", "Wald P")
  wald.print <- wald.print[, c(2,1,3:5)]
  wald.print <- rndRup(wald.print)
  return(wald.print)
}
convert_mr_results_BSE <- function(data){
  data_wide <- data[, -c(1:3)]
  ivw <- subset(data_wide, method == "Inverse variance weighted") %>% .[, c("exposure", "b", "se", "pval")]
  ivw$adj_pval <- p.adjust(ivw$pval, method = "fdr")
  names(ivw)[2:ncol(ivw)] <- paste0("ivw.", names(ivw)[2:ncol(ivw)])
  
  simp.med <- subset(data_wide, method == "Simple median") %>% .[, c("exposure", "b", "se", "pval")]
  names(simp.med)[2:ncol(simp.med)] <- paste0("simp.med.", names(simp.med)[2:ncol(simp.med)])
  
  weight.med <- subset(data_wide, method == "Weighted median") %>% .[, c("exposure", "b", "se", "pval")]
  names(weight.med)[2:ncol(weight.med)] <- paste0("weight.med.", names(weight.med)[2:ncol(weight.med)])
  
  egger <- subset(data_wide, method == "MR Egger") %>% .[, c("exposure", "b", "se", "pval")]
  names(egger)[2:ncol(egger)] <- paste0("egger.", names(egger)[2:ncol(egger)])
  
  merged <- merge(data_wide[, c("exposure", "nsnp")] %>% unique(), ivw, by = "exposure", all.x = T)
  merged <- merge(merged, egger, by = "exposure", all.x = T, all.y = T)
  merged <- merge(merged, simp.med, by = "exposure", all.x = T, all.y = T)
  merged <- merge(merged, weight.med, by = "exposure", all.x = T, all.y = T)
  merged <- merge(harmonized.df[, c("exposure", "SeqId.use")] %>% unique(), merged, by = "exposure")
  merged$SeqId.use <- gsub("_", "-", merged$SeqId.use)
  
  key <- readxl::read_xlsx("aric.prot.real.archive/scd.assoc/6_manuscript/ver6/BthRcSx.Meta.Results/Meta.BthRcSx.AllMods_FullResults.xlsx") %>% as.data.frame() %>% .[, c("SeqId", "Target")]
  merged <- merge(key, merged, by.y = "SeqId.use", by.x = "SeqId")
  merged <- merged[order(merged$ivw.adj_pval),]
  
  merged$ivw.CI.lower <- ci.lower.noexp(merged, "ivw.b", "ivw.se")
  merged$ivw.CI.upper <- ci.upper.noexp(merged, "ivw.b", "ivw.se")
  
  merged$egger.CI.lower <- ci.lower.noexp(merged, "egger.b", "egger.se")
  merged$egger.CI.upper <- ci.upper.noexp(merged, "egger.b", "egger.se")
  
  merged$simp.med.CI.lower <- ci.lower.noexp(merged, "simp.med.b", "simp.med.se")
  merged$simp.med.CI.upper <- ci.upper.noexp(merged, "simp.med.b", "simp.med.se")
  
  merged$weight.med.CI.lower <- ci.lower.noexp(merged, "weight.med.b", "weight.med.se")
  merged$weight.med.CI.upper <- ci.upper.noexp(merged, "weight.med.b", "weight.med.se")
  
  names(merged) <- gsub(".b", ".Beta", names(merged))
  names(merged) <- gsub(".se", ".SE", names(merged))
  names(merged) <- gsub("pval", "P", names(merged))
  names(merged) <- gsub("_P", ".P", names(merged))
  
  merged.formatted <- rndRup(merged)
  merged.formatted.print <- merged.formatted[, c("Target", "SeqId", "nsnp", "ivw.Beta", "ivw.SE", "ivw.P", "ivw.adj.P",
                                                 "egger.Beta", "egger.SE", "egger.P",
                                                 "simp.med.Beta", "simp.med.SE", "simp.med.P",
                                                 "weight.med.Beta", "weight.med.SE", "weight.med.P")]
  merged.formatted.print <- merged.formatted.print[merged.formatted.print$ivw.Beta!= "NA",]
  names(merged.formatted.print) <- c("Target", "SeqId", "nSNP", "IVW Beta", "IVW SE", "IVW P", "IVW Adj. P",
                                     "MR-Egger Beta", "MR-Egger SE", "MR-Egger P",
                                     "Simple Median Beta", "Simple Median SE", "Simple Median P",
                                     "Median Weighted Beta", "Median Weighted SE", "Median Weighted P")
  return(merged.formatted.print)
}

####Load data####
load("aric.prot.real.archive/scd.assoc/7_MR/analysis/Pietzner.2021/topmed.panel/run07/ARIC.TM.p5e2.maf001.rsq01.kb1e3_MR.Analysis.rds")

####Plot SVEP1 result####
#SVEP1_11109_56 - probe that passed 0.05 FDR
#SVEP1*_11178_21 - other probe

## Plot 56 probe
svep1.56_harmonized.df <- subset(harmonized.df, exposure == "SVEP1_11109_56" & mr_keep == TRUE)

## Flip betas
index <- svep1.56_harmonized.df$beta.exposure < 0
svep1.56_harmonized.df$beta.exposure[index] <- svep1.56_harmonized.df$beta.exposure[index] * -1
svep1.56_harmonized.df$beta.outcome[index] <- svep1.56_harmonized.df$beta.outcome[index] * -1

svep1.56_mr_results <- subset(mr_results, exposure == "SVEP1_11109_56")
svep1.56_mr_results$a <- 0
temp <- mr_egger_regression(svep1.56_harmonized.df$beta.exposure, svep1.56_harmonized.df$beta.outcome, svep1.56_harmonized.df$se.exposure, svep1.56_harmonized.df$se.outcome, default_parameters())
svep1.56_mr_results$a[svep1.56_mr_results$method == "MR Egger"] <- temp$b_i
svep1.56_mr_results <- generate_odds_ratios(svep1.56_mr_results)

svep1.56_plot <-
  ggplot(data=svep1.56_harmonized.df, aes(x=beta.exposure, y=beta.outcome)) +
  geom_errorbar(aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
  geom_errorbarh(aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
  geom_point(size = 0.8) +
  geom_abline(data=subset(svep1.56_mr_results, method == "Inverse variance weighted"),
              aes(intercept=a, slope=b, colour=method), linewidth = 0.2, colour = "royalblue") +
  geom_abline(data=subset(svep1.56_mr_results, method == "Inverse variance weighted"),
              aes(intercept=a, slope=lo_ci, colour=method), linewidth = 0.2, colour = "royalblue", linetype = "dashed") +
  geom_abline(data=subset(svep1.56_mr_results, method == "Inverse variance weighted"),
              aes(intercept=a, slope=up_ci, colour=method), linewidth = 0.2, colour = "royalblue", linetype = "dashed") +
  geom_hline(yintercept=0, linewidth = 0.1) +
  labs(colour="MR Test", x="SNP Effect on SVEP1 (in SD units)", y="SNP Effect on SCD (log odds of risk)") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) + theme_bw() +
  theme(legend.position="none", 
        plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), units="line"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title=element_text(size=8))

tiff(file = "aric.prot.real.archive/scd.assoc/6_manuscript/ver17/MR/SVEP1_56_scatterplot_HiRES.tiff", res = 1600, units = "in", height = 2.4, width = 2.4)
svep1.56_plot
dev.off()
#Panel F (ver19)

## Plot 21 probe
svep1.21_harmonized.df <- subset(harmonized.df, exposure == "SVEP1*_11178_21" & mr_keep == TRUE)

## Flip betas
index <- svep1.21_harmonized.df$beta.exposure < 0
svep1.21_harmonized.df$beta.exposure[index] <- svep1.21_harmonized.df$beta.exposure[index] * -1
svep1.21_harmonized.df$beta.outcome[index] <- svep1.21_harmonized.df$beta.outcome[index] * -1

svep1.21_mr_results <- subset(mr_results, exposure == "SVEP1*_11178_21")
svep1.21_mr_results$a <- 0
temp <- mr_egger_regression(svep1.21_harmonized.df$beta.exposure, svep1.21_harmonized.df$beta.outcome, svep1.21_harmonized.df$se.exposure, svep1.21_harmonized.df$se.outcome, default_parameters())
svep1.21_mr_results$a[svep1.21_mr_results$method == "MR Egger"] <- temp$b_i
svep1.21_mr_results <- generate_odds_ratios(svep1.21_mr_results)

svep1.21_plot <-
  ggplot(data=svep1.21_harmonized.df, aes(x=beta.exposure, y=beta.outcome)) +
  geom_errorbar(aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
  geom_errorbarh(aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
  geom_point(size = 1.5) +
  geom_abline(data=subset(svep1.21_mr_results, method == "Inverse variance weighted"),
              aes(intercept=a, slope=b, colour=method), size = 0.3, colour = "royalblue") +
  geom_abline(data=subset(svep1.21_mr_results, method == "Inverse variance weighted"),
              aes(intercept=a, slope=lo_ci, colour=method), size = 0.3, colour = "royalblue", linetype = "dashed") +
  geom_abline(data=subset(svep1.21_mr_results, method == "Inverse variance weighted"),
              aes(intercept=a, slope=up_ci, colour=method), size = 0.3, colour = "royalblue", linetype = "dashed") +
  geom_hline(yintercept=0, size = 0.2) +
  labs(colour="MR Test", x="SNP Effect on SVEP1* (in SD units)", y="SNP Effect on SCD (log odds of risk)") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) + theme_bw() +
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

tiff(file = "aric.prot.real.archive/scd.assoc/6_manuscript/ver17/MR/SVEP1_21_scatterplot.tiff", res = 1200, units = "in", height = 4.5, width = 4.5)
svep1.21_plot
dev.off()

####Make supplemental MR tables####
#Note that function corrects FDR P value - only use max number of IVW analyses, doesn't include Wald ratio proteins

p05 <- convert_mr_results(mr_results_edit)
p05.wald <- convert_wald_results(mr_results_edit)

mr_results_edit_p1e3 <- readRDS("aric.prot.real.archive/scd.assoc/7_MR/analysis/Pietzner.2021/topmed.panel/run05/ARIC.TM.p1e3.maf001.rsq01.kb1e3_mr.results.edit.RDS")
p1e3 <- convert_mr_results(mr_results_edit_p1e3)
p1e3.wald <- convert_wald_results(mr_results_edit_p1e3)

mr_results_edit_p1e5 <- readRDS("aric.prot.real.archive/scd.assoc/7_MR/analysis/Pietzner.2021/topmed.panel/run17/ARIC.TM.p1e5.maf001.rsq01.kb1e3_mr.results.edit.RDS")
p1e5 <- convert_mr_results(mr_results_edit_p1e5)
p1e5.wald <- convert_wald_results(mr_results_edit_p1e5)

mr_results_edit_p1e11 <- readRDS("aric.prot.real.archive/scd.assoc/7_MR/analysis/Pietzner.2021/topmed.panel/run16/ARIC.TM.p1e11.maf001.rsq01.kb1e3_mr.results.edit.RDS")
p1e11 <- convert_mr_results(mr_results_edit_p1e11)
p1e11.wald <- convert_wald_results(mr_results_edit_p1e11)

sheets <- list("p05_results" = p05, 
               "p1e3_results" = p1e3, 
               "p1e5_results" = p1e5,
               "p1e11_results" = p1e11)
write_xlsx(sheets, "aric.prot.real.archive/scd.assoc/6_manuscript/ver17/MR/p005.p1e3.p1e5.p1e11.Full.MR.Results.xlsx")

p05$P <- "0.05"
p1e3$P <- "1e-3"
p1e5$P <- "1e-5"
p1e11$P <- "1e-11"
allresults <- rbind(p05, p1e3, p1e5, p1e11)
allresults_names <- names(allresults)
allresults_names[21] <- "pQTL P"
allresults <- as.data.frame(lapply(allresults, stringr::str_remove_all, "NA"))
names(allresults) <- allresults_names
allresults$`IVW Beta [95% CI]` <- paste0(allresults[,"IVW Beta"], " [", allresults[,"IVW Lower 95% CI"], "-", allresults[,"IVW Upper 95% CI"], "]")
allresults$`MR-Egger Beta [95% CI]` <- paste0(allresults[,"MR-Egger Beta"], " [", allresults[,"MR-Egger Lower 95% CI"], "-", allresults[,"MR-Egger Upper 95% CI"], "]")
allresults$`Simple Median Beta [95% CI]` <- paste0(allresults[,"Simple Median Beta"], " [", allresults[,"Simple Median Lower 95% CI"], "-", allresults[,"Simple Median Upper 95% CI"], "]")
allresults$`Median Weighted Beta [95% CI]` <- paste0(allresults[,"Median Weighted Beta"], " [", allresults[,"Median Weighted Lower 95% CI"], "-", allresults[,"Median Weighted Upper 95% CI"], "]")
write_xlsx(allresults[, c("Target", "SeqId", "nSNP", "IVW Beta [95% CI]", "IVW P", "IVW Adj. P", 
                          "MR-Egger Beta [95% CI]", "MR-Egger P",
                          "Simple Median Beta [95% CI]", "Simple Median P", "Median Weighted Beta [95% CI]", "Median Weighted P", "pQTL P")], 
           "aric.prot.real.archive/scd.assoc/6_manuscript/ver17/MR/p05.p1e3.p1e5.p1e11.Full.MR.Results_SuppTable.xlsx")

####Make supplemental table with Beta, SE, and P####

p05_BSE <- convert_mr_results_BSE(mr_results_edit)
p1e3_BSE <- convert_mr_results_BSE(mr_results_edit_p1e3)
p1e5_BSE <- convert_mr_results_BSE(mr_results_edit_p1e5)
p1e11_BSE <- convert_mr_results_BSE(mr_results_edit_p1e11)

p05_BSE$P <- "0.05"
p1e3_BSE$P <- "1.0e-03"
p1e5_BSE$P <- "1.0e-05"
p1e11_BSE$P <- "1.0e-11"

allresults_BSE <- rbind(p05_BSE, p1e3_BSE, p1e5_BSE, p1e11_BSE)
allresults_BSE_names <- names(allresults_BSE)
allresults_BSE_names[17] <- "pQTL P"
allresults_BSE <- as.data.frame(lapply(allresults_BSE, stringr::str_remove_all, "NA"))
names(allresults_BSE) <- allresults_BSE_names

write_xlsx(allresults_BSE, "aric.prot.real.archive/scd.assoc/6_manuscript/ver17/MR/p05.p1e3.p1e5.p1e11.Full.MR.Results.BSE_SuppTable.xlsx")

####Make supplemental table for SVEP1 SNPs####
getEffect <- function(data, cutoff){
  sub <- subset(data, exposure == "SVEP1_11109_56" & mr_keep == TRUE) %>% .[, c("SNP", "beta.exposure", "beta.outcome", "se.outcome")]
  sub$weight <- (1/sub$se.outcome^2)
  names(sub)[2:5] <- paste0(names(sub)[2:5], "_p", cutoff)
  return(sub)
}

svep1.56_snps_p05 <- getEffect(harmonized.df, "05")

harmonized.df_p1e3 <- readRDS("aric.prot.real.archive/scd.assoc/7_MR/analysis/Pietzner.2021/topmed.panel/run05/ARIC.TM.p1e3.maf001.rsq01.kb1e3_harmonized.df.RDS")
svep1.56_snps_p1e3 <- getEffect(harmonized.df_p1e3, "1e3")

harmonized.df_p1e5 <- readRDS("aric.prot.real.archive/scd.assoc/7_MR/analysis/Pietzner.2021/topmed.panel/run17/ARIC.TM.p1e5.maf001.rsq01.kb1e3_harmonized.df.RDS")
svep1.56_snps_p1e5 <- getEffect(harmonized.df_p1e5, "1e5")

harmonized.df_p1e11 <- readRDS("aric.prot.real.archive/scd.assoc/7_MR/analysis/Pietzner.2021/topmed.panel/run16/ARIC.TM.p1e11.maf001.rsq01.kb1e3_harmonized.df.RDS")
svep1.56_snps_p1e11 <- getEffect(harmonized.df_p1e11, "1e11")

svep1.56_all_snps <- merge(svep1.56_snps_p05, svep1.56_snps_p1e3, by = "SNP", all.x = TRUE, all.y = TRUE)
svep1.56_all_snps <- merge(svep1.56_all_snps, svep1.56_snps_p1e5, by = "SNP", all.x = TRUE, all.y = TRUE)
svep1.56_all_snps <- merge(svep1.56_all_snps, svep1.56_snps_p1e11, by = "SNP", all.x = TRUE, all.y = TRUE)

svep1.56_all_snps_final <- svep1.56_all_snps[, -c(4,8,12,16)]
svep1.56_all_snps_final[, c(2,3,5,6,8,9,11,12)] <- lapply(svep1.56_all_snps_final[, c(2,3,5,6,8,9,11,12)] ,sprintf,fmt="%.3f")
svep1.56_all_snps_final[, c(4,7,10,13)] <- lapply(svep1.56_all_snps_final[, c(4,7,10,13)], sprintf,fmt="%.3f")
svep1.56_all_snps_final <- as.data.frame(lapply(svep1.56_all_snps_final, stringr::str_remove_all, "NA"))
pvalues <- c("0.05", "1e-3", "1e-5", "1e-11")
col.names <- c("SVEP1 Effect", "SCA Effect", "Weight")
names(svep1.56_all_snps_final) <- c("SNP", apply(expand.grid(col.names, pvalues), 1, paste, collapse=" P="))
write_xlsx(svep1.56_all_snps_final, "aric.prot.real.archive/scd.assoc/6_manuscript/ver17/MR/SVEP1_11109-56_SNP.Effect.Weights.xlsx")

####Make supplemental table for SVEP1* SNPs####
getEffect2 <- function(data, cutoff){
  sub <- subset(data, exposure == "SVEP1*_11178_21" & mr_keep == TRUE) %>% .[, c("SNP", "beta.exposure", "beta.outcome", "se.outcome")]
  sub$weight <- (1/sub$se.outcome^2)
  names(sub)[2:5] <- paste0(names(sub)[2:5], "_p", cutoff)
  return(sub)
}

svep1.21_snps_p05 <- getEffect(harmonized.df, "05")

harmonized.df_p1e3 <- readRDS("aric.prot.real.archive/scd.assoc/7_MR/analysis/Pietzner.2021/topmed.panel/run05/ARIC.TM.p1e3.maf001.rsq01.kb1e3_harmonized.df.RDS")
svep1.21_snps_p1e3 <- getEffect(harmonized.df_p1e3, "1e3")

harmonized.df_p1e5 <- readRDS("aric.prot.real.archive/scd.assoc/7_MR/analysis/Pietzner.2021/topmed.panel/run17/ARIC.TM.p1e5.maf001.rsq01.kb1e3_harmonized.df.RDS")
svep1.21_snps_p1e5 <- getEffect(harmonized.df_p1e5, "1e5")

harmonized.df_p1e11 <- readRDS("aric.prot.real.archive/scd.assoc/7_MR/analysis/Pietzner.2021/topmed.panel/run16/ARIC.TM.p1e11.maf001.rsq01.kb1e3_harmonized.df.RDS")
svep1.21_snps_p1e11 <- getEffect(harmonized.df_p1e11, "1e11")

svep1.21_all_snps <- merge(svep1.21_snps_p05, svep1.21_snps_p1e3, by = "SNP", all.x = TRUE, all.y = TRUE)
svep1.21_all_snps <- merge(svep1.21_all_snps, svep1.21_snps_p1e5, by = "SNP", all.x = TRUE, all.y = TRUE)
svep1.21_all_snps <- merge(svep1.21_all_snps, svep1.21_snps_p1e11, by = "SNP", all.x = TRUE, all.y = TRUE)

svep1.21_all_snps_final <- svep1.21_all_snps[, -c(4,8,12,16)]
svep1.21_all_snps_final[, c(2,3,5,6,8,9,11,12)] <- lapply(svep1.21_all_snps_final[, c(2,3,5,6,8,9,11,12)] ,sprintf,fmt="%.3f")
svep1.21_all_snps_final[, c(4,7,10,13)] <- lapply(svep1.21_all_snps_final[, c(4,7,10,13)], sprintf,fmt="%.3f")
svep1.21_all_snps_final <- as.data.frame(lapply(svep1.21_all_snps_final, stringr::str_remove_all, "NA"))
pvalues <- c("0.05", "1e-3", "1e-5", "1e-11")
col.names <- c("SVEP1 Effect", "SCA Effect", "Weight")
names(svep1.21_all_snps_final) <- c("SNP", apply(expand.grid(col.names, pvalues), 1, paste, collapse=" P="))
write_xlsx(svep1.21_all_snps_final, "aric.prot.real.archive/scd.assoc/6_manuscript/ver17/MR/SVEP1*_11178-21_SNP.Effect.Weights.xlsx")

####Save spot###
save.image("aric.prot.real.archive/scd.assoc/6_manuscript/ver17/MR/MR.Results.ForManuscript.RData")
writeLines(capture.output(sessionInfo()), "aric.prot.real.archive/scd.assoc/6_manuscript/ver17/MR/sessionInfo.txt")
