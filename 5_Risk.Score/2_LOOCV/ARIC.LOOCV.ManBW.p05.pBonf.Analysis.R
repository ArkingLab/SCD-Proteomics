####Description###
## Run PtRS association models for ARIC LOOCV

####Load packages####
library(magrittr)
library(survival)
library(ggplot2)
library(dplyr)

####Load data####
load("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/ThuyVy.SCA.Assoc.ARIC.V2Prot_BthRcSx.RData")

allresults <- readRDS("aric.prot.real.archive/scd.assoc/9_other/combine.data/All.SCAxProteomics.Results_10072022.RDS")

annotation <- subset(allresults, cohort == "aric" & model == "M3" & strat == "all") %>% unique()

scores_p05 <- readRDS("aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/LOOCV.p05.ManBW/ARIC.LOOCV.ManBW.p05.Proteins.SumScore.RDS")

scores_pBonf <- readRDS("aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/LOOCV.pBonf.ManBW/ARIC.LOOCV.ManBW.pBonf.Proteins.SumScore.RDS")

####Merge risk score####
scores_p05$sum_scaled <- scale(scores_p05$sum)
scores_pBonf$sum_scaled <- scale(scores_pBonf$sum)
names(scores_p05)[2:3] <- c("sum_p05", "sum_scaled_p05")
names(scores_pBonf)[2:3] <- c("sum_pBonf", "sum_scaled_pBonf")
scores <- merge(scores_p05, scores_pBonf, by = "pid")
pheno <- merge(final.alldata, scores, by = "pid")

####Run model####
#rskfct_covar <- "sex + race + age + bmi + height + center + smoke + eGFR_ckdepi + chd_chs + hf + diabetes + sbp + hnt_med + hdl + total_chol + tri + stroke + afib"

## p05
coxph(Surv(fuday_age_censored, scd_age_filt) ~ sum_scaled_p05 + sex + race + age + bmi + height + center + smoke + eGFR_ckdepi + chd_chs + hf + diabetes + sbp + hnt_med + hdl + total_chol + tri + stroke + afib, data = pheno) %>% 
  summary() %>% coef()

## pBonf
coxph(Surv(fuday_age_censored, scd_age_filt) ~ sum_scaled_pBonf + sex + race + age + bmi + height + center + smoke + eGFR_ckdepi + chd_chs + hf + diabetes + sbp + hnt_med + hdl + total_chol + tri + stroke + afib, data = pheno) %>% 
  summary() %>% coef()

####Create quintiles - scaled####
q_p05 <- quantile(pheno$sum_scaled_p05, probs = seq(0,1,0.2))
q_pBonf <- quantile(pheno$sum_scaled_pBonf, probs = seq(0,1,0.2))

pheno$quintile_p05 <- cut(pheno$sum_scaled_p05,breaks=q_p05,labels=c('1','2','3','4','5'), include.lowest=TRUE)
pheno$quintile_p05[is.na(pheno$quintile_p05)] <- 1
pheno$quintile_pBonf <- cut(pheno$sum_scaled_pBonf,breaks=q_pBonf,labels=c('1','2','3','4','5'), include.lowest=TRUE)
pheno$quintile_pBonf[is.na(pheno$quintile_pBonf)] <- 1

coxph(Surv(fuday_age_censored, scd_age_filt) ~ as.factor(quintile_p05) + sex + race + age + bmi + height + center + smoke + eGFR_ckdepi + chd_chs + hf + diabetes + sbp + hnt_med + hdl + total_chol + tri + stroke + afib, data = pheno) %>% summary() %>% coef() %>% .[1:4,]
coxph(Surv(fuday_age_censored, scd_age_filt) ~ as.factor(quintile_pBonf) + sex + race + age + bmi + height + center + smoke + eGFR_ckdepi + chd_chs + hf + diabetes + sbp + hnt_med + hdl + total_chol + tri + stroke + afib, data = pheno) %>% summary() %>% coef() %>% .[1:4,]

####Write out table for main results - p05####
nonquint_p05 <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sum_scaled_p05 + sex + race + age + bmi + height + center + smoke + eGFR_ckdepi + chd_chs + hf + diabetes + sbp + hnt_med + hdl + total_chol + tri + stroke + afib, data = pheno)
quint_p05 <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ as.factor(quintile_p05) + sex + race + age + bmi + height + center + smoke + eGFR_ckdepi + chd_chs + hf + diabetes + sbp + hnt_med + hdl + total_chol + tri + stroke + afib, data = pheno)

nonquint_p05_summary<- summary(nonquint_p05)$conf.int
quint_p05_summary <- summary(quint_p05)$conf.int

results.table_p05 <- data.frame(score=c("Sum", paste0("Q", 1:5)),
                            hr = c(nonquint_p05_summary[1], 1, quint_p05_summary[1:4]),
                            lower.ci = c(nonquint_p05_summary[1,3], NA, quint_p05_summary[1:4,3]),
                            upper.ci = c(nonquint_p05_summary[1,4], NA, quint_p05_summary[1:4,4]),
                            P=c(coef(summary(nonquint_p05))[1,5], NA, coef(summary(quint_p05))[1:4,5]))

results.table_p05.print <- results.table_p05
results.table_p05.print[,2:4] <- round(results.table_p05.print[,2:4], digits = 2)
results.table_p05.print[,5] <- signif(results.table_p05.print[,5], digits = 2)
results.table_p05.print$hr.ci <- paste0(results.table_p05.print$hr, " [", results.table_p05.print$lower.ci, "-", results.table_p05.print$upper.ci, "]")
results.table_p05.print$hr.ci[2] <- 1
results.table_p05.print.final <- results.table_p05.print[, c("score", "hr.ci", "P")]
names(results.table_p05.print.final) <- c("Variable", "HR [95% CI]", "P")
writexl::write_xlsx(results.table_p05.print.final, "aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/ARIC.LOOCV.Risk.Score.Results.ManBW.p05_Main.Table.xlsx")

####Write out table for main results - pBonf####
nonquint_pBonf <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sum_scaled_pBonf + sex + race + age + bmi + height + center + smoke + eGFR_ckdepi + chd_chs + hf + diabetes + sbp + hnt_med + hdl + total_chol + tri + stroke + afib, data = pheno)
quint_pBonf <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ as.factor(quintile_pBonf) + sex + race + age + bmi + height + center + smoke + eGFR_ckdepi + chd_chs + hf + diabetes + sbp + hnt_med + hdl + total_chol + tri + stroke + afib, data = pheno)

nonquint_pBonf_summary<- summary(nonquint_pBonf)$conf.int
quint_pBonf_summary <- summary(quint_pBonf)$conf.int

results.table_pBonf <- data.frame(score=c("Sum", paste0("Q", 1:5)),
                                hr = c(nonquint_pBonf_summary[1], 1, quint_pBonf_summary[1:4]),
                                lower.ci = c(nonquint_pBonf_summary[1,3], NA, quint_pBonf_summary[1:4,3]),
                                upper.ci = c(nonquint_pBonf_summary[1,4], NA, quint_pBonf_summary[1:4,4]),
                                P=c(coef(summary(nonquint_pBonf))[1,5], NA, coef(summary(quint_pBonf))[1:4,5]))

results.table_pBonf.print <- results.table_pBonf
results.table_pBonf.print[,2:4] <- round(results.table_pBonf.print[,2:4], digits = 2)
results.table_pBonf.print[,5] <- signif(results.table_pBonf.print[,5], digits = 2)
results.table_pBonf.print$hr.ci <- paste0(results.table_pBonf.print$hr, " [", results.table_pBonf.print$lower.ci, "-", results.table_pBonf.print$upper.ci, "]")
results.table_pBonf.print$hr.ci[2] <- 1
results.table_pBonf.print.final <- results.table_pBonf.print[, c("score", "hr.ci", "P")]
names(results.table_pBonf.print.final) <- c("Variable", "HR [95% CI]", "P")
writexl::write_xlsx(results.table_pBonf.print.final, "aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/ARIC.LOOCV.Risk.Score.Results.ManBW.pBonf_Main.Table.xlsx")

####Densityplot for score####
tiff("aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/ARIC.Sum.Scaled.Prot.Risk.Score.ManBW.p05_Density.tiff", units = "in", width = 5, height = 5, res = 600)
ggplot(pheno, aes(x=sum_scaled_p05)) + geom_density(fill = "royalblue", alpha = 0.5, colour = NA) + 
  geom_rug(alpha = 0.5) + theme_classic(base_size = 13) + ylab("Density") + xlab("Protein Risk Score")
dev.off()

tiff("aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/ARIC.Sum.Scaled.Prot.Risk.Score.ManBW.pBonf_Density.tiff", units = "in", width = 5, height = 5, res = 600)
ggplot(pheno, aes(x=sum_scaled_pBonf)) + geom_density(fill = "royalblue", alpha = 0.5, colour = NA) + 
  geom_rug(alpha = 0.5) + theme_classic(base_size = 13) + ylab("Density") + xlab("Protein Risk Score")
dev.off()

####Save spot####
save.image("aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/ARIC.LOOCV.ManBW.p05.pBonf.Analysis.RData")
