####Description####
## Calculate C-stat and NRI for clinical only risk factor SCD model (base), base + PtRS model
## Archive version removes printouts and unused code

####Load packages####
library(survival)
library(magrittr)
library(boot)
# library(nricens)
library(survIDINRI)
library(dplyr)
# library(CsChange)
# library(survC1)

####Load functions####
CI <- function(beta, se){
  up <- beta + (1.96*se)
  down <- beta - (1.96*se)
  print("UP:")
  print(up)
  print("LOW:")
  print(down)
}

####Load data####
pheno <- readRDS("aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/ARIC.LOOCV.Pheno.with.Risk.Score.RDS")

## Set seed for calculations
set.seed(123)

####Max C-stat####

## Define models
base.mod <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + race + age + bmi + height + 
                    smoke + chd_chs + hf + diabetes + sbp + 
                    hnt_med + hdl + total_chol + tri + stroke + 
                    afib, data = pheno)

new.mod <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + race + age + bmi + height + 
                   smoke + chd_chs + hf + diabetes + sbp + 
                   hnt_med + hdl + total_chol + tri + stroke + 
                   afib + sum_scaled_p05, data = pheno)

## Calculate C-stat
harrell <- concordance(base.mod, new.mod, timewt = "n")
uno <- concordance(base.mod, new.mod, timewt = "n/G2")

## Harrell change CI
#Harrell discards censored observations, so metric depends on censoring distribution
#Use simple resampling, from manual: The simplest is case resampling which simply resamples with replacement from the observations. 

#Get CI for change via bootstrapping
func <- function(data, ind, type){
  mod1 <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + race + age + bmi + height + 
                  smoke + chd_chs + hf + diabetes + sbp + 
                  hnt_med + hdl + total_chol + tri + stroke + 
                  afib, data = data[ind,])
  mod2 <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + race + age + bmi + height + 
                  smoke + chd_chs + hf + diabetes + sbp + 
                  hnt_med + hdl + total_chol + tri + stroke + 
                  afib + sum_scaled_p05, data = data[ind,])
  if(type == "h"){
    cstat <- concordance(mod1, mod2, timewt = "n")
  } else if(type == "u"){
    cstat <- concordance(mod1, mod2, timewt = "n/G2")
  }
  cstat.diff <- cstat$concordance[["mod2"]]-cstat$concordance[["mod1"]]
  cstat.diff
}

res <- censboot(pheno[, c("pid", "fuday_age_censored", "scd_age_filt", "sex", "race", 
                          "age", "bmi", "height", "smoke", "chd_chs", "hf", 
                          "diabetes", "sbp", "hnt_med", "hdl", "total_chol", 
                          "tri", "stroke", "afib", "sum_scaled_p05")], 
                func, R = 1000, index = c(3,2), type = "h")
res.ci <- boot.ci(res, type = c("norm", "perc"))

## Uno change CI
#Uno is censoring independent
res.uno <- censboot(pheno[, c("pid", "fuday_age_censored", "scd_age_filt", "sex", "race", 
                          "age", "bmi", "height", "smoke", "chd_chs", "hf", 
                          "diabetes", "sbp", "hnt_med", "hdl", "total_chol", 
                          "tri", "stroke", "afib", "sum_scaled_p05")], 
                func, R = 1000, index = c(3,2), type = "u")

res.uno.ci <- boot.ci(res.uno, type = c("norm", "perc"))

####5YR C-stat####

## Harrell
harrell5YR <- concordance(base.mod, new.mod, timewt = "n", ymin = 0, ymax = 365.25*5)

#Get CI for change via bootstrapping
func5YR <- function(data, ind, type){
  mod1 <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + race + age + bmi + height + 
                  smoke + chd_chs + hf + diabetes + sbp + 
                  hnt_med + hdl + total_chol + tri + stroke + 
                  afib, data = data[ind,])
  mod2 <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + race + age + bmi + height + 
                  smoke + chd_chs + hf + diabetes + sbp + 
                  hnt_med + hdl + total_chol + tri + stroke + 
                  afib + sum_scaled_p05, data = data[ind,])
  if(type == "h"){
    cstat <- concordance(mod1, mod2, timewt = "n", ymin = 0, ymax = 365.25*5)
  } else if (type == "u"){
    cstat <- concordance(mod1, mod2, timewt = "n/G2", ymin = 0, ymax = 365.25*5)
  }
  cstat.diff <- cstat$concordance[["mod2"]]-cstat$concordance[["mod1"]]
  cstat.diff
}

res5YR <- censboot(pheno[, c("pid", "fuday_age_censored", "scd_age_filt", "sex", "race", 
                             "age", "bmi", "height", "smoke", "chd_chs", "hf", 
                             "diabetes", "sbp", "hnt_med", "hdl", "total_chol", 
                             "tri", "stroke", "afib", "sum_scaled_p05")], 
                   func5YR, R = 1000, index = c(3,2), type = "h")
res5YR.ci <- boot.ci(res5YR, type = c("norm", "perc"))

## Uno
uno5YR <- concordance(base.mod, new.mod, timewt = "n/G2", ymin = 0, ymax = 365.25*5)

res5YR.uno <- censboot(pheno[, c("pid", "fuday_age_censored", "scd_age_filt", "sex", "race", 
                             "age", "bmi", "height", "smoke", "chd_chs", "hf", 
                             "diabetes", "sbp", "hnt_med", "hdl", "total_chol", 
                             "tri", "stroke", "afib", "sum_scaled_p05")], 
                   func5YR, R = 1000, index = c(3,2), type = "u")
res5YR.uno.ci <- boot.ci(res5YR.uno, type = c("norm", "perc"))

####10-YR C-stat####

## Harrell
harrell10YR <- concordance(base.mod, new.mod, timewt = "n", ymin = 0, ymax = 365.25*10)

#Get CI for change via bootstrapping
func10YR <- function(data, ind, type){
  mod1 <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + race + age + bmi + height + 
                  smoke + chd_chs + hf + diabetes + sbp +
                  hnt_med + hdl + total_chol + tri + stroke + 
                  afib, data = data[ind,])
  mod2 <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + race + age + bmi + height + 
                  smoke + chd_chs + hf + diabetes + sbp + 
                  hnt_med + hdl + total_chol + tri + stroke + 
                  afib + sum_scaled_p05, data = data[ind,])
  if(type == "h"){
    cstat <- concordance(mod1, mod2, timewt = "n", ymin = 0, ymax = 365.25*10)
  } else if (type == "u"){
    cstat <- concordance(mod1, mod2, timewt = "n/G2", ymin = 0, ymax = 365.25*10)
  }
  cstat.diff <- cstat$concordance[["mod2"]]-cstat$concordance[["mod1"]]
  cstat.diff
}
res10YR <- censboot(pheno[, c("pid", "fuday_age_censored", "scd_age_filt", "sex", "race", 
                              "age", "bmi", "height", "smoke", "chd_chs", "hf", 
                              "diabetes", "sbp", "hnt_med", "hdl", "total_chol", 
                              "tri", "stroke", "afib", "sum_scaled_p05")], 
                    func10YR, R = 1000, index = c(3,2), type = "h")
res10YR.ci <- boot.ci(res10YR, type = c("norm", "perc"))

## Uno
uno10YR <- concordance(base.mod, new.mod, timewt = "n/G2", ymin = 0, ymax = 365.25*10)

res10YR.uno <- censboot(pheno[, c("pid", "fuday_age_censored", "scd_age_filt", "sex", "race", 
                                 "age", "bmi", "height", "smoke", "chd_chs", "hf", 
                                 "diabetes", "sbp", "hnt_med", "hdl", "total_chol", 
                                 "tri", "stroke", "afib", "sum_scaled_p05")], 
                       func10YR, R = 1000, index = c(3,2), type = "u")
res10YR.uno.ci <- boot.ci(res10YR.uno, type = c("norm", "perc"))

####Continuous NRI from nricens####
#Note that NRI only considers a single time point

base.modx <- coxph(Surv(as.numeric(fuday_age_censored), scd_age_filt) ~ sex + race + age + bmi + height +
                     smoke + chd_chs + hf + diabetes + sbp +
                     hnt_med + hdl + total_chol + tri + stroke +
                     afib, data = pheno, x=TRUE)

new.modx <- coxph(Surv(as.numeric(fuday_age_censored), scd_age_filt) ~ sex + race + age + bmi + height +
                    smoke + chd_chs + hf + diabetes + sbp +
                    hnt_med + hdl + total_chol + tri + stroke +
                    afib + sum_scaled_p05, data = pheno, x=TRUE)

## 5YR cfNRI
nri5YR <- nricens(mdl.std = base.modx, mdl.new = new.modx, updown = "diff",
                  cut = 0, t0=365.25*5, niter = 1000)

nri5YR.ipw <- nricens(mdl.std = base.modx, mdl.new = new.modx, updown = "diff",
                  cut = 0, t0=365.25*5, niter = 1000, point.method = "ipw")

## 10YR cfNRI
nri10YR <- nricens(mdl.std = base.modx, mdl.new = new.modx, updown = "diff",
                   cut = 0, t0=365.25*10, niter = 1000)

####Continuous NRI from survIDINRI####

## Numeric only dataframe
pheno.num <- pheno[, c("pid", "fuday_age_censored", "scd_age_filt", "sex", "race", 
                       "age", "bmi", "height", "smoke", "chd_chs", "hf", 
                       "diabetes", "sbp", "hnt_med", "hdl", "total_chol", 
                       "tri", "stroke", "afib", "sum_scaled_p05")]
pheno.num$sex <- as.numeric(pheno.num$sex)
pheno.num$race <- as.numeric(pheno.num$race)
pheno.num$chd_chs <- as.numeric(pheno.num$chd_chs)
pheno.num$hf <- as.numeric(pheno.num$hf)
pheno.num$diabetes <- as.numeric(pheno.num$diabetes)
pheno.num$hnt_med <- as.numeric(pheno.num$hnt_med)
pheno.num$stroke <- as.numeric(pheno.num$stroke)
pheno.num$afib <- as.numeric(pheno.num$afib)
pheno.num$smokeC <- ifelse(pheno.num$smoke == "C", 1, 0)
pheno.num$smoke <- ifelse(pheno.num$smoke == "F", 1, 0)

## 5YR
nri5YR.survIDI <- IDI.INF(indata = pheno.num[,c("fuday_age_censored","scd_age_filt")],
                          covs0 = pheno.num[,c("sex", "race", "age", "bmi", "height", "smoke",
                                               "smokeC", "chd_chs", "hf", "diabetes", "sbp", 
                                               "hnt_med", "hdl", "total_chol", "tri", "stroke", 
                                               "afib")] %>% as.matrix(),
                          covs1 = pheno.num[,c("sex", "race", "age", "bmi", "height", "smoke",
                                               "smokeC", "chd_chs", "hf", "diabetes", "sbp", 
                                               "hnt_med", "hdl", "total_chol", "tri", "stroke", 
                                               "afib", "sum_scaled_p05") %>% as.matrix()], 
                          t0 = 365.25*5, npert = 1000, seed1 = 123, alpha = 0.05)
IDI.INF.OUT(nri5YR.survIDI)
#M2 is NRI

## 10YR
nri10YR.survIDI <- IDI.INF(indata = pheno.num[,c("fuday_age_censored","scd_age_filt")],
                          covs0 = pheno.num[,c("sex", "race", "age", "bmi", "height", "smoke",
                                               "smokeC", "chd_chs", "hf", "diabetes", "sbp", 
                                               "hnt_med", "hdl", "total_chol", "tri", "stroke", 
                                               "afib")] %>% as.matrix(),
                          covs1 = pheno.num[,c("sex", "race", "age", "bmi", "height", "smoke",
                                               "smokeC", "chd_chs", "hf", "diabetes", "sbp", 
                                               "hnt_med", "hdl", "total_chol", "tri", "stroke", 
                                               "afib", "sum_scaled_p05") %>% as.matrix()], 
                          t0 = 365.25*10, npert = 1000, seed1 = 123, alpha = 0.05)
IDI.INF.OUT(nri10YR.survIDI)
#M2 is NRI

## USE RESULTS FROM PACKAGE: survIDINRI (widely used, previously used in lab)

####C-stats results table####

## Harrell
harrell.results <- data.frame(model = c("base", "risk"), type = "harrell", outcome = "Max", 
                              cstat = coef(harrell) %>% as.data.frame(), 
                              se = c(harrell$var %>% sqrt() %>% .[1,1], harrell$var %>% sqrt() %>% .[2,2]),
                              change = coef(harrell)[2] - coef(harrell)[1],
                              change.norm.lci = res.ci$normal[2],
                              change.norm.uci = res.ci$normal[3],
                              change.perc.lci = res.ci$percent[4],
                              change.perc.uci = res.ci$percent[5]) %>% setnames(., ".", "cstat")

harrell.5YR.results <- data.frame(model = c("base", "risk"), type = "harrell", outcome = "5YR", 
                                  cstat = coef(harrell5YR) %>% as.data.frame(), 
                                  se = c(harrell5YR$var %>% sqrt() %>% .[1,1], harrell5YR$var %>% sqrt() %>% .[2,2]),
                                  change = coef(harrell5YR)[2] - coef(harrell5YR)[1],
                                  change.norm.lci = res5YR.ci$normal[2],
                                  change.norm.uci = res5YR.ci$normal[3],
                                  change.perc.lci = res5YR.ci$percent[4],
                                  change.perc.uci = res5YR.ci$percent[5]) %>% setnames(., ".", "cstat")

harrell.10YR.results <- data.frame(model = c("base", "risk"), type = "harrell", outcome = "10YR", 
                                   cstat = coef(harrell10YR) %>% as.data.frame(), 
                                   se = c(harrell10YR$var %>% sqrt() %>% .[1,1], harrell10YR$var %>% sqrt() %>% .[2,2]),
                                   change = coef(harrell10YR)[2] - coef(harrell10YR)[1],
                                   change.norm.lci = res10YR.ci$normal[2],
                                   change.norm.uci = res10YR.ci$normal[3],
                                   change.perc.lci = res10YR.ci$percent[4],
                                   change.perc.uci = res10YR.ci$percent[5]) %>% setnames(., ".", "cstat")

## Uno
uno.results <- data.frame(model = c("base", "risk"), type = "uno", outcome = "Max", 
                          cstat = coef(uno) %>% as.data.frame(), 
                          se = c(uno$var %>% sqrt() %>% .[1,1], uno$var %>% sqrt() %>% .[2,2]),
                          change = coef(uno)[2] - coef(uno)[1],
                          change.norm.lci = res.uno.ci$normal[2],
                          change.norm.uci = res.uno.ci$normal[3],
                          change.perc.lci = res.uno.ci$percent[4],
                          change.perc.uci = res.uno.ci$percent[5]) %>% setnames(., ".", "cstat")

uno.5YR.results <- data.frame(model = c("base", "risk"), type = "uno", outcome = "5YR", 
                              cstat = coef(uno5YR) %>% as.data.frame(), 
                              se = c(uno5YR$var %>% sqrt() %>% .[1,1], uno5YR$var %>% sqrt() %>% .[2,2]),
                              change = coef(uno5YR)[2] - coef(uno5YR)[1],
                              change.norm.lci = res5YR.uno.ci$normal[2],
                              change.norm.uci = res5YR.uno.ci$normal[3],
                              change.perc.lci = res5YR.uno.ci$percent[4],
                              change.perc.uci = res5YR.uno.ci$percent[5]) %>% setnames(., ".", "cstat")

uno.10YR.results <- data.frame(model = c("base", "risk"), type = "uno", outcome = "10YR", 
                               cstat = coef(uno10YR) %>% as.data.frame(), 
                               se = c(uno10YR$var %>% sqrt() %>% .[1,1], uno10YR$var %>% sqrt() %>% .[2,2]),
                               change = coef(uno10YR)[2] - coef(uno10YR)[1],
                               change.norm.lci = res10YR.uno.ci$normal[2],
                               change.norm.uci = res10YR.uno.ci$normal[3],
                               change.perc.lci = res10YR.uno.ci$percent[4],
                               change.perc.uci = res10YR.uno.ci$percent[5]) %>% setnames(., ".", "cstat")

cstats.results <- rbind(harrell.results, harrell.5YR.results, harrell.10YR.results, uno.results, uno.5YR.results, uno.10YR.results)
cstats.results$upper <- cstats.results$cstat+(1.96*cstats.results$se)
cstats.results$lower <- cstats.results$cstat-(1.96*cstats.results$se)
cstats.results <- cstats.results[, c("model", "type", "outcome", "cstat", "se", "lower", "upper", "change", "change.norm.lci", "change.norm.uci", "change.perc.lci", "change.perc.uci")]

####NRI results table####

nri5YR.results <- rbind(nri5YR.survIDI$m1, nri5YR.survIDI$m2, nri5YR.survIDI$m3) %>% as.data.frame() %>% mutate(outcome = "5YR") %>% 
  setNames(., c("Est.", "Lower", "Upper", "P", "Outcome"))
nri10YR.results <- rbind(nri10YR.survIDI$m1, nri10YR.survIDI$m2, nri10YR.survIDI$m3) %>% as.data.frame() %>% mutate(outcome = "10YR") %>%
  setNames(., c("Est.", "Lower", "Upper", "P", "Outcome"))

nri.results <- rbind(nri5YR.results, nri10YR.results)
nri.results <- nri.results[,c("Outcome", "Est.", "Lower", "Upper", "P")]
nri.results.print <- nri.results[, c("Est.", "Lower", "Upper", "P", "Outcome")]
nri.results.print[,1:3] <- signif(nri.results.print[,1:3], digits = 3)

####Save data####
write.csv(cstats.results, "ARIC.Prot.Risk.Score.Cindex_Results.csv", row.names = F, col.names = T, quote = F, sep = "\t")
write.csv(nri.results.print, "ARIC.Prot.Risk.Score.NRI_Results.csv", row.names = F, col.names = T, quote = F, sep = "\t")
save.image("ARIC.Prot.Risk.Score.C-stat.NRI.Analysis.RData")
