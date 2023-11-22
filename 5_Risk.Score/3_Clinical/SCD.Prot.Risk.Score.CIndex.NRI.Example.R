####Description####
## Calculates Harrell's and Uno's C-statistic, continuous NRI for SCD clinical model and model with protein risk score
## Script sent to CHS for analysis

####Load packages####
library(survival)
library(boot)
library(survIDINRI)
library(magrittr)
library(data.table)
library(dplyr)

####Load functions####
#CI for C-stat change, with bootstrapping
#IMPORTANT: variable names in models need to be changed to match dataset----

func <- function(data, ind, type){
  mod1 <- coxph(Surv(futime, scd) ~ sex + race + age + bmi + height + 
                  smoke + chd + hf + diabetes + sbp + 
                  hnt_med + hdl + total_chol + tri + stroke + 
                  afib, data = data[ind,])
  mod2 <- coxph(Surv(futime, scd) ~ sex + race + age + bmi + height + 
                  smoke + chd + hf + diabetes + sbp + 
                  hnt_med + hdl + total_chol + tri + stroke + 
                  afib + risk_score, data = data[ind,])
  if(type == "h"){
    cstat <- concordance(mod1, mod2, timewt = "n")
  } else if(type == "u"){
    cstat <- concordance(mod1, mod2, timewt = "n/G2")
  }
  cstat.diff <- cstat$concordance[["mod2"]]-cstat$concordance[["mod1"]]
  cstat.diff
}

func5YR <- function(data, ind, type){
  mod1 <- coxph(Surv(futime, scd) ~ sex + race + age + bmi + height + 
                  smoke + chd + hf + diabetes + sbp + 
                  hnt_med + hdl + total_chol + tri + stroke + 
                  afib, data = data[ind,])
  mod2 <- coxph(Surv(futime, scd) ~ sex + race + age + bmi + height + 
                  smoke + chd + hf + diabetes + sbp + 
                  hnt_med + hdl + total_chol + tri + stroke + 
                  afib + risk_score, data = data[ind,])
  if(type == "h"){
    cstat <- concordance(mod1, mod2, timewt = "n", ymin = 0, ymax = 365.25*5)
  } else if (type == "u"){
    cstat <- concordance(mod1, mod2, timewt = "n/G2", ymin = 0, ymax = 365.25*5)
  }
  cstat.diff <- cstat$concordance[["mod2"]]-cstat$concordance[["mod1"]]
  cstat.diff
}

func10YR <- function(data, ind, type){
  mod1 <- coxph(Surv(futime, scd) ~ sex + race + age + bmi + height + 
                  smoke + chd + hf + diabetes + sbp +
                  hnt_med + hdl + total_chol + tri + stroke + 
                  afib, data = data[ind,])
  mod2 <- coxph(Surv(futime, scd) ~ sex + race + age + bmi + height + 
                  smoke + chd + hf + diabetes + sbp + 
                  hnt_med + hdl + total_chol + tri + stroke + 
                  afib + risk_score, data = data[ind,])
  if(type == "h"){
    cstat <- concordance(mod1, mod2, timewt = "n", ymin = 0, ymax = 365.25*10)
  } else if (type == "u"){
    cstat <- concordance(mod1, mod2, timewt = "n/G2", ymin = 0, ymax = 365.25*10)
  }
  cstat.diff <- cstat$concordance[["mod2"]]-cstat$concordance[["mod1"]]
  cstat.diff
}

####Load data####
#pheno: every row is individual, columns are pheno, proteins, etc.

## Set seed for calculations
set.seed(350)

####Define models####
#Base model contains: sex, race, age, bmi, height, smoking status, CHD, HF, diabetes, sbp, hypertension medication, HDL, total cholesterol, triglycerides, stroke, AF
#Risk score model is base + risk score (with NPPB)
#IMPORTANT: variable names in models need to be changed to match dataset----

pheno <- alldata[, c("ID", "futime", "scd", "sex", "race", 
                     "age", "bmi", "height", "smoke", "chd", "hf", 
                     "diabetes", "sbp", "hnt_med", "hdl", "total_chol", 
                     "tri", "stroke", "afib", "risk_score")]

base.mod <- coxph(Surv(futime, scd) ~ sex + race + age + bmi + height + 
                    smoke + chd + hf + diabetes + sbp + 
                    hnt_med + hdl + total_chol + tri + stroke + 
                    afib, data = pheno)

new.mod <- coxph(Surv(futime, scd) ~ sex + race + age + bmi + height + 
                   smoke + chd + hf + diabetes + sbp + 
                   hnt_med + hdl + total_chol + tri + stroke + 
                   afib + risk_score, data = pheno)

####Harrell C-stat####

## Max outcome
harrell <- concordance(base.mod, new.mod, timewt = "n")
res <- censboot(pheno, func, R = 1000, index = c(3,2), type = "h")
res.ci <- boot.ci(res, type = c("norm", "perc"))

## 5YR outcome
harrell5YR <- concordance(base.mod, new.mod, timewt = "n", ymin = 0, ymax = 365.25*5)
res5YR <- censboot(pheno, func5YR, R = 1000, index = c(3,2), type = "h")
res5YR.ci <- boot.ci(res5YR, type = c("norm", "perc"))

## 10YR outcome
harrell10YR <- concordance(base.mod, new.mod, timewt = "n", ymin = 0, ymax = 365.25*10)
res10YR <- censboot(pheno, func10YR, R = 1000, index = c(3,2), type = "h")
res10YR.ci <- boot.ci(res10YR, type = c("norm", "perc"))

## Make results table
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

####Uno C-stat####

## Max outcome
uno <- concordance(base.mod, new.mod, timewt = "n/G2")
res.uno <- censboot(pheno, func, R = 1000, index = c(3,2), type = "u")
res.uno.ci <- boot.ci(res.uno, type = c("norm", "perc"))

## 5YR outcome
uno5YR <- concordance(base.mod, new.mod, timewt = "n/G2", ymin = 0, ymax = 365.25*5)
res5YR.uno <- censboot(pheno, func5YR, R = 1000, index = c(3,2), type = "u")
res5YR.uno.ci <- boot.ci(res5YR.uno, type = c("norm", "perc"))

## 10YR outcome
uno10YR <- concordance(base.mod, new.mod, timewt = "n/G2", ymin = 0, ymax = 365.25*10)
res10YR.uno <- censboot(pheno, func10YR, R = 1000, index = c(3,2), type = "u")
res10YR.uno.ci <- boot.ci(res10YR.uno, type = c("norm", "perc"))

## Make results table
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

## Make final table
cstats.results <- rbind(harrell.results, harrell.5YR.results, harrell.10YR.results, uno.results, uno.5YR.results, uno.10YR.results)

####Continuous NRI####

## Add dummy variables, change factors and categorical to numeric
#IMPORTANT: variable names in models need to be changed to match dataset----
pheno.num <- pheno
pheno.num$sex <- as.numeric(pheno.num$sex)
pheno.num$race <- as.numeric(pheno.num$race)
pheno.num$chd <- as.numeric(pheno.num$chd)
pheno.num$hf <- as.numeric(pheno.num$hf)
pheno.num$diabetes <- as.numeric(pheno.num$diabetes)
pheno.num$hnt_med <- as.numeric(pheno.num$hnt_med)
pheno.num$stroke <- as.numeric(pheno.num$stroke)
pheno.num$afib <- as.numeric(pheno.num$afib)
pheno.num$smokeC <- ifelse(pheno.num$smoke == "C", 1, 0)
pheno.num$smoke <- ifelse(pheno.num$smoke == "F", 1, 0)

## 5YR
#IMPORTANT: variable names in models need to be changed to match dataset---
nri5YR <- IDI.INF(indata = pheno.num[,c("futime","scd")],
                          covs0 = pheno.num[,c("sex", "race", "age", "bmi", "height", "smoke",
                                               "smokeC", "chd", "hf", "diabetes", "sbp", 
                                               "hnt_med", "hdl", "total_chol", "tri", "stroke", 
                                               "afib")] %>% as.matrix(),
                          covs1 = pheno.num[,c("sex", "race", "age", "bmi", "height", "smoke",
                                               "smokeC", "chd", "hf", "diabetes", "sbp", 
                                               "hnt_med", "hdl", "total_chol", "tri", "stroke", 
                                               "afib", "risk_score") %>% as.matrix()], 
                          t0 = 365.25*5, npert = 1000, seed1 = 123, alpha = 0.05)

nri5YR.results <- rbind(nri5YR$m1, nri5YR$m2, nri5YR$m3) %>% as.data.frame() %>% mutate(outcome = "5YR") %>% 
  setNames(., c("Est", "Lower", "Upper", "P", "Outcome"))

## 10YR
#IMPORTANT: variable names in models need to be changed to match dataset----
nri10YR <- IDI.INF(indata = pheno.num[,c("futime","scd")],
                           covs0 = pheno.num[,c("sex", "race", "age", "bmi", "height", "smoke",
                                                "smokeC", "chd", "hf", "diabetes", "sbp", 
                                                "hnt_med", "hdl", "total_chol", "tri", "stroke", 
                                                "afib")] %>% as.matrix(),
                           covs1 = pheno.num[,c("sex", "race", "age", "bmi", "height", "smoke",
                                                "smokeC", "chd", "hf", "diabetes", "sbp", 
                                                "hnt_med", "hdl", "total_chol", "tri", "stroke", 
                                                "afib", "risk_score") %>% as.matrix()], 
                           t0 = 365.25*10, npert = 1000, seed1 = 123, alpha = 0.05)

nri10YR.results <- rbind(nri10YR$m1, nri10YR$m2, nri10YR$m3) %>% as.data.frame() %>% mutate(outcome = "10YR") %>% 
  setNames(., c("Est", "Lower", "Upper", "P", "Outcome"))

nri.results <- rbind(nri5YR.results, nri10YR.results)
nri.results <- nri.results[,c("Outcome", "Est", "Lower", "Upper", "P")]

####Write out final results####
#Fill in date
date <- ""
write.csv(cstats.results, file = paste0("CHS.SCD.Prot.Cindex.Results_", date, ".csv"), row.names = F, quote = F)
write.csv(nri.results, file = paste0("CHS.SCD.Prot.NRI.Results_", date, ".csv"), row.names = F, quote = F)
