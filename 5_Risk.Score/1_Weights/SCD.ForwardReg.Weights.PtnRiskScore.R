####Description####
## Run forwards regression on ARIC-only significant hits under Min (Primary) and Dz+RskFct (Full) models to decrease proteins used for PtRS

####Load packages####
library(magrittr)
library(survival)
library(stringr)
# library(randomForestSRC)
# library(corrplot)

####Load data####
load("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/ThuyVy.SCA.Assoc.ARIC.V2Prot_BthRcSx.RData")
allresults <- readRDS("aric.prot.real.archive/scd.assoc/9_other/combine.data/All.SCAxProteomics.Results_10072022.RDS")
annotation <- subset(allresults, cohort == "aric" & model == "M3" & strat == "all") %>% unique()

####Get ARIC Min, Dz+RskFct Bonf proteins####
aric.min.bonf.seqid <- subset(allresults, cohort == "aric" & model == "M2" & P<(0.05/4955) & strat == "all")$SeqId.use
aric.rskfct.bonf <- subset(allresults, cohort == "aric" & model == "M3" & P<(0.05/4955) & strat == "all" & SeqId.use %in% aric.min.bonf.seqid) %>% .[order(.$P),]
aric.rskfct.bonf.seqid <- subset(allresults, cohort == "aric" & model == "M3" & P<(0.05/4955) & strat == "all" & SeqId.use %in% aric.min.bonf.seqid) %>% .[order(.$P),"SeqId.use"]

## Convert to ARIC format of SeqId
aric.rskfct.bonf.seqid <- gsub("-", "_", aric.rskfct.bonf.seqid) %>% paste0("SeqId_", .)

## Make dataset with just the proteins
columns <- which(names(final.alldata) %in% aric.rskfct.bonf.seqid)
row.names(final.alldata) <- final.alldata$pid
probes <- final.alldata[, columns]

## Check how many of the ARIC set overlaps with the meta set
meta.min.bonf.seqid <- subset(allresults, cohort == "meta" & model == "M2" & P<(0.05/4955) & strat == "all")$SeqId.use
meta.rskfct.bonf.seqid <- subset(allresults, cohort == "meta" & model == "M3" & P<(0.05/4955) & strat == "all" & SeqId.use %in% meta.min.bonf.seqid)$SeqId.use %>% gsub("-", "_", .) %>% paste0("SeqId", "_", .)
table(meta.rskfct.bonf.seqid %in% aric.rskfct.bonf.seqid)
# FALSE  TRUE
# 7    17
aric.meta.overlap <- meta.rskfct.bonf.seqid %>% gsub("SeqId_", "", .) %>% gsub("_", "-",.)
aric.meta.overlap <- subset(annotation, SeqId.use %in% aric.meta.overlap)

####Subset for model####
rskfct <- c("sex", "age", "eGFR_ckdepi", "smoke", "bmi", "height", 
            "race", "center", "hf", "chd_chs", "diabetes", "sbp", 
            "hnt_med", "hdl", "total_chol", "tri", "stroke", "afib")
scd <- c("scd_age_filt", "fuday_age_censored")
rskfct.columns <- which(names(final.alldata) %in% rskfct)
scd.columns <- which(names(final.alldata) %in% scd)
pheno <- final.alldata[, c(scd.columns, rskfct.columns)]
#pheno and probes

proteins <- merge(pheno, probes, by = "row.names")
names(proteins)[1] <- "pid"
row.names(proteins) <- proteins$pid
proteins$pid <- NULL
save.image("Base.RData")

####Run forward regression####

## Functions

#Gets protein P values from model
get.model.P <- function(coxph.mod, protein){
  model.summary <- coef(summary(coxph.mod)) %>% as.data.frame()
  model.summary$term <- row.names(model.summary)
  terms.keep <- grep("SeqId", model.summary$term)
  protein.terms <- model.summary[terms.keep,] %>% setNames(., c("beta", "hr", "se", "z", "p", "term"))
  P <- subset(protein.terms, term == protein)$p
  return(P)
}

#Gets proteins from model
get.model.proteins <- function(coxph.mod){
  model.summary <- coef(summary(coxph.mod)) %>% as.data.frame()
  model.summary$term <- row.names(model.summary)
  terms.keep <- grep("SeqId", model.summary$term)
  protein.terms <- model.summary[terms.keep,] %>% setNames(., c("beta", "hr", "se", "z", "p", "term"))
  proteins <- protein.terms$term
  return(proteins)
}

#Outputs protein with most significant P, passing alpha
get.top.protein <- function(basemodel, seqids, df, alpha){
  models <- list()
  significance <- as.data.frame(matrix(nrow = length(seqids), ncol = 2))
  names(significance) <- c("SeqId", "P")
  significance$SeqId <- seqids
  #newdata <- merge(pheno, proteins[, seqids], by = "row.names")
  for(i in 1:length(seqids)){
    current.protein <- seqids[i]
    #print(current.protein)
    x <- paste0(as.formula(basemodel), " + ", current.protein)[3]
    y <- "Surv(fuday_age_censored, scd_age_filt)"
    update <- as.formula(paste0(y, "~", x))
    if(i == 1){
      print("Model to determine next protein to add:")
      print(update)
     }
    model <- coxph(update, data = df)
    models[[i]] <- model
    models[[i]] <- setNames(models[i], current.protein)
    significance$P[which(significance$SeqId %in% current.protein)] <- get.model.P(model, current.protein)
  }
  models[[length(seqids)+1]] <- significance
  significance <- significance[order(significance$P),]
  if(significance$P[1] < alpha){
  protein.add <- significance$SeqId[1]
  return(protein.add)
  } else {
    print("Protein > alpha")
  }
  #return(models)
}

#Actual function that does the selection
runFW <- function(base, seqid, df, alpha){
  first <- get.top.protein(base, seqid, df, alpha)
  new.model.x <- paste0(as.formula(base), " + ", first)[3]
  y <- "Surv(fuday_age_censored, scd_age_filt)"
  base <- as.formula(paste0(y, "~", new.model.x))
  seqid <- seqid[-which(seqid %in% first)]
  print("Forward selection model:")
  print(base)
  update <- coxph(base, data = df)
  return(update)
}

#Remove seqid from full list
update.seqid <- function(full,tbr){
  remain <- full[-which(full %in% tbr)]
  return(remain)
}

## Run actual regression
base.listed.model <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi + smoke + bmi +
                             height + center + hf + chd_chs + diabetes + sbp + 
                             hnt_med + hdl + total_chol + tri + stroke + afib, data = proteins)
basemod <- base.listed.model
seqid_use <- aric.rskfct.bonf.seqid
y <- "Surv(fuday_age_censored, scd_age_filt)"

models_keep <- list()
for(i in 1:length(seqid_use)){
  first <- runFW(basemod, seqid_use, proteins, 0.05)
  added <- get.model.proteins(first)
  seqid_use <- update.seqid(seqid_use, added)
  models_keep[[i]] <- first
  basemod <- first
}
#models_keep[[8]] was last one model to run before erroring out, meaning no additional proteins have P < 0.05

## Run manual starting from new base model with 8 proteins added

#Gets model object returned
get.top.model <- function(basemodel, seqids, df){
  models <- list()
  significance <- as.data.frame(matrix(nrow = length(seqids), ncol = 2))
  names(significance) <- c("SeqId", "P")
  significance$SeqId <- seqids
  #newdata <- merge(pheno, proteins[, seqids], by = "row.names")
  for(i in 1:length(seqids)){
    current.protein <- seqids[i]
    #print(current.protein)
    x <- paste0(as.formula(basemodel), " + ", current.protein)[3]
    y <- "Surv(fuday_age_censored, scd_age_filt)"
    update <- as.formula(paste0(y, "~", x))
    if(i == 1){
      print("Model to determine next protein to add:")
      print(update)
    }
    model <- coxph(update, data = df)
    models[[i]] <- model
    models[[i]] <- setNames(models[i], current.protein)
    significance$P[which(significance$SeqId %in% current.protein)] <- get.model.P(model, current.protein)
  }
  significance <- significance[order(significance$P),]
  models[[length(seqids)+1]] <- significance
  # if(significance$P[1] < alpha){
  #   protein.add <- significance$SeqId[1]
  #   return(protein.add)
  # } else {
  #   print("Protein > alpha")
  # }
  return(models)
}

base2.auto <- models_keep[[8]]
base2.manual.mod <- 
  coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi +
          smoke + bmi + height + center + hf + chd_chs + diabetes +
          sbp + hnt_med + hdl + total_chol + tri + stroke + afib +
          SeqId_7655_11 + SeqId_4496_60 + SeqId_11089_7 + SeqId_11388_75 +
          SeqId_6586_19 + SeqId_8969_49 + SeqId_4124_24 + SeqId_16785_45, data = proteins)

seqid.tbr <- c("SeqId_7655_11", "SeqId_4496_60", "SeqId_11089_7", "SeqId_11388_75",
                   "SeqId_6586_19", "SeqId_8969_49", "SeqId_4124_24", "SeqId_16785_45")
seqid.base2 <- aric.rskfct.bonf.seqid[-which(aric.rskfct.bonf.seqid %in% seqid.tbr)]
base2.next <- get.top.model(base2.auto, seqid.base2, proteins)
#Manually checked all models as well - no additionally added protein is P <0.05

base2.proteins <- proteins[, seqid.tbr]
M <- cor(base2.proteins)
corrplot::corrplot(M, method = "number")
#Not much correlation, which is to be expected since we did this to remove collinearity
#Highest cor is 0.50: MMP12 and HE4

####Do all manual calculation to double check####
step1 <- get.top.model(base.listed.model, aric.rskfct.bonf.seqid, proteins)
base1 <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi +
                   smoke + bmi + height + center + hf + chd_chs + diabetes +
                   sbp + hnt_med + hdl + total_chol + tri + stroke + afib +
                   SeqId_7655_11, data = proteins)
add1 <- "SeqId_7655_11"
seqid1 <- aric.rskfct.bonf.seqid[-which(aric.rskfct.bonf.seqid %in% add1)]

step2 <- get.top.model(base1, seqid1, proteins)
base2 <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi +
                 smoke + bmi + height + center + hf + chd_chs + diabetes +
                 sbp + hnt_med + hdl + total_chol + tri + stroke + afib +
                 SeqId_7655_11 + SeqId_4496_60, data = proteins)
add2 <- c(add1, "SeqId_4496_60")
seqid2 <- aric.rskfct.bonf.seqid[-which(aric.rskfct.bonf.seqid %in% add2)]

step3 <- get.top.model(base2, seqid2, proteins)
base3 <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi +
                 smoke + bmi + height + center + hf + chd_chs + diabetes +
                 sbp + hnt_med + hdl + total_chol + tri + stroke + afib +
                 SeqId_7655_11 + SeqId_4496_60 + SeqId_11089_7, data = proteins)
add3 <- c(add2, "SeqId_11089_7")
seqid3 <- aric.rskfct.bonf.seqid[-which(aric.rskfct.bonf.seqid %in% add3)]

step4 <- get.top.model(base3, seqid3, proteins)
base4 <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi +
                 smoke + bmi + height + center + hf + chd_chs + diabetes +
                 sbp + hnt_med + hdl + total_chol + tri + stroke + afib +
                 SeqId_7655_11 + SeqId_4496_60 + SeqId_11089_7 + SeqId_11388_75, data = proteins)
add4 <- c(add3, "SeqId_11388_75")
seqid4 <- aric.rskfct.bonf.seqid[-which(aric.rskfct.bonf.seqid %in% add4)]

step5 <- get.top.model(base4, seqid4, proteins)
base5 <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi +
                 smoke + bmi + height + center + hf + chd_chs + diabetes +
                 sbp + hnt_med + hdl + total_chol + tri + stroke + afib +
                 SeqId_7655_11 + SeqId_4496_60 + SeqId_11089_7 + SeqId_11388_75 + SeqId_6586_19, data = proteins)
add5 <- c(add4, "SeqId_6586_19")
seqid5 <- aric.rskfct.bonf.seqid[-which(aric.rskfct.bonf.seqid %in% add5)]

step6 <- get.top.model(base5, seqid5, proteins)
base6 <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi +
                 smoke + bmi + height + center + hf + chd_chs + diabetes +
                 sbp + hnt_med + hdl + total_chol + tri + stroke + afib +
                 SeqId_7655_11 + SeqId_4496_60 + SeqId_11089_7 + SeqId_11388_75 + SeqId_6586_19 +
                 SeqId_8969_49, data = proteins)
add6 <- c(add5, "SeqId_8969_49")
seqid6 <- aric.rskfct.bonf.seqid[-which(aric.rskfct.bonf.seqid %in% add6)]

step7 <- get.top.model(base6, seqid6, proteins)
base7 <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi +
                 smoke + bmi + height + center + hf + chd_chs + diabetes +
                 sbp + hnt_med + hdl + total_chol + tri + stroke + afib +
                 SeqId_7655_11 + SeqId_4496_60 + SeqId_11089_7 + SeqId_11388_75 + SeqId_6586_19 +
                 SeqId_8969_49 + SeqId_4124_24, data = proteins)
add7 <- c(add6, "SeqId_4124_24")
seqid7 <- aric.rskfct.bonf.seqid[-which(aric.rskfct.bonf.seqid %in% add7)]

step8 <- get.top.model(base7, seqid7, proteins)
base8 <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi +
                 smoke + bmi + height + center + hf + chd_chs + diabetes +
                 sbp + hnt_med + hdl + total_chol + tri + stroke + afib +
                 SeqId_7655_11 + SeqId_4496_60 + SeqId_11089_7 + SeqId_11388_75 + SeqId_6586_19 +
                 SeqId_8969_49 + SeqId_4124_24 + SeqId_16785_45, data = proteins)
add8 <- c(add7, "SeqId_16785_45")
seqid8 <- aric.rskfct.bonf.seqid[-which(aric.rskfct.bonf.seqid %in% add8)]

step9 <- get.top.model(base8, seqid8, proteins)
#None with P < 0.05
#final model is base8

identical(base8, base2.manual.mod)
#TRUE
#Manual and automated gives same answer
p05.ManFW.model <- base2.manual.mod

####Run forward regression - Bonf cutoff####
base.listed.model <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi + smoke + bmi +
                             height + center + hf + chd_chs + diabetes + sbp + 
                             hnt_med + hdl + total_chol + tri + stroke + afib, data = proteins)
basemod <- base.listed.model
seqid_use <- aric.rskfct.bonf.seqid
y <- "Surv(fuday_age_censored, scd_age_filt)"

models_keep_pBonf <- list()
for(i in 1:length(seqid_use)){
  first <- runFW(basemod, seqid_use, proteins, 0.05/33)
  added <- get.model.proteins(first)
  seqid_use <- update.seqid(seqid_use, added)
  models_keep_pBonf[[i]] <- first
  basemod <- first
}
#stopped at: models_keep_pBonf[[6]]
base2.auto.pBonf <- models_keep_pBonf[[6]]
base2.auto.pBonf.man <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi +
                                smoke + bmi + height + center + hf + chd_chs + diabetes +
                                sbp + hnt_med + hdl + total_chol + tri + stroke + afib +
                                SeqId_7655_11 + SeqId_4496_60 + SeqId_11089_7 + SeqId_11388_75 +
                                SeqId_6586_19 + SeqId_8969_49, data  = proteins)

####Conclusion####
## These are the same proteins as determined by backwards regression, so use the weights from there

save.image("aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/forward/SCD.ForwardReg.Weights.PtnRiskScore.RData")
