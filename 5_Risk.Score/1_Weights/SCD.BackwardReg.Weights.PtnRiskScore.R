####Description####
## Run backwards regression on ARIC-only significant hits under Min (Primary) and Dz+RskFct (Full) models to decrease proteins used for PtRS

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

####Manual backwards regression####
full.og <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ ., data = proteins)
full <- full.og

models <- list()
for(i in 1:33){
  print(paste0("Running model number ", i))
  model.summary <- coef(summary(full)) %>% as.data.frame()
  model.summary$term <- row.names(model.summary)
  terms.keep <- grep("SeqId", model.summary$term)
  protein.terms <- model.summary[terms.keep,] %>% setNames(., c("beta", "hr", "se", "z", "p", "term"))
  remove <- protein.terms[order(protein.terms$p),] %>% row.names(.) %>% .[length(.)]
  if(i == 1){
    col.tbr <- which(names(proteins) %in% remove)
    prot.data.temp <- proteins[, -col.tbr]
    print(paste0("Number of columns in data: ", ncol(prot.data.temp)))
  } else {
    col.tbr <- which(names(prot.data) %in% remove)
    prot.data.temp <- prot.data[, -col.tbr]
    print(paste0("Number of columns in data: ", ncol(prot.data.temp)))
  }
  updated.model <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ ., data = prot.data.temp)
  old.model <- full
  models[[i]] <- updated.model
  full <- updated.model
  prot.data <- prot.data.temp
}

#Model #25 is using p<0.05

#Bonferroni cut-off: 0.05/33 (0.001515152)
#Model #27 is using (0.05/33) cut-off

p05.manbw.model <- models[[25]]
bonf.manbw.model <- models[[27]]

####Automated backwards regression, based on AIC####
full.listed.model <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi + smoke + bmi +
                             height + center + hf + chd_chs + diabetes + sbp +
                             hnt_med + hdl + total_chol + tri + stroke + afib + 
                             SeqId_10866_60 + SeqId_11089_7 + SeqId_11109_56 + 
                             SeqId_11178_21 + SeqId_11212_7 + SeqId_11388_75 + 
                             SeqId_13133_73 + SeqId_15533_97 + SeqId_16322_10 + 
                             SeqId_16785_45 + SeqId_17456_53 + SeqId_18242_8 + 
                             SeqId_18871_24 + SeqId_19183_164 + SeqId_2602_2 + 
                             SeqId_2677_1 + SeqId_2944_66 + SeqId_2948_58 + 
                             SeqId_3457_57 + SeqId_3485_28 + SeqId_3710_49 + 
                             SeqId_4124_24 + SeqId_4374_45 + SeqId_4496_60 + 
                             SeqId_5008_51 + SeqId_6077_63 + SeqId_6366_38 + 
                             SeqId_6586_19 + SeqId_7551_33 + SeqId_7655_11 + 
                             SeqId_8304_50 + SeqId_8480_29 + SeqId_8969_49, data = proteins)

autobw.step <- step(full.listed.model,
                    scope = list(upper = ~ sex + age + race + eGFR_ckdepi + smoke + bmi +
                                   height + center + hf + chd_chs + diabetes + sbp +
                                   hnt_med + hdl + total_chol + tri + stroke + afib + 
                                   SeqId_10866_60 + SeqId_11089_7 + SeqId_11109_56 + 
                                   SeqId_11178_21 + SeqId_11212_7 + SeqId_11388_75 + 
                                   SeqId_13133_73 + SeqId_15533_97 + SeqId_16322_10 + 
                                   SeqId_16785_45 + SeqId_17456_53 + SeqId_18242_8 + 
                                   SeqId_18871_24 + SeqId_19183_164 + SeqId_2602_2 + 
                                   SeqId_2677_1 + SeqId_2944_66 + SeqId_2948_58 + 
                                   SeqId_3457_57 + SeqId_3485_28 + SeqId_3710_49 + 
                                   SeqId_4124_24 + SeqId_4374_45 + SeqId_4496_60 + 
                                   SeqId_5008_51 + SeqId_6077_63 + SeqId_6366_38 + 
                                   SeqId_6586_19 + SeqId_7551_33 + SeqId_7655_11 + 
                                   SeqId_8304_50 + SeqId_8480_29 + SeqId_8969_49,
                                 lower = ~ sex + age + race + eGFR_ckdepi + smoke + bmi +
                                   height + center + hf + chd_chs + diabetes + sbp +
                                   hnt_med + hdl + total_chol + tri + stroke + afib), 
                       direction = "backward")

final.autobw.model <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi + smoke + bmi +
                              height + center + hf + chd_chs + diabetes + sbp +
                              hnt_med + hdl + total_chol + tri + stroke + afib + 
                              SeqId_11089_7 + SeqId_11178_21 + SeqId_11212_7 +
                              SeqId_11388_75 + SeqId_16785_45 + SeqId_18242_8 + 
                              SeqId_19183_164 + SeqId_2944_66 + SeqId_4124_24 + 
                              SeqId_4496_60 + SeqId_5008_51 + SeqId_6077_63 + 
                              SeqId_6586_19 + SeqId_7655_11 + SeqId_8969_49,
                            data = proteins)
#15 proteins

####Automated forwards regression, based on AIC####
base.listed.model <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi + smoke + bmi +
                             height + center + hf + chd_chs + diabetes + sbp + 
                             hnt_med + hdl + total_chol + tri + stroke + afib, data = proteins)
autofw.step <- step(base.listed.model, 
                     scope = list(upper = ~ sex + age + race + eGFR_ckdepi + smoke + bmi +
                                    height + center + hf + chd_chs + diabetes + sbp +
                                    hnt_med + hdl + total_chol + tri + stroke + afib + 
                                    SeqId_10866_60 + SeqId_11089_7 + SeqId_11109_56 + 
                                    SeqId_11178_21 + SeqId_11212_7 + SeqId_11388_75 + 
                                    SeqId_13133_73 + SeqId_15533_97 + SeqId_16322_10 + 
                                    SeqId_16785_45 + SeqId_17456_53 + SeqId_18242_8 + 
                                    SeqId_18871_24 + SeqId_19183_164 + SeqId_2602_2 + 
                                    SeqId_2677_1 + SeqId_2944_66 + SeqId_2948_58 + 
                                    SeqId_3457_57 + SeqId_3485_28 + SeqId_3710_49 + 
                                    SeqId_4124_24 + SeqId_4374_45 + SeqId_4496_60 + 
                                    SeqId_5008_51 + SeqId_6077_63 + SeqId_6366_38 + 
                                    SeqId_6586_19 + SeqId_7551_33 + SeqId_7655_11 + 
                                    SeqId_8304_50 + SeqId_8480_29 + SeqId_8969_49, 
                                  lower = ~ sex + age + race + eGFR_ckdepi + smoke + bmi +
                                    height + center + hf + chd_chs + diabetes + sbp +
                                    hnt_med + hdl + total_chol + tri + stroke + afib), 
                     direction = "forward")

final.autofw.model <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi + smoke + bmi +
                              height + center + hf + chd_chs + diabetes + sbp +
                              hnt_med + hdl + total_chol + tri + stroke + afib + 
                              SeqId_7655_11 + SeqId_4496_60 + SeqId_11089_7 +
                              SeqId_11388_75 + SeqId_8969_49 + SeqId_6586_19 + 
                              SeqId_4124_24 + SeqId_16785_45 + SeqId_11178_21 + 
                              SeqId_5008_51 + SeqId_2944_66 + SeqId_2948_58,
                            data = proteins)
#12 proteins

####Automated stepwise regression, based on AIC####
autoboth.fwstart.step <- step(base.listed.model, 
                              scope = list(upper = ~ sex + age + race + eGFR_ckdepi + smoke + bmi +
                                             height + center + hf + chd_chs + diabetes + sbp +
                                             hnt_med + hdl + total_chol + tri + stroke + afib + 
                                             SeqId_10866_60 + SeqId_11089_7 + SeqId_11109_56 + 
                                             SeqId_11178_21 + SeqId_11212_7 + SeqId_11388_75 + 
                                             SeqId_13133_73 + SeqId_15533_97 + SeqId_16322_10 + 
                                             SeqId_16785_45 + SeqId_17456_53 + SeqId_18242_8 + 
                                             SeqId_18871_24 + SeqId_19183_164 + SeqId_2602_2 + 
                                             SeqId_2677_1 + SeqId_2944_66 + SeqId_2948_58 + 
                                             SeqId_3457_57 + SeqId_3485_28 + SeqId_3710_49 + 
                                             SeqId_4124_24 + SeqId_4374_45 + SeqId_4496_60 + 
                                             SeqId_5008_51 + SeqId_6077_63 + SeqId_6366_38 + 
                                             SeqId_6586_19 + SeqId_7551_33 + SeqId_7655_11 + 
                                             SeqId_8304_50 + SeqId_8480_29 + SeqId_8969_49, 
                                           lower = ~ sex + age + race + eGFR_ckdepi + smoke + bmi +
                                             height + center + hf + chd_chs + diabetes + sbp +
                                             hnt_med + hdl + total_chol + tri + stroke + afib), 
                              direction = "both")

final.autoboth.fwstart.model <- coxph(formula = Surv(fuday_age_censored, scd_age_filt) ~ sex +
                                        age + race + eGFR_ckdepi + smoke + bmi + height + center +
                                        hf + chd_chs + diabetes + sbp + hnt_med + total_chol + tri +
                                        stroke + afib + 
                                        SeqId_7655_11 + SeqId_4496_60 + SeqId_11089_7 +
                                        SeqId_11388_75 + SeqId_8969_49 + SeqId_6586_19 + 
                                        SeqId_4124_24 + SeqId_16785_45 + SeqId_11178_21 + 
                                        SeqId_5008_51 + SeqId_2944_66 + SeqId_2948_58, data = proteins)
#12 proteins - same as forward


####Compare models####
#bw:
# SeqId_11089_7 + SeqId_11178_21 + SeqId_11212_7 +
#   SeqId_11388_75 + SeqId_16785_45 + SeqId_18242_8 + 
#   SeqId_19183_164 + SeqId_2944_66 + SeqId_4124_24 + 
#   SeqId_4496_60 + SeqId_5008_51 + SeqId_6077_63 + 
#   SeqId_6586_19 + SeqId_7655_11 + SeqId_8969_49

#fw:
# SeqId_7655_11 + SeqId_4496_60 + SeqId_11089_7 +
#   SeqId_11388_75 + SeqId_8969_49 + SeqId_6586_19 + 
#   SeqId_4124_24 + SeqId_16785_45 + SeqId_11178_21 + 
#   SeqId_5008_51 + SeqId_2944_66 + SeqId_2948_58

#Common: SeqId_7655_11 + SeqId_4496_60 + SeqId_11089_7 + 
#   SeqId_11388_75 + SeqId_8969_49 + SeqId_6586_19 + 
#   SeqId_4124_24 + SeqId_16785_45 + SeqId_11178_21 + 
#   SeqId_5008_51 + SeqId_2944_66

#Only in fw: SeqId_2948_58
#Only in bw: SeqId_11212_7, SeqId_18242_8,  SeqId_19183_164, SeqId_4496_60, SeqId_6077_63

####Try MASS package####
library(MASS)
#Scope explanation: The set of models searched is determined by the scope argument. The right-hand-side of its lower component is always included in the model, 
#and right-hand-side of the model is included in the upper component. If scope is a single formula, it specifies the upper component, and the lower model is empty. 
#If scope is missing, the initial model is used as the upper model.

stepAIC(base.listed.model, scope = list(upper = ~ sex + age + race + eGFR_ckdepi + smoke + bmi + 
                                               height + center + hf + chd_chs + diabetes + 
                                               sbp + hnt_med + hdl + total_chol + tri + stroke + afib + 
                                               SeqId_10866_60 + SeqId_11089_7 + SeqId_11109_56 + 
                                               SeqId_11178_21 + SeqId_11212_7 + SeqId_11388_75 + 
                                               SeqId_13133_73 + SeqId_15533_97 + SeqId_16322_10 + 
                                               SeqId_16785_45 + SeqId_17456_53 + SeqId_18242_8 + 
                                               SeqId_18871_24 + SeqId_19183_164 + SeqId_2602_2 + 
                                               SeqId_2677_1 + SeqId_2944_66 + SeqId_2948_58 + 
                                               SeqId_3457_57 + SeqId_3485_28 + SeqId_3710_49 + 
                                               SeqId_4124_24 + SeqId_4374_45 + SeqId_4496_60 + 
                                               SeqId_5008_51 + SeqId_6077_63 + SeqId_6366_38 + 
                                               SeqId_6586_19 + SeqId_7551_33 + SeqId_7655_11 + 
                                               SeqId_8304_50 + SeqId_8480_29 + SeqId_8969_49,
                                        lower = ~ sex + age + race + eGFR_ckdepi + smoke + bmi + 
                                               height + center + hf + chd_chs + diabetes + 
                                               sbp + hnt_med + hdl + total_chol + tri + stroke + afib), 
        direction = "both", k = 2)
#gives same result as stepwise using step

#BIC
#set k='log(n)', where n is the number of events for event models
autoboth.bic.step <- 
  stepAIC(base.listed.model, scope = list(upper = ~ sex + age + race + eGFR_ckdepi + smoke + bmi + 
                                          height + center + hf + chd_chs + diabetes + 
                                          sbp + hnt_med + hdl + total_chol + tri + stroke + afib + 
                                          SeqId_10866_60 + SeqId_11089_7 + SeqId_11109_56 + 
                                          SeqId_11178_21 + SeqId_11212_7 + SeqId_11388_75 + 
                                          SeqId_13133_73 + SeqId_15533_97 + SeqId_16322_10 + 
                                          SeqId_16785_45 + SeqId_17456_53 + SeqId_18242_8 + 
                                          SeqId_18871_24 + SeqId_19183_164 + SeqId_2602_2 + 
                                          SeqId_2677_1 + SeqId_2944_66 + SeqId_2948_58 + 
                                          SeqId_3457_57 + SeqId_3485_28 + SeqId_3710_49 + 
                                          SeqId_4124_24 + SeqId_4374_45 + SeqId_4496_60 + 
                                          SeqId_5008_51 + SeqId_6077_63 + SeqId_6366_38 + 
                                          SeqId_6586_19 + SeqId_7551_33 + SeqId_7655_11 + 
                                          SeqId_8304_50 + SeqId_8480_29 + SeqId_8969_49,
                                        lower = ~ sex + age + race + eGFR_ckdepi + smoke + bmi + 
                                          height + center + hf + chd_chs + diabetes + 
                                          sbp + hnt_med + hdl + total_chol + tri + stroke + afib), 
        direction = "both", k=log(353))

final.autoboth.bic.step.model <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi +
                                         smoke + bmi + height + center + hf + chd_chs + diabetes +
                                         sbp + hnt_med + hdl + total_chol + tri + stroke + afib +
                                         SeqId_7655_11 + SeqId_4496_60 + SeqId_11089_7 + SeqId_11388_75 +
                                         SeqId_8969_49 + SeqId_6586_19, data = proteins)

####Use rms package####
library(rms)
full.listed.model_cph <- cph(Surv(fuday_age_censored, scd_age_filt) ~ sex + age + race + eGFR_ckdepi + smoke + bmi +
                             height + center + hf + chd_chs + diabetes + sbp +
                             hnt_med + hdl + total_chol + tri + stroke + afib + 
                             SeqId_10866_60 + SeqId_11089_7 + SeqId_11109_56 + 
                             SeqId_11178_21 + SeqId_11212_7 + SeqId_11388_75 + 
                             SeqId_13133_73 + SeqId_15533_97 + SeqId_16322_10 + 
                             SeqId_16785_45 + SeqId_17456_53 + SeqId_18242_8 + 
                             SeqId_18871_24 + SeqId_19183_164 + SeqId_2602_2 + 
                             SeqId_2677_1 + SeqId_2944_66 + SeqId_2948_58 + 
                             SeqId_3457_57 + SeqId_3485_28 + SeqId_3710_49 + 
                             SeqId_4124_24 + SeqId_4374_45 + SeqId_4496_60 + 
                             SeqId_5008_51 + SeqId_6077_63 + SeqId_6366_38 + 
                             SeqId_6586_19 + SeqId_7551_33 + SeqId_7655_11 + 
                             SeqId_8304_50 + SeqId_8480_29 + SeqId_8969_49, data = proteins)
autobw.p05.step <- fastbw(full.listed.model_cph, rule = "p", type = "individual", sls = 0.05, force = c(1:21))
#same as manual BW P results

autobw.pBonf.step <- fastbw(full.listed.model_cph, rule = "p", type = "individual", sls = 0.05/33, force = c(1:21))
#same as manual BW Bonf results

####Protein correlation####
library(corrplot)
M <- cor(probes)
corrplot(M)
#Not as correlated as I thought...

####Get weights####
## Same proteins were found for forward regression

## P<0.05
p05.manbw.seqid <- coef(summary(p05.manbw.model)) %>% as.data.frame() %>% row.names() %>% .[22:29]
p05.manbw.data <- final.alldata[, c("scd_age_filt", "fuday_age_censored", rskfct, p05.manbw.seqid)]
#8 proteins
p05.manbw.weights <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ ., data = p05.manbw.data)
#26 covariates
p05.manbw.weights.send <- coef(summary(p05.manbw.weights)) %>% as.data.frame() %>% .[22:29,]
p05.manbw.weights.send$SeqId <- gsub("SeqId_", "", rownames(p05.manbw.weights.send)) %>% gsub("_", "-", .)
p05.manbw.weights.send.real <- p05.manbw.weights.send[, c("SeqId", "coef")] %>% setNames(., c("SeqId", "Beta"))
p05.manbw.weights.send.real$Set <- "ManBW.Setp05"

## P<(0.05/33)
bonf.manbw.seqid <- coef(summary(bonf.manbw.model)) %>% as.data.frame() %>% row.names() %>% .[22:27]
bonf.manbw.data <- final.alldata[, c("scd_age_filt", "fuday_age_censored", rskfct, bonf.manbw.seqid)]
#8 proteins
bonf.manbw.weights <- coxph(Surv(fuday_age_censored, scd_age_filt) ~ ., data = bonf.manbw.data)
#26 covariates
bonf.manbw.weights.send <- coef(summary(bonf.manbw.weights)) %>% as.data.frame() %>% .[22:27,]
bonf.manbw.weights.send$SeqId <- gsub("SeqId_", "", rownames(bonf.manbw.weights.send)) %>% gsub("_", "-", .)
bonf.manbw.weights.send.real <- bonf.manbw.weights.send[, c("SeqId", "coef")] %>% setNames(., c("SeqId", "Beta"))
bonf.manbw.weights.send.real$Set <- "ManBW.SetpBonf"

####Save output####
all.weights <- rbind(p05.manbw.weights.send.real, bonf.manbw.weights.send.real)
all.weights <- merge(all.weights, annotation[, c("SeqId.use", "target_use")], by.x = "SeqId", by.y = "SeqId.use", sort = T)
all.weights <- all.weights[order(all.weights$Set),]
all.weights <- all.weights[, c("SeqId", "Beta", "target_use", "Set")] %>% setNames(., c("SeqId", "Beta", "Protein", "Analysis"))
all.weights$Protein <- make.names(all.weights$Protein)

write.table(all.weights, file = "aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/backward/ARIC.RiskScore.Weights_11292022.txt", row.names = F, col.names = T, quote = F)
save.image("aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/backward/SCD.BackwardReg.Weights.PtnRiskScore.RData")
