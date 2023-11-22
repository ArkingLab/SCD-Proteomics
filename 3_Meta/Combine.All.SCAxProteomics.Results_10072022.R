####Description####
## Merge main datasets from SCA x Proteomics analyses
## Only ARIC V2 results are merged
## Min = Primary model
## RskFct = Full model

####Load packages####
library(magrittr)
library(dplyr)
library(stringr)

####Meta, ARIC, CHS - BthRcSx####
aric.chs.meta.results <- readRDS("aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/Meta.Max.BthRcSx.Processed.RDS")

columns <- c("SeqId.use", "target_use", "Beta", "SE", "Z", "P", "CI.upper", "CI.lower", "HR", "HR.CI.upper", "HR.CI.lower", "model", "cohort", "PASSED_PLOT", "model_use", "cohort_use")

aric.chs.meta.results.trimmed <- aric.chs.meta.results[, columns]
aric.chs.meta.results.trimmed$strat <- "all"
aric.chs.meta.results.trimmed$strat_use <- "All"

####Meta, ARIC, CHS - CHD stratification####
aric.chs.meta.yeschd.nochd.results <- readRDS("aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/chd/Meta.Max.BthRcSx.yesCHD.noCHD.Processed.RDS")

aric.chs.meta.yeschd.nochd.results.trimmed <- aric.chs.meta.yeschd.nochd.results[aric.chs.meta.yeschd.nochd.results$chd.stats != "both", c(columns, "chd.stats")]
names(aric.chs.meta.yeschd.nochd.results.trimmed)[17] <- "strat"
aric.chs.meta.yeschd.nochd.results.trimmed$strat_use <- ifelse(aric.chs.meta.yeschd.nochd.results.trimmed$strat == "yesCHD", "CHD", "No CHD")

####Meta, ARIC, CHS - Race stratification####
aric.chs.meta.white.black.results <- readRDS("aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/race/Meta.Max.BthRcSx.White.Black.Processed.RDS")

aric.chs.meta.white.black.results.trimmed <- aric.chs.meta.white.black.results[aric.chs.meta.white.black.results$race != "both", c(columns, "race")]
names(aric.chs.meta.white.black.results.trimmed)[17] <- "strat"
aric.chs.meta.white.black.results.trimmed$strat_use <- ifelse(aric.chs.meta.white.black.results.trimmed$strat == "white", "White", "Black")

####Meta, ARIC, CHS - Sex stratification####
aric.chs.meta.fem.male.results <- readRDS("aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/sex/Meta.Max.BthRcSx.Fem.Male.Processed.RDS")

aric.chs.meta.fem.male.results.trimmed <- aric.chs.meta.fem.male.results[aric.chs.meta.fem.male.results$sex != "both", c(columns, "sex")]
names(aric.chs.meta.fem.male.results.trimmed)[17] <- "strat"
aric.chs.meta.fem.male.results.trimmed$strat_use <- ifelse(aric.chs.meta.fem.male.results.trimmed$strat == "female", "Female", "Male")

####ARIC yesCHD and noCHD - Dz and RskFct models####
aric.yeschd.results <- readRDS("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/postStrat/heart.disease.only/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/rds.files/ARIC.SCA.V2.Final_Norm_Results.RDS")

aric.yeschd.results$SeqId.use <- gsub("SeqId_", "", aric.yeschd.results$seqid_in_sample)
aric.yeschd.results$SeqId.use <- gsub("_", "-", aric.yeschd.results$SeqId.use)
aric.yeschd.results.trimmed <- 
  rbind(aric.yeschd.results[, c("SeqId.use", "M1.Beta", "M1.SE", "M1.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M1", model_use = "Disease", strat = "yesCHD", strat_use = "CHD") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
      aric.yeschd.results[, c("SeqId.use", "M3.Beta", "M3.SE", "M3.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M3", model_use = "Risk Factor", strat = "yesCHD", strat_use = "CHD") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")))

aric.nochd.results <- readRDS("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/postStrat/no.heart.disease/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/rds.files/ARIC.SCA.V2.Final_Norm_Results.RDS")

aric.nochd.results$SeqId.use <- gsub("SeqId_", "", aric.nochd.results$seqid_in_sample)
aric.nochd.results$SeqId.use <- gsub("_", "-", aric.nochd.results$SeqId.use)
aric.nochd.results.trimmed <- 
  rbind(aric.nochd.results[, c("SeqId.use", "M1.Beta", "M1.SE", "M1.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M1", model_use = "Disease", strat = "noCHD", strat_use = "No CHD") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.nochd.results[, c("SeqId.use", "M3.Beta", "M3.SE", "M3.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M3", model_use = "Risk Factor", strat = "noCHD", strat_use = "No CHD") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")))

####ARIC female and male - Dz and RskFct models####
aric.female.results <- readRDS("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/postStrat/female/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/rds.files/ARIC.SCA.V2.Final_Norm_Results.RDS")

aric.female.results$SeqId.use <- gsub("SeqId_", "", aric.female.results$seqid_in_sample)
aric.female.results$SeqId.use <- gsub("_", "-", aric.female.results$SeqId.use)
aric.female.results.trimmed <- 
  rbind(aric.female.results[, c("SeqId.use", "M1.Beta", "M1.SE", "M1.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M1", model_use = "Disease", strat = "female", strat_use = "Female") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.female.results[, c("SeqId.use", "M3.Beta", "M3.SE", "M3.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M3", model_use = "Risk Factor", strat = "female", strat_use = "Female") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")))

aric.male.results <- readRDS("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/postStrat/male/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/rds.files/ARIC.SCA.V2.Final_Norm_Results.RDS")

aric.male.results$SeqId.use <- gsub("SeqId_", "", aric.male.results$seqid_in_sample)
aric.male.results$SeqId.use <- gsub("_", "-", aric.male.results$SeqId.use)
aric.male.results.trimmed <- 
  rbind(aric.male.results[, c("SeqId.use", "M1.Beta", "M1.SE", "M1.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M1", model_use = "Disease", strat = "male", strat_use = "Male") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.male.results[, c("SeqId.use", "M3.Beta", "M3.SE", "M3.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M3", model_use = "Risk Factor", strat = "male", strat_use = "Male") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")))

####ARIC Black and White - Dz and RskFct models####
aric.white.results <- readRDS("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/postStrat/EA/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/rds.files/ARIC.SCA.V2.Final_Norm_Results.RDS")

aric.white.results$SeqId.use <- gsub("SeqId_", "", aric.white.results$seqid_in_sample)
aric.white.results$SeqId.use <- gsub("_", "-", aric.white.results$SeqId.use)
aric.white.results.trimmed <- 
  rbind(aric.white.results[, c("SeqId.use", "M1.Beta", "M1.SE", "M1.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M1", model_use = "Disease", strat = "white", strat_use = "White") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.white.results[, c("SeqId.use", "M3.Beta", "M3.SE", "M3.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M3", model_use = "Risk Factor", strat = "white", strat_use = "White") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")))

aric.black.results <- readRDS("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/postStrat/AA/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/rds.files/ARIC.SCA.V2.Final_Norm_Results.RDS")

aric.black.results$SeqId.use <- gsub("SeqId_", "", aric.black.results$seqid_in_sample)
aric.black.results$SeqId.use <- gsub("_", "-", aric.black.results$SeqId.use)
aric.black.results.trimmed <- 
  rbind(aric.black.results[, c("SeqId.use", "M1.Beta", "M1.SE", "M1.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M1", model_use = "Disease", strat = "black", strat_use = "Black") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.black.results[, c("SeqId.use", "M3.Beta", "M3.SE", "M3.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M3", model_use = "Risk Factor", strat = "black", strat_use = "Black") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")))

####ARIC Age stratification####
## Stratify by median age - 57
aric.lowerAge.results <- readRDS("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/postStrat/low.age/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/rds.files/ARIC.SCA.V2.Final_Norm_Results.RDS")
  
aric.lowerAge.results$SeqId.use <- gsub("SeqId_", "", aric.lowerAge.results$seqid_in_sample)
aric.lowerAge.results$SeqId.use <- gsub("_", "-", aric.lowerAge.results$SeqId.use)
aric.lowerAge.results.trimmed <- 
  rbind(aric.lowerAge.results[, c("SeqId.use", "M2.Beta", "M2.SE", "M2.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M2", model_use = "Minimum", strat = "lowerAge", strat_use = "Young") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.lowerAge.results[, c("SeqId.use", "M1.Beta", "M1.SE", "M1.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M1", model_use = "Disease", strat = "lowerAge", strat_use = "Young") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.lowerAge.results[, c("SeqId.use", "M3.Beta", "M3.SE", "M3.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M3", model_use = "Risk Factor", strat = "lowerAge", strat_use = "Young") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")))

aric.higherAge.results <- readRDS("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/postStrat/high.age/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/rds.files/ARIC.SCA.V2.Final_Norm_Results.RDS")

aric.higherAge.results$SeqId.use <- gsub("SeqId_", "", aric.higherAge.results$seqid_in_sample)
aric.higherAge.results$SeqId.use <- gsub("_", "-", aric.higherAge.results$SeqId.use)
aric.higherAge.results.trimmed <- 
  rbind(aric.higherAge.results[, c("SeqId.use", "M2.Beta", "M2.SE", "M2.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M2", model_use = "Minimum", strat = "higherAge", strat_use = "Old") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.higherAge.results[, c("SeqId.use", "M1.Beta", "M1.SE", "M1.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M1", model_use = "Disease", strat = "higherAge", strat_use = "Old") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.higherAge.results[, c("SeqId.use", "M3.Beta", "M3.SE", "M3.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M3", model_use = "Risk Factor", strat = "higherAge", strat_use = "Old") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")))

####ARIC female pre vs. post-menopause results####
aric.meno.results <- readRDS("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/postStrat/female/menopause/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/rds.files/ARIC.SCA.V2.Final_Norm_Results.RDS")

aric.meno.results$SeqId.use <- gsub("SeqId_", "", aric.meno.results$seqid_in_sample)
aric.meno.results$SeqId.use <- gsub("_", "-", aric.meno.results$SeqId.use)
aric.meno.results.trimmed <- 
  rbind(aric.meno.results[, c("SeqId.use", "M2.Beta", "M2.SE", "M2.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M2", model_use = "Minimum", strat = "meno", strat_use = "Post-menopause") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.meno.results[, c("SeqId.use", "M1.Beta", "M1.SE", "M1.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M1", model_use = "Disease", strat = "meno", strat_use = "Post-menopause") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.meno.results[, c("SeqId.use", "M3.Beta", "M3.SE", "M3.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M3", model_use = "Risk Factor", strat = "meno", strat_use = "Post-menopause") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")))

aric.premeno.results <- readRDS("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/postStrat/female/premenopause/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/rds.files/ARIC.SCA.V2.Final_Norm_Results.RDS")

aric.premeno.results$SeqId.use <- gsub("SeqId_", "", aric.premeno.results$seqid_in_sample)
aric.premeno.results$SeqId.use <- gsub("_", "-", aric.premeno.results$SeqId.use)
aric.premeno.results.trimmed <- 
  rbind(aric.premeno.results[, c("SeqId.use", "M2.Beta", "M2.SE", "M2.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M2", model_use = "Minimum", strat = "premeno", strat_use = "Pre-menopause") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.premeno.results[, c("SeqId.use", "M1.Beta", "M1.SE", "M1.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M1", model_use = "Disease", strat = "premeno", strat_use = "Pre-menopause") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.premeno.results[, c("SeqId.use", "M3.Beta", "M3.SE", "M3.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M3", model_use = "Risk Factor", strat = "premeno", strat_use = "Pre-menopause") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")))

####ARIC 5 YR models####
aric.five.results <- readRDS("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/rds.files/ARIC.SCA.V2.Final_Five_Results.RDS")

aric.five.results$SeqId.use <- gsub("SeqId_", "", aric.five.results$seqid_in_sample)
aric.five.results$SeqId.use <- gsub("_", "-", aric.five.results$SeqId.use)
aric.five.results.trimmed <- 
  rbind(aric.five.results[, c("SeqId.use", "M25.Beta", "M25.SE", "M25.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M2", model_use = "Minimum", strat = "five", strat_use = "5 YR") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.five.results[, c("SeqId.use", "M15.Beta", "M15.SE", "M15.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M1", model_use = "Disease", strat = "five", strat_use = "5 YR") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.five.results[, c("SeqId.use", "M35.Beta", "M35.SE", "M35.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M3", model_use = "Risk Factor", strat = "five", strat_use = "5 YR") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")))

####ARIC 10 YR models####
aric.ten.results <- readRDS("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/rds.files/ARIC.SCA.V2.Final_Ten_Results.RDS")

aric.ten.results$SeqId.use <- gsub("SeqId_", "", aric.ten.results$seqid_in_sample)
aric.ten.results$SeqId.use <- gsub("_", "-", aric.ten.results$SeqId.use)
aric.ten.results.trimmed <- 
  rbind(aric.ten.results[, c("SeqId.use", "M210.Beta", "M210.SE", "M210.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M2", model_use = "Minimum", strat = "ten", strat_use = "10 YR") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.ten.results[, c("SeqId.use", "M110.Beta", "M110.SE", "M110.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M1", model_use = "Disease", strat = "ten", strat_use = "10 YR") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")),
        aric.ten.results[, c("SeqId.use", "M310.Beta", "M310.SE", "M310.P")] %>% mutate(cohort = "aric", cohort_use = "ARIC", model = "M3", model_use = "Risk Factor", strat = "ten", strat_use = "10 YR") %>% set_colnames(c("SeqId.use", "Beta", "SE", "P", "cohort", "cohort_use", "model", "model_use", "strat", "strat_use")))

####Combine ARIC Dz and RskFct models####
aric.other.mods <- rbind(aric.black.results.trimmed, aric.white.results.trimmed, aric.yeschd.results.trimmed, aric.nochd.results.trimmed, aric.female.results.trimmed, aric.male.results.trimmed, aric.premeno.results.trimmed, aric.meno.results.trimmed, aric.lowerAge.results.trimmed, aric.higherAge.results.trimmed, aric.five.results.trimmed, aric.ten.results.trimmed)
aric.other.mods$HR <- exp(aric.other.mods$Beta)
aric.other.mods$Z <- aric.other.mods$Beta/aric.other.mods$SE
aric.other.mods$CI.upper <- (aric.other.mods$Beta) + (aric.other.mods$SE*1.96)
aric.other.mods$CI.lower <- (aric.other.mods$Beta) - (aric.other.mods$SE*1.96)
aric.other.mods$HR.CI.upper <- exp((aric.other.mods$Beta) + (aric.other.mods$SE*1.96))
aric.other.mods$HR.CI.lower <- exp((aric.other.mods$Beta) - (aric.other.mods$SE*1.96))
aric.other.mods$PASSED_PLOT <- ifelse(aric.other.mods$P < (0.05/4955), "*", "")
aric.other.mods <- merge(aric.chs.meta.results.trimmed[, c("SeqId.use", "target_use")] %>% unique(), aric.other.mods, by = "SeqId.use")
aric.other.mods <- aric.other.mods[, names(aric.chs.meta.results.trimmed)]

####Combine all data####
alldata <- rbind(aric.chs.meta.results.trimmed, aric.chs.meta.white.black.results.trimmed, aric.chs.meta.fem.male.results.trimmed, aric.chs.meta.yeschd.nochd.results.trimmed, aric.other.mods)
alldata <- merge(alldata, aric.chs.meta.results[, c("SeqId.use", "order_meta.M2.p.beta", "order_meta.M2.p", "order_meta.M2.p.bonf", "order_rand")] %>% unique(), by = "SeqId.use")
alldata <- merge(alldata, aric.chs.meta.results[, c(1,15:55)] %>% unique(), by = "SeqId.use")
alldata$model.cohort <- paste0(alldata$model, ".", alldata$cohort)
alldata$model.cohort_use <- paste0(alldata$cohort_use, "-", alldata$model_use)
alldata$model.cohort.strat <- paste0(alldata$model, ".", alldata$cohort, ".", alldata$strat)
alldata$model.cohort.strat_use <- paste0(alldata$cohort_use, "-", alldata$model_use, "-", alldata$strat_use)

####Figure out heatmap color limit####
meta.min.bonf.seqid <- subset(alldata, cohort == "meta" & strat == "all" & model == "M2")
meta.min.bonf.seqid <- subset(meta.min.bonf.seqid, P<(0.05/4955))$SeqId.use
limits.df <- subset(alldata, SeqId.use %in% meta.min.bonf.seqid)
heatmap.color.limits <- range(limits.df$HR)
saveRDS(heatmap.color.limits, "Heatmap.Color.Limits_10072022.RDS")

####Fix fibulin 5 spelling####
alldata$target_use <- ifelse(alldata$target_use == "fibulin 5", "Fibulin 5", alldata$target_use)

####Change RskFct model naming####
alldata$model_use <- as.character(alldata$model_use)
alldata$model_use <- ifelse(alldata$model_use == "Risk Factor", "Disease + Risk Factor", alldata$model_use)
alldata$model.cohort_use <- paste0(alldata$cohort_use, "-", alldata$model_use)
alldata$model.cohort.strat_use <- paste0(alldata$cohort_use, "-", alldata$strat_use, "-", alldata$model_use)
alldata$model.abrv_use <- NA
alldata$model.abrv_use <- ifelse(alldata$model == "M2", "Min", alldata$model.abrv_use)
alldata$model.abrv_use <- ifelse(alldata$model == "M1", "Dz", alldata$model.abrv_use)
alldata$model.abrv_use <- ifelse(alldata$model == "M3", "Dz+RskFct", alldata$model.abrv_use)
alldata$model.abrv_use <- factor(alldata$model.abrv_use, ordered = T, levels = c("Min", "Dz", "Dz+RskFct"))
alldata$model_use <- factor(alldata$model_use, ordered = T, levels = c("Minimum", "Disease", "Disease + Risk Factor"))

####Add stars to duplicated targets####
dup.targets <- subset(alldata, cohort == "meta" & model == "M2") %>% .[, c("SeqId.use", "target_use")] %>% unique()
dup.targets <- subset(dup.targets, duplicated(target_use) | duplicated(target_use, fromLast = T))
dup.targets <- dup.targets[order(dup.targets$target_use),]
dup.targets.freq <- as.data.frame(table(dup.targets$target_use))
doubles <- subset(dup.targets.freq, Freq == 2)
triples <- subset(dup.targets.freq, Freq == 3)
triples.fixed <- subset(dup.targets, target_use %in% triples$Var1)
triples.fixed$target_use <- c("APOE", "APOE*", "APOE**", "DLK1", "DLK1*", "DLK1**")
doubles.fixed <- subset(dup.targets, target_use %in% doubles$Var1)
doubles.fixed$order <- 1:nrow(doubles.fixed)
for(i in 1:nrow(doubles.fixed)){
  if(doubles.fixed$order[i]%%2==0){
    doubles.fixed$target_use[i] <- paste0(doubles.fixed$target_use[i], "*")
  }
}
doubles.fixed$order <- NULL
dup.targets.fixed <- rbind(doubles.fixed, triples.fixed)
names(dup.targets.fixed)[2] <- "target_use2"
alldata <- merge(alldata, dup.targets.fixed, by = "SeqId.use", all.x = T)
alldata$target_use <- ifelse(!is.na(alldata$target_use2), alldata$target_use2, alldata$target_use)
alldata$target_use2 <- NULL

####Fix lowercase first letter to upper case####
capped.targets <- data.table::fread("capped.targets.txt", header = F)
names(capped.targets) <- c("SeqId.use", "target_use2")
alldata <- merge(alldata, capped.targets, by = "SeqId.use", all.x = T)
alldata$target_use <- ifelse(!is.na(alldata$target_use2), alldata$target_use2, alldata$target_use)
alldata$target_use2 <- NULL
table(duplicated(subset(alldata, cohort == "meta" & model == "M2" & strat == "all")$target_use))
#FALSE
#4955

####Save spot####
saveRDS(alldata[, c("SeqId.use", "target_use")] %>% unique(), file = "aric.prot.real.archive/scd.assoc/9_other/combine.data/Unique.SeqId.Target.RDS")
saveRDS(alldata, file = "aric.prot.real.archive/scd.assoc/9_other/combine.data/All.SCAxProteomics.Results_10072022.RDS")
save.image("aric.prot.real.archive/scd.assoc/9_other/combine.data/Combine.All.SCAxProteomics.Results_10072022.RData")