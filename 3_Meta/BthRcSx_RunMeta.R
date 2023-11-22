####Description####
## Run meta analysis for BthRcSx min, dz, and rskfct models
## Min = Primary model
## RskFct = Full model

####Load packages####
library(magrittr)
library(metafor)
library(dplyr)

####Load data####

## CHS data
chs.min.results <- read.csv("aric.prot.real.archive/scd.assoc/4_CHS.results/10032022/2022_10_03_SCD_proteomics_results.csv")
chs.dz.results <- read.csv("aric.prot.real.archive/scd.assoc/4_CHS.results/10032022/2022_10_03_SCD_HD_proteomics_results.csv")
chs.rskfct.results <- read.csv("aric.prot.real.archive/scd.assoc/4_CHS.results/10032022/2022_10_03_SCD_HD_RF_proteomics_results.csv")

names(chs.min.results) <- paste0("CHS.", names(chs.min.results))
names(chs.dz.results) <- paste0("CHS.", names(chs.dz.results))
names(chs.rskfct.results) <- paste0("CHS.", names(chs.rskfct.results))

names(chs.min.results)[10:13] <- c("CHS.M2.Beta", "CHS.M2.SE", "CHS.M2.P", "CHS.M2.N")
names(chs.dz.results)[10:13] <- c("CHS.M1.Beta", "CHS.M1.SE", "CHS.M1.P", "CHS.M1.N")
names(chs.rskfct.results)[10:13] <- c("CHS.M3.Beta", "CHS.M3.SE", "CHS.M3.P", "CHS.M3.N")

chs.norm.results <- merge(chs.min.results, chs.dz.results[, c("CHS.SeqId", "CHS.M1.Beta", "CHS.M1.SE", "CHS.M1.P", "CHS.M1.N")], by = "CHS.SeqId")
chs.norm.results <- merge(chs.norm.results, chs.rskfct.results[, c("CHS.SeqId", "CHS.M3.Beta", "CHS.M3.SE", "CHS.M3.P", "CHS.M3.N")], by = "CHS.SeqId")

## ARIC data
aric.norm.results <- readRDS("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/rds.files/ARIC.SCA.V2.Final_Norm_Results.RDS")

aric.norm.results$ARIC.SeqId <- gsub("SeqId_","",aric.norm.results$seqid_in_sample)
aric.norm.results$ARIC.SeqId <- gsub("_","-",aric.norm.results$ARIC.SeqId)

names(aric.norm.results)[1:47] <- paste0("ARIC.", names(aric.norm.results)[1:47])

aric.chs.norm.results <- merge(aric.norm.results, chs.norm.results, by.y = "CHS.SeqId", by.x = "ARIC.SeqId")
names(aric.chs.norm.results)[1] <- "SeqId.use"

####Run meta####
run.meta.feinv <- function(data, model){
  aric.beta <- paste0("ARIC.", model, ".Beta")
  aric.se <- paste0("ARIC.", model, ".SE")
  chs.beta <- paste0("CHS.", model, ".Beta")
  chs.se <- paste0("CHS.", model, ".SE")
  
  meta.results <- data[, c("SeqId.use", "ARIC.target", "CHS.Target")]
  meta.results$Meta.Beta <- NA
  meta.results$Meta.SE <- NA
  meta.results$Meta.Z <- NA
  meta.results$Meta.P <- NA
  meta.results$Meta.CI.upper <- NA
  meta.results$Meta.CI.lower <- NA
  meta.results$Meta.HR <- NA
  
  meta.list <- list()
  for(i in 1:nrow(data)){
    if(i %% 1000 == 0){
      print(i)
    }
    yi <- c(data[i, aric.beta], data[i, chs.beta])
    sei <- c(data[i, aric.se], data[i, chs.se])
    res = rma(yi = yi, sei = sei, method="FE")
    meta.results$Meta.Beta[i] <- res$beta[1]
    meta.results$Meta.SE[i] <- res$se[1]
    meta.results$Meta.Z[i] <- res$zval[1]
    meta.results$Meta.P[i] <- res$pval[1]
    meta.results$Meta.CI.upper[i] <- res$ci.lb[1]
    meta.results$Meta.CI.lower[i] <- res$ci.ub[1]
    meta.results$Meta.HR[i] <- exp(res$beta[1])
    meta.results$Meta.HR.CI.upper[i] <- exp(res$beta[1] + (1.96*res$se[1]))
    meta.results$Meta.HR.CI.lower[i] <- exp(res$beta[1] - (1.96*res$se[1]))
    
    meta.list[[i]] <- res
  }
  meta.results <- meta.results[order(meta.results$Meta.P),]
  meta.results$Meta.PASSED <- ifelse(meta.results$Meta.P < (0.05/4955), "YES", "NO")
  meta.results$model <- model
  meta.list[[4956]] <- meta.results
  
  names(meta.list)[1:4955] <- aric.chs.norm.results$SeqId
  names(meta.list)[4956] <- "meta.results"
  
  return(meta.list)
}

## Min
min.meta.list <- run.meta.feinv(aric.chs.norm.results, "M2")
min.meta.results <- min.meta.list[[4956]] %>% as.data.frame()

## Dz
dz.meta.list <- run.meta.feinv(aric.chs.norm.results, "M1")
dz.meta.results <- dz.meta.list[[4956]] %>% as.data.frame()

## RskFct
rskfct.meta.list <- run.meta.feinv(aric.chs.norm.results, "M3")
rskfct.meta.results <- rskfct.meta.list[[4956]] %>% as.data.frame()

####Combine data for forestplot####

## Min Meta
min.all.results <- min.meta.results[, -c(2,3)]
min.all.results$cohort <- "meta"
min.all.results$model.cohort <- "M2.meta"
names(min.all.results)[2:11] <- gsub("Meta.", "", names(min.all.results)[2:11])

## Dz Meta
dz.all.results <- dz.meta.results[, -c(2,3)]
dz.all.results$cohort <- "meta"
dz.all.results$model.cohort <- "M1.meta"
names(dz.all.results)[2:11] <- gsub("Meta.", "", names(dz.all.results)[2:11])

## RskFct Meta
rskfct.all.results <- rskfct.meta.results[, -c(2,3)]
rskfct.all.results$cohort <- "meta"
rskfct.all.results$model.cohort <- "M3.meta"
names(rskfct.all.results)[2:11] <- gsub("Meta.", "", names(rskfct.all.results)[2:11])

meta.formatted <- as.data.frame(rbind(min.all.results, dz.all.results, rskfct.all.results))

## ARIC
load("aric.prot.real.archive/scd.assoc/6_manuscript/ver1/Figures.Tables/05162022/ARIC.Max.Risk/ARIC.V2.Max.Risk.Ez.Extract.rds")
rm(v2.norm.results)
aric.formatted <- v2.norm.supp.heatmap.df[, c("seqid_in_sample", "model", "P", "PASSED", "Beta", "SE", "HR")]
aric.formatted$Z <- (aric.formatted$Beta/aric.formatted$SE)
aric.formatted$CI.upper <- (aric.formatted$Beta + (1.96*aric.formatted$SE))
aric.formatted$CI.lower <- (aric.formatted$Beta - (1.96*aric.formatted$SE))
aric.formatted$HR.CI.upper <- exp(aric.formatted$Beta + (1.96*aric.formatted$SE))
aric.formatted$HR.CI.lower <- exp(aric.formatted$Beta - (1.96*aric.formatted$SE))
aric.formatted$cohort <- "aric"
aric.formatted$model.cohort <- paste0(aric.formatted$model, ".aric")
aric.formatted$SeqId.use <- gsub("SeqId_","",aric.formatted$seqid_in_sample)
aric.formatted$SeqId.use <- gsub("_","-",aric.formatted$SeqId.use)

columns <- names(min.all.results)
aric.formatted <- aric.formatted[, columns]

## CHS
chs.formatted <- rbind(read.csv("aric.prot.real.archive/scd.assoc/4_CHS.results/10032022/2022_10_03_SCD_proteomics_results.csv") %>% mutate(model = "M2", cohort = "chs", model.cohort = "M2.chs"),
                       read.csv("aric.prot.real.archive/scd.assoc/4_CHS.results/10032022/2022_10_03_SCD_HD_proteomics_results.csv") %>% mutate(model = "M1", cohort = "chs", model.cohort = "M1.chs"),
                       read.csv("aric.prot.real.archive/scd.assoc/4_CHS.results/10032022/2022_10_03_SCD_HD_RF_proteomics_results.csv") %>% mutate(model = "M3", cohort = "chs", model.cohort = "M3.chs")) %>% as.data.frame()

chs.formatted$Z <- (chs.formatted$beta/chs.formatted$SE)
chs.formatted$CI.upper <- (chs.formatted$beta + (1.96*chs.formatted$SE))
chs.formatted$CI.lower <- (chs.formatted$beta - (1.96*chs.formatted$SE))
chs.formatted$HR.CI.upper <- exp(chs.formatted$beta + (1.96*chs.formatted$SE))
chs.formatted$HR.CI.lower <- exp(chs.formatted$beta - (1.96*chs.formatted$SE))
names(chs.formatted)[10] <- "Beta"
chs.formatted$HR <- exp(chs.formatted$Beta)
names(chs.formatted)[1] <- "SeqId.use"
chs.formatted$PASSED <- ifelse(chs.formatted$P < (0.05/4955), "YES", "NO")
chs.formatted <- chs.formatted[, columns]

aric.chs.meta.results <- rbind(aric.formatted, chs.formatted, meta.formatted) %>% as.data.frame()
aric.chs.meta.results <- merge(aric.chs.meta.results, aric.chs.norm.results[, c(1,3,51,13,14,19:50,52:56)], by = "SeqId.use", sort = F)

## Reference list for proteins
prot.ref.list <- aric.chs.norm.results[, c("SeqId.use", "ARIC.target", "CHS.Target")]
prot.ref.list$RefNum <- 1:nrow(prot.ref.list)

find.ptn.target <- function(search, input){
  if(input == "seqid"){
    subset(prot.ref.list, SeqId.use == search)
  } else if(input == "ptn"){
    subset(prot.ref.list, ARIC.target == search)
  }
}

####Save data for export####
saveRDS(aric.chs.meta.results, file = "aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/BthRcSx.Meta.Results.RDS")
save(prot.ref.list, find.ptn.target, file = "aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/PtnRefList.rds")
save.image("aric.prot.real.archive/scd.assoc/5_meta/with.10032022.CHS/BthRcSx_RunMeta.RData")
