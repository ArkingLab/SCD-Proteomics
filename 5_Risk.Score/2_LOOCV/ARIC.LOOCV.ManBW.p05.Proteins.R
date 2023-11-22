####Description####
## Run leave-one-out cross-validation ARIC PtRS analysis
## End result will generate a PtRS for each ARIC individual

####Load packages####
library(magrittr)
library(survival)
#library(rms)

####Load data####
load("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/ThuyVy.SCA.Assoc.ARIC.V2Prot_BthRcSx.RData")

weights <- read.table("aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/backward/ARIC.RiskScore.Weights_11292022.txt", header = T)
p05.seqid <- subset(weights, Analysis == "ManBW.Setp05")$SeqId %>% paste0("SeqId_", .) %>% gsub("-", "_",.)

####Run LOOCV####
pheno <- final.alldata[, c(1:106,5062:5072)]
pheno <- merge(pheno, final.alldata[, c("pid", p05.seqid)], by = "pid")

rskfct_covar <- "sex + race + age + bmi + height + center + smoke + eGFR_ckdepi +
    chd_chs + hf + diabetes + sbp + hnt_med + hdl + total_chol +
    tri + stroke + afib"

x <- paste0(p05.seqid, collapse = " + ")
x <- paste0(x, " + ", rskfct_covar)
y <- "Surv(fuday_age_censored, scd_age_filt)"
form <- as.formula(paste(y, "~", x))

#Run for all individuals
scores <- data.frame(pid = pheno$pid)
scores$sum <- NA
model.list <- list()
for(i in 1:nrow(pheno)){
  if(i == 1 | i%%250==0){
    print(paste0("Working on ", i, " pid"))
  }
  pid <- pheno$pid[i]
  train <- pheno[-i, ]
  model <- coxph(form, data = train)
  model.list[[i]] <- model
  weights <- coef(model)[1:length(p05.seqid)] %>% as.data.frame()
  names(weights)[1] <- "beta"
  proteins <- t(pheno[pheno$pid == pid, p05.seqid]) %>% as.data.frame()
  names(proteins) <- pid
  merged <- merge(proteins, weights, by = "row.names")
  merged$sum <- merged[,pid] * merged$beta
  sumScore <- sum(merged$sum)
  scores$sum[which(scores$pid == pid)] <- sumScore
}
saveRDS(scores, "aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/LOOCV.p05.ManBW/ARIC.LOOCV.ManBW.p05.Proteins.SumScore.RDS")
saveRDS(model.list, "aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/LOOCV.p05.ManBW/ARIC.LOOCV.ManBW.p05.Proteins.Models.RDS")
save.image("aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/LOOCV.p05.ManBW/ARIC.LOOCV.ManBW.p05.Proteins.LOOCV.RData")
