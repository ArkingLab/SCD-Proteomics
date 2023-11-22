####Description#####
## Create figures and data for PtRS clinical measures for research letter
## Output .txt is base data for Panel I (ver19)

####Load packages####
library(magrittr)

####Load data####
aric_cstat <- read.csv("aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/clinical.metrics/ARIC.Prot.Risk.Score.Cindex_Results.csv")
aric_nri <- read.csv("aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/clinical.metrics/ARIC.Prot.Risk.Score.NRI_Results.csv")

chs_cstat <- read.csv("aric.prot.real.archive/scd.assoc/4_CHS.results/10032022/risk.score_12052022/clinical.measures/CHS.SCD.Prot.Cindex.Results_2023_03_06.csv")
chs_nri <- read.csv("aric.prot.real.archive/scd.assoc/4_CHS.results/10032022/risk.score_12052022/clinical.measures/CHS.SCD.Prot.NRI.Results_2023_03_06.csv")

####Edit ARIC data####
aric_cstat$cohort <- "aric"
aric_nri$cohort <- "aric"
chs_cstat$cohort <- "chs"
chs_nri$cohort <- "chs"

aric_cstat_max <- aric_cstat[1:2,]
aric_cstat_max[, 4:7] <- round(aric_cstat_max[, 4:7], digits = 2)
aric_cstat_max[, 8:12] <- round(aric_cstat_max[, 8:12], digits = 3)
aric_cstat_max$combo <- paste0(aric_cstat_max$cstat, " [", aric_cstat_max$lower, "-", aric_cstat_max$upper, "]")

aric_nri_max <- aric_nri
aric_nri_max[,1:3] <- round(aric_nri_max [,1:3], digits = 2)

####Edit CHS data####
chs_cstat$upper <- chs_cstat$cstat+(1.96*chs_cstat$se)
chs_cstat$lower <- chs_cstat$cstat-(1.96*chs_cstat$se)
chs_cstat_max <- chs_cstat[1:2,]
chs_cstat_max <- chs_cstat_max[, names(aric_cstat_max)[1:13]]
chs_cstat_max[, 4:7] <- round(chs_cstat_max[, 4:7], digits = 2)
chs_cstat_max[, 8] <- round(chs_cstat_max[, 8], digits = 4)
chs_cstat_max[, 9:12] <- round(chs_cstat_max[, 9:12], digits = 3)
chs_cstat_max$combo <- paste0(chs_cstat_max$cstat, " [", chs_cstat_max$lower, "-", chs_cstat_max$upper, "]")

chs_nri_max <- chs_nri
chs_nri_max[,2:4] <- round(chs_nri_max[,2:4], digits = 2)
chs_nri_max$Value <- c(paste0("M", 1:3))

####Combine data####
final.table <- as.data.frame(matrix(nrow=2, ncol=6))
names(final.table) <- c("Cohort", "Base_Cstat", "PtRS_Cstat", "Change_Cstat", "5YR_NRI", "10YR_NRI")
final.table$Cohort <- c("ARIC", "CHS")

final.table$Base_Cstat[1] <- aric_cstat_max$combo[1]
final.table$PtRS_Cstat[1] <- aric_cstat_max$combo[2]
final.table$Change_Cstat[1] <- paste0(aric_cstat_max$change[1], " [", aric_cstat_max$change.norm.lci[1], "-", aric_cstat_max$change.norm.uci[1], "]")
final.table$`5YR_NRI`[1] <- paste0(aric_nri_max$`Est.`[1], " [", aric_nri_max$Lower[1], "-", aric_nri_max$Upper[1], "]")
final.table$`10YR_NRI`[1] <- paste0(aric_nri_max$`Est.`[2], " [", aric_nri_max$Lower[2], "-", aric_nri_max$Upper[2], "]")

final.table$Base_Cstat[2] <- chs_cstat_max$combo[1]
final.table$PtRS_Cstat[2] <- chs_cstat_max$combo[2]
final.table$Change_Cstat[2] <- paste0(chs_cstat_max$change[1], " [", chs_cstat_max$change.norm.lci[1], "-", chs_cstat_max$change.norm.uci[1], "]")
final.table$`5YR_NRI`[2] <- paste0(chs_nri_max$`Est`[2], " [", chs_nri_max$Lower[2], "-", chs_nri_max$Upper[2], "]")
final.table$`10YR_NRI`[2] <- paste0(chs_nri_max$`Est`[5], " [", chs_nri_max$Lower[5], "-", chs_nri_max$Upper[5], "]")

####Save output####
write.table(final.table, "aric.prot.real.archive/scd.assoc/6_manuscript/ver13/Risk.Score/ARIC.CHS.Protein.Risk.Score.NPPB_Clinical.Measures.txt", col.names = T,
            row.names = F, sep="\t", quote = F)
#Panel I (ver19)

save.image("aric.prot.real.archive/scd.assoc/6_manuscript/ver13/Risk.Score/ARIC.CHS.Protein.Risk.Score.NPPB_Clinical.Measures_ForManuscript.RData")
