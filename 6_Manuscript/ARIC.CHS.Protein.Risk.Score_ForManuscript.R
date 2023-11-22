####Description####
## Creates figures and data of PtRS for research letter
## Panel H (ver19) is generated

####Load packages####
library(magrittr)
library(readxl)

####Edit ARIC data#####
aric <- read_xlsx("aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/ARIC.LOOCV.Risk.Score.Results.ManBW.p05_Main.Table.xlsx") %>% as.data.frame()
aric.num <- read.table("aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/ARIC.LOOCV.ManBW.p05.Quintile.Numbers.txt", header = T)
aric.num.final <- merge(subset(aric.num, SCD == "Cases"), subset(aric.num, SCD != "Cases"), by = "Quintile")
aric.num.final <- aric.num.final[, c("Quintile", "Freq.x", "Freq.y")] %>% setNames(., c("Quintile", "Cases", "Non-cases")) %>% as.data.frame()
aric.num.final[6,] <- c("Linear", "353", "11,048")
aric.num.final <- aric.num.final[c(6,1:5),]
aric <- as.data.frame(cbind(aric.num.final, aric))

#Add in 0
aric[6,5] <- "6.20 [3.92-9.83]"

####Edit CHS data####
chs <- read.csv("aric.prot.real.archive/scd.assoc/4_CHS.results/10032022/risk.score_12052022/Risk_Score_Results_2022_12_02.csv") %>% as.data.frame()
chs$Quintile <- rep(c(aric$Quintile),2)
chs <- chs[, -c(8:10)]

chs$Beta <- as.numeric(chs$Beta)
chs$SE <- as.numeric(chs$SE)
chs$P <- as.numeric(chs$P) %>% signif(., digits = 2) 

chs$HR <- exp(chs$Beta) %>% round(., digits = 3) %>% sprintf(fmt = "%.2f",.)
chs$lower.ci <- exp(chs$Beta - 1.96*chs$SE) %>% round(., digits = 3) %>% sprintf(fmt = "%.2f",.)
chs$upper.ci <- exp(chs$Beta + 1.96*chs$SE) %>% round(., digits = 3) %>% sprintf(fmt = "%.2f",.)
chs$`HR [95% CI]` <- paste0(chs$HR, " [", chs$lower.ci, "-", chs$upper.ci, "]")

names(chs)[3] <- "Non-cases"
names(chs)[7] <- "NPPB"

####Edit ARIC data w/o NPPB####
aric.no.nppb <- read_xlsx("aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/no.NPPB/ARIC.LOOCV.Risk.Score.Results.ManBW.p05.noNPPB_Main.Table.xlsx")
aric.num.no.nppb <- read.table("aric.prot.real.archive/scd.assoc/8_risk.score/stepwise/no.NPPB/ARIC.LOOCV.ManBW.p05.no.NPPB.Quintile.Numbers.txt", header = T)
aric.num.no.nppb.final <- merge(subset(aric.num.no.nppb, SCD == "Cases"), subset(aric.num.no.nppb, SCD != "Cases"), by = "Quintile")
aric.num.no.nppb.final <- aric.num.no.nppb.final[, c("Quintile", "Freq.x", "Freq.y")] %>% setNames(., c("Quintile", "Cases", "Non-cases")) %>% as.data.frame()
aric.num.no.nppb.final[6,] <- c("Linear", "353", "11,048")
aric.num.no.nppb.final <- aric.num.no.nppb.final[c(6,1:5),]
aric.no.nppb <- as.data.frame(cbind(aric.num.no.nppb.final, aric.no.nppb))
aric.no.nppb$NPPB <- 'N'
aric$NPPB <- 'Y'

aric <- rbind(aric, aric.no.nppb) %>% as.data.frame()
aric$Cohort <- "ARIC"
aric$Variable <- NULL
chs$Cohort <- "CHS"

####Combine cohorts####
columns <- names(aric)
final.table <- rbind(aric, chs[, columns]) %>% as.data.frame()

final.table$Quintile <- rep(c("Linear", paste0("Q", 1:5)), 4)

####Save and write output####
write.table(subset(final.table, NPPB == "Y"), 
            sep = "\t", quote = F, row.names = F, col.names = T, 
            file = "aric.prot.real.archive/scd.assoc/6_manuscript/ver13/Risk.Score/ARIC.CHS.Protein.Risk.Score.NPPB_Final.Results.txt")
#Base data used for Panel H (ver19)

write.table(subset(final.table, NPPB == "N"), 
            sep = "\t", quote = F, row.names = F, col.names = T, 
            file = "aric.prot.real.archive/scd.assoc/6_manuscript/ver13/Risk.Score/ARIC.CHS.Protein.Risk.Score.no.NPPB_Final.Results.txt")

save.image("aric.prot.real.archive/scd.assoc/6_manuscript/ver13/Risk.Score/ARIC.CHS.Protein.Risk.Score_ForManuscript.RData")
