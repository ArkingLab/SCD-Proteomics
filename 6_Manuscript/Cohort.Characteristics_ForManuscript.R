####Description####
## Generates summary statistics of cohort characteristics for research letter
## Output .txt is used as base data for Panel B (ver19)

####Load packages####
library(magrittr)

####ARIC Visit 2 Table####
meansdcalc <- function(data, variable){
  mean_num <- mean(data[,variable]) %>% round(., digits = 2) %>% sprintf("%.1f", .)
  sd_num <- sd(data[,variable]) %>% round(., digits = 2) %>% sprintf("%.1f", .)
  out <- paste0(mean_num, " (", sd_num, ")")
  return(out)
}
npercalc <- function(data, variable, condition){
  filt <- data[data[,variable] == condition,]
  n <- nrow(filt)
  total <- nrow(data)
  nper <- ((n/total) %>% round(., digits = 3))*100
  out <- paste0(n, " (", nper, "%)")
  return(out)
}

v2_pheno <- readRDS("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit2/combined/exclude.egfr/with.afib/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/rds.files/ARIC.SCA.V2.Final_Pheno.RDS")
v2_cases <- subset(v2_pheno, scd_age_filt == 1)
v2_controls <- subset(v2_pheno, scd_age_filt == 0)
visit2 <- as.data.frame(matrix(nrow = 19, ncol = 4))
names(visit2) <- c("trait", "all", "controls", "cases")
visit2$trait <- c("n", "cases", "age", "black", "female", "height", "bmi", "curr_smoke", "egfr", "hf", "chd", "stroke", "diabetes", "afib", "sbp", "hnt_med", "total_chol", "tri", "hdl")

visit2$all <- c(nrow(v2_pheno),
                npercalc(v2_pheno, "scd_age_filt", "1"),
                meansdcalc(v2_pheno, "age"),
                npercalc(v2_pheno, "race", "Black"),
                npercalc(v2_pheno, "sex", "F"),
                meansdcalc(v2_pheno, "height"),
                meansdcalc(v2_pheno, "bmi"),
                npercalc(v2_pheno, "smoke", "C"),
                meansdcalc(v2_pheno, "eGFR_ckdepi"),
                npercalc(v2_pheno, "hf", "Y"),
                npercalc(v2_pheno, "chd_chs", "Y"),
                npercalc(v2_pheno, "stroke", "Y"),
                npercalc(v2_pheno, "diabetes", "Y"),
                npercalc(v2_pheno, "afib", "Y"),
                meansdcalc(v2_pheno, "sbp"),
                npercalc(v2_pheno, "hnt_med", "Y"),
                meansdcalc(v2_pheno, "total_chol"),
                meansdcalc(v2_pheno, "tri"),
                meansdcalc(v2_pheno, "hdl"))

visit2$cases<- c(nrow(v2_cases),
                 npercalc(v2_cases, "scd_age_filt", "1"),
                 meansdcalc(v2_cases, "age"),
                 npercalc(v2_cases, "race", "Black"),
                 npercalc(v2_cases, "sex", "F"),
                 meansdcalc(v2_cases, "height"),
                 meansdcalc(v2_cases, "bmi"),
                 npercalc(v2_cases, "smoke", "C"),
                 meansdcalc(v2_cases, "eGFR_ckdepi"),
                 npercalc(v2_cases, "hf", "Y"),
                 npercalc(v2_cases, "chd_chs", "Y"),
                 npercalc(v2_cases, "stroke", "Y"),
                 npercalc(v2_cases, "diabetes", "Y"),
                 npercalc(v2_cases, "afib", "Y"),
                 meansdcalc(v2_cases, "sbp"),
                 npercalc(v2_cases, "hnt_med", "Y"),
                 meansdcalc(v2_cases, "total_chol"),
                 meansdcalc(v2_cases, "tri"),
                 meansdcalc(v2_cases, "hdl"))

visit2$controls<- c(nrow(v2_controls),
                    npercalc(v2_controls, "scd_age_filt", "0"),
                    meansdcalc(v2_controls, "age"),
                    npercalc(v2_controls, "race", "Black"),
                    npercalc(v2_controls, "sex", "F"),
                    meansdcalc(v2_controls, "height"),
                    meansdcalc(v2_controls, "bmi"),
                    npercalc(v2_controls, "smoke", "C"),
                    meansdcalc(v2_controls, "eGFR_ckdepi"),
                    npercalc(v2_controls, "hf", "Y"),
                    npercalc(v2_controls, "chd_chs", "Y"),
                    npercalc(v2_controls, "stroke", "Y"),
                    npercalc(v2_controls, "diabetes", "Y"),
                    npercalc(v2_controls, "afib", "Y"),
                    meansdcalc(v2_controls, "sbp"),
                    npercalc(v2_controls, "hnt_med", "Y"),
                    meansdcalc(v2_controls, "total_chol"),
                    meansdcalc(v2_controls, "tri"),
                    meansdcalc(v2_controls, "hdl"))

####ARIC Visit 3 Table####
v3_pheno <- readRDS("aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit3/combined/exclude.egfr/with.afib/association.stroke.covar.af.rskfct.chs.harm.hrt.angio/rds.files/ARIC.SCA.V3.Final_Pheno.RDS")

v3_cases <- subset(v3_pheno, scd_age_filt == 1)
v3_controls <- subset(v3_pheno, scd_age_filt == 0)
visit3 <- as.data.frame(matrix(nrow = 19, ncol = 4))
names(visit3) <- c("trait", "all", "controls", "cases")
visit3$trait <- c("n", "cases", "age", "black", "female", "height", "bmi", "curr_smoke", "egfr", "hf", "chd", "stroke", "diabetes", "afib" ,"sbp", "hnt_med", "total_chol", "tri", "hdl")

visit3$all <- c(nrow(v3_pheno),
                npercalc(v3_pheno, "scd_age_filt", "1"),
                meansdcalc(v3_pheno, "age"),
                npercalc(v3_pheno, "race", "Black"),
                npercalc(v3_pheno, "sex", "F"),
                meansdcalc(v3_pheno, "height"),
                meansdcalc(v3_pheno, "bmi"),
                npercalc(v3_pheno, "smoke", "C"),
                meansdcalc(v3_pheno, "eGFR_ckdepi"),
                npercalc(v3_pheno, "hf", "Y"),
                npercalc(v3_pheno, "chd_chs", "Y"),
                npercalc(v3_pheno, "stroke", "Y"),
                npercalc(v3_pheno, "diabetes", "Y"),
                npercalc(v3_pheno, "afib", "Y"),
                meansdcalc(v3_pheno, "sbp"),
                npercalc(v3_pheno, "hnt_med", "Y"),
                meansdcalc(v3_pheno, "total_chol"),
                meansdcalc(v3_pheno, "tri"),
                meansdcalc(v3_pheno, "hdl"))

visit3$cases <- c(nrow(v3_cases),
                  npercalc(v3_cases, "scd_age_filt", "1"),
                  meansdcalc(v3_cases, "age"),
                  npercalc(v3_cases, "race", "Black"),
                  npercalc(v3_cases, "sex", "F"),
                  meansdcalc(v3_cases, "height"),
                  meansdcalc(v3_cases, "bmi"),
                  npercalc(v3_cases, "smoke", "C"),
                  meansdcalc(v3_cases, "eGFR_ckdepi"),
                  npercalc(v3_cases, "hf", "Y"),
                  npercalc(v3_cases, "chd_chs", "Y"),
                  npercalc(v3_cases, "stroke", "Y"),
                  npercalc(v3_cases, "diabetes", "Y"),
                  npercalc(v3_cases, "afib", "Y"),
                  meansdcalc(v3_cases, "sbp"),
                  npercalc(v3_cases, "hnt_med", "Y"),
                  meansdcalc(v3_cases, "total_chol"),
                  meansdcalc(v3_cases, "tri"),
                  meansdcalc(v3_cases, "hdl"))

visit3$controls <- c(nrow(v3_controls),
                     npercalc(v3_controls, "scd_age_filt", "0"),
                     meansdcalc(v3_controls, "age"),
                     npercalc(v3_controls, "race", "Black"),
                     npercalc(v3_controls, "sex", "F"),
                     meansdcalc(v3_controls, "height"),
                     meansdcalc(v3_controls, "bmi"),
                     npercalc(v3_controls, "smoke", "C"),
                     meansdcalc(v3_controls, "eGFR_ckdepi"),
                     npercalc(v3_controls, "hf", "Y"),
                     npercalc(v3_controls, "chd_chs", "Y"),
                     npercalc(v3_controls, "stroke", "Y"),
                     npercalc(v3_controls, "diabetes", "Y"),
                     npercalc(v3_controls, "afib", "Y"),
                     meansdcalc(v3_controls, "sbp"),
                     npercalc(v3_controls, "hnt_med", "Y"),
                     meansdcalc(v3_controls, "total_chol"),
                     meansdcalc(v3_controls, "tri"),
                     meansdcalc(v3_controls, "hdl"))

####CHS Year 5 Table####
chs.stats <- read.csv("aric.prot.real.archive/scd.assoc/4_CHS.results/10032022/Table_1_2022_10_04.csv", header = T)
names(chs.stats)[2] <- "All"

chs.table <- visit2
chs.table$all <- NA
chs.table$controls <- NA
chs.table$cases <- NA
chs.table <- chs.table[, c("trait","all", "cases", "controls")]

chs.nper <- function(sub,total){
  percent <- (sub %>% as.numeric()/total %>% as.numeric())*100
  percent <- round(percent, digits = 1)
  n <- sub
  output <- paste0(n, " (", percent, "%)")
  return(output)
}

#n
chs.table[1,2:4] <- chs.stats[1,c(2:4)]
#cases
chs.table[2,2:4] <- c(chs.nper(chs.stats$Cases[1], chs.stats$All[1]),NA,NA)
#age
chs.table[3,2:4] <- chs.stats[16,2:4]
#black
chs.table[4,2:4] <- c(chs.nper(chs.stats$All[5], chs.stats$All[1]),
                      chs.nper(chs.stats$Cases[5], chs.stats$Cases[1]),
                      chs.nper(chs.stats$Controls[5], chs.stats$Controls[1]))
#female
chs.table[5,2:4] <- c(chs.nper(chs.stats$All[4], chs.stats$All[1]),
                      chs.nper(chs.stats$Cases[4], chs.stats$Cases[1]),
                      chs.nper(chs.stats$Controls[4], chs.stats$Controls[1]))
#height
chs.table[6,2:4] <- chs.stats[17,2:4]
#bmi
chs.table[7,2:4] <- chs.stats[18,2:4]
#current smoke
chs.table[8,2:4] <- c(chs.nper(chs.stats$All[9], chs.stats$All[1]),
                      chs.nper(chs.stats$Cases[9], chs.stats$Cases[1]),
                      chs.nper(chs.stats$Controls[9], chs.stats$Controls[1]))
#egfr
chs.table[9,2:4] <- chs.stats[19,2:4]
#hf
chs.table[10,2:4] <- c(chs.nper(chs.stats$All[15], chs.stats$All[1]),
                      chs.nper(chs.stats$Cases[15], chs.stats$Cases[1]),
                      chs.nper(chs.stats$Controls[15], chs.stats$Controls[1]))
#chd
chs.table[11,2:4] <- c(chs.nper(chs.stats$All[10], chs.stats$All[1]),
                      chs.nper(chs.stats$Cases[10], chs.stats$Cases[1]),
                      chs.nper(chs.stats$Controls[10], chs.stats$Controls[1]))
#stroke
chs.table[12,2:4] <- c(chs.nper(chs.stats$All[13], chs.stats$All[1]),
                       chs.nper(chs.stats$Cases[13], chs.stats$Cases[1]),
                       chs.nper(chs.stats$Controls[13], chs.stats$Controls[1]))
#diabetes
chs.table[13,2:4] <- c(chs.nper(chs.stats$All[12], chs.stats$All[1]),
                       chs.nper(chs.stats$Cases[12], chs.stats$Cases[1]),
                       chs.nper(chs.stats$Controls[12], chs.stats$Controls[1]))
#afib
chs.table[14,2:4] <- c(chs.nper(chs.stats$All[14], chs.stats$All[1]),
                       chs.nper(chs.stats$Cases[14], chs.stats$Cases[1]),
                       chs.nper(chs.stats$Controls[14], chs.stats$Controls[1]))
#sbp
chs.table[15,2:4] <- chs.stats[20,2:4]
#hnt_med
chs.table[16,2:4] <- c(chs.nper(chs.stats$All[11], chs.stats$All[1]),
                       chs.nper(chs.stats$Cases[11], chs.stats$Cases[1]),
                       chs.nper(chs.stats$Controls[11], chs.stats$Controls[1]))
#total_chol
chs.table[17,2:4] <- chs.stats[21,2:4]
#tri
chs.table[18,2:4] <- chs.stats[22,2:4]
#hdl
chs.table[19,2:4] <- chs.stats[23,2:4]

####Convert to SI units####
#Formula source: https://www.ncbi.nlm.nih.gov/books/NBK83505/

## Means
convert.means <- chs.table[17:19,]
convert.means$all[c(1,3)] <- gsub("[(].*", "", convert.means$all[c(1,3)]) %>% as.numeric(.)/38.67
convert.means$cases[c(1,3)] <- gsub("[(].*", "", convert.means$cases[c(1,3)]) %>% as.numeric(.)/38.67
convert.means$controls[c(1,3)] <- gsub("[(].*", "", convert.means$controls[c(1,3)]) %>% as.numeric(.)/38.67

convert.means$all[2] <- gsub("[(].*", "", convert.means$all[2]) %>% as.numeric(.)/88.57
convert.means$cases[2] <- gsub("[(].*", "", convert.means$cases[2]) %>% as.numeric(.)/88.57
convert.means$controls[2] <- gsub("[(].*", "", convert.means$controls[2]) %>% as.numeric(.)/88.57

convert.means$all <- as.numeric(convert.means$all)
convert.means$cases <- as.numeric(convert.means$cases)
convert.means$controls <- as.numeric(convert.means$controls)

convert.means[,c(2:4)] <- round(convert.means[,c(2:4)], digits = 1)

## SE
convert.se <- chs.table[17:19,]
convert.se$all[c(1,3)] <- gsub(".*[(]", "", convert.se$all[c(1,3)]) %>% gsub(")", "", .) %>% as.numeric(.)/38.67
convert.se$cases[c(1,3)] <- gsub(".*[(]", "", convert.se$cases[c(1,3)]) %>% gsub(")", "", .)%>% as.numeric(.)/38.67
convert.se$controls[c(1,3)] <- gsub(".*[(]", "", convert.se$controls[c(1,3)]) %>% gsub(")", "", .) %>% as.numeric(.)/38.67

convert.se$all[2] <- gsub(".*[(]", "", convert.se$all[2]) %>% gsub(")", "", .) %>% as.numeric(.)/88.57
convert.se$cases[2] <- gsub(".*[(]", "", convert.se$cases[2]) %>% gsub(")", "", .) %>% as.numeric(.)/88.57
convert.se$controls[2] <- gsub(".*[(]", "", convert.se$controls[2]) %>% gsub(")", "", .) %>% as.numeric(.)/88.57

convert.se$all <- as.numeric(convert.se$all)
convert.se$cases <- as.numeric(convert.se$cases)
convert.se$controls <- as.numeric(convert.se$controls)

convert.se[,c(2:4)] <- round(convert.se[,c(2:4)], digits = 1)

## Paste together
chs.table[17:19,2] <- paste0(convert.means$all, " (", convert.se$all, ")")
chs.table[17:19,3] <- paste0(convert.means$cases, " (", convert.se$cases, ")")
chs.table[17:19,4] <- paste0(convert.means$controls, " (", convert.se$controls, ")")

## Reorganize
chs.table <- chs.table[, c("trait", "all", "controls", "cases")]

####Final table####
stats.table <- merge(visit2, visit3, by = "trait", sort = "F")
stats.table <- merge(stats.table, chs.table, by = "trait", sort = "F")
names(stats.table) <- c("trait", "v2_all", "v2_controls", "v2_cases", "v3_all", "v3_controls", "v3_cases", "chs_all", "chs_controls", "chs_cases")
write.table(stats.table, file = "aric.prot.real.archive/scd.assoc/6_manuscript/ver13/Cohort.Characteristics/Cohort.Characteristics.txt", row.names = F, col.names = T, sep = "\t", quote = F)
#Base data for Panel B in ver19

save.image("aric.prot.real.archive/scd.assoc/6_manuscript/ver13/Cohort.Characteristics/Cohort.Characteristics.Source.RData")