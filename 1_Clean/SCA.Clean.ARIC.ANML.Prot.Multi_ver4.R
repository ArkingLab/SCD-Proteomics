####Description####
## Script used to prepare ARIC ANML V2 or V3 proteomics data to run association with SCA. Script will:
## 1) Merge phenotype and proteomics data [option to analyze V2 or V3]
## 2) Apply exclusions [option to change exclusion criteria based on mi_hf, cvd, diabetes, eGFR < 30]
## 3) Run PCA to remove outliers
## 4) Make final dataset for CoxPH models

## Changes from ver2 to ver3: added option to run association or not, added commands to make postStrat folder, fixed issue where cvd wasn't marked correctly in exam 3 phenotype, 
##                    added option to create/run stratified datasets from the combined dataset, added commands to make sex stratified combined race datasets, automate appendix creation

## Changes from ver3 to ver4: took out stroke as a covariate, added in AF phenotype

####Options####
## race_strat - White, Black or Both; automatic from file path
## sex_strat - M, F, or Both; automatic from file path
## exclude_only - name of column in phenotype file, separated by + [eGFR_ckdepi+cvd+mi_hf]
## density_on - Yes or No
## egfr_crit - eGFR exclusion criteria
## appendix - text to append to output files
## assoc_cont - Yes or No [No will save .rds files but not run association script]
## postStrat - Yes or No [No will not create statified .rds files from the combined (sex_strat & race_strat = Both) datasets]
## visit - 2 or 3; automatic from file path

####Load packages####
library(ggplot2)
library(magrittr)
library(gridExtra)
library(RColorBrewer)
library(data.table)
library(readxl)
#library(ggfortify)
#library(lattice)
#library(cowplot)
print("Finish loading packages")

####Load functions####
num.prot <- function(data){
  grep("Seq", names(data)) %>% length()
}

## From CC, color PCs by covariates
panel.cor<-function(x,y,digits=2,prefix="",cex.cor,...){
  usr<-par("usr");on.exit(par(usr))
  par(usr=c(0,1,0,1))
  r<-abs(cor(x,y, use = "complete.obs"))
  txt<-format(c(r,0.123456789),digits=digits)[1]
  txt<-paste(prefix,txt,sep="")
  if(missing(cex.cor))cex.cor<-0.8/strwidth(txt)
  text(0.5,0.5,txt,cex=cex.cor*r)
}

plot.pairs <- function(pca.data, samp.data, samp.id.col, interest, type, title){
  
  #' @param pca.data - prcomp output
  #' @param samp.data - sample/subject dataframe
  #' @param samp.id.col - name of column with sample ID
  #' @param interest - variable to color pairs plot by
  #' @param type - class of variable of interest
  #' @param title - title to add to plot
  
  ## Merge data into one large frame to plot by
  pca.matrix <- pca.data$x %>% as.data.frame()
  pca.matrix[,samp.id.col] <- row.names(pca.matrix)
  samp.sub <- samp.data[, c(samp.id.col, interest)]
  merged <- merge(samp.sub, pca.matrix, by.x = samp.id.col, by.y = samp.id.col)
  row.names(merged) <- merged[, samp.id.col]
  
  if(any(is.na(merged[, interest]))){
    merged <- subset(merged, !is.na(merged[, interest]))
  }
  
  if(type == "numeric"){
    merged[, interest] <- as.numeric(merged[, interest])
    ## Assign color
    #color_pal <- colorRampPalette(brewer.pal(9, "Blues"))
    color_pal <- colorRampPalette(c("red","blue"))
    color_unique <- color_pal((length(unique(merged[, interest]))))
    unique_values <- unique(merged[, interest]) %>% sort()
    color <- merged[, interest]
    color <- ifelse(color == unique_values[1], color_unique[1], color)
    for (i in 2:length(unique_values)){
      color <- ifelse(color == unique_values[i], color_unique[i], color)
    }
    
    ## Plot
    pairs(~merged$PC1+merged$PC2+merged$PC3+merged$PC4+merged$PC5+merged$PC6+merged$PC7+merged$PC8+merged$PC9+merged$PC10+merged[, interest], labels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", interest), col = color, upper.panel = panel.cor, pch = 20, main = paste("n=", nrow(merged), title, sep = " "), oma = c(11, 2, 5, 2))
  }
  else if(type == "factor"){
    #merged[, interest] <- as.factor(merged[, interest])
    ## Assign colors
    color <- merged[,interest]
    color <- as.character(color)
    var_levels <- unique(merged[,interest])
    palette <- brewer.pal(8, "Dark2")
    color <- ifelse(color == var_levels[1], palette[1], color)
    if (length(var_levels) <= 8) {
      for (i in 2:length(var_levels)){
        color <- ifelse(color == var_levels[i], palette[i], color)
      }
      ## Plot
      pairs(~merged$PC1+merged$PC2+merged$PC3+merged$PC4+merged$PC5+merged$PC6+merged$PC7+merged$PC8+merged$PC9+merged$PC10+as.factor(merged[, interest]), labels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", interest), col = color, upper.panel = panel.cor, pch = 20, main = paste(nrow(merged), title, sep = " "), oma = c(11, 2, 5, 2))
      legend("bottom", fill = palette[1:length(var_levels)], legend = var_levels, box.lty = 0, horiz = TRUE, cex = 0.8)  
    }
    else if(length(var_levels>8)){
      print("More than 8 levels in factor! Need to manually plot")
    }    
  }
  #table(color)
}

find.PCA.outliers <- function(data, sd_num){
  #' @param data - prcomp object
  #' @param sd_num - number of SDs outliers should be discovered by
  
  sample_outliers=c()
  all_outliers=c()
  
  for(i in 1:10){
    a<-subset(rownames(data$x), data$x[,i] > (mean(data$x[,i])+sd_num*sd(data$x[,i])))
    b<-subset(rownames(data$x), data$x[,i] < (mean(data$x[,i])-sd_num*sd(data$x[,i])))
    out<-c(a,b)
    sample_outliers <- c(sample_outliers,out)
    all_outliers=c(all_outliers,sample_outliers)
    sample_outliers=c()
  }
  
  outlier<-unique(all_outliers)
  return(outlier)
}
print("Finish loading functions")

####Arguments####
args <- commandArgs(trailingOnly=TRUE)
exclude_only <- args[1]
density_on <- args[2]
egfr_crit <- args[3] %>% as.numeric()
#appendix <- args[4]

## Add option to run association or just make dataset
assoc_cont <- args[4]

## Add option to stratify from combined dataset
postStrat <- args[5]

dir_out <- paste0(getwd(), "/")
if(grepl("combined", dir_out)){
  race_strat <- "Both"
  race_app <- "BthRce"
}else if(grepl("AA", dir_out)){
  race_strat <- "Black"
  race_app <- "AA"
} else if(grepl("EA", dir_out)){
  race_strat <- "White"
  race_app <- "EA"
}
if(grepl("female", dir_out)){
  sex_strat <- "F"
  sex_app <- "Fem"
}else if(grepl("male", dir_out)){
  sex_strat <- "M"
  sex_app <- "Male"
} else {
  sex_strat <- "Both"
  sex_app <- "BthSx"
  if(race_strat == "Both" & sex_strat == "Both"){
    sex_app <- "Sx"
  }
}
if(grepl("visit2", dir_out)){
  visit <- 2
}else if(grepl("visit3", dir_out)){
  visit <- 3
}
appendix <- paste0(race_app, sex_app)

print(paste0("Files will be saved in ", dir_out))
print(paste0("SCA Analysis to run: ", race_strat, " race(s), ", sex_strat, " sex(es), ", "exclude ", exclude_only, " criteria, ", "Visit ", visit))

####Load data####

## Protein measurements
load("aric.prot.real.archive/ecg.assoc/1_data.cleaning/ANML/ARIC.V2V3V5.ANML.Prot.NoFlag2.rds")
if(visit == 2){
  rm(v3data.nof2, v5data.nof2)
}else if(visit == 3){
  rm(v5data.nof2, v2data.nof2)
}
print("Finish loading ANML protein data")

## Edit protein annotation
if(visit == 2){
  full.annot <- fread("aric.anml.prot/Annotation/soma_visit_2_annot_ANML_SMP_updated.txt") %>% as.data.frame()
} else if(visit == 3){
  full.annot <- fread("aric.anml.prot/Annotation/soma_visit_3_annot_ANML_SMP_updated.txt") %>% as.data.frame()
}
abbrv.annot <- read_xlsx("aric.anml.prot/Annotation/Protein Dictionary and Abbreviated Annotation-soma visit 2&3&5.xlsx", sheet = "Annotation") %>% as.data.frame()
annot.final <- merge(abbrv.annot, annot.nof2[, c("seqid_in_sample", "multi.apt")], by = "seqid_in_sample")

## Phenotype file used for ECG omics for covariates
ecg_pheno <- fread("ecg.omics.archive/ARIC.ECG/pheno/02152022_Pheno/ARIC_ECG_Omics_TVD_02152022_Phenotype.txt") %>% as.data.frame()
dz <- ecg_pheno[, c("aric_iid", "exam", "visit_mi", "visit_hf", "visit_stroke", "visit_tia", "visit_chd")]
names(ecg_pheno)[ncol(ecg_pheno)] <- "ID"
ecg_pheno_visit <- subset(ecg_pheno, exam == visit)

## Add in separated CHD, stroke, MI variables
exam1 <- subset(dz, exam == 1)
miv1 <- subset(exam1, visit_mi == 1)$aric_iid
hfv1 <- subset(exam1, visit_hf == 1)$aric_iid
strokev1 <- subset(exam1, visit_stroke == 1)$aric_iid
tiav1 <- subset(exam1, visit_tia == 1)$aric_iid
chdv1 <- subset(exam1, visit_chd == 1)$aric_iid

exam2 <- subset(dz, exam == 2)
exam2$mi <- ifelse(exam2$visit_mi == 1 | exam2$aric_iid %in% miv1, "Y", "N")
exam2$hf <- ifelse(exam2$visit_hf == 1 | exam2$aric_iid %in% hfv1, "Y", "N")
exam2$stroke <- ifelse(exam2$visit_stroke == 1 | exam2$aric_iid %in% strokev1, "Y", "N")
exam2$tia <- NA
exam2$tia <- ifelse(exam2$visit_tia == 1 & !is.na(exam2$visit_tia), "Y", exam2$tia)
exam2$tia <- ifelse(exam2$visit_tia == 0 & !is.na(exam2$visit_tia), "N", exam2$tia)
exam2$tia <- ifelse(exam2$aric_iid %in% tiav1, 1, exam2$tia)
exam2$chd <- ifelse(exam2$visit_chd == 1 | exam2$aric_iid %in% chdv1, "Y", "N")

miv12 <- c(subset(exam2, mi == "Y")$aric_iid, miv1) %>% unique()
hfv12 <- c(subset(exam2, hf == "Y")$aric_iid, hfv1) %>% unique()
strokev12 <- c(subset(exam2, stroke == "Y")$aric_iid, strokev1) %>% unique()
tiav12 <- c(subset(exam2, tia == "Y")$aric_iid, tiav1) %>% unique()
chdv12 <- c(subset(exam2, chd == "Y")$aric_iid, chdv1) %>% unique()

exam3 <- subset(dz, exam == 3)
exam3$mi <- ifelse(exam3$visit_mi == 1 | exam3$aric_iid %in% miv12, "Y", "N")
exam3$hf <- ifelse(exam3$visit_hf == 1 | exam3$aric_iid %in% hfv12, "Y", "N")
exam3$stroke <- ifelse(exam3$visit_stroke == 1 | exam3$aric_iid %in% strokev12, "Y", "N")
exam3$tia <- NA
exam3$tia <- ifelse(exam3$visit_tia == 1 & !is.na(exam3$visit_tia), "Y", exam3$tia)
exam3$tia <- ifelse(exam3$visit_tia == 0 & !is.na(exam3$visit_tia), "N", exam3$tia)
exam3$tia <- ifelse(exam3$aric_iid %in% tiav12, 1, exam3$tia)
exam3$chd <- ifelse(exam3$visit_chd == 1 | exam3$aric_iid %in% chdv12, "Y", "N")

dz_extra <- as.data.frame(rbind(exam2[, c("aric_iid", "exam", "mi", "hf", "stroke", "tia", "chd")], exam3[, c("aric_iid", "exam", "mi", "hf", "stroke", "tia", "chd")]))
dz_extra_exam <- dz_extra[dz_extra$exam == visit,]
dz_extra_exam$exam <- NULL
names(dz_extra_exam)[1] <- "ID"
ecg_pheno_visit <- merge(ecg_pheno_visit, dz_extra_exam, by = "ID")
ecg_pheno_visit$chd_nomi <- ifelse(ecg_pheno_visit$mi == "N" & ecg_pheno_visit$chd == "Y", "Y", "N")

## AF phenotype
af_pheno <- fread("aric.prot.real.archive/af.assoc/1_data.cleaning/ARIC_V1-5_AF.txt") %>% as.data.frame()
if(visit == 2){
  af_pheno <- af_pheno[, c("ID", "prevafv2")]
}else if(visit == 3){
  af_pheno <- af_pheno[, c("ID", "prevafv3")]
}
names(af_pheno) <- c("ID", "afib")
ecg_pheno_visit <- merge(ecg_pheno_visit, af_pheno, by = "ID")

## SCA phenotype
if(visit == 2){
  scd_pheno <- fread("aric.prot.real.archive/scd.assoc/1_data.cleaning/ARIC_SCD_Pheno_For_Proteomics_From_V2_Age_5YR_10YR_Censored.txt") %>% as.data.frame()
  names(scd_pheno) <- gsub("v2_", "", names(scd_pheno))
  names(scd_pheno) <- gsub("v2", "", names(scd_pheno))
} else if(visit == 3){
  scd_pheno <- fread("aric.prot.real.archive/scd.assoc/1_data.cleaning/ARIC_SCD_Pheno_For_Proteomics_From_V3_Age_5YR_10YR_Censored.txt") %>% as.data.frame()
  names(scd_pheno) <- gsub("v3_", "", names(scd_pheno))
  names(scd_pheno) <- gsub("v3", "", names(scd_pheno))
}

####Edit data####

## Merge phenotypes
pheno <- merge(ecg_pheno_visit, scd_pheno, by = "ID")

## Merge to protein
if(visit == 2){
  alldata <- merge(pheno, v2data.nof2, by.x = "ID", by.y = "SampleId")
} else if(visit == 3){
  alldata <- merge(pheno, v3data.nof2, by.x = "ID", by.y = "SampleId")
}

## Remove non-whites and non-blacks
alldata <- alldata[alldata$race == "Black" | alldata$race == "White",]

## Remove mismatch center & race
#In previous iteration, most were removed due to Geneva ID exclusion
center.mismatch <- c(alldata[alldata$race == "Black" & alldata$center == "M","ID"], alldata[alldata$race == "Black" & alldata$center == "W","ID"])
alldata <- alldata[!(alldata$ID %in% center.mismatch),]

## Remove those without Geneva IDs
#Geneva ID is submitted dbGAP ID, also ID submitted for genotyping
#Worried that there is an issue with these samples and that's why they weren't used for genotyping so remove
alldata <- alldata[!(is.na(alldata$subject_id)),]

## Remove samples with missing covariates
#covariates <- c("height", "bmi", "age", "sex", "race", "smoke", "eGFR_ckdepi", "cvd", "diabetes", "center", "sbp", "hnt_med", "total_chol", "hdl", "tri", "stroke", "mi", "hf", "chd", "chd_nomi")
#covariates <- c("height", "bmi", "age", "sex", "race", "smoke", "eGFR_ckdepi", "diabetes", "center", "sbp", "hnt_med", "total_chol", "hdl", "tri", "stroke", "hf", "chd")
covariates <- c("height", "bmi", "age", "sex", "race", "smoke", "eGFR_ckdepi", "cvd", "diabetes", "center", "sbp", "hnt_med", "total_chol", "hdl", "tri", "mi", "hf", "chd", "afib")
check.covar <- alldata[, c("ID", covariates)] %>% unique()
check.covar$sum <- rowSums(is.na(check.covar[, 2:ncol(check.covar)]))
missing.covar <- check.covar[check.covar$sum != 0, "ID"]
alldata <- alldata[!(alldata$ID %in% missing.covar),]
#Other lipids: ldl

## Remove smoking F or C
alldata <- alldata[!(alldata$smoke == "F or C"),]

## Remove missing SCD
alldata <- alldata[!(is.na(alldata$scd_age_filt)),]
print("Finish merging phenotype and protein data, excluding those missing covariates")
print(paste0("Before exclusions: ", nrow(alldata), " people and ", num.prot(alldata), " proteins"))

####Remove cohort overlap####
overlap <- read.table("aric.prot.real.archive/ARIC.JHS.CHS.MESA.Overlap.txt", header = T)
## Only remove CHS [cohort with SCD phenotype]
chs.tbr <- overlap[overlap$other.cohort == "CHS", "pid"]
alldata <- alldata[!(alldata$ID %in% chs.tbr),]
print("Finish removing overlapping samples")
print(paste0("Before exclusions, after overlap: ", nrow(alldata), " people and ", num.prot(alldata), " proteins"))

####Phenotype exclusions####
alldata.og <- alldata

exclude_vars <- strsplit(exclude_only, "[+]") %>% as.data.frame()
print("Exclusion variables:")
print(exclude_vars[,1])
for(i in 1:nrow(exclude_vars)){
  var <- exclude_vars[i,1]
  if(var == "eGFR_ckdepi"){
    print(paste0("eGFR exclusion criteria: ", egfr_crit))
    alldata$eGFR_ckdepi_exclude <- ifelse(alldata$eGFR_ckdepi < egfr_crit, 1, 0)
    print(var)
    print(table(alldata[,"eGFR_ckdepi_exclude"]))
  } else if(var == "none"){
    break()
    print("No exclusions will be applied")
  } else if(var != "eGFR_ckdepi"){
    new.col <- paste0(var, "_exclude")
    alldata[,new.col] <- ifelse(alldata[,var] == "Y", 1, 0)
    print(var)
    print(table(alldata[,new.col]))
  }
}
if(exclude_vars[1,1] != "none"){
  exclude_cols_names <- names(alldata)[grep("exclude", names(alldata))]
  exclude_cols <- grep("exclude", names(alldata))

  if(length(exclude_cols) > 1){
    alldata$exclude <- rowSums(alldata[, exclude_cols])
  } else if(length(exclude_cols) == 1){
    alldata$exclude <- alldata[, exclude_cols]
  }
  alldata <- alldata[alldata$exclude == 0,]
  
  col_tbr <- grep("_exclude", names(alldata))
  alldata <- alldata[,-col_tbr]
  #for(i in 1:length(exclude_cols)){
  #  alldata[,exclude_cols[i]] <- NULL
  #}
}

####Sex and/or race filtering####
if(race_strat == "Both"){
    print("Running both EA and AA")
    print(paste0("Before exclusion, total: ", alldata.og %>% nrow(), "; cases: ", alldata.og[alldata.og$scd_age_filt == 1,] %>% nrow(), "; ctrls: ", alldata.og[alldata.og$scd_age_filt == 0,] %>% nrow()))
} else if(race_strat == "White" | race_strat == "Black"){
    alldata <- subset(alldata, race == race_strat)
    print(paste0("Running ", race_strat, " only"))
    print(paste0("Before exclusion, total: ", alldata.og[alldata.og$race == race_strat,] %>% nrow(), "; cases: ", alldata.og[alldata.og$scd_age_filt == 1 & alldata.og$race == race_strat,] %>% nrow(), "; ctrls: ", alldata.og[alldata.og$scd_age_filt == 0 & alldata.og$race == race_strat,] %>% nrow()))
}

if(sex_strat == "Both"){
  print("Running both sexes")
  print(paste0("Before exclusion, total: ", alldata.og %>% nrow(), "; cases: ", alldata.og[alldata.og$scd_age_filt == 1,] %>% nrow(), "; ctrls: ", alldata.og[alldata.og$scd_age_filt == 0,] %>% nrow()))
} else if(sex_strat == "M" | sex_strat == "F"){
  alldata <- subset(alldata, sex == sex_strat)
  print(paste0("Running ", sex_strat, " only"))
  print(paste0("Before exclusion, total: ", alldata.og[alldata.og$sex == sex_strat,] %>% nrow(), "; cases: ", alldata.og[alldata.og$scd_age_filt == 1 & alldata.og$sex == sex_strat,] %>% nrow(), "; ctrls: ", alldata.og[alldata.og$scd_age_filt == 0 & alldata.og$sex == sex_strat,] %>% nrow()))
}

print("Finish applying phenotype, sex (if applicable), race (if applicable) exclusions")
print(paste0("After all exclusions, after overlap: ", nrow(alldata), " people and ", num.prot(alldata), " proteins; ", alldata[alldata$scd_age_filt == 1,] %>% nrow(), " cases and ", alldata[alldata$scd_age_filt == 0,] %>% nrow(), " controls"))

####Density plots####
start=grep("Seq", names(alldata)) %>% head() %>% .[1]
last=grep("Seq", names(alldata)) %>% tail() %>% .[6]

plot_list = list()
for (i in start:last){
  prot.name <- subset(annot.final, seqid_in_sample == names(alldata)[i])$target
  seq <- subset(annot.final, seqid_in_sample == names(alldata)[i])$seqid_in_sample %>% gsub("SeqId_", "", .)
  p = ggplot(alldata, aes_string(x=names(alldata)[i])) + geom_density(fill="lightblue") + theme_bw() + xlab(prot.name) + ylab(seq)
  plot_list[[i]] = p
}

if(density_on == "Yes"){
  pdf(file = paste0(dir_out, "Protein_Density.pdf"), width = 20, height = 15)
  plot_num <- seq(start, length(plot_list), 25)
  last_plot_num <- plot_num[length(plot_num)]
  plot_num <- plot_num[-length(plot_num)]
  ith <- plot_num[seq(1, length(plot_num), 14)]
  for (z in plot_num){
    if(z %in% ith){
      print(paste0("Plot #", (z-start)+1))
    }
    do.call(grid.arrange, list(grobs = plot_list[z:(z+24)], ncol=5, top=paste0("ARIC n=", nrow(alldata), ", p=", num.prot(alldata))))
  }
  do.call(grid.arrange, list(grobs = plot_list[last_plot_num:length(plot_list)], ncol=2, top=paste0("ARIC n=", nrow(alldata), ", p=", num.prot(alldata))))
  dev.off()
} else {
  print("Skipping density plot")
}

## Standardize
alldata[, c(start:last)] <- scale(alldata[, start:last])
save(alldata, annot.final, full.annot, file = paste0(dir_out, "PrePCA.rds"))

####PCA - Iteration 1####
names(alldata)[grep("^ID$", names(alldata))] <- "pid"
row.names(alldata) <- alldata$pid
print("Running PCA iteration #1")
#Data is already scaled and centered
pca.it1 <- prcomp(alldata[, start:last], scale = F, center = F)
pca.it1.outliers <- find.PCA.outliers(pca.it1, 4)
print("Finish PCA iteration #1")
print(paste0("Number of outliers: ", length(pca.it1.outliers)))

## Skree plot and color by covariates
pdf(file = paste0(dir_out, "PCA_It1.pdf"), width = 10, height = 10)
par(cex = 1.1)
plot(summary(pca.it1)$importance[2,1:10]*100, ylab="% Variance Explained", xlab="PC", type="b", main = paste0("PCA for n=", nrow(alldata), "; p=", num.prot(alldata)), pch = 19)
axis(1, at = c(1:10))
if(sex_strat == "Both"){
  plot.pairs(pca.it1, alldata, "pid", "sex", "factor", "; p=4955")
} else if(sex_strat == "F" | sex_strat == "M"){
  print("No sex plot")
}
plot.pairs(pca.it1, alldata, "pid", "center", "factor", "; p=4955")
if(race_strat == "Both"){
  plot.pairs(pca.it1, alldata, "pid", "race", "factor", "; p=4955")
} else if(race_strat == "White" | race_strat == "Black"){
  print("No race plot")
}
plot.pairs(pca.it1, alldata, "pid", "age", "numeric", "; p=4955")
plot.pairs(pca.it1, alldata, "pid", "eGFR_ckdepi", "numeric", "; p=4955")
plot.pairs(pca.it1, alldata, "pid", "smoke", "factor", "; p=4955")
plot.pairs(pca.it1, alldata, "pid", "bmi", "numeric", "; p=4955")
plot.pairs(pca.it1, alldata, "pid", "height", "numeric", "; p=4955")
plot.pairs(pca.it1, alldata, "pid", "sbp", "numeric", "; p=4955")
plot.pairs(pca.it1, alldata, "pid", "dbp", "numeric", "; p=4955")
plot.pairs(pca.it1, alldata, "pid", "NormScale_20", "numeric", "; p=4955")
plot.pairs(pca.it1, alldata, "pid", "hdl", "numeric", "; p=4955")
plot.pairs(pca.it1, alldata, "pid", "ldl", "numeric", "; p=4955")
plot.pairs(pca.it1, alldata, "pid", "total_chol", "numeric", "; p=4955")
plot.pairs(pca.it1, alldata, "pid", "tri", "numeric", "; p=4955")
dev.off()
alldata.it1 <- alldata[-which(alldata$pid %in% pca.it1.outliers),]
alldata.it1 <- alldata.it1[, -grep("PC", names(alldata.it1))]
#alldata.it1 <- merge(alldata.it1, pca.it1$x[, 1:10] %>% as.data.frame(), by.x = "pid", by.y = "row.names")
save(alldata.it1, pca.it1, annot.final, full.annot, file = paste0(dir_out, "PostPCA_It1.rds"))

####PCA - Subsequent iterations####
row.names(alldata.it1) <- alldata.it1$pid
start.new <- grep("Seq", names(alldata.it1)) %>% head() %>% .[1]
last.new <- grep("Seq", names(alldata.it1)) %>% tail() %>% .[6]

## Iteration 2
print("Running PCA iteration #2")
pca.it2 <- prcomp(alldata.it1[, c(start.new:last.new)], scale = F, center = F)
pca.it2.outliers <- find.PCA.outliers(pca.it2, 4)
print("Finish PCA iteration #2")
print(paste0("Number of outliers: ", length(pca.it2.outliers)))

if(length(pca.it2.outliers)<=3){
  print("No more PCA, move on to generating final datasets. Only 2 iterations of PCA done.")
  final.pca <- pca.it2
  final.pca.dataset <- alldata.it1
  
  # ## Skree plot and color by covariates
  # pdf(file = paste0(dir_out, "Final_PCA.pdf"), width = 10, height = 10)
  # par(cex = 1.1)
  # plot(summary(final.pca)$importance[2,1:10]*100, ylab="% Variance Explained", xlab="PC", type="b", main = paste0("PCA for ", nrow(final.pca.dataset), " samples; p=4955"), pch = 19)
  # axis(1, at = c(1:10))
  # if(sex_strat == "Both"){
  #   plot.pairs(final.pca, final.pca.dataset, "pid", "sex", "factor", "; p=4955")
  # } else if(sex_strat == "F" | sex_strat == "M"){
  #   print("No sex plot")
  # }
  # plot.pairs(final.pca, final.pca.dataset, "pid", "center", "factor", "; p=4955")
  # if(race_strat == "Both"){
  #   plot.pairs(final.pca, final.pca.dataset, "pid", "race", "factor", "; p=4955")
  # } else if(race_strat == "White" | race_strat == "Black"){
  #   print("No race plot")
  # }
  # plot.pairs(final.pca, final.pca.dataset, "pid", "age", "numeric", "; p=4955")
  # plot.pairs(final.pca, final.pca.dataset, "pid", "eGFR_ckdepi", "numeric", "; p=4955")
  # plot.pairs(final.pca, final.pca.dataset, "pid", "smoke", "factor", "; p=4955")
  # plot.pairs(final.pca, final.pca.dataset, "pid", "bmi", "numeric", "; p=4955")
  # plot.pairs(final.pca, final.pca.dataset, "pid", "height", "numeric", "; p=4955")
  # plot.pairs(final.pca, final.pca.dataset, "pid", "sbp", "numeric", "; p=4955")
  # plot.pairs(final.pca, final.pca.dataset, "pid", "dbp", "numeric", "; p=4955")
  # plot.pairs(final.pca, final.pca.dataset, "pid", "NormScale_20", "numeric", "; p=4955")
  # plot.pairs(final.pca, final.pca.dataset, "pid", "ldl", "numeric", "; p=4955")
  # plot.pairs(final.pca, final.pca.dataset, "pid", "total_chol", "numeric", "; p=4955")
  # plot.pairs(final.pca, final.pca.dataset, "pid", "tri", "numeric", "; p=4955")
  # dev.off()
  
  #final.alldata <- merge(final.pca.dataset, final.pca$x[, 1:10] %>% as.data.frame(), by.x = "pid", by.y = "row.names")
  #save(final.alldata, annot.final, file = paste0(dir_out, "Final.rds"))
  print("Final data made after PCA iteration #2.")
  
} else if(length(pca.it2.outliers)>3){
  print("Running PCA iteration #3")
  alldata.it2 <- alldata.it1[-which(alldata.it1$pid %in% pca.it2.outliers),]
  #save(alldata.it2, pca.it2, annot.final, file = paste0(dir_out, "PostPCA_It2.rds"))
  row.names(alldata.it2) <- alldata.it2$pid
  ## Iteration 3
  pca.it3 <- prcomp(alldata.it2[, c(start.new:last.new)], scale = F, center = F)
  pca.it3.outliers <- find.PCA.outliers(pca.it3, 4)
  print(paste0("Number of outliers: ", length(pca.it3.outliers)))
  print("Finish PCA iteration #3")
}

if(!exists("pca.it3.outliers")){
  print("Ran less than 3 PCA iterations.")
} else if(exists("pca.it3.outliers") & length(pca.it3.outliers)<=3){
  print("No more PCA, move on to generating final datasets. Only 3 iterations of PCA done.")
  final.pca <- pca.it3
  final.pca.dataset <- alldata.it2
  print("Final data made after PCA iteration #3")
  
} else if(exists("pca.it3.outliers") & length(pca.it3.outliers)>3){
  print("Running PCA iteration #4")
  alldata.it3 <- alldata.it2[-which(alldata.it2$pid %in% pca.it3.outliers),]
  row.names(alldata.it3) <- alldata.it3$pid
  ## Iteration 4
  pca.it4 <- prcomp(alldata.it3[, c(start.new:last.new)], scale = F, center = F)
  pca.it4.outliers <- find.PCA.outliers(pca.it4, 4)
  print(paste0("Number of outliers: ", length(pca.it4.outliers)))
  print("Finish PCA iteration #4")
}

if(!exists("pca.it4.outliers")){
  print("Ran less than 4 PCA iterations.")
} else if(exists("pca.it4.outliers") & length(pca.it4.outliers)<=3){
  print("No more PCA, move on to generating final datasets. Only 4 iterations of PCA done.")
  final.pca <- pca.it4
  final.pca.dataset <- alldata.it3
  #final.alldata <- merge(final.pca.dataset, final.pca$x[, 1:10] %>% as.data.frame(), by.x = "pid", by.y = "row.names")
} else if(length(pca.it4.outliers)>3){
  print("Running PCA iteration #5")
  alldata.it4 <- alldata.it3[-which(alldata.it3$pid %in% pca.it4.outliers),]
  #Iteration 5
  pca.it5 <- prcomp(alldata.it4[, c(start.new:last.new)], scale = F, center = F)
  pca.it5.outliers <- find.PCA.outliers(pca.it5, 4)
  print(paste0("Number of outliers: ", length(pca.it5.outliers)))
  print("Finish PCA iteration #5")
  print("No more PCA, move on to generating final datasets. Only 5 iterations of PCA done.")
  final.pca <- pca.it5
  final.pca.dataset <- alldata.it4
  #final.alldata <- merge(final.pca.dataset, final.pca$x[, 1:10] %>% as.data.frame(), by.x = "pid", by.y = "row.names")
}

## Skree plot and color by covariates only for 4th and 5th iteration
#if(exists("pca.it4.outliers") | exists("pca.it5.outliers")){
  pdf(file = paste0(dir_out, "Final_PCA.pdf"), width = 10, height = 10)
  par(cex = 1.1)
  plot(summary(final.pca)$importance[2,1:10]*100, ylab="% Variance Explained", xlab="PC", type="b", main = paste0("PCA for ", nrow(final.pca.dataset), " samples; p=4955"), pch = 19)
  axis(1, at = c(1:10))
  if(sex_strat == "Both"){
    plot.pairs(final.pca, final.pca.dataset, "pid", "sex", "factor", "; p=4955")
  } else if(sex_strat == "F" | sex_strat == "M"){
    print("No sex plot")
  }
  plot.pairs(final.pca, final.pca.dataset, "pid", "center", "factor", "; p=4955")
  if(race_strat == "Both"){
    plot.pairs(final.pca, final.pca.dataset, "pid", "race", "factor", "; p=4955")
  } else if(race_strat == "White" | race_strat == "Black"){
    print("No race plot")
  }
  plot.pairs(final.pca, final.pca.dataset, "pid", "age", "numeric", "; p=4955")
  plot.pairs(final.pca, final.pca.dataset, "pid", "eGFR_ckdepi", "numeric", "; p=4955")
  plot.pairs(final.pca, final.pca.dataset, "pid", "smoke", "factor", "; p=4955")
  plot.pairs(final.pca, final.pca.dataset, "pid", "bmi", "numeric", "; p=4955")
  plot.pairs(final.pca, final.pca.dataset, "pid", "height", "numeric", "; p=4955")
  plot.pairs(final.pca, final.pca.dataset, "pid", "sbp", "numeric", "; p=4955")
  plot.pairs(final.pca, final.pca.dataset, "pid", "dbp", "numeric", "; p=4955")
  plot.pairs(final.pca, final.pca.dataset, "pid", "NormScale_20", "numeric", "; p=4955")
  plot.pairs(final.pca, final.pca.dataset, "pid", "hdl", "numeric", "; p=4955")
  plot.pairs(final.pca, final.pca.dataset, "pid", "ldl", "numeric", "; p=4955")
  plot.pairs(final.pca, final.pca.dataset, "pid", "total_chol", "numeric", "; p=4955")
  plot.pairs(final.pca, final.pca.dataset, "pid", "tri", "numeric", "; p=4955")
  dev.off()
#} else if (!(exists("pca.it4.outliers")) | !(exists("pca.it5.outliers"))){
#  print("4th and/or 5th PCA iteration(s) not run")
#}

####Save final data####
print("Saving final dataset")
final.alldata <- merge(final.pca.dataset, final.pca$x[, 1:10] %>% as.data.frame(), by.x = "pid", by.y = "row.names")
print(paste0("Final numbers: ", nrow(final.alldata), " people and ", num.prot(final.alldata), " proteins; ", final.alldata[final.alldata$scd_age_filt == 1,] %>% nrow(), " cases and ", final.alldata[final.alldata$scd_age_filt == 0,] %>% nrow(), " controls"))
save(final.alldata, annot.final, full.annot, final.pca, final.pca.dataset, file = paste0(dir_out, "Final.rds"))

if(file.exists(paste0(dir_out, "Final.rds")) & assoc_cont == "Yes"){
  print("Final output made, run associations.")
  system("mkdir ./association")
  setwd("./association")
  system(paste0("export appendix=", appendix, "; qsub -v appendix aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/scripts/runSCAAssocMulti_ver3.sh"))
} else if(file.exists(paste0(dir_out, "Final.rds")) & assoc_cont == "No"){
  print("Final output made.")
}

if(postStrat == "Yes" & assoc_cont == "Yes"){
  #setwd(dir_out)
  #system("mkdir ./postStrat")
  #setwd(paste0(dir_out, "postStrat"))
  #template <- final.alldata
  
  if(sex_strat == "Both" & race_strat == "Both"){
    setwd(dir_out)
    system("mkdir ./postStrat")
    template <- final.alldata
    
    final.alldata <- template[template$sex == "F",]
    system(paste0("mkdir ", dir_out, "/postStrat/female"))
    save(final.alldata, annot.final, full.annot, file = paste0(dir_out, "postStrat/female/Final.rds"))
    system(paste0("mkdir ", dir_out, "postStrat/female/association"))
    setwd(paste0(dir_out, "postStrat/female/association"))
    appendix.ps = "BthRceFem.FmCombPS"
    system(paste0("export appendix=", appendix.ps, "; qsub -v appendix aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/scripts/runSCAAssocMulti_ver3.sh"))
    
    setwd(dir_out)
    final.alldata <- template[template$sex == "M",]
    system(paste0("mkdir ", dir_out, "/postStrat/male"))
    save(final.alldata, annot.final, full.annot, file = paste0(dir_out, "postStrat/male/Final.rds"))
    system(paste0("mkdir ", dir_out, "postStrat/male/association"))
    setwd(paste0(dir_out, "postStrat/male/association"))
    appendix.ps = "BthRceMale.FmCombPS"
    system(paste0("export appendix=", appendix.ps, "; qsub -v appendix aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/scripts/runSCAAssocMulti_ver3.sh"))
    
    setwd(dir_out)
    final.alldata <- template[template$race == "White",]
    system(paste0("mkdir ", dir_out, "/postStrat/EA"))
    save(final.alldata, annot.final, full.annot, file = paste0(dir_out, "postStrat/EA/Final.rds"))
    system(paste0("mkdir ", dir_out, "postStrat/EA/association"))
    setwd(paste0(dir_out, "postStrat/EA/association"))
    appendix.ps = "EABthSx.FmCombPS"
    system(paste0("export appendix=", appendix.ps, "; qsub -v appendix aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/scripts/runSCAAssocMulti_ver3.sh"))
    
    setwd(dir_out)
    final.alldata <- template[template$race == "Black",]
    system(paste0("mkdir ", dir_out, "/postStrat/AA"))
    save(final.alldata, annot.final, full.annot, file = paste0(dir_out, "postStrat/AA/Final.rds"))
    system(paste0("mkdir ", dir_out, "postStrat/AA/association"))
    setwd(paste0(dir_out, "postStrat/AA/association"))
    appendix.ps = "AABthSx.FmCombPS"
    system(paste0("export appendix=", appendix.ps, "; qsub -v appendix aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/scripts/runSCAAssocMulti_ver3.sh"))
    
    setwd(dir_out)
    final.alldata <- template[template$race == "White" & template$sex == "F",]
    system(paste0("mkdir ", dir_out, "/postStrat/EA/female"))
    save(final.alldata, annot.final, full.annot, file = paste0(dir_out, "postStrat/EA/female/Final.rds"))
    system(paste0("mkdir ", dir_out, "postStrat/EA/female/association"))
    setwd(paste0(dir_out, "postStrat/EA/female/association"))
    appendix.ps = "EAFem.FmCombPS"
    system(paste0("export appendix=", appendix.ps, "; qsub -v appendix aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/scripts/runSCAAssocMulti_ver3.sh"))
    
    setwd(dir_out)
    final.alldata <- template[template$race == "Black" & template$sex == "F",]
    system(paste0("mkdir ", dir_out, "/postStrat/AA/female"))
    save(final.alldata, annot.final, full.annot, file = paste0(dir_out, "postStrat/AA/female/Final.rds"))
    system(paste0("mkdir ", dir_out, "postStrat/AA/female/association"))
    setwd(paste0(dir_out, "postStrat/AA/female/association"))
    appendix.ps = "AAFem.FmCombPS"
    system(paste0("export appendix=", appendix.ps, "; qsub -v appendix aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/scripts/runSCAAssocMulti_ver3.sh"))
    
    setwd(dir_out)
    final.alldata <- template[template$race == "White" & template$sex == "M",]
    system(paste0("mkdir ", dir_out, "/postStrat/EA/male"))
    save(final.alldata, annot.final, full.annot, file = paste0(dir_out, "postStrat/EA/male/Final.rds"))
    system(paste0("mkdir ", dir_out, "postStrat/EA/male/association"))
    setwd(paste0(dir_out, "postStrat/EA/male/association"))
    appendix.ps = "EAMale.FmCombPS"
    system(paste0("export appendix=", appendix.ps, "; qsub -v appendix aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/scripts/runSCAAssocMulti_ver3.sh"))
    
    setwd(dir_out)
    final.alldata <- template[template$race == "Black" & template$sex == "M",]
    system(paste0("mkdir ", dir_out, "/postStrat/AA/male"))
    save(final.alldata, annot.final, full.annot, file = paste0(dir_out, "postStrat/AA/male/Final.rds"))
    system(paste0("mkdir ", dir_out, "postStrat/AA/male/association"))
    setwd(paste0(dir_out, "postStrat/AA/male/association"))
    appendix.ps = "AAMale.FmCombPS"
    system(paste0("export appendix=", appendix.ps, "; qsub -v appendix aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/scripts/runSCAAssocMulti_ver3.sh"))
    
  } else if(sex_strat == "Both" & race_strat != "Both"){
    if(race_strat == "White"){
      race_dir = "EA"
    } else {
      race_dir = "AA"
    }
    
    setwd(dir_out)
    system("mkdir ./postStrat")
    #setwd(paste0(dir_out, "postStrat"))
    template <- final.alldata
    
    final.alldata <- template[template$sex == "F",]
    system(paste0("mkdir ", dir_out, "postStrat/female"))
    save(final.alldata, annot.final, full.annot, file = paste0(dir_out, "postStrat/female/Final.rds"))
    system(paste0("mkdir ", dir_out, "postStrat/female/association"))
    setwd(paste0(dir_out, "postStrat/female/association"))
    appendix.ps = paste0(race_dir, "Fem.PS")
    system(paste0("export appendix=", appendix.ps, "; qsub -v appendix aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/scripts/runSCAAssocMulti_ver3.sh"))
    
    setwd(dir_out)
    final.alldata <- template[template$sex == "M",]
    system(paste0("mkdir ", dir_out, "postStrat/male"))
    save(final.alldata, annot.final, full.annot, file = paste0(dir_out, "postStrat/male/Final.rds"))
    system(paste0("mkdir ", dir_out, "postStrat/male/association"))
    setwd(paste0(dir_out, "postStrat/male/association"))
    appendix.ps = paste0(race_dir, "Male.PS")
    system(paste0("export appendix=", appendix.ps, "; qsub -v appendix aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/scripts/runSCAAssocMulti_ver3.sh"))
  }
} else if(postStrat == "Yes" & assoc_cont == "No"){
  print("Final output made.")
  
  if(sex_strat == "Both" & race_strat == "Both"){
    setwd(dir_out)
    system("mkdir ./postStrat")
    template <- final.alldata
    
    final.alldata <- template[template$sex == "F",]
    system("mkdir ./postStrat/female")
    save(final.alldata, annot.final, full.annot, appendix, file = paste0(dir_out, "postStrat/female/Final.rds"))
    
    final.alldata <- template[template$sex == "M",]
    system("mkdir ./postStrat/male")
    save(final.alldata, annot.final, full.annot, appendix, file = paste0(dir_out, "postStrat/male/Final.rds"))
    
    final.alldata <- template[template$race == "White",]
    system("mkdir ./postStrat/EA")
    save(final.alldata, annot.final, full.annot, appendix, file = paste0(dir_out, "postStrat/EA/Final.rds"))

    final.alldata <- template[template$race == "Black",]
    system("mkdir ./postStrat/AA")
    save(final.alldata, annot.final, full.annot, appendix, file = paste0(dir_out, "postStrat/AA/Final.rds"))
   
    final.alldata <- template[template$race == "White" & template$sex == "F",]
    system("mkdir ./postStrat/EA/female")
    save(final.alldata, annot.final, full.annot, appendix, file = paste0(dir_out, "postStrat/EA/female/Final.rds"))
   
    final.alldata <- template[template$race == "Black" & template$sex == "F",]
    system("mkdir ./postStrat/AA/female")
    save(final.alldata, annot.final, full.annot, appendix, file = paste0(dir_out, "postStrat/AA/female/Final.rds"))
  
    final.alldata <- template[template$race == "White" & template$sex == "M",]
    system("mkdir ./postStrat/EA/male")
    save(final.alldata, annot.final, full.annot, appendix, file = paste0(dir_out, "postStrat/EA/male/Final.rds"))

    final.alldata <- template[template$race == "Black" & template$sex == "M",]
    system("mkdir ./postStrat/AA/male")
    save(final.alldata, annot.final, full.annot, appendix, file = paste0(dir_out, "postStrat/AA/male/Final.rds"))
   
  } else if(sex_strat == "Both" & race_strat != "Both"){
    if(race_strat == "White"){
      race_dir = "EA"
    } else {
      race_dir = "AA"
    }
    setwd(dir_out)
    system("mkdir ./postStrat")
    template <- final.alldata
    
    final.alldata <- template[template$sex == "F",]
    system("mkdir ./postStrat/female")
    save(final.alldata, annot.final, full.annot, appendix, file = paste0(dir_out, "postStrat/female/Final.rds"))

    final.alldata <- template[template$sex == "M",]
    system("mkdir ./postStrat/male")
    save(final.alldata, annot.final, full.annot, appendix, file = paste0(dir_out, "postStrat/male/Final.rds"))
  }
}
