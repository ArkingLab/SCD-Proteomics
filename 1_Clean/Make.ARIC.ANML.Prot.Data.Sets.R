####Description####
## This script will merge proteins and phenotype information, generating datasets to move forward with
## This script is for the new ANML normalization

####Load packages####
library(magrittr)
library(data.table)
library(lattice)
library(ggplot2)
library(gridExtra)
library(plyr)
library(cowplot)
library(readxl)

####Load functions####
lookc <- function(data){
  data[1:10, 1:10]
}

isna <- function(data, column){
  table(is.na(data[,column]))
}

find.omic.callrate <- function(data, omic.name, type){
  #' @param data - dataframe with omic measurement in column, rows are individuals
  #' @param omic.name - character vector of omic measurement (i.e. protein, metabolite names, etc.)
  #' @param type - type of omics (i.e. protein, metabolite, etc.)
  
  all=list()
  for (i in 1:length(omic.name)){
    omic=as.character(omic.name[i])
    testing=data[which(names(data)==omic)]
    all[[omic]]=c(sum(is.na(testing)), min(testing, na.rm=TRUE), max(testing, na.rm=TRUE), mean(testing[,1], na.rm=TRUE), median(testing[,1],na.rm=TRUE))
  }
  all.miss.list=as.data.frame(all)
  names(all.miss.list)=omic.name
  all.miss.list=as.data.frame(all.miss.list)
  all.miss.list.t=as.data.frame(t(all.miss.list))
  all.miss.list.t[,type]=rownames(all.miss.list.t)
  colnames(all.miss.list.t)=c("N_Missing","Min","Max","Mean","Median",type)
  
  return(all.miss.list.t)
}

num.prot <- function(data){
  grep("Seq", names(data)) %>% length()
}

col.prot <- function(data){
  print("Col num of first ptn:")
  start = grep("Seq", names(data)) %>% head() %>% .[1]
  print(start)
  print("Col num of last ptn:")
  last=grep("Seq", names(data)) %>% tail() %>% .[6]
  print(last)
}

####Arguments####
dir_out <- "aric.prot.real.archive/ecg.assoc/1_data.cleaning/ANML/"

####Load data####
## Fluorescence measurement, log2 transformed already, SMP dataset (normalized to pooled, healthy controls), ANML normalization
v3data=fread("aric.anml.prot/Data/Visit_3/soma_visit_3_log2_ANML_SMP_updated.txt") %>% as.data.frame()
v5data=fread("aric.anml.prot/Data/Visit_5/soma_visit_5_log2_ANML_SMP_updated.txt") %>% as.data.frame()
v2data=fread("aric.anml.prot/Data/Visit_2/soma_visit_2_log2_ANML_SMP_updated.txt") %>% as.data.frame()

row.names(v3data) <- v3data$SampleId
row.names(v5data) <- v5data$SampleId
row.names(v2data) <- v2data$SampleId

v3proteins <- names(v3data)[grep("Seq", names(v3data))]
v5proteins <- names(v5data)[grep("Seq", names(v5data))]
v2proteins <- names(v2data)[grep("Seq", names(v2data))]
#Same number and proteins observed at all visits

## Re-arrange so protein data is last
v3data=v3data[, c(1:34, 5319:5331, grep("Seq", names(v3data)))]
v5data=v5data[, c(1:33, 5318:5329, grep("Seq", names(v5data)))]
v2data=v2data[, c(1:34, 5319:5332, grep("Seq", names(v2data)))]

## Sample and apatamer flagging, per visit (row = ARIC person, column = f_SeqID (for person in row X, column = 1 if the protein was flagged for this person); N_Sample_Flag (total # of proteins that were flagged for each person)
v3.sample.flag <- fread("aric.anml.prot/Data/Visit_3/soma_visit_3_sample_flag_ANML.txt") %>% as.data.frame()
v5.sample.flag <- fread("aric.anml.prot/Data/Visit_5/soma_visit_5_sample_flag_ANML.txt") %>% as.data.frame()
v2.sample.flag <- fread("aric.anml.prot/Data/Visit_2/soma_visit_2_sample_flag_ANML.txt") %>% as.data.frame()

## Protein annotation file
annot <- read_xlsx("aric.anml.prot/Annotation/Protein Dictionary and Abbreviated Annotation-soma visit 2&3&5.xlsx") %>% as.data.frame()

## Proteins with multiple aptamers
multi_apt <- fread("aric.static.prot.archive/Targets with Multiple Aptamers.txt") %>% as.data.frame()
row.names(multi_apt) <- multi_apt$seqid_in_sample
  
## Add multi_apt column into annot file
annot$multi.apt <- NA
annot$multi.apt <- ifelse(annot$seqid_in_sample %in% multi_apt$seqid_in_sample, 1, 0)
row.names(annot) <- annot$seqid_in_sample

####Apply flag2 filter####
flag2_proteins <- subset(annot, flag2==1)$seqid_in_sample
v3data.nof2 <- v3data[,!(names(v3data) %in% flag2_proteins)]
v5data.nof2 <- v5data[,!(names(v5data) %in% flag2_proteins)]
v2data.nof2 <- v2data[,!(names(v2data) %in% flag2_proteins)]

####Multiple aptamer correlation####
## Look at correlations between multiple aptamers after flag2=1 aptamers are removed
annot_multi <- subset(annot, flag2 == 0) %>% subset(., multi.apt == 1)
multi_apt_unique_uniprot_id <- unique(annot_multi$uniprot_id)
multi_count <- ddply(annot_multi,.(uniprot_id),nrow)
#Examine aptamers with duplicates only, for ease of plotting
multi_count_double <- subset(multi_count, V1 == 2)

## Visit 3
plot_list = list()
for (i in 1:length(multi_count_double$uniprot_id)){
  ptn <- multi_count_double$uniprot_id[i]
  print(paste0("Making plot for ", ptn, ": #", i))
  seq_ids <- annot_multi$seqid_in_sample[grep(ptn, annot_multi$uniprot_id)]
  seq_idA <- seq_ids[1]
  seq_idB <- seq_ids[2]
  full_nameA <- annot_multi[annot_multi$seqid_in_sample == seq_idA, "uniprot.full.name"]
  full_nameB <- annot_multi[annot_multi$seqid_in_sample == seq_idB, "uniprot.full.name"]
  p = ggplot(v3data, aes_string(x=seq_idA,y=seq_idB)) + geom_point() + 
    #geom_smooth(method=lm, se=FALSE) + 
    #geom_abline(intercept=0, slope = 1, col = "blue", size = 1.7) + 
    theme_bw() + ylab(paste0(seq_idB, "-", full_nameB)) + xlab(paste0(seq_idA, "-", full_nameA)) + ggtitle(ptn, subtitle = paste0("r = ", round(cor(v3data[,seq_idA], v3data[,seq_idB]),3))) + theme(text = element_text(size=10))
  plot_list[[i]] = p
  print(paste0("Finish plot for ", ptn))
}
pdf(file = paste0(dir_out, "ARIC_ANML_Prot_V3_Multiple_Aptamers_Doubles_Cor.pdf"), width = 20, height = 13)
for (z in seq(1, length(plot_list)-2, 5)){
  grid.arrange(grobs=plot_list[z:(z+4)], 
               ncol=2)
}
print(plot_list[206])
print(plot_list[207])
dev.off()

## Visit 5
plot_list = list()
for (i in 1:length(multi_count_double$uniprot_id)){
  ptn <- multi_count_double$uniprot_id[i]
  print(paste0("Making plot for ", ptn, ": #", i))
  seq_ids <- annot_multi$seqid_in_sample[grep(ptn, annot_multi$uniprot_id)]
  seq_idA <- seq_ids[1]
  seq_idB <- seq_ids[2]
  full_nameA <- annot_multi[annot_multi$seqid_in_sample == seq_idA, "uniprot.full.name"]
  full_nameB <- annot_multi[annot_multi$seqid_in_sample == seq_idB, "uniprot.full.name"]
  p = ggplot(v5data, aes_string(x=seq_idA,y=seq_idB)) + geom_point() + 
    #geom_smooth(method=lm, se=FALSE) + 
    #geom_abline(intercept=0, slope = 1, col = "blue", size = 1.7) + 
    theme_bw() + ylab(paste0(seq_idB, "-", full_nameB)) + xlab(paste0(seq_idA, "-", full_nameA)) + ggtitle(ptn, subtitle = paste0("r = ", round(cor(v5data[,seq_idA], v5data[,seq_idB]),3))) + theme(text = element_text(size=10))
  plot_list[[i]] = p
  print(paste0("Finish plot for ", ptn))
}
pdf(file = paste0(dir_out, "ARIC_ANML_Prot_V5_Multiple_Aptamers_Doubles_Cor.pdf"), width = 20, height = 13)
for (z in seq(1, length(plot_list)-2, 5)){
  grid.arrange(grobs=plot_list[z:(z+4)], 
               ncol=2)
}
print(plot_list[206])
print(plot_list[207])
dev.off()

## Visit 2
plot_list = list()
for (i in 1:length(multi_count_double$uniprot_id)){
  ptn <- multi_count_double$uniprot_id[i]
  print(paste0("Making plot for ", ptn, ": #", i))
  seq_ids <- annot_multi$seqid_in_sample[grep(ptn, annot_multi$uniprot_id)]
  seq_idA <- seq_ids[1]
  seq_idB <- seq_ids[2]
  full_nameA <- annot_multi[annot_multi$seqid_in_sample == seq_idA, "uniprot.full.name"]
  full_nameB <- annot_multi[annot_multi$seqid_in_sample == seq_idB, "uniprot.full.name"]
  p = ggplot(v2data, aes_string(x=seq_idA,y=seq_idB)) + geom_point() + 
    #geom_smooth(method=lm, se=FALSE) + 
    #geom_abline(intercept=0, slope = 1, col = "blue", size = 1.7) + 
    theme_bw() + ylab(paste0(seq_idB, "-", full_nameB)) + xlab(paste0(seq_idA, "-", full_nameA)) + ggtitle(ptn, subtitle = paste0("r = ", round(cor(v2data[,seq_idA], v2data[,seq_idB]),3))) + theme(text = element_text(size=10))
  plot_list[[i]] = p
  print(paste0("Finish plot for ", ptn))
}
pdf(file = paste0(dir_out, "ARIC_ANML_Prot_V2_Multiple_Aptamers_Doubles_Cor.pdf"), width = 20, height = 13)
for (z in seq(1, length(plot_list)-2, 5)){
  grid.arrange(grobs=plot_list[z:(z+4)], 
               ncol=2)
}
print(plot_list[206])
print(plot_list[207])
dev.off()

## Remove plots, otherwise .RData file is very large - DO LATER
rm(plot_list, p)

####Protein call rate, flag2 subset####

## Check protein call rate, looking at proteins in V3 and in all samples
v3data.nof2.pMissing <- find.omic.callrate(v3data.nof2, as.character(names(v3data.nof2)[grep("Seq", names(v3data.nof2))]), "Protein")
v3data.nof2.pMissing$callrate <- 1-(v3data.nof2.pMissing$N_Missing/nrow(v3data.nof2))

## Check protein call rate, looking at proteins in V5 and in all samples
v5data.nof2.pMissing <- find.omic.callrate(v5data.nof2, as.character(names(v5data.nof2)[grep("Seq", names(v5data.nof2))]), "Protein")
v5data.nof2.pMissing$callrate <- 1-(v5data.nof2.pMissing$N_Missing/nrow(v5data.nof2))

## Check protein call rate, looking at proteins in V5 and in all samples
v2data.nof2.pMissing <- find.omic.callrate(v2data.nof2, as.character(names(v2data.nof2)[grep("Seq", names(v2data.nof2))]), "Protein")
v2data.nof2.pMissing$callrate <- 1-(v2data.nof2.pMissing$N_Missing/nrow(v2data.nof2))

plot_grid(densityplot(v3data.nof2.pMissing$callrate, xlab="ANML Protein Call Rate", main="V3 Samples, Proteins After flag2 Exclusion"),
          densityplot(v5data.nof2.pMissing$callrate, xlab="ANML Protein Call Rate", main="V5 Samples, Proteins After flag2 Exclusion"),
          densityplot(v2data.nof2.pMissing$callrate, xlab="ANML Protein Call Rate", main="V2 Samples, Proteins After flag2 Exclusion"))
#No missing values for any individual

## After multiple discussions at TOPMed meeting, will not remove aptamers based on individual flags so sections are removed ##
## See non-ANML data for this section's code ##

####Apply flag3 filter####
flag3_proteins <- subset(annot, flag3==1)$seqid_in_sample
v3data.nof3 <- v3data[,!(names(v3data) %in% flag3_proteins)]
v5data.nof3 <- v5data[,!(names(v5data) %in% flag3_proteins)]
v2data.nof3 <- v2data[,!(names(v2data) %in% flag3_proteins)]

####Make datasets####
## Will make two datasets: (1) Remove proteins with Flag2=1; (2) Remove proteins with Flag3=1
annot.nof2 <- subset(annot, flag2 == 0)
annot.nof3 <- subset(annot, flag3 == 0)

## Flag2=1 removed
save(annot.nof2, v3data.nof2, v5data.nof2, v2data.nof2, multi_apt, file = paste0(dir_out, "ARIC.V2V3V5.ANML.Prot.NoFlag2.rds"))

## Flag3=1 removed
save(annot.nof3, v3data.nof3, v5data.nof3, v2data.nof3, multi_apt, file = paste0(dir_out, "ARIC.V2V3V5.ANML.Prot.NoFlag3.rds"))

save.image(file = paste0(dir_out, "Make.ARIC.ANML.Prot.Data.Sets.RData"))
