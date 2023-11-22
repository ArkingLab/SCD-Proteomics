####Description####
## Format SCA GWAS statistics for MR

####Load packages####
library(TwoSampleMR)
library(magrittr)

####Load data####
outcome <- readRDS("aric.prot.real.archive/scd.assoc/7_MR/SCA.Sum.Stats/SCA.2018.GWAS_SumStats_hg38.updated.RDS")
#Effect is for Allele1

## Remove SNPs without positions
outcome.filtered <- subset(outcome, POS != "")

## Duplicated SNPs still exist
table(duplicated(outcome.filtered$CHR.POS.MIN.MAJ))
#1752
duplicates <- outcome.filtered[duplicated(outcome.filtered$CHR.POS.MIN.MAJ)|duplicated(outcome.filtered$CHR.POS.MIN.MAJ, fromLast = TRUE),] %>% .[order(.$CHR.POS.MIN.MAJ),]
#3504 SNPs
table(duplicated(duplicates$MarkerName), duplicated(duplicates$CHR.POS.MIN.MAJ))
#FALSE TRUE
#FALSE  1752 1752
#rsIDs are unique but duplicated POS, Allele; SNPs have different statistics...
#Remove these SNPs...easiest solution
outcome.filtered <- subset(outcome.filtered, !(CHR.POS.MIN.MAJ %in% unique(duplicates$CHR.POS.MIN.MAJ)))

####Reformat outcome for MR####
outcome.rf <- format_data(outcome.filtered, type = "outcome", header = TRUE, phenotype_col = "Phenotype", 
                          snp_col = "CHR.POS.MIN.MAJ", beta_col = "Effect", se_col = "StdErr", 
                          eaf_col = "Freq1", effect_allele_col = "Allele1", 
                          other_allele_col = "Allele2", pval_col = "P.value", 
                          chr_col = "CHR", pos_col = "POS")
#No errors

## Add rsids
outcome.rf$mergeid <- toupper(outcome.rf$SNP) %>% gsub("CHR", "chr", .)
outcome.rf <- merge(outcome.rf, outcome.filtered[, c("CHR.POS.MIN.MAJ", "MarkerName")], by.x = "mergeid", by.y = "CHR.POS.MIN.MAJ")
outcome.rf$mergeid <- NULL
names(outcome.rf)[ncol(outcome.rf)] <- "rsid"

####Save output####
saveRDS(outcome.rf, "aric.prot.real.archive/scd.assoc/7_MR/SCA.Sum.Stats/SCA.2018.GWAS_SumStats_hg38_MR.formatted.RDS")