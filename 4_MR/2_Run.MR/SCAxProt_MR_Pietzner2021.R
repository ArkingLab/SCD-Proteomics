####Description####
## Runs MR using TwoSampleMR package
## pQTL database is from Pietzner 2021 paper
## Options to choose pQTL P, MAF, various clump settings
## Does NOT use --clump-best function
## All exposures are analyzed in the same script
## Only exposures that past pQTL p=0.05 are included in pQTL input
## Modified from Ferkingstad 2021 script

## Changes from version 1: updated meta data to version with 10/03/2022 CHS data

####Load packages####
library(TwoSampleMR)
library(magrittr)
library(data.table)
library(readxl)
library(ggplot2)
library(writexl)
library(dplyr)
library(gridExtra)
library(kableExtra)

# library(cowplot)
# library(biomaRt)
# library(ieugwasr)

####Load functions####
mr_scatter_plot_vy <- function(mr_results, dat, auto.axis, minYval, maxYval, maxXval, minXval)
{
  # dat <- subset(dat, paste(id.outcome, id.exposure) %in% paste(mr_results$id.outcome, mr_results$id.exposure))
  requireNamespace("ggplot2", quietly=TRUE)
  requireNamespace("plyr", quietly=TRUE)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(nrow(d) < 2 | sum(d$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
    mrres$a <- 0
    if("MR Egger" %in% mrres$method)
    {
      temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {
      temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    ## Set plot limits
    if(auto.axis == "Y"){
      inkme <- subset(dat, mr_keep == "TRUE")
      minY <- min(inkme$beta.outcome, (inkme$beta.outcome-inkme$se.outcome)) %>% 
        plyr::round_any(., 0.1, floor)
      maxY <- max(inkme$beta.outcome, (inkme$beta.outcome+inkme$se.outcome)) %>% 
        plyr::round_any(., 0.1, ceiling) + 0.1
      maxX <- max(inkme$beta.exposure, (inkme$beta.exposure+inkme$se.exposure)) %>% 
        plyr::round_any(., 0.1, ceiling)
      #minX <- min(inkme$beta.exposure, (inkme$beta.exposure-inkme$se.exposure)) %>% 
        #plyr::round_any(., 0.1, ceiling) + 0.15
      minX <- -0.1
      # print(paste0("minY=", minY))
      # print(paste0("maxY=", maxY))
      # print(paste0("minX=", minX))
      # print(paste0("maxX=", maxX))
    } else if(auto.axis == "N"){
      minY <- minYval
      maxY <- maxYval
      maxX <- maxXval
      minX <- minXval
    }
    
    ggplot2::ggplot(data=d, ggplot2::aes(x=beta.exposure, y=beta.outcome)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::geom_abline(data=mrres, ggplot2::aes(intercept=a, slope=b, colour=method), show.legend=TRUE, size = 0.3) +
      #ggplot2::geom_segment(x = 0, y = 0, yend = maxY, xend = subset(mrres, exposure == ))
      ggplot2::scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
      ggplot2::labs(colour="MR Test", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
      #ggplot2::coord_cartesian(xlim = c(minX,maxX), ylim = c(minY,maxY), expand = FALSE) +
      #ggplot2::scale_x_continuous(breaks = seq(0, maxX, by = 0.2)) + 
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(minX,maxX)) +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(minY,maxY)) +
      #ggplot2::scale_y_continuous(breaks = seq(minY, maxY, by = 0.2)) +
      ggplot2::theme(legend.position="top", legend.direction="vertical", legend.title=element_blank()) +
      ggplot2::guides(colour=ggplot2::guide_legend(ncol=2))
  })
  mrres
}

blank_plot <- function(message)
{
  requireNamespace("ggplot2", quietly=TRUE)
  ggplot2::ggplot(data.frame(a=0,b=0,n=message)) + ggplot2::geom_text(ggplot2::aes(x=a,y=b,label=n)) + ggplot2::labs(x=NULL,y=NULL) + ggplot2::theme(axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank())
}

####Arguments####
args <- commandArgs(trailingOnly=TRUE)

## pQTL filters
pqtl_p <- args[1] %>% as.numeric()
pqtl_maf <- args[2]
#Default setting from what paper defined pQTL as
pqtl_distTSS <- 500000
pqtl_hetdf <- args[3]
pqtl_hetpval <- args[4]
pqtl_hetIsq <- args[5]

## Outcome filters
sca_hetdf <- args[6]
sca_hetpval <- args[7]
sca_hetIsq <- args[8]

## Clumping filters
r2 <- args[9] %>% as.numeric()
kb <- args[10] %>% as.numeric()
panel <- args[11]
appendix <- args[12]

####Test arguments####
## pQTL filters
# pqtl_p <- 1.004e-11
# pqtl_maf <- 0.01
# pqtl_distTSS <- 500000
# pqtl_hetdf <- 2
# pqtl_hetpval <- 0.05
# pqtl_hetIsq <- 75

## Outcome filters
# sca_hetdf <- 2
# sca_hetpval <- 0.05
# sca_hetIsq <- 75

## Clumping filters
# r2 <- 0.1
# kb <- 1000
# panel <- "arkinglab/resources/1000G/mrcieu.plink/EUR"
# panel <- "arkinglab/active/projects/scd.meta/analyses/scd.meta.ver2/ARIC.b35.b37.liftover/aric.f3v2.imputed.b37"

write(paste0("####MR parameters####\npQTL P value: ", pqtl_p, "\npQTL MAF: ", pqtl_maf, "\npQTL distTSS: ", pqtl_distTSS, "\nR2 threshold: ", r2, "\nclumping kb: ", kb, "\nclumping panel: ", panel, "\nHetDf: ", pqtl_hetdf, "\nHetPval: ", pqtl_hetpval, "\nHetISq: ", pqtl_hetIsq, "\nSCA HetDf: ", sca_hetdf, "\nSCA HetPval: ", sca_hetpval, "\nSCA HetISq: ", sca_hetIsq, "\n"), file = "output.txt", append = T)

####Read in pQTL data####
## Coordinates are available in hg19 and hg38
## Start and End gene body coordinates are hg19 only
write("Reading in pQTL data\n", file = "output.txt", append = T)

pqtl_file <- "aric.prot.real.archive/scd.assoc/7_MR/pQTL.databases/Pietzner.2021/p_05/All.Exposures_Filtered.RDS"
pqtl <- readRDS(pqtl_file)

####Filter pQTL data####
expo <- pqtl

# ## Additional P value filter beyond 0.05
# print("Filtering pQTL based on P value")
# if(pqtl_p != 0.05){
#   expo <- pqtl[(pqtl$Pval < pqtl_p), ]
# }

## MAF filter
write("Filtering pQTL based on MAF\n", file = "output.txt", append = T)
if(pqtl_maf != "none"){
  pqtl_maf <- as.numeric(pqtl_maf)
  expo <- expo[(expo$MAF > pqtl_maf),]
}

## Apply exposure meta-analysis filters
write("Filtering pQTL based on pQTL meta-analysis metrics\n", file = "output.txt", append = T)
if(pqtl_hetdf == 3){
  expo <- expo[(expo$HetDf != 0),]
} else if(pqtl_hetdf == 1 | pqtl_hetdf == 2){
  pqtl_hetdf <- as.numeric(pqtl_hetdf)
  expo <- expo[(expo$HetDf == pqtl_hetdf),]
} else if(pqtl_hetdf == "none"){
  write("No HetDf filter\n", file = "output.txt", append = T)
}

if(pqtl_hetpval != "none"){
  pqtl_hetpval <- as.numeric(pqtl_hetpval)
  expo <- expo[(expo$HetPVal > pqtl_hetpval),]
} else if(pqtl_hetpval == "none"){
  write("No HetPval filter\n", file = "output.txt", append = T)
}

if(pqtl_hetIsq != "none"){
  pqtl_hetIsq <- as.numeric(pqtl_hetIsq)
  expo <- expo[(expo$HetISq < pqtl_hetIsq),]
} else if(pqtl_hetIsq  == "none"){
  write("No HetIsq filter\n", file = "output.txt", append = T)
}

## Keep cis-pQTLs only within +/- 500kb of gene body
# write("Keep cis-pQTLS within a certain distance only\n", file = "output.txt", append = T)
# expo <- expo[(expo$distTSS < pqtl_distTSS),]
# if(nrow(expo) == 0){
#   write("No cis-pQTLS found, script will stop\n", file = "output.txt", append = T)
#   stop()
# }
write("Keep cis-pQTLS only within 500kb of gene body\n", file = "output.txt", append = T)
expo <- expo[expo$cisStatus == "cis",]

####Read in summary statistics####
write("Reading in GWAS summary statistics\n", file = "output.txt", append = T)
#Already formatted for MR, duplicate snps + those without positions already removed
outcome.rf <- readRDS("aric.prot.real.archive/scd.assoc/7_MR/SCA.Sum.Stats/SCA.2018.GWAS_SumStats_hg38_MR.formatted.RDS")

## Filter outcome by MAF
outcome.rf$oaf.outcome <- 1-outcome.rf$eaf.outcome
outcome.rf$MAF <- ifelse(outcome.rf$eaf.outcome < outcome.rf$oaf.outcome, outcome.rf$eaf.outcome, outcome.rf$oaf.outcome)
outcome.rf <- subset(outcome.rf, MAF > pqtl_maf)
outcome.rf$oaf.outcome <- NULL
outcome.rf$MAF <- NULL
names(outcome.rf)[ncol(outcome.rf)] <- "rsid.outcome"

## Add in SCA meta-analysis statistics
sca <- fread("arkinglab/static/papers/Ashar.Mitchell.SCA.2018/GWAS/data/METAANALYSIS_woHSDS.SUDS_Jun20151.tbl.gz") %>% as.data.frame()
outcome.rf <- merge(outcome.rf, sca[, c(1,5:7,11:15)], by.x = "rsid.outcome", by.y = "MarkerName")
og <- outcome.rf

## Apply outcome meta-analysis filters
write("Filtering SCA SNPs based on meta-analysis metrics\n", file = "output.txt", append = T)
if(sca_hetdf == "none"){
  write("No sCA HetDf filter\n", file = "output.txt", append = T)
} else if(sca_hetdf == 9){
  outcome.rf <- outcome.rf[(outcome.rf$HetDf != 0),]
} else if(grepl("g", sca_hetdf)){
  sca_hetdf_num <- gsub("g", "", sca_hetdf) %>% as.numeric()
  outcome.rf <- outcome.rf[(outcome.rf$HetDf>=sca_hetdf_num),]
} else if(grepl("only", pqtl_hetdf)){
  sca_hetdf_num <- gsub("only", "", sca_hetdf) %>% as.numeric()
  outcome.rf <- outcome.rf[(outcome.rf$HetDf==sca_hetdf_num),]
}

if(sca_hetpval != "none"){
  sca_hetpval <- as.numeric(sca_hetpval)
  outcome.rf <- outcome.rf[(outcome.rf$HetPVal > sca_hetpval),]
} else if(sca_hetpval == "none"){
  write("No SCA HetPval filter\n", file = "output.txt", append = T)
}

if(sca_hetIsq != "none"){
  sca_hetIsq <- as.numeric(sca_hetIsq)
  outcome.rf <- outcome.rf[(outcome.rf$HetISq < sca_hetIsq),]
} else if(sca_hetIsq  == "none"){
  write("No SCA HetIsq filter\n", file = "output.txt", append = T)
}

names(outcome.rf)[15:22] <- paste0(names(outcome.rf)[15:22], ".outcome")

####Reformat exposure####
write("Reformatting exposure\n", file = "output.txt", append = T)
#expo$CHR.POS_hg38.MIN.MAJ <- paste0("chr", expo$chr, ":", expo$pos_hg38, ":", expo$Amin, ":", expo$Amaj)
#expo$Chrom <- gsub("chr", "", expo$Chrom)
#expo$Phenotype2 <- paste0(expo$target_use, "_", expo$SeqId.use)
#expo$Phenotype2 <- gsub(" ", ".", expo$Phenotype2)
# expo$Phenotype <- ifelse(expo$target_use == "SVEP1*", paste0("SVEP1", "_", expo$SeqId.use), expo$Phenotype)
# expo$Phenotype <- ifelse(expo$target_use == "GOLM1*", paste0("GOLM1", "_", expo$SeqId.use), expo$Phenotype)
expo <- expo[!is.na(expo$pos_hg38),]
expo.rf <- format_data(expo, type = "exposure", header = TRUE, phenotype_col = "Phenotype", 
                       snp_col = "uniqueid_hg38", beta_col = "Effect", 
                       se_col = "StdErr", eaf_col = "Freq1", 
                       effect_allele_col = "Allele1", other_allele_col = "Allele2", 
                       pval_col = "Pvalue", chr_col = "chr", pos_col = "pos_hg38", 
                       min_pval = 1e-400, samplesize_col = "TotalSampleSize", gene_col = "geneid")

expo.rf$mergeid <- gsub("_.*", "", expo.rf$SNP) %>% toupper()
expo.rf$mergeid <- gsub("CHR", "chr", expo.rf$mergeid)
expo.rf$mergeid <- paste0(expo.rf$mergeid, "_", expo.rf$exposure)
expo.rf$mergeid <- gsub(" ", ".", expo.rf$mergeid)
expo.rf <- merge(expo.rf, expo[, c("SeqId.use","uniqueid_hg38","rsid", "CHR.POS_hg38.MIN.MAJ","FreqSE", "MinFreq", "MaxFreq", "Freq2", "MAF", "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal", "pos_hg19", "strand_prot", "distTSS", "prot_start_hg19", "prot_end_hg19", "distStart", "distEnd", "cisStart", "cisEnd")], by.x = "mergeid", by.y = "uniqueid_hg38")

names(expo.rf)[18] <- "rsid.exposure"
names(expo.rf)[19:38] <- paste0(names(expo.rf)[19:38], ".exposure")
#expo.rf$SNP <- expo.rf$CHR.POS_hg38.MIN.MAJ
#expo.rf$SNP <- tolower(expo.rf$SNP)
names(expo.rf)[1] <- "uniqueid"

#table(expo.rf$rsid.exposure %in% outcome.rf$rsid.outcome)
#FALSE  TRUE
#78993 29976
expo.rf$SNP <- expo.rf$rsid.exposure
outcome.rf$SNP <- outcome.rf$rsid.outcome
#More SNPs with rsid matching

####Finding overlap of SNPs####
write("Finding overlapping SNPs\n", file = "output.txt", append = T)
# overlap <- subset(expo.rf, SNP %in% outcome.rf$SNP)
# overlap <- merge(expo.rf, outcome.rf[, c("SNP", "rsid")], by = "SNP")
# names(overlap)[ncol(overlap)] <- "snpid"
# overlap$P <- overlap$pval.exposure
# #overlap <- merge(overlap, outcome.rf[, c("SNP", "rsid")], by = "SNP")

#action = 3: Correct strand for non-palindromic SNPs, and drop all palindromic SNPs from the analysis (more conservative). If a single value is passed then this action is applied to all outcomes. But multiple values can be supplied as a vector, each element relating to a different outcome.
overlap <- harmonise_data(exposure_dat=expo.rf, outcome_dat=outcome.rf, action = 3)
#overlap2 <- harmonise_data(exposure_dat=expo.rf, outcome_dat=outcome.rf, action = 3)
#More SNPs with rsid

# ## Check rsid overlap with 1KG panel
#bim <- fread("arkinglab/resources/1000G/mrcieu.plink/EUR.bim")
# table(bim$V2 %in% outcome.rf$rsid)
# #FALSE    TRUE
# #6083346 2466810
# 
#table(bim$V2 %in% expo.rf$rsid)
#FALSE    TRUE
#8454736   95420
overlap$exposure <- gsub(" ", ".", overlap$exposure)

####Clump SNPs####
write("Clumping SNPs\n", file = "output.txt", append = T)

system("mkdir ./clumping.input")
system("mkdir ./clumping.output")

# ## Convert P value
# pqtl_p.Z <- qnorm(1-(pqtl_p/2))
# max_zscore <- 116.7863
# scaling_factor <- max_zscore/37.02
# pqtl_p.Z.pseudo <- pqtl_p.Z/scaling_factor
# pqtl_p.converted <- 2*pnorm(abs(pqtl_p.Z.pseudo), lower.tail = F)

## Write out SNPs for clumping
for(i in 1:length(unique(overlap$exposure))){
  gene <- unique(overlap$exposure)[i]
  assign(paste0(gene, "_clump.input"), subset(overlap, exposure == gene) %>% 
           .[, c("SNP", "pval.exposure")])
  write.table(get(paste0(gene, "_clump.input")), 
              file = paste0("./clumping.input/", gene, "_clump.input.txt"), 
              row.names = F, quote = F, col.names = T, sep = "\t")}
  
## Clump
toclump <- list.files(path = "./clumping.input", pattern = "_clump.input.txt")
for(i in 1:length(toclump)){
  filetoclump <- toclump[i]
  exposure <- gsub("_clump.input.txt", "", filetoclump)
  system(paste0("arkinglab/software/src/plink/plink-1.90-beta-20180410/plink --bfile ", panel, " --clump ./clumping.input/", filetoclump, " --clump-r2 ", r2, " --clump-kb ", kb, " --clump-p1 ", pqtl_p, " --maf ", pqtl_maf, " --clump-snp-field SNP --clump-field pval.exposure --out ./clumping.output/", exposure))}

## Read in clumped results
doneclumped <- list.files(path = "./clumping.output", pattern = ".clumped")
for(i in 1:length(doneclumped)){
  filedoneclumped <- doneclumped[i]
  exposure <- gsub("[.]clumped", "", filedoneclumped)
  clumped.df <- fread(paste0("./clumping.output/", filedoneclumped)) %>% .[, c("SNP", "BP")]
  #names(clumped.df)[1] <- "rsid.outcome"
  afterclump <- merge(clumped.df, get(paste0(exposure, "_clump.input")), by = "SNP")
  afterclump$rsid.outcome.expo.uniqueid <- paste0(afterclump$SNP, "_", exposure)
  assign(paste0(exposure, ".clumped.final"), afterclump)}

all.afterclump <- bind_rows(mget(ls(pattern = ".clumped.final")))
names(all.afterclump)[1] <- "SNP"
harmonized.df <- overlap %>% mutate(rsid.outcome.expo.uniqueid = paste0(overlap$SNP, "_", overlap$exposure)) %>% 
  merge(., all.afterclump[, c("rsid.outcome.expo.uniqueid")] %>% as.data.frame(), by = "rsid.outcome.expo.uniqueid")

####Run MR####
included_snps <- subset(harmonized.df, mr_keep == TRUE)
if(nrow(included_snps) == 0){
  write("No SNPs found after harmonizing exposure and outcome, clumping: script will stop\n", file = "output.txt", append = T)
  stop()
} else {
  write("SNPs found after harmonizing exposure and outcome, clumping: analysis will continue\n", file = "output.txt", append = T)
}
write("Running MR...\n", file = "output.txt", append = T)
mr_results <- mr(harmonized.df, method_list = c("mr_wald_ratio", "mr_egger_regression", "mr_simple_median", "mr_weighted_median", "mr_ivw"))

####IVW only####
mr_results_ivw_only <- subset_on_method(mr_results, single_snp_method = "Wald ratio", 
                                         multi_snp_method = "Inverse variance weighted")

####Single SNP tests####
mr_results_single <- mr_singlesnp(harmonized.df)
 
####Tests for robustness####
mr_pleio_results <- mr_pleiotropy_test(harmonized.df)
mr_het_results <- mr_heterogeneity(harmonized.df)

####Check SNP correlation####
write("Checking SNP correlation...\n", file = "output.txt", append = T)
system("mkdir included.snps")

## Write out final SNP list
for(i in 1:length(unique(included_snps$exposure))){
  gene <- unique(included_snps$exposure)[i]
  subset(included_snps, exposure == gene) %>% .[, "rsid.outcome"] %>% 
  write.table(., file = paste0("./included.snps/", gene, ".MR.Final.SNPs.txt"), 
              row.names = F, quote = F, col.names = F, sep = "\t")}

## Calculate R2
system("mkdir snps.R2")
finalsnps <- list.files(path = "./included.snps", pattern = ".MR.Final.SNPs.txt")
for(i in 1:length(finalsnps)){
  filetor2 <- finalsnps[i]
  exposure <- gsub(".MR.Final.SNPs.txt", "", filetor2)
  system(paste0("arkinglab/software/src/plink/plink-1.90-beta-20180410/plink --bfile ", panel, " --extract ./included.snps/", exposure, ".MR.Final.SNPs.txt", " --r2 'with-freqs' --ld-window-r2 ", r2, " --threads 1 --out ./snps.R2/", exposure, ".MR.Final.SNPs"))}

## Read in R2 matrix
doner2 <- list.files(path = "./snps.R2", pattern = ".MR.Final.SNPs.ld")
for(i in 1:length(doner2)){
  r2toread <- doner2[i]
  exposure <- gsub(".MR.Final.SNPs.ld", "", r2toread)
  snps.ld.corr <- fread(paste0("./snps.R2/", r2toread), header = T)
  if(nrow(snps.ld.corr) >= 1){
    snps.ld.corr$Phenotype <- exposure
    snps.ld.corr <- as.data.frame(snps.ld.corr)
    assign(paste0(exposure, ".R2ofSNPs.final"), snps.ld.corr)}}
all.R2 <- bind_rows(mget(ls(pattern = ".R2ofSNPs.final")))
if(nrow(all.R2)>=1){
  write("There are SNPs with R2 past the R2 threshold:", file = "output.txt", append = T)
  density <- ggplot(all.R2, aes(x = R2)) + geom_density(fill = "lightblue") + theme_bw()
  ggsave(paste0(appendix, "_MR.Results.SNPS.with.R2.Over.the.Limit_densityplot.tiff"), plot = density, 
         units = "in", height = 5, width = 5)
  write.table(table(all.R2$Phenotype), file = "output.txt", append = T, row.names = F, quote = F, col.names = F)
  write_xlsx(all.R2, "SNPs.with.R2.Over.the.Limit.xlsx")
} else {
  write("\nNo SNPs past R2 threshold.", file = "output.txt", append = T)
}

####Correct P values####
## Use fdr
n.expo <- length(unique(mr_results$exposure))
mr_results_ivw_only$adj_pval <- p.adjust(mr_results_ivw_only$pval, method = "fdr")

mr_results_edit <- mr_results
mr_results_edit$adj_pval <- NA
mr_results_edit <- subset(mr_results_edit, !(method == "Inverse variance weighted" | method == "Wald ratio"))
mr_results_edit <- rbind(mr_results_edit, mr_results_ivw_only)
mr_results_edit <- mr_results_edit[order(mr_results_edit$exposure),]

####Plots####

## Find limits for scatterplots
# minY <- min(included_snps$beta.outcome) %>% plyr::round_any(., 0.2, ceiling)
# maxY <- max(abs(included_snps$beta.outcome)) %>% plyr::round_any(., 0.2, ceiling)
# maxX <- max(abs(included_snps$beta.exposure), beta.outcome+se.outcome) %>% plyr::round_any(., 0.2, ceiling)

write("\nGenerating plots", file = "output.txt", append = T)
## With all methods
scatterplots.all.methods <- mr_scatter_plot_vy(mr_results_edit, harmonized.df, auto.axis = "Y")
plot_list <- list()
for(i in 1:length(scatterplots.all.methods)){
  plot <- scatterplots.all.methods[i]
  plot_list[i] <- plot
}
ggsave(paste0(appendix, "_All.Methods_Scatterplots.pdf"), marrangeGrob(plot_list, nrow=4, ncol=3), width = 12, height = 13, units = "in")

## With all methods, natural axes
scatterplots.all.methods.nat <- mr_scatter_plot(mr_results_edit, harmonized.df)
plot_list <- list()
for(i in 1:length(scatterplots.all.methods.nat)){
  plot <- scatterplots.all.methods.nat[i]
  plot_list[i] <- plot
}
ggsave(paste0(appendix, "_All.Methods_Scatterplots.NatAxes.pdf"), marrangeGrob(plot_list, nrow=4, ncol=3), width = 12, height = 13, units = "in")

## Just IVW
scatterplots.IVW <- mr_scatter_plot_vy(mr_results_ivw_only, harmonized.df, auto.axis = "Y")
plot_list <- list()
for(i in 1:length(scatterplots.IVW)){
  plot <- scatterplots.IVW[i]
  plot_list[i] <- plot
}
ggsave(paste0(appendix, "_IVW_Scatterplots.pdf"), marrangeGrob(plot_list, nrow=4, ncol=3), width = 12, height = 13, units = "in")

## Just IVW, natural axes
scatterplots.IVW.nat <- mr_scatter_plot(mr_results_ivw_only, harmonized.df)
plot_list <- list()
for(i in 1:length(scatterplots.IVW.nat)){
  plot <- scatterplots.IVW.nat[i]
  plot_list[i] <- plot
}
ggsave(paste0(appendix, "_IVW_Scatterplots.NatAxes.pdf"), marrangeGrob(plot_list, nrow=4, ncol=3), width = 12, height = 13, units = "in")

####Plots of FDR significant results only####
sig.fdr <- subset(mr_results_edit, method == "Inverse variance weighted" & adj_pval < 0.05)
if(nrow(sig.fdr)>=1){
  ## All
  scatterplots.all.methods.fdr <- mr_scatter_plot(mr_results_edit, subset(harmonized.df, exposure %in% sig.fdr$exposure))
  plot_list <- list()
  for(i in 1:length(scatterplots.all.methods.fdr)){
    plot <- scatterplots.all.methods.fdr[i]
    plot_list[i] <- plot
  }
  ggsave(paste0(appendix, "_All.Methods.FDR_Scatterplots.pdf"), marrangeGrob(plot_list, nrow=4, ncol=3), width = 12, height = 13, units = "in")
  
  ## IVW only
  scatterplots.IVW.fdr <- mr_scatter_plot(mr_results_ivw_only, subset(harmonized.df, exposure %in% sig.fdr$exposure))
  plot_list <- list()
  for(i in 1:length(scatterplots.IVW.fdr)){
    plot <- scatterplots.IVW.fdr[i]
    plot_list[i] <- plot
  }
  ggsave(paste0(appendix, "_IVW.FDR_Scatterplots.pdf"), marrangeGrob(plot_list, nrow=4, ncol=3), width = 12, height = 13, units = "in")
} else {
  write("\nNo exposure passed FDR<0.05.", file = "output.txt", append = T)
}

####Summarize results####
write("\nSummarizing results", file = "output.txt", append = T)
all.assoc.results <- readRDS("aric.prot.real.archive/scd.assoc/9_other/combine.data/All.SCAxProteomics.Results_10072022.RDS")
meta.results.min <- subset(all.assoc.results, cohort == "meta" & strat == "all" & model == "M2")

ivw.print <- mr_results_ivw_only
ivw.print <- merge(ivw.print, harmonized.df[, c("id.exposure", "SeqId.use")] %>% unique(), by = "id.exposure")
ivw.print$SeqId.use <- gsub("_", "-", ivw.print$SeqId.use)

ivw.print <- merge(ivw.print, meta.results.min[, c("SeqId.use", "Beta")], by = "SeqId.use")
names(ivw.print)[ncol(ivw.print)] <- "SCAxProt_Beta"
ivw.print$Beta.Dir <- NA
ivw.print$Beta.Dir <- ifelse(ivw.print$b < 0 & ivw.print$SCAxProt_Beta < 0, "--", ivw.print$Beta.Dir)
ivw.print$Beta.Dir <- ifelse(ivw.print$b > 0 & ivw.print$SCAxProt_Beta > 0, "++", ivw.print$Beta.Dir)
ivw.print$Beta.Dir <- ifelse(ivw.print$b < 0 & ivw.print$SCAxProt_Beta > 0, "-+", ivw.print$Beta.Dir)
ivw.print$Beta.Dir <- ifelse(ivw.print$b > 0 & ivw.print$SCAxProt_Beta < 0, "+-", ivw.print$Beta.Dir)
ivw.print <- ivw.print[order(ivw.print$pval),]
ivw.print$method <- ifelse(ivw.print$method == "Inverse variance weighted", "IVW", ivw.print$method)
#ivw.print <- merge(ivw.print, expo[, c("SeqId.use", "target_use")] %>% unique(), by = "SeqId.use")
ivw.print <- merge(ivw.print, meta.results.min[, c("SeqId.use", "target_use")] %>% unique(), by = "SeqId.use")
ivw.print <- ivw.print[order(ivw.print$adj_pval),]
ivw.print.final <- ivw.print %>% mutate(b = round(b, digits = 3), 
                                        se = round(se, digits = 3), 
                                        pval = signif(pval, digits = 3) %>% format(., scientific=T),
                                        adj_pval = signif(adj_pval, digits = 3) %>% format(., scientific=T),
                                        SCAxProt_Beta = round(SCAxProt_Beta, digits = 3))
ivw.print.final$id.exposure <- NULL
ivw.print.final$id.outcome <- NULL
ivw.print.final$outcome <- NULL
row.names(ivw.print.final) <- NULL
names(ivw.print.final)[ncol(ivw.print.final)] <- "Protein"
ivw.print.final <- ivw.print.final[,c(1:2,11,3:10)]

ivw.print.final %>% 
  kbl() %>% kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>% kable_classic(full_width = F) %>% save_kable(file = paste0(appendix, "_IVW.Results.Table.html"), self_contained = T)

sheets <- list("IVW.only" = ivw.print.final, "All.Methods" = mr_results_edit, "Pleio.Results" = mr_pleio_results, "Het.Results" = mr_het_results)
write_xlsx(sheets, paste0(appendix, "_MR.Results.Table.xlsx"))

if(nrow(mr_pleio_results %>% subset(., pval < 0.05)) >=1){
  mr_pleio_results %>% subset(., pval < 0.05) %>% 
    write.table(., file = paste0(appendix, "_Sig.Pleio.Results.txt"), 
                col.names = T, row.names = F, quote = F, sep = "\t")
  write("\nExposures with significant pleiotropy:", file = "output.txt", append = T)
  write(subset(mr_pleio_results, pval < 0.05)$exposure %>% unique(), file = "output.txt", append = T)
} else {
  write("\nNo exposure with significant pleiotropy", file = "output.txt", append = T)
}

if(nrow(mr_het_results %>% subset(., Q_pval < 0.05)) >=1){
  mr_het_results %>% subset(., Q_pval < 0.05) %>% 
    write.table(., file = paste0(appendix, "_Sig.Het.Results.txt"), 
                col.names = T, row.names = F, quote = F, sep = "\t")
  write("\nExposures with significant heterogeniety:", file = "output.txt", append = T)
  write(subset(mr_het_results, Q_pval < 0.05)$exposure %>% unique(), file = "output.txt", append = T)
} else {
  write("\nNo exposure with significant heterogeniety", file = "output.txt", append = T)
}

fdr.order <- ivw.print.final$exposure %>% as.data.frame()
fdr.order$order <- 1:nrow(fdr.order)
names(fdr.order)[1] <- "exposure"
mr_results_edit_print <- merge(mr_results_edit, fdr.order, by = "exposure") %>% .[order(.$order),]
mr_results_edit_print$method <- factor(mr_results_edit_print$method, ordered = T, levels = c("Inverse variance weighted", "MR Egger", "Simple median", "Weighted median"))
mr_results_edit_print <- mr_results_edit_print[order(mr_results_edit_print$order, mr_results_edit_print$method),]
mr_results_edit_print$id.exposure <- NULL
mr_results_edit_print$id.outcome <- NULL
mr_results_edit_print$outcome <- NULL
mr_results_edit_print$order <- NULL
mr_results_edit_print <- mr_results_edit_print %>% dplyr::mutate(b = round(b, digits = 3),
                                       se = round(se, digits = 3),
                                       pval = signif(pval, digits = 3) %>%
                                         format(., scientific = TRUE),
                                       adj_pval = signif(adj_pval, digits = 3) %>%
                                         format(., scientific = TRUE))
row.names(mr_results_edit_print) <- NULL
mr_results_edit_print %>% kbl() %>% kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  kable_classic(full_width = F) %>%
  save_kable(file = paste0(appendix, "_MR.Results.Table.html"), self_contained = T)

####Exposure without SNPs####
all.exposures <- unique(expo$Phenotype) %>% gsub(" ", ".", .)
included.exposures <- unique(mr_results$exposure)
not.included.exposures <- all.exposures[-(which(all.exposures %in% included.exposures))]
not.included.exposures <- subset(expo[, c("Phenotype", "target_use")] %>% unique(), Phenotype %in% not.included.exposures)
if(length(not.included.exposures)>=1){
  write("\nThere are exposures not included in the final MR analysis:", file = "output.txt", append = T)
  write.table(not.included.exposures[,1], file = "output.txt", 
              append = T, row.names = F, col.names = F, quote = F)
  write.table(not.included.exposures, file = paste0(appendix, "_MR.Results.Exposures.NotIncluded.txt"), 
              row.names = F, col.names = F, quote = F)
}

####Print out information####
if(nrow(sig.fdr)>=1){
  write("\nNumber of FDR proteins:", file = "output.txt", append = TRUE)
  num <- subset(mr_results_ivw_only, adj_pval < 0.05) %>% nrow()
  write(num, file = "output.txt", append = TRUE)
}
write("\nNumber of exposures after clumping:", file = "output.txt", append = TRUE)
num <- length(unique(harmonized.df$exposure))
write(num, file = "output.txt", append = TRUE)

write("\nNumber of exposures for MR:", file = "output.txt", append = TRUE)
num <- subset(harmonized.df, mr_keep == TRUE) %>% .$exposure %>% unique() %>% length()
write(num, file = "output.txt", append = TRUE)

####Save Spot####
write("\nSaving data...", file = "output.txt", append = T)
save(expo.rf, outcome.rf, mr_results, overlap, all.afterclump, harmonized.df, included_snps,
     mr_results_ivw_only, mr_results_single, mr_pleio_results, mr_het_results, all.R2,
     pqtl_p, pqtl_maf, pqtl_distTSS, r2, kb, panel, appendix, mr_results_edit, mr_results_edit_print, ivw.print.final,
     file = paste0(appendix, "_MR.Analysis.rds"))
#save.image(paste0(appendix, "_Full.MR.Analysis.RData"))
if(file.exists(paste0(appendix, "_MR.Analysis.rds"))){
  write("\nMR analysis finished. Good job. Now go play some videogames!\n", file = "output.txt", append = T)
}
