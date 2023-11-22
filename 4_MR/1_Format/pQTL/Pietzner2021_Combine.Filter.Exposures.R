####Description####
## Combine all protein exposures into an .RDS file to use for MR, after filtering out SNPs based on P-value cutoff (0.05)
## Convert coordinates to hg38, filter out aptamers with missing gene distance information

## filter.pqtl.array.sh and filter.pqtl.array_round2.sh were used to filter SNPs, based on files2filt.txt and files2filt_round2.txt files (file names unchanged from Pietzner 2021 download site)
## Meta.Min.Bonf.Prot.List.txt made with Make.Meta.Min.Bonf.Prot.Key_10072022.R to match to protein target names

####Load packages####
library(magrittr)
library(data.table)
library(dplyr)
library(readxl)
library(writexl)

####Combine all SNPs into one file####
key <- fread("aric.prot.real.archive/scd.assoc/7_MR/pQTL.databases/Meta.Min.Bonf.Prot.List.txt")
files <- list.files(pattern = "_05.txt")
files <- as.data.frame(files)
files$SeqId <- gsub("_Fenland_MA_auto_chrX_filtered_1pc_05.txt", "", files$files) %>% gsub("res_invn_X", "", .)
names(files)[1] <- "file"
key$SeqId <- gsub("-", "_", key$SeqId.use)
key <- merge(key, files, by = "SeqId", all.x = TRUE)
#EPHA10 has no file
key$gene <- key$ARIC.entrezgenesymbol
key$gene <- ifelse(key$ARIC.entrezgenesymbol == "C4A|C4B", "C4A.C4B", key$gene)
key.og <- key
key <- subset(key, !is.na(file))

#names(key) <- c("chr", "gene", "filename")
gene_list <- list()
for(i in 1:nrow(key)){
  target <- key$target_use[i]
  file <- key$file[i]
  chrom <- gsub("chr","",key$chr[i])
  seqid <- key$SeqId[i]
  geneid <- key$gene[i]
  #seqid.file <- gsub("_[A-Z].*","",key$file[i])
  #print(paste0("Working on ", file))
  if(i==1 | i%%20==0){
    print(paste0("Working on ", i, ": ",file))
  }
  if(file.exists((file))){
    data <- fread(file)
    data$target_use <- target
    data$Phenotype <- paste0(target, "_", seqid)
    data$SeqId.use <- seqid
    data$geneid <- geneid
    #Keep only SNPs on the same chromosome
    data <- subset(data, chr == chrom)
    assign(paste0(geneid, ".", seqid, ".pqtl"), data)
    check <- get(paste0(geneid, ".", seqid, ".pqtl"))
    gene_list[i] <- geneid
    if(!exists("check")){
      print(paste0("pQTL file for ", geneid, "-", seqid," does not exist. Loop will stop."))
    }
  } else {
    print("One file unexpectedly does not exist. Loop will quit.")
  }
}
df <- bind_rows(mget(ls(pattern = ".pqtl")))
df.og <- df

####Get hg38 position####
df$pos_hg19 <- gsub("[A-Z]", "", df$MarkerName) %>% gsub("_", "", .) %>% gsub(".*:", "",.)
bed <- data.frame(chr = paste0("chr", df$chr), pos = df$pos, pos1 = df$pos + 1) %>% unique()
bed$pos <- format(bed$pos, scientific = FALSE)
bed$pos1 <- format(bed$pos1, scientific = FALSE)
bed$pos <- gsub(" ", "", bed$pos)
bed$pos1 <- gsub(" ", "", bed$pos1)
write.table(bed, file = "aric.prot.real.archive/scd.assoc/7_MR/pQTL.databases/Pietzner.2021/p_05/liftover/Pietzner.p05.input.bed", row.names = F, col.names = F, quote = F, sep = "\t")
#Use run.Liftover.sh to get hg38 coordinates

output <- fread("aric.prot.real.archive/scd.assoc/7_MR/pQTL.databases/Pietzner.2021/p_05/liftover/Pietzner.p05.Liftover.Output.bed") 
unlifted <- read.table("aric.prot.real.archive/scd.assoc/7_MR/pQTL.databases/Pietzner.2021/p_05/liftover/Pietzner.p05.LiftOver.Unlifted.bed")

bed$round1 <- ifelse(bed$pos %in% unlifted$V2, "unlifted", "lifted")
df.lifted.round1 <- subset(df, !(pos_hg19 %in% unlifted$V2))

bed2 <- data.frame(chr = paste0("chr", df.lifted.round1$chr), pos = df.lifted.round1$pos, pos1 = df.lifted.round1$pos + 1) %>% unique()
bed2$pos <- format(bed2$pos, scientific = FALSE)
bed2$pos1 <- format(bed2$pos1, scientific = FALSE)
bed2$pos <- gsub(" ", "", bed2$pos)
bed2$pos1 <- gsub(" ", "", bed2$pos1)
write.table(bed2, file = "aric.prot.real.archive/scd.assoc/7_MR/pQTL.databases/Pietzner.2021/p_05/liftover/Pietzner.p05.input.round2.bed", row.names = F, col.names = F, quote = F, sep = "\t")
#Re-run Liftover using run.Liftover2.sh

output2 <- fread("aric.prot.real.archive/scd.assoc/7_MR/pQTL.databases/Pietzner.2021/p_05/liftover/Pietzner.p05.Liftover.Output.round2.bed") 
unlifted2 <- read.table("aric.prot.real.archive/scd.assoc/7_MR/pQTL.databases/Pietzner.2021/p_05/liftover/Pietzner.p05.LiftOver.Unlifted.round2.bed")
#No unlifted SNPs in second round

bed3 <- cbind(bed2, output2[, 2])
names(bed3) <- c("chr", "pos_hg19", "pos_hg19+", "pos_hg38")
bed3$chrpos <- paste0(bed3$chr, ":", bed3$pos_hg19)
df$chrpos <- paste0("chr", df$chr, ":", df$pos_hg19)

df2 <- merge(df, bed3[, c("chrpos", "pos_hg38")], by = "chrpos", all.x = TRUE)
#1940 SNPs with no hg38 position - kept in dataset but pos_hg38 column will be NA

####Add in other allele frequency#####
df2$Freq2 <- 1-df2$Freq1
df2$MAF <- ifelse(df2$Freq1<df2$Freq2,df2$Freq1,df2$Freq2)

####Add in TSS####
df3 <- merge(df2, key[, c("SeqId","strand_prot", "tss_prot")], by.y = "SeqId", by.x = "SeqId.use")
df3$distTSS <- abs(as.numeric(df3$tss_prot) - df3$pos_hg38)

####Create ID to merge to summary statistics####
df3$Amin <- ifelse(df3$Freq1<df3$Freq2, df3$Allele1, df3$Allele2) %>% toupper()
df3$Amaj <- ifelse(df3$Freq1>df3$Freq2, df3$Allele1, df3$Allele2) %>% toupper()
df3$Allele1 <- toupper(df3$Allele1)
df3$Allele2 <- toupper(df3$Allele2)

df3$CHR.POS_hg19.MIN.MAJ <- paste0(df3$chrpos, ":", df3$Amin, ":", df3$Amaj)
df3$CHR.POS_hg38.MIN.MAJ <- paste0("chr", df3$chr,":",df3$pos_hg38,":", df3$Amin, ":", df3$Amaj)
df3$uniqueid_hg19 <- paste0(df3$CHR.POS_hg19.MIN.MAJ, "_", df3$Phenotype)
df3$uniqueid_hg19 <- gsub(" ", ".", df3$uniqueid_hg19)
df3$uniqueid_hg38 <- paste0(df3$CHR.POS_hg38.MIN.MAJ, "_", df3$Phenotype)
df3$uniqueid_hg38 <- gsub(" ", ".", df3$uniqueid_hg38)

####Remove INDELs####
df3 <- subset(df3, Allele1 != "D")
#Need to know allele...D/I only labeled as such so can't use

####Create pseudo P value for clumping####
min(df3$Pvalue)
#2.98e-308
#P value still able to be printed, will not make a pseudo P value

####Add in gene body start and end####
suppt3 <- read_xlsx("aric.prot.real.archive/scd.assoc/7_MR/pQTL.databases/Pietzner.2021/Pietzner.Supp.T3.xlsx") %>% as.data.frame()
ref <- fread("aric.prot.real.archive/scd.assoc/7_MR/pQTL.databases/Meta.Min.Bonf.Prot.List.txt")

## Merge SeqIds
seqids <- df3[, c("SeqId.use", "target_use")] %>% unique()
seqids <- ref %>% mutate(SeqId.use = gsub("-", "_", SeqId.use)) %>% .[, c("SeqId.use", "ARIC.uniprot_id")] %>% merge(., seqids, by = "SeqId.use")
#128/129
suppt3$SeqId.use <- gsub("SeqId_", "", suppt3$SomaScan.ID)
merged <- merge(seqids, suppt3, by = "SeqId.use", all.x = T)
write_xlsx(merged, "aric.prot.real.archive/7_MR/pQTL.databases/Pietzner.2021/p_05/Pietzner.2021.Edited.SuppT3.Merged.xlsx")

## Get start, end for NA SeqIds
#For SeqIds without start and end pos, manually found them on Biomart
updated.locations <- read_xlsx("aric.prot.real.archive/7_MR/pQTL.databases/Pietzner.2021/p_05/Pietzner.2021.Edited.SuppT3.Merged.xlsx", sheet = "Sheet2") %>% as.data.frame()
names(updated.locations) <- paste0(names(updated.locations), "_new")
updated.merged <- merge(merged, updated.locations, by.x = "SeqId.use", by.y = "SeqId.use_new", all.x = T)
updated.merged$Chr <- ifelse(is.na(updated.merged$Chr), updated.merged$Chr_new, updated.merged$Chr)
updated.merged$start <- ifelse(is.na(updated.merged$start), updated.merged$start_new, updated.merged$start)
updated.merged$end <- ifelse(is.na(updated.merged$end), updated.merged$end_new, updated.merged$end)
final.merged <- updated.merged[, -c(25:29)]
final.short <- updated.merged[, c("SeqId.use", "Chr", "start", "end")]
names(final.short) <- c("SeqId.use", "prot_chr", "prot_start", "prot_end")
pqtl.updated <- merge(df3, final.short, by = "SeqId.use")
names(pqtl.updated)[39:40] <- c("prot_start_hg19", "prot_end_hg19")

####Remove C4A/C4B, KIAA0040####
## No idea how to calculate distance since not sure which is the correct target gene
pqtl.updated.filt <- subset(pqtl.updated, geneid != "C4A.C4B")

## Ambiguous start position
pqtl.updated.filt <- subset(pqtl.updated.filt, geneid != "KIAA0040")
## This leaves 126 aptamers for analysis

####Get distance from gene body####
pqtl.updated.filt$prot_start_hg19 <- as.numeric(pqtl.updated.filt$prot_start_hg19)
pqtl.updated.filt$prot_end_hg19 <- as.numeric(pqtl.updated.filt$prot_end_hg19)
pqtl.updated.filt$pos_hg19 <- as.numeric(pqtl.updated.filt$pos_hg19)
pqtl.updated.filt$distStart <- abs(pqtl.updated.filt$prot_start_hg19 - pqtl.updated.filt$pos_hg19)
pqtl.updated.filt$distEnd <- abs(pqtl.updated.filt$prot_end_hg19 - pqtl.updated.filt$pos_hg19)
pqtl.updated.filt$cisStart <- ifelse(pqtl.updated.filt$distStart<500000, "cis", "trans")
pqtl.updated.filt$cisEnd <- ifelse(pqtl.updated.filt$distEnd<500000, "cis", "trans")
pqtl.updated.filt$cisStatus <- ifelse(pqtl.updated.filt$cisStart == "cis" | pqtl.updated.filt$cisEnd == "cis", "cis", "trans")

####Write out final file####
printout <- pqtl.updated.filt[,c("SeqId.use", "rsid", "Allele1", "Allele2", 
                              "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Freq2", "MAF",
                              "Effect", "StdErr", "Pvalue", "Direction", "HetISq",
                              "HetChiSq", "HetDf", "HetPVal", "TotalSampleSize",
                              "chr", "target_use", "Phenotype", "geneid",
                              "pos_hg19", "pos_hg38", "strand_prot", "distTSS",
                              "Amin", "Amaj", "CHR.POS_hg19.MIN.MAJ", "CHR.POS_hg38.MIN.MAJ",
                              "uniqueid_hg19", "uniqueid_hg38", "prot_start_hg19",
                              "prot_end_hg19", "distStart", "distEnd", 
                              "cisStart", "cisEnd", "cisStatus")]
saveRDS(printout, "aric.prot.real.archive/scd.assoc/7_MR/pQTL.databases/Pietzner.2021/p_05/All.Exposures_Filtered.RDS")

####Save data####
save.image("aric.prot.real.archive/scd.assoc/7_MR/pQTL.databases/Pietzner.2021/p_05/Pietzner2021_Combine.Filter.Exposure.RData")
