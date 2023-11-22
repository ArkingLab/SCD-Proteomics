####Description####
## Add hg38 coordinates to 2018 SCA GWAS

####Load packages####
library(data.table)
library(rsnps)
library(magrittr)

####Read in summary statistics####
zz <- gzfile("arkinglab/static/papers/Ashar.Mitchell.SCA.2018/GWAS/data/METAANALYSIS_woHSDS.SUDS_Jun20151.tbl.gz")  
outcome <- read.table(zz, header=T)

####Update rsID for GWAS####
#rsids from GWAS was converted to hg38 position and dbsnp build 151 rsids so POS in chrpos is from hg38

chrpos <- fread("aric.prot.real.archive/scd.assoc/2_analyses/MR/analysis/SCA.GWAS.rsids_chrpos_ver2.txt", fill = TRUE, header =F)
names(chrpos) <- c("rsid", "chr.pos.ref.alt")
chrpos$CHROM <- gsub(":.*", "", chrpos$chr.pos.ref.alt) %>% gsub("chr", "", .)
chrpos$POS <- sub(":[^:]+$", "", chrpos$chr.pos.ref.alt) %>% sub(":[^:]+$", "", .) %>% gsub(".*:", "", .)
chrpos$REF <-  sub(":[^:]+$", "", chrpos$chr.pos.ref.alt) %>% gsub(".*:", "", .)
chrpos$ALT <-  sub(".*:", "", chrpos$chr.pos.ref.alt)
names(chrpos) <- c("MarkerName", "CHR.POS.REF.ALT", "CHR", "POS", "REF", "ALT")

####Merge file####
outcome <- merge(outcome, chrpos, by = "MarkerName")

outcome$Phenotype <- "SCA"
outcome$Freq2 <- 1-outcome$Freq1
outcome$MAF <- ifelse(outcome$Freq1 < outcome$Freq2, outcome$Freq1, outcome$Freq2)
outcome$Amin <- ifelse(outcome$Freq1 < outcome$Freq2, outcome$Allele1, outcome$Allele2)
outcome$Amaj <- ifelse(outcome$Freq1 > outcome$Freq2, outcome$Allele1, outcome$Allele2)

####Edit SNPs without POS####
outcome.missing.POS <- subset(outcome, POS == "")

## Manually edit VG SNPs, remove those where rsID can't be found
outcome.missing.POS$MarkerName[which(outcome.missing.POS$MarkerName == "VG01S1077")] <- "rs77931234"
outcome.missing.POS$MarkerName[which(outcome.missing.POS$MarkerName == "VG07S29628")] <- "rs78655421"
outcome.missing.POS.noVG <- outcome.missing.POS[-grep("VG", outcome.missing.POS$MarkerName),]

## Use dbSNP137 
dbSNP137 <- fread("arkinglab/resources/dbSNP/137/dbSNP137.chrpos.txt", skip = 58)

outcome.missing.POS.noVG.merged <- merge(outcome.missing.POS.noVG, dbSNP137, by.x = "MarkerName", by.y = "ID", all.x = TRUE)

## Liftover SNPs found
for.liftover <- subset(outcome.missing.POS.noVG.merged, !is.na(POS.y))
for.liftover$submit <- paste0("chr", for.liftover$`#CHROM`, ":", for.liftover$POS.y, "-", for.liftover$POS.y)
write.table(for.liftover$submit, file = "for.liftover.txt", row.names = F, col.names = F, sep = "\t", quote = F)

#Liftover done manually on website
liftover.results <- fread("hglft_genome_cd77_8ab340_liftover.results.bed", header = F)
names(liftover.results) <- "chr.hg38.pos"
liftover.fail <- fread("liftover.fail.txt", header = F)
liftover.fail <- subset(liftover.fail, V1 != "#Deleted in new")
#65121/65375 SNPs lifted over
#254 SNPs unable to be lifted over
#No way to merge back to original file, so remove the failed liftovers and re-do

liftover.convertible <- subset(for.liftover, !(submit %in% liftover.fail$V1))
write.table(liftover.convertible$submit, file = "for.liftover2.txt", row.names = F, col.names = F, sep = "\t", quote = F)
liftover.results2 <- fread("hglft_genome_12372_8ae6a0_liftover2.results.bed", header = F)

final.liftover <- cbind(liftover.convertible, liftover.results) %>% as.data.frame()
final.liftover$POS.hg38 <- gsub("-.*", "", final.liftover$chr.hg38.pos) %>% gsub(".*:", "", .)

#final.liftover contains the 65121 SNPs with chr, and hg38 Position

####Get chr, pos for remining SNPs####
## ~2500 SNPs still missing
still.missing <- subset(outcome.missing.POS.noVG.merged, is.na(outcome.missing.POS.noVG.merged$POS.y))
#Looks like these are merged to other IDs

# ptm <- proc.time()
# dbsnp_info <- ncbi_snp_query(weird$MarkerName[1:100])
# proc.time() - ptm

chunk_size <- 100
chunks <- ggplot2::cut_interval(1:nrow(still.missing), length=chunk_size, labels=FALSE)
still.missing$chunks <- chunks
datalist <- list()
for (i in unique(chunks)){
  print(paste0("Working on chunk ", i))
  query <- still.missing[which(still.missing$chunks==i), "MarkerName"]
  res <- ncbi_snp_query(query)
  datalist[[i]] <- res
}
#Chunk 3 has an issue, chunk 1 and 2 are done

chunks2 <- unique(chunks)[-c(1:3)]
for (i in unique(chunks2)){
  print(paste0("Working on chunk ", i))
  query <- still.missing[which(still.missing$chunks==i), "MarkerName"]
  res <- ncbi_snp_query(query)
  datalist[[i]] <- res
}
#Chunk 7,8 has an issue, chunks 4-6 are done

chunks3 <- unique(chunks2)[-c(1:5)]
for (i in unique(chunks3)){
  print(paste0("Working on chunk ", i))
  query <- still.missing[which(still.missing$chunks==i), "MarkerName"]
  res <- ncbi_snp_query(query)
  datalist[[i]] <- res
}
#Chunk 10 has an issue, chunk 9 is done
#start at chunk 11

chunks3 <- unique(chunks2)[-c(1:7)]
for (i in unique(chunks3)){
  print(paste0("Working on chunk ", i))
  query <- still.missing[which(still.missing$chunks==i), "MarkerName"]
  res <- ncbi_snp_query(query)
  datalist[[i]] <- res
}
#Chunk 13 has an issue, chunks 11,12 are done
#start at chunk 14

chunks3 <- unique(chunks2)[-c(1:10)]
for (i in unique(chunks3)){
  print(paste0("Working on chunk ", i))
  query <- still.missing[which(still.missing$chunks==i), "MarkerName"]
  res <- ncbi_snp_query(query)
  datalist[[i]] <- res
}
#Chunk 15,16,17,18 has an issue, chunks 14 are done
#start at chunk 19

chunks3 <- unique(chunks2)[-c(1:15)]
for (i in unique(chunks3)){
  print(paste0("Working on chunk ", i))
  query <- still.missing[which(still.missing$chunks==i), "MarkerName"]
  res <- ncbi_snp_query(query)
  datalist[[i]] <- res
}
#Chunk 24,25,26,27 has an issue

#chunks that worked: 1,2,4,5,6,9,11,12,14,19-23
#chunks that didn't work: 3,7,8,10,13,15-18,24-27

big_data <- do.call(rbind, datalist)

snps.found <- merge(still.missing, big_data[, c("query", "rsid", "chromosome", "bp")], by.x = "MarkerName", by.y = "query", all.x = T)
#1377/2697 SNPs found
names(snps.found)[31] <- "pos_hg38"
#will lose my mind if i lose this so save this please
saveRDS(big_data, "found.snps.RDS")

still.missing2 <- subset(snps.found, is.na(chromosome))
#1320 SNPs still missing
snps.found <- subset(snps.found, !is.na(chromosome))
#snps.found contains the 1377 SNPs where chromosome, pos hg38 were found

####Do last round of finding missing children####
#saveRDS(still.missing2, file = "still.missing2.RDS")
find.the.children <- as.data.frame(matrix(nrow = nrow(still.missing2), ncol = 4))
names(find.the.children) <- c("query", "rsid", "chromosome", "pos_hg38")
find.the.children$query <- still.missing2$MarkerName

for (i in 1:nrow(find.the.children)){
  query=find.the.children$query[i]
  if(i == 1){
    print("Begin SNP conversion...")
  }
  if(i%%50 == 0){ 
    print(paste0("Working on ", query, ", SNP #", i))}
  if(i == 200){
    saveRDS(find.the.children, "find.the.children.RDS")
  }
  skip_to_next <- FALSE
  tryCatch({
    res=ncbi_snp_query(query)
    find.the.children$chromosome[i] <- res$chromosome
    find.the.children$rsid[i] <- res$rsid
    find.the.children$pos_hg38[i] <- res$bp
  }, error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }
}

table(is.na(find.the.children$rsid))
#107 SNPs still unable to be converted - these children will perish

####Merge all children####
filled.snps <- subset(find.the.children, !is.na(rsid))
filled.snps <- final.liftover[, c(1,1,26,30)] %>% setNames(., c("query", "rsid", "chromosome", "pos_hg38")) %>% rbind(., filled.snps)
filled.snps <- snps.found[, c("MarkerName", "rsid", "chromosome", "pos_hg38")] %>% setNames(., c("query", "rsid", "chromosome", "pos_hg38")) %>% rbind(., filled.snps)
subset(filled.snps, duplicated(query))
#None
subset(filled.snps, duplicated(rsid))
#There are duplicated new rsids

####Fill final summary stat data frame####
outcome.merged <- merge(outcome, filled.snps, by.x = "MarkerName", by.y = "query", all.x = TRUE)
outcome.merged$POS <- ifelse(outcome.merged$POS=="", outcome.merged$pos_hg38, outcome.merged$POS)
outcome.merged$CHR <- ifelse(outcome.merged$CHR=="", outcome.merged$chromosome, outcome.merged$CHR)
outcome.merged$CHR <- ifelse(outcome.merged$CHR=="X", "23", outcome.merged$CHR)
outcome.merged$CHR.POS.MIN.MAJ <- paste0("chr", outcome.merged$CHR, ":", outcome.merged$POS, ":", outcome.merged$Amin, ":", outcome.merged$Amaj)

####Save spot####
final.outcome <- outcome.merged[, -c(16,19,20,26:28)]
saveRDS(final.outcome, "aric.prot.real.archive/scd.assoc/7_MR/SCA.Sum.Stats/SCA.2018.GWAS_SumStats_hg38.updated.RDS")
save.image("aric.prot.real.archive/scd.assoc/7_MR/SCA.Sum.Stats/Add.hg38.POS.to.SCA.GWAS_updated.RData")
