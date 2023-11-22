####Description####
## Script runs COXPH model to examine association of ARIC proteomics with SCA. Script can be used for either visits and will run the following models:
## Mod 1 (Dz): [race+sex+]age+bmi+height+center+smoke+eGFR_ckdepi+[heart disease options]
## Mod 2 (Min): [race+sex+]age+bmi+height+center+smoke+eGFR_ckdepi (THIS IS THE PRIMARY MODEL)
## Mod 3 (RskFct): [race+sex+]age+bmi+height+center+smoke+eGFR_ckdepi+[heart disease options]+diabetes+sbp+hnt_med+hdl+total_chol+trigycerides+stroke+AF (THIS IS THE FULL MODEL)

## Version 8 uses the CHD definition that has been harmonized with CHS. This new variable is called chd_chs. It differs from Version 7 because those who have angioplasty on their legs are not counted as having CHD.
## The ARIC CHD variable is called chd.

## This version also saves all model results into a list for further downstream analysis.

####Load packages####
library(magrittr)
library(survival)
library(writexl)
library(kableExtra)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(data.table)
library(stringr)

####Load functions####
drug_cleaner <- function(data, drug_list){
  drugs <- drug_list
  drugs_tbr <- list()
  for(i in 1:length(drugs)){
    drug_test <- drugs[i]
    if(dim(table(as.character(data[, drug_test]))) == 1){
      drugs_tbr[i] <- drug_test
    } else {
      drugs_tbr[i] <- NA
    }
  }
  drugs_tbr <- drugs_tbr %>% as.character() %>% .[. != "NA"]
  drugs <- drugs[!(drugs %in% drugs_tbr)]
  drugs <- str_c(drugs,collapse="+")
  return(drugs)
}

make.qqplot <- function(df, beta.name, se.name, p.name, title){
  z.sq <- (df[,beta.name] / df[,se.name])^2
  lambda<- (median(z.sq,na.rm=TRUE) / 0.456) %>% round(., digits = 3)
  pvals <- df[, p.name]
  observed <- sort(pvals)
  lobs <- -(log10(observed))
  expected <- c(1:length(observed)) 
  lexp <- -(log10(expected / (length(expected)+1)))
  s<-paste("lambda=",lambda,sep="")
  plot(c(0,11), c(0,11), col="red", lwd=3, type="l", xlab="Expected (-logP)",
       ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l", sub=s, main=title)
  points(lexp, lobs, pch=23, cex=.5, col="black", bg="black")
}

run.sca.coxph.list <- function(outcome, model, covariates){
  
  return_list <- list()
  
  ## Make results table
  results <- annot.final
  results$beta <- NA
  results$SE <- NA
  results$P <- NA
  results$n_total <- nrow(final.alldata)
  results$mod <- model
  results$outcome <- outcome
  
  ## Outcome choice
  if(outcome == "normal"){
    y = "Surv(fuday_age_censored, scd_age_filt)" 
    results$n_cases <- subset(final.alldata, scd_age_filt == 1) %>% nrow()
    results$n_ctrls <- subset(final.alldata, scd_age_filt == 0) %>% nrow()
  } else if(outcome == "five"){
    y = "Surv(FIVEYR_SCD_fuyear, FIVEYR_SCD)"
    results$n_cases <- subset(final.alldata, FIVEYR_SCD == 1) %>% nrow()
    results$n_ctrls <- subset(final.alldata, FIVEYR_SCD == 0) %>% nrow()
  } else if (outcome == "ten"){
    y = "Surv(TENYR_SCD_fuyear, TENYR_SCD)"
    results$n_cases <- subset(final.alldata, TENYR_SCD == 1) %>% nrow()
    results$n_ctrls <- subset(final.alldata, TENYR_SCD == 0) %>% nrow()
  }
  
  ## Run model
  print(paste0("Running ", outcome, " ", model))
  for(i in 1:nrow(results)){
    protein <- results$seqid_in_sample[i]
    prot.run <- grep(protein, names(final.alldata))
    x = paste0(protein, " + ", covariates)
    form = as.formula(paste(y, "~", x))
    if(i == 1){
      print(form)
    }
    if(i%%500 == 0){
      print(paste0("Working on protein # ", i))
    }
    test <- coxph(form, data = final.alldata)
    return_list[[i]] <- test
    results$beta[i] <- summary(test)$coefficients[1,1]
    results$SE[i] <- summary(test)$coefficients[1,3]
    results$P[i] <- summary(test)$coefficients[1,5]
  }
  print(paste0("Finish running ", outcome, " ", model))
  first.prot <- grep("Seq", names(final.alldata))[1]
  last.prot <- grep("Seq", names(final.alldata)) %>% tail() %>% .[6]
  results$PASSED <- ifelse(results$P < (0.05/(ncol(final.alldata[, first.prot:last.prot]))), "YES", "NO")
  results <- results[order(results$P),]
  i <- nrow(results)
  return_list[[i+1]] <- results
  return(return_list)
}

mrg.mkr <- function(data,model,outcome){
  final <- data[, c("seqid_in_sample", "beta", "SE", "P")]
  names(final)[2] <- "Beta"
  if(outcome == 0){
    names(final)[2:4] <- paste0("M", model, ".", names(final)[2:4])
  } else {
    names(final)[2:4] <- paste0("M", model, outcome, ".", names(final)[2:4])
  }
  return(final)
}

prettyPresOnly <- function(results){
  present <- results[, c("seqid_in_sample", "target", "beta", "SE", "P", "flag3", "multi.apt", "PASSED", "n_total", "n_cases", "n_ctrls", "outcome", "mod")]
  names(present) <- c("SeqId", "Protein", "Beta", "SE", "P", "Best", "Multi", "PassBonf", "n", "Cases", "Ctrls", "Outcome", "Model")
  present$Best <- ifelse(present$Best == 0, "Y", "N")
  present$Multi <- ifelse(present$Multi == 0, "N", "Y")
  present$Beta <- round(present$Beta, digits = 3)
  present$SE <- round(present$SE, digits = 3)
  present$P <- signif(present$P, digits = 3)
  return(present)
}

prettyPresMany <- function(results){
  results[,grep("Beta", names(results))] <- round(results[,grep("Beta", names(results))], digits = 3)
  results[,grep("SE$", names(results))] <- round(results[,grep("SE$", names(results))], digits = 3)
  results[,grep("P$", names(results))] <- signif(results[,grep("P$", names(results))], digits = 3)
  names(results)[grep("seqid", names(results))] <- "SeqId"
  names(results)[grep("target", names(results))] <- "Protein"
  names(results)[grep("flag3", names(results))] <- "Best"
  results$Best <- ifelse(results$Best == 0, "Y", "N")
  names(results)[grep("multi.apt", names(results))] <- "Multi"
  results$Multi <- ifelse(results$Multi == 0, "N", "Y")
  names(results)[grep("n_total", names(results))] <- "n"
  names(results)[grep("cases", names(results))] <- "Cases"
  names(results)[grep("ctrls", names(results))] <- "Ctrls"
  names(results)[grep("outcome", names(results))] <- "Out"
  return(results)
}

kableModCombo <- function(data, savename, final_dz_covar_print, exclude_vars){
  names(data) <- gsub("M[0-9][0-9][0-9][.]", "", names(data))
  names(data) <- gsub("M[0-9][0-9][.]", "", names(data))
  names(data) <- gsub("M[0-9][.]", "", names(data))
  data[,1] <- gsub("SeqId_", "", data[,1])
  exclusions <- paste0(exclude_vars, collapse = ",")
  if(length(final_dz_covar_print) == 0){
    final_dz_covar_print <- "none"
    exclusions <- "all dz"
  }
  dzO_name <- paste0("Dz:", final_dz_covar_print, "\n[exclude:", exclusions, "]")
  data %>% kable(., row.names = F, align = c('l', 'l', rep('r', 15)), digits = 50) %>% 
    kable_styling(c("striped", "hover", "condensed"), full_width = F, font_size = 13, fixed_thead = TRUE) %>% 
    column_spec(5, color = ifelse(data[,5] < bonf, "red", "#73797D"), border_right = T) %>%
    column_spec(8, color = ifelse(data[,8] < bonf, "red", "#73797D"), border_right = T) %>%
    column_spec(11, color = ifelse(data[,11] < bonf, "red", "#73797D"), border_right = T) %>%
    # column_spec(14, color = ifelse(data[,14] < bonf, "red", "#73797D"), border_right = T) %>%
    column_spec(2, width = "9em") %>% 
    #add_header_above(c(" " = 2, "Min" = 3, setNames(3,dzO_name), "RskFct" = 3, "AF" = 3, " " = 6), line = TRUE, font_size = 13) %>%
    add_header_above(c(" " = 2, "Min" = 3, 
                       setNames(3,dzO_name),
                       "RskFct" = 3, " " = 6), 
                     line = TRUE, font_size = 13) %>% 
    save_kable(file = paste0(dir_out, savename, "_", appendix, ".html"), self_contained = T)
}

kableOutComboSimp <- function(df, model_num, savename, final_dz_covar_print, exclude_vars){
  cols <- grep(paste0("^", model_num), names(df))
  best <- grep("Best", names(df))
  multi <- grep("Multi", names(df))
  n <- grep("^n$", names(df))
  rear <- c(best, multi, n)
  data <- df[, c(1,2,cols,rear)]
  data <- data[order(data[,5]),]
  names(data) <- gsub("M[0-9][0-9][0-9][.]", "", names(data))
  names(data) <- gsub("M[0-9][0-9][.]", "", names(data))
  names(data) <- gsub("M[0-9][.]", "", names(data))
  
  model <- gsub(paste0("ARIC.SCA.V", visit, ".ALL_"), "", savename)
  exclusions <- paste0(exclude_vars, collapse = ",")
  if(length(final_dz_covar_print) == 0){
    final_dz_covar_print <- "none"
    exclusions <- "all dz"
  }
  if(model == "Dz" | model == "RskFct"){
    model <- paste0(model, ":", final_dz_covar_print, "\n[exclude:", exclusions, "]")
  } else {
    model <- paste0(model, ": \n[exclude:", exclusions, "]")
  }
  
  data <- data[1:150,]
  data[,1] <- gsub("SeqId_", "", data[,1])
  if(test == "Yes"){
    data <- data[1:3,]
  }
  data %>% kable(., row.names = F, align = c('l', 'l', rep('c', 12)), digits = 50) %>% 
    kable_styling(c("striped", "hover", "condensed"), full_width = F, font_size = 13, fixed_thead = TRUE) %>% 
    column_spec(5, color = ifelse(data[,5] < bonf, "red", "#73797D"), border_right = T) %>% 
    column_spec(8, color = ifelse(data[,8] < bonf, "red", "#73797D"), border_right = T) %>% 
    column_spec(11, color = ifelse(data[,11] < bonf, "red", "#73797D"), border_right = T) %>%
    column_spec(2, width = "10em") %>% 
    add_header_above(c(setNames(2, model), "Norm" = 3, "10 YR" = 3, "5 YR" = 3, " " = 3), line = TRUE, font_size = 12) %>% 
    save_kable(file = paste0(dir_out, savename, "_", appendix, ".html"), self_contained = T)
}

firaga <- function(df, title){
  ## Colors
  colors <- brewer.pal(9, "Set1")
  colors <- colors[-6]
  colors2 <- brewer.pal(12, "Set3")
  a <- "Min"
  b <- "Dz"
  c <- "RskFct"
  data <- df
  
  ## Label assignment
  data$pcol <- NA
  data$pcol <- ifelse(data[,5] < bonf & data[,8] >= bonf & data[,11] >= bonf, paste0(a), data$pcol)
  data$pcol <- ifelse(data[,5] >= bonf & data[,8] < bonf & data[,11] >= bonf, paste0(b), data$pcol)
  data$pcol <- ifelse(data[,5] >= bonf & data[,8] >= bonf & data[,11] < bonf, paste0(c), data$pcol)
  data$pcol <- ifelse(data[,5] < bonf & data[,8] < bonf & data[,11] >= bonf, paste0(a,"+",b), data$pcol)
  data$pcol <- ifelse(data[,5] >= bonf & data[,8] < bonf & data[,11] < bonf, paste0(b,"+",c), data$pcol)
  data$pcol <- ifelse(data[,5] < bonf & data[,8] >= bonf & data[,11] < bonf, paste0(a,"+",c), data$pcol)
  data$pcol <- ifelse(data[,5] < bonf & data[,8] < bonf & data[,11] < bonf, "All", data$pcol)
  data$pcol <- ifelse(data[,5] >= bonf & data[,8] >= bonf & data[,11] >= bonf, "None", data$pcol)
  
  data$label <- ifelse(data$pcol != "None", data$Protein, NA)
  data$label <- ifelse(nchar(data$label)>10, NA, data$label)
  data$label <- ifelse(data$SeqId == "SeqId_7655_11", "N-terminal pro-BNP", data$label)
  
  ## Color assignment
  data$col.assign <- NA
  data$col.assign <- ifelse(data$pcol == paste0(a), colors[1], data$col.assign)
  data$col.assign <- ifelse(data$pcol == paste0(b), colors[5], data$col.assign)
  data$col.assign <- ifelse(data$pcol == paste0(c), colors[4], data$col.assign)
  data$col.assign <- ifelse(data$pcol == paste0(a,"+",b), colors[2], data$col.assign)
  data$col.assign <- ifelse(data$pcol == paste0(b,"+",c), colors[6], data$col.assign)
  data$col.assign <- ifelse(data$pcol == paste0(a,"+",c), colors[7], data$col.assign)
  data$col.assign <- ifelse(data$pcol == "All", colors[3], data$col.assign)
  data$col.assign <- ifelse(data$pcol == "None", colors2[9], data$col.assign)
  
  ## Count
  data <- data %>% add_count(pcol) %>% mutate(znn = paste0(pcol, ' (', nn, ')')) 
  uniq.col.assign <- unique(data$col.assign)
  uniq.znn <- unique(data$znn)
  names(uniq.col.assign) <- uniq.znn
  
  ## Plot
  bonf <- -log10(0.05/4955)
  ggplot(data, aes(x=data[,3], y=(as.numeric(data[,5]) %>% -log10(.)), col=znn, label=label)) + 
    geom_point(size = 0.8, alpha = 0.8) + geom_line(y=bonf, linetype = "dashed", color = "black", size = 0.1) + 
    geom_text_repel(size = 1.8, max.overlaps = 8, col = "black", min.segment.length = 0.1, hjust = 0.5, nudge_y=0.3, nudge_x = 0.01, segment.size = 0.1) + 
    ylab("-log10(Min P)") + xlab("Min Beta") + 
    scale_color_manual(values = uniq.col.assign, 
                       breaks = c(names(uniq.col.assign)[grep("None [(]", names(uniq.col.assign))], 
                                  names(uniq.col.assign)[grep("^Min [(]", names(uniq.col.assign))], 
                                  names(uniq.col.assign)[grep("Dz [(]", names(uniq.col.assign))], 
                                  names(uniq.col.assign)[grep("^RskFct [(]", names(uniq.col.assign))], 
                                  names(uniq.col.assign)[grep("Min[+]Dz [(]", names(uniq.col.assign))], 
                                  names(uniq.col.assign)[grep("Min[+]RskFct [(]", names(uniq.col.assign))], 
                                  names(uniq.col.assign)[grep("Dz[+]RskFct [(]", names(uniq.col.assign))], 
                                  names(uniq.col.assign)[grep("All [(]", names(uniq.col.assign))])) + 
    ggtitle(paste0("ARIC SCA ", title)) + theme_classic() + theme(legend.title = element_blank())
}

####Arguments####
args <- commandArgs(trailingOnly=TRUE)
#appendix <- args[1]
input <- args[2]
test <- args[3]

dir_out <- paste0(getwd(), "/")
if(grepl("combined", dir_out)){
    race_strat <- "Both"
  } else if(grepl("AA", dir_out)){
    race_strat <- "Black"
  } else if(grepl("EA", dir_out)){
    race_strat <- "White"
}
if(grepl("female", dir_out)){
    sex_strat <- "F"
  } else if(grepl("male", dir_out)){
    sex_strat <- "M"
  } else {
    sex_strat <- "Both"
}
if(grepl("visit2", dir_out)){
  visit <- 2
} else if(grepl("visit3", dir_out)){
  visit <- 3
}

if(grepl("combined", dir_out) & grepl("postStrat", dir_out) & grepl("EA", dir_out)){
  race_strat <- "White"
} else if(grepl("combined", dir_out) & grepl("postStrat", dir_out) & grepl("AA", dir_out)){
  race_strat <- "Black"
}

if(race_strat == "Both" & sex_strat == "Both"){
  print.covar <- "sex + race + age + bmi + height + center + smoke + eGFR_ckdepi"
} else if(race_strat == "Both" & sex_strat != "Both"){
  print.covar <- "race + age + bmi + height + center + smoke + eGFR_ckdepi"
} else if(race_strat != "Both" & sex_strat == "Both"){
  print.covar <- "sex + age + bmi + height + center + smoke + eGFR_ckdepi"
} else if(race_strat != "Both" & sex_strat != "Both"){
  print.covar <- "age + bmi + height + center + smoke + eGFR_ckdepi"
}
print(paste0("Race stratification: ", race_strat, "; Sex stratification: ", sex_strat, "; Visit: ", visit, "; Minimum covariates: ", print.covar))

####Load data####
load(input)
appendix <- args[1]
print(paste0("Appendix: ", appendix))
row.names(final.alldata) <- final.alldata$pid

## Fix 5 and 10 YR SCA risk column name
names(final.alldata)[grep("5YR", names(final.alldata))] <- gsub("5YR", "FIVEYR", names(final.alldata)[grep("5YR", names(final.alldata))])
names(final.alldata)[grep("10YR", names(final.alldata))] <- gsub("10YR", "TENYR", names(final.alldata)[grep("10YR", names(final.alldata))])

## Read in CHS harmonized CHD variable
new.chd.var <- fread("aric.prot.real.archive/scd.assoc/1_data.cleaning/ARIC_V1-3_CHD_CHS.Harmonized.txt") %>% as.data.frame()
new.chd.var <- new.chd.var[new.chd.var$exam == visit, c("ID", "chd_chs")]
names(new.chd.var)[1] <- "pid"
final.alldata <- merge(new.chd.var, final.alldata, by = "pid")

## Change covariate class
final.alldata$sex <- as.factor(final.alldata$sex)
final.alldata$age <- as.numeric(final.alldata$age)
final.alldata$bmi <- as.numeric(final.alldata$bmi)
final.alldata$race <- as.factor(final.alldata$race)
final.alldata$center <- as.factor(final.alldata$center)
final.alldata$eGFR_ckdepi <- as.numeric(final.alldata$eGFR_ckdepi)

final.alldata$smoke <- factor(final.alldata$smoke, levels = c("N", "C", "F"))
print("Smoking:")
print(table(final.alldata$smoke))

final.alldata$sbp <- as.numeric(final.alldata$sbp)

final.alldata$hnt_med <- as.factor(final.alldata$hnt_med)
print("Hypertension medication:")
print(table(final.alldata$hnt_med))

final.alldata$diabetes <- as.factor(final.alldata$diabetes)
print("Diabetes:")
print(table(final.alldata$diabetes))

final.alldata$height <- as.numeric(final.alldata$height)
final.alldata$weight <- as.numeric(final.alldata$weight)
final.alldata$hdl <- as.numeric(final.alldata$hdl)
final.alldata$total_chol <- as.numeric(final.alldata$total_chol)
final.alldata$tri <- as.numeric(final.alldata$tri)

# final.alldata$cvd <- as.factor(final.alldata$cvd)
# print("CVD:")
# print(table(final.alldata$cvd))

final.alldata$chd_chs <- as.factor(final.alldata$chd_chs)
print("CHD (CHS harmonized):")
print(table(final.alldata$chd_chs))

final.alldata$stroke <- as.factor(final.alldata$stroke)
print("Stroke:")
print(table(final.alldata$stroke))

final.alldata$hf <- as.factor(final.alldata$hf)
print("Heart failure:")
print(table(final.alldata$hf))

# final.alldata$mi <- as.factor(final.alldata$mi)
# print("MI:")
# print(table(final.alldata$mi))

final.alldata$afib <- ifelse(final.alldata$afib == 0, "N", "Y")
final.alldata$afib <- as.factor(final.alldata$afib)
print("AF:")
print(table(final.alldata$afib))

first.prot <- grep("Seq", names(final.alldata))[1]
last.prot <- grep("Seq", names(final.alldata)) %>% tail() %>% .[6]
print("Finish loading and editing data")

####Covariates####

## Disease covariates
#dir_out <- "aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/visit3/combined/exclude.egfr.chd.hf.strk.dia/"
exclude_vars <- gsub(".*exclude.", "", dir_out)
exclude_vars <- gsub("/.*", "", exclude_vars) %>% strsplit(., "[.]") %>% as.data.frame()
names(exclude_vars) <- "exclusions"
#all_dz_covar <- c("egfr", "chd", "hf", "stroke", "diabetes")
#all_dz_covar <- c("egfr", "chd", "hf")
all_dz_covar <- c("egfr", "chd_chs", "hf")
exclude_vars <- exclude_vars$exclusions
exclude_vars_print <- exclude_vars

# if(!is.na(match("dia", exclude_vars))){
#   exclude_vars[grep("dia", exclude_vars)] <- "diabetes"
# }
# if(!is.na(match("strk", exclude_vars))){
#   exclude_vars[grep("strk", exclude_vars)] <- "stroke"
# }
if(!is.na(match("chd", exclude_vars))){
   exclude_vars[grep("chd", exclude_vars)] <- "chd_chs"
}

final_dz_covar <- all_dz_covar[-which(all_dz_covar %in% exclude_vars)]
if(length(final_dz_covar) != 0){
    final_dz_covar <- drug_cleaner(final.alldata, final_dz_covar)
}
final_dz_covar_print <- final_dz_covar
final_dz_covar <- paste0("+", final_dz_covar)
if(final_dz_covar == "+"){
    final_dz_covar <- ""
}
#else{
  #final_dz_covar_print <- gsub("stroke", "strk", final_dz_covar_print)
  #final_dz_covar_print <- gsub("diabetes", "dia", final_dz_covar_print)
#}                                

## Final model covariates
mod1_covar <- paste0(print.covar, final_dz_covar)
mod2_covar <- print.covar
mod3_covar <- paste0(print.covar, final_dz_covar, "+diabetes+sbp+hnt_med+hdl+total_chol+tri+stroke+afib")
# mod4_covar <- paste0(print.covar, final_dz_covar, "+afib")

print(paste0("Exclusions: ", paste0(exclude_vars_print, collapse = ",")))
print(paste0("Dz covariates: ", mod1_covar))

## To test
if(test == "Yes"){
final.alldata.og <- final.alldata
annot.final.og <- annot.final
annot.final <- annot.final[1:10,]
final.alldata <- final.alldata[,1:115]
}

####Normal COXPH####
mod1_norm_results_list <- run.sca.coxph.list(outcome = "normal", model = "Dz", covariates = mod1_covar)
mod1_norm_results <- mod1_norm_results_list[[length(mod1_norm_results_list)]] %>% as.data.frame()

mod2_norm_results_list <- run.sca.coxph.list(outcome = "normal", model = "Min", covariates = mod2_covar)
mod2_norm_results <- mod2_norm_results_list[[length(mod2_norm_results_list)]] %>% as.data.frame()

mod3_norm_results_list <- run.sca.coxph.list(outcome = "normal", model = "RskFct", covariates = mod3_covar)
mod3_norm_results <- mod3_norm_results_list[[length(mod3_norm_results_list)]] %>% as.data.frame()

## Merge results
norm_noresults_main_head <- mod1_norm_results[, c("seqid_in_sample", "target")]
norm_noresults_main_tail <- mod1_norm_results[, c("seqid_in_sample", "flag3", "multi.apt", "n_total", "mod", "outcome", "n_cases", "n_ctrls")]
norm_noresults_other <- mod1_norm_results[, c("seqid_in_sample", "uniprot_id", "uniprot.full.name", "flag1",
                                              "flag2", "SomaCVBA_v2", "SomaCVBA_v3", "SomaCVBA_v5", "ARICCVBA_v2", "ARICCVBA_v3",
                                              "ARICCVBA_v5", "ARICpcorr_v2", "ARICpcorr_v3", "ARICpcorr_v5", "SOMAtotal.cv.plasma", "Plateflag_any_v2",
                                              "Plateflag_any_v3", "Plateflag_any_v5", "Plateflag_n_v2", "Plateflag_n_v3", "Plateflag_n_v5", "targetfullname",
                                              "apparent.kd..m.", "status", "characterization.info", "entrez.gene.id", "entrezgeneid", "entrezgenesymbol", "organism", "type")]
mod1_norm_results_base <- mod1_norm_results
names(mod1_norm_results_base)[grep("beta", names(mod1_norm_results_base))] <- "M1.Beta"
names(mod1_norm_results_base)[grep("^SE$", names(mod1_norm_results_base))] <- "M1.SE"
names(mod1_norm_results_base)[grep("^P$", names(mod1_norm_results_base))] <- "M1.P"
mod1_norm_merge <- mod1_norm_results_base[, c("seqid_in_sample", "M1.Beta", "M1.SE", "M1.P")]
mod2_norm_merge <- mrg.mkr(mod2_norm_results,2,0)
mod3_norm_merge <- mrg.mkr(mod3_norm_results,3,0)

all_norm_models <- Reduce(function(x,y) merge(x,y,by="seqid_in_sample",all=TRUE), list(mod2_norm_merge, mod1_norm_merge, mod3_norm_merge))
all_norm_results <- merge(norm_noresults_main_head, all_norm_models, by = "seqid_in_sample")
all_norm_results <- merge(all_norm_results, norm_noresults_main_tail, by = "seqid_in_sample")
all_norm_results <- merge(all_norm_results, norm_noresults_other, by = "seqid_in_sample")
all_norm_results$mod <- NULL
all_norm_results$PASSED <- NULL

save(mod1_norm_results_list, mod2_norm_results_list, mod3_norm_results_list, mod1_norm_results, mod2_norm_results, mod3_norm_results, 
     all_norm_results, final.alldata, file = paste0(dir_out, "ARIC.SCA.V", visit, ".Final_NormMods_", appendix, ".rds"))

tmp <- all_norm_results
tmp$exam <- visit
assign(paste0("aricv", visit, "_all_norm_results"), tmp)

####5YR COXPH####
mod1_five_results_list <- run.sca.coxph.list(outcome = "five", model = "Dz", covariates = mod1_covar)
mod1_five_results <- mod1_five_results_list[[length(mod1_five_results_list)]] %>% as.data.frame()

mod2_five_results_list <- run.sca.coxph.list(outcome = "five", model = "Min", covariates = mod2_covar)
mod2_five_results <- mod2_five_results_list[[length(mod2_five_results_list)]] %>% as.data.frame()

mod3_five_results_list <- run.sca.coxph.list(outcome = "five", model = "RskFct", covariates = mod3_covar)
mod3_five_results <- mod3_five_results_list[[length(mod3_five_results_list)]] %>% as.data.frame()

## Merge results
five_noresults_main_head <- mod1_five_results[, c("seqid_in_sample", "target")]
five_noresults_main_tail <- mod1_five_results[, c("seqid_in_sample", "flag3", "multi.apt", "n_total", "mod", "outcome", "n_cases", "n_ctrls")]
five_noresults_other <- mod1_five_results[, c("seqid_in_sample", "uniprot_id", "uniprot.full.name", "flag1",
                                              "flag2", "SomaCVBA_v2", "SomaCVBA_v3", "SomaCVBA_v5", "ARICCVBA_v2", "ARICCVBA_v3",
                                              "ARICCVBA_v5", "ARICpcorr_v2", "ARICpcorr_v3", "ARICpcorr_v5", "SOMAtotal.cv.plasma", "Plateflag_any_v2",
                                              "Plateflag_any_v3", "Plateflag_any_v5", "Plateflag_n_v2", "Plateflag_n_v3", "Plateflag_n_v5", "targetfullname",
                                              "apparent.kd..m.", "status", "characterization.info", "entrez.gene.id", "entrezgeneid", "entrezgenesymbol", "organism", "type")]
mod1_five_results_base <- mod1_five_results
names(mod1_five_results_base)[grep("beta", names(mod1_five_results_base))] <- "M15.Beta"
names(mod1_five_results_base)[grep("^SE$", names(mod1_five_results_base))] <- "M15.SE"
names(mod1_five_results_base)[grep("^P$", names(mod1_five_results_base))] <- "M15.P"
mod1_five_merge <- mod1_five_results_base[, c("seqid_in_sample", "M15.Beta", "M15.SE", "M15.P")]
mod2_five_merge <- mrg.mkr(mod2_five_results,2,5)
mod3_five_merge <- mrg.mkr(mod3_five_results,3,5)
all_five_models <- Reduce(function(x,y) merge(x,y,by="seqid_in_sample",all=TRUE), list(mod2_five_merge, mod1_five_merge, mod3_five_merge))
all_five_results <- merge(five_noresults_main_head, all_five_models, by = "seqid_in_sample")
all_five_results <- merge(all_five_results, five_noresults_main_tail, by = "seqid_in_sample")
all_five_results <- merge(all_five_results, five_noresults_other, by = "seqid_in_sample")
all_five_results$mod <- NULL
all_five_results$PASSED <- NULL
save(mod1_five_results_list, mod2_five_results_list, mod3_five_results_list, mod1_five_results, mod2_five_results, mod3_five_results, 
     all_five_results, final.alldata, file = paste0(dir_out, "ARIC.SCA.V", visit, ".Final_5YRMods_", appendix, ".rds"))

tmp <- all_five_results
tmp$exam <- visit
assign(paste0("aricv", visit, "_all_five_results"), tmp)

####10YR COXPH####
mod1_ten_results_list <- run.sca.coxph.list(outcome = "ten", model = "Dz", covariates = mod1_covar)
mod1_ten_results <- mod1_ten_results_list[[length(mod1_ten_results_list)]] %>% as.data.frame()

mod2_ten_results_list <- run.sca.coxph.list(outcome = "ten", model = "Min", covariates = mod2_covar)
mod2_ten_results <- mod2_ten_results_list[[length(mod2_ten_results_list)]] %>% as.data.frame()

mod3_ten_results_list <- run.sca.coxph.list(outcome = "ten", model = "RskFct", covariates = mod3_covar)
mod3_ten_results <- mod3_ten_results_list[[length(mod3_ten_results_list)]] %>% as.data.frame()


## Merge results
ten_noresults_main_head <- mod1_ten_results[, c("seqid_in_sample", "target")]
ten_noresults_main_tail <- mod1_ten_results[, c("seqid_in_sample", "flag3", "multi.apt", "n_total", "mod", "outcome", "n_cases", "n_ctrls")]
ten_noresults_other <- mod1_ten_results[, c("seqid_in_sample", "uniprot_id", "uniprot.full.name", "flag1",
                                            "flag2", "SomaCVBA_v2", "SomaCVBA_v3", "SomaCVBA_v5", "ARICCVBA_v2", "ARICCVBA_v3",
                                            "ARICCVBA_v5", "ARICpcorr_v2", "ARICpcorr_v3", "ARICpcorr_v5", "SOMAtotal.cv.plasma", "Plateflag_any_v2",
                                            "Plateflag_any_v3", "Plateflag_any_v5", "Plateflag_n_v2", "Plateflag_n_v3", "Plateflag_n_v5", "targetfullname",
                                            "apparent.kd..m.", "status", "characterization.info", "entrez.gene.id", "entrezgeneid", "entrezgenesymbol", "organism", "type")]
mod1_ten_results_base <- mod1_ten_results
names(mod1_ten_results_base)[grep("beta", names(mod1_ten_results_base))] <- "M110.Beta"
names(mod1_ten_results_base)[grep("^SE$", names(mod1_ten_results_base))] <- "M110.SE"
names(mod1_ten_results_base)[grep("^P$", names(mod1_ten_results_base))] <- "M110.P"
mod1_ten_merge <- mod1_ten_results_base[, c("seqid_in_sample", "M110.Beta", "M110.SE", "M110.P")]
mod2_ten_merge <- mrg.mkr(mod2_ten_results,2,10)
mod3_ten_merge <- mrg.mkr(mod3_ten_results,3,10)

all_ten_models <- Reduce(function(x,y) merge(x,y,by="seqid_in_sample",all=TRUE), list(mod2_ten_merge, mod1_ten_merge, mod3_ten_merge))
all_ten_results <- merge(ten_noresults_main_head, all_ten_models, by = "seqid_in_sample")
all_ten_results <- merge(all_ten_results, ten_noresults_main_tail, by = "seqid_in_sample")
all_ten_results <- merge(all_ten_results, ten_noresults_other, by = "seqid_in_sample")
all_ten_results$mod <- NULL
all_ten_results$PASSED <- NULL

save(mod1_ten_results_list, mod2_ten_results_list, mod3_ten_results_list, mod1_ten_results, mod2_ten_results, 
     mod3_ten_results, all_ten_results, final.alldata, file = paste0(dir_out, "ARIC.SCA.V", visit, ".Final_10YRMods_", appendix, ".rds"))

tmp <- all_five_results
tmp$exam <- visit
assign(paste0("aricv", visit, "_all_ten_results"), tmp)

## Save output for visit comparison
visit.files <- ls()[grep("aric", ls())]
saveRDS(get(visit.files[grep("five_results$", visit.files)]), file = paste0(dir_out, "ARIC.SCA.V", visit, ".Final_Five_Results.RDS"))
saveRDS(get(visit.files[grep("norm_results$", visit.files)]), file = paste0(dir_out, "ARIC.SCA.V", visit, ".Final_Norm_Results.RDS"))
saveRDS(get(visit.files[grep("ten_results$", visit.files)]), file = paste0(dir_out, "ARIC.SCA.V", visit, ".Final_Ten_Results.RDS"))
saveRDS(final.alldata, file = paste0(dir_out, "ARIC.SCA.V", visit, ".Final_Pheno.RDS"))

print("All models complete.")

####QQ Plots####
pdf(file = paste0(dir_out, "ARIC.SCA.V", visit, ".Final_QQ_", appendix, ".pdf"), width = 10, height = 8)
par(mfrow = c(3, 3))
make.qqplot(mod2_norm_results, "beta", "SE", "P", paste0("ARIC SCA Norm Min [n=", nrow(final.alldata), "; \ncases=", subset(final.alldata, scd_age_filt == 1) %>% nrow(), "]"))
make.qqplot(mod1_norm_results, "beta", "SE", "P", paste0("Norm Dz"))
make.qqplot(mod3_norm_results, "beta", "SE", "P", paste0("Norm RskFct"))
# make.qqplot(mod4_norm_results, "beta", "SE", "P", paste0("Norm AF"))

make.qqplot(mod2_ten_results, "beta", "SE", "P", paste0("ARIC SCA 10YR Min [n=", nrow(final.alldata), "; \ncases=", subset(final.alldata, TENYR_SCD == 1) %>% nrow(), "]"))
make.qqplot(mod1_ten_results, "beta", "SE", "P", paste0("10YR Dz"))
make.qqplot(mod3_ten_results, "beta", "SE", "P", paste0("10YR RskFct"))
# make.qqplot(mod4_ten_results, "beta", "SE", "P", paste0("10YR AF"))

make.qqplot(mod2_five_results, "beta", "SE", "P", paste0("ARIC SCA 5YR Min [n=", nrow(final.alldata), "; \ncases=", subset(final.alldata, FIVEYR_SCD == 1) %>% nrow(), "]"))
make.qqplot(mod1_five_results, "beta", "SE", "P", paste0("5YR Dz"))
make.qqplot(mod3_five_results, "beta", "SE", "P", paste0("5YR RskFct"))
# make.qqplot(mod4_five_results, "beta", "SE", "P", paste0("5YR AF"))
dev.off()
print("QQ plots made.")

####Merge outcomes####
all_models <- Reduce(function(x,y) merge(x,y,by="seqid_in_sample",all=TRUE), list(all_norm_models, all_ten_models, all_five_models))
all_results <- merge(norm_noresults_main_head, all_models, by = "seqid_in_sample")
all_results <- merge(all_results, norm_noresults_main_tail[, c(c("seqid_in_sample", "flag3", "multi.apt", "n_total"))], by = "seqid_in_sample")
all_results <- merge(all_results, norm_noresults_other, by = "seqid_in_sample")
all_results <- all_results[order(all_results$M2.P),]

all_norm_results <- all_norm_results[order(all_norm_results$M2.P),]
all_five_results <- all_five_results[order(all_five_results$M25.P),]
all_ten_results <- all_ten_results[order(all_ten_results$M210.P),]
print("All outcome results merged, moving to saving files.")

####Save Excel####
sheets <- list("Min_Norm" = mod2_norm_results, "Dz_Norm" = mod1_norm_results, "RskFct_Norm" = mod3_norm_results,
               "Min_10YR" = mod2_ten_results,  "Dz_10YR" = mod1_ten_results, "RskFct_10YR" = mod3_ten_results,
               "Min_5YR" = mod2_five_results, "Dz_5YR" = mod1_five_results, "RskFct_5YR" = mod3_five_results,
               "All_Norm" = all_norm_results, "All_10YR" = all_ten_results, "All_5YR" = all_five_results,
               "All" = all_results)
write_xlsx(sheets, paste0(dir_out, "ARIC.SCA.V", visit, ".Final_Results_", appendix, ".xlsx"))

####Save Kable####
bonf <- 0.05/(nrow(annot.final))
all_norm_results_present <- all_norm_results[, c("seqid_in_sample", "target", names(all_norm_results)[grep("^M", names(all_norm_results))], "flag3", "multi.apt", "n_total", "outcome", "n_cases", "n_ctrls")]
all_norm_results_present <- prettyPresMany(all_norm_results_present)
if(test == "Yes"){
  kableModCombo(all_norm_results_present[1:3,], paste0("ARIC.SCA.V", visit, ".ALL_NORM"), final_dz_covar_print, exclude_vars_print)
} else {
  kableModCombo(all_norm_results_present[1:150,], paste0("ARIC.SCA.V", visit, ".ALL_NORM"), final_dz_covar_print, exclude_vars_print)
}

all_five_results_present <- all_five_results[, c("seqid_in_sample", "target", names(all_five_results)[grep("^M", names(all_five_results))], "flag3", "multi.apt", "n_total", "outcome", "n_cases", "n_ctrls")]
all_five_results_present <- prettyPresMany(all_five_results_present)
if(test == "Yes"){
  kableModCombo(all_five_results_present[1:3,],  paste0("ARIC.SCA.V", visit, ".ALL_FIVE"), final_dz_covar_print, exclude_vars_print)
} else {
  kableModCombo(all_five_results_present[1:150,],  paste0("ARIC.SCA.V", visit, ".ALL_FIVE"), final_dz_covar_print, exclude_vars_print) 
}

all_ten_results_present <- all_ten_results[, c("seqid_in_sample", "target", names(all_ten_results)[grep("^M", names(all_ten_results))], "flag3", "multi.apt", "n_total", "outcome", "n_cases", "n_ctrls")]
all_ten_results_present <- prettyPresMany(all_ten_results_present)
if(test == "Yes"){
  kableModCombo(all_ten_results_present[1:3,],  paste0("ARIC.SCA.V", visit, ".ALL_TEN"), final_dz_covar_print, exclude_vars_print)
} else {
  kableModCombo(all_ten_results_present[1:150,],  paste0("ARIC.SCA.V", visit, ".ALL_TEN"), final_dz_covar_print, exclude_vars_print)
}

all_results_present <- all_results[, c("seqid_in_sample", "target", names(all_results)[grep("^M", names(all_results))], "flag3", "multi.apt", "n_total")]
all_results_present  <- prettyPresMany(all_results_present)

kableOutComboSimp(all_results_present, "M1", paste0("ARIC.SCA.V", visit, ".ALL_Dz"), final_dz_covar_print, exclude_vars_print)
kableOutComboSimp(all_results_present, "M2", paste0("ARIC.SCA.V", visit, ".ALL_Min"), final_dz_covar_print, exclude_vars_print)
kableOutComboSimp(all_results_present, "M3", paste0("ARIC.SCA.V", visit, ".ALL_RskFct"), final_dz_covar_print, exclude_vars_print)

####Volcano Plots####
pdf(file = paste0(dir_out, "ARIC.SCA.V", visit, ".Final_Vcno_", appendix, ".pdf"), width = 13, height = 4)
plot_grid(firaga(all_norm_results_present, "Norm"), firaga(all_ten_results_present, "10YR"), firaga(all_five_results_present, "5YR"), nrow = 1, ncol = 3)
dev.off()

####Save data####
rm(list = ls(pattern = "_results_list"))
save.image(paste0(dir_out, "/SCA.Assoc.ARIC.V", visit, "Prot_", appendix, ".RData"))
if(file.exists(paste0(dir_out, "/SCA.Assoc.ARIC.V", visit, "Prot_", appendix, ".RData"))){
  print("RData file saved. Hooray!")
} else {
  print("No RData file found.")
}
system("mkdir rds.files")
system("mv *.rds rds.files")
system("mv *.RDS rds.files")