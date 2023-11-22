# SCD-Proteomics

## Directories

**1_Clean:** filter SomaScan files received from ARIC-central and generate output for merging to SCD phenotype data for further exclusions in preparation for associations

**2_Association:** perform aptamer-SCD associations in ARIC

**3_Meta:** run meta-analysis of ARIC and CHS results, combine all meta-analyses into a single .rds file for downstream analyses, generate files used for manuscript figure and table creation. _No stratified results were included in the manuscript (only combined meta-analysis using both races and sexes)

**4_MR:** format pQTL and SCD summary statistics then run MR

**5_Risk.Score:** protein-based risk score (PtRS) creation, using weights from ARIC-only data followed by leave-one-out cross-validation analysis in ARIC and assessing the clinical utility of the score

**6_Manuscript:** create figures and base .txt files for tables utilized in multi-panel figure
