#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -R y
#$ -l mem_free=3G,h_vmem=3G
##$ -m e
##$ -M 

####Arguments####
## pQTL filters
pqtl_p=1.004e-11
pqtl_maf=0.01
pqtl_hetdf=2
pqtl_hetpval=0.05
pqtl_hetIsq=75

## Outcome filters
sca_hetdf=2
sca_hetpval=0.05
sca_hetIsq=75

## Clumping filters
rsq=0.1
kb=1000
##1KG##
#panel=arkinglab/resources/1000G/mrcieu.plink/EUR
##ARIC##
panel=arkinglab/active/projects/scd.meta/analyses/scd.meta.ver2/ARIC.b35.b37.liftover/aric.f3v2.imputed.b37
##ARIC TOPMed##
#panel=ARIC_Data/GWAS/TOPMed/EA/plink/ARICEA_TOPMedimputed_maf005.imp30.rsid

## File naming
appendix=ARIC.p1e11.maf001.rsq01.kb1e3

####Run script####
Rscript aric.prot.real.archive/scd.assoc/7_MR/scripts/SCAxProt_MR_Pietzner2021.R $pqtl_p $pqtl_maf $pqtl_hetdf $pqtl_hetpval $pqtl_hetIsq $sca_hetdf $sca_hetpval $sca_hetIsq $rsq $kb $panel $appendix
