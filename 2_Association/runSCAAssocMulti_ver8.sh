##export appendix=BthRcSx test=No; qsub -v appendix -v test ./runSCAAssocMulti_ver8.sh
#!/bin/bash
##$ -N Mod1
#$ -cwd
#$ -l mem_free=15G,h_vmem=15G
#$ -m e
#$ -M 

echo $appendix
input_data=$(realpath ../Final.rds)
Rscript aric.prot.real.archive/scd.assoc/3b_analyses_ANML_ver2/scripts/SCA.Assoc.ARIC.ANML.Prot.Multi_ver8.R $appendix $input_data $test

