##To Run: export exclude=eGFR_ckdepi density=No egfrcrit=30 assoc=No postStrat=Yes; qsub -v exclude -v density -v egfrcrit -v assoc -v postStrat -l mem_free=15G,h_vmem=15G ./runSCACleanMulti_ver4.sh
#exclude format: eGFR_ckdepi+chd+hf+stroke

#!/bin/bash
##$ -N
#$ -cwd
##$ -R y
##$ -pe local 5
##$ -l mem_free=4G,h_vmem=4G
##$ -e aric.prot.real.archive/scd.assoc/2_analyses/combined
##$ -o aric.prot.real.archive/scd.assoc/2_analyses/combined
#$ -m e
#$ -M 

Rscript aric.prot.real.archive/scd.assoc/2_analyses/combined/exclude.egfr/with.afib/SCA.Clean.ARIC.ANML.Prot.Multi_ver4.R $exclude $density $egfrcrit $assoc $postStrat

