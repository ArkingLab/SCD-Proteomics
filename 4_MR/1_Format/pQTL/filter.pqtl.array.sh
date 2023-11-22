#qsub -t 1:102 filter.pqtl.array.sh

#!/bin/bash
##$ -N Mod1
#$ -cwd
#$ -l mem_free=1.5G,h_vmem=1.5G
#$ -o aric.prot.real.archive/7_MR/pQTL.databases/Pietzner.2021/p_05/outfiles
#$ -e aric.prot.real.archive/7_MR/pQTL.databases/Pietzner.2021/p_05/outfiles
##$ -m e
##$ -M 

i=$SGE_TASK_ID
f=$(awk -v a="$i" 'NR==a' files2filt.txt)
b=$(basename "$f" .txt.gz)

#FILES="../*.txt.gz"
#for f in $FILES
#do
echo "Processing $f file..."
c="${b}_temp.txt"
zcat "$f" | awk '{ if ($11 < 0.05) { print } }' > $c
#echo $b
b2="${b}_05.txt"
#echo "$b2"
cat header.txt $c > $b2
rm $c
echo "Done $b2 file..."
#done
