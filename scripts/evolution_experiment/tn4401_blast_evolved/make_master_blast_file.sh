#!/bin/sh
cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/evolution_experiment/tn4401_blast_evolved
# Curate results  
first_isolate=$(ls *_clean.tsv | head -n1)
## 0. Create file
touch tn4401_blast_master_evolved.txt
# 1. Add first line of first isolate
cat $first_isolate | head -n1 >> tn4401_blast_master_evolved.txt
for isolate in `ls *_clean.tsv`
do 
echo $isolate
tail -n +2 $isolate  >> tn4401_blast_master_evolved.txt
done
