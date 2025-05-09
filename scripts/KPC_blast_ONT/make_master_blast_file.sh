#!/bin/sh
cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/KPC_blast_ONT
# Curate results  
first_isolate=$(ls *_clean.tsv | head -n1)
## 0. Create file
touch KPC_blast_master_ONT.txt
# 1. Add first line of first isolate
cat $first_isolate | head -n1 >> KPC_blast_master_ONT.txt
for isolate in `ls *_clean.tsv`
do 
echo $isolate
tail -n +2 $isolate  >> KPC_blast_master_ONT.txt
done
