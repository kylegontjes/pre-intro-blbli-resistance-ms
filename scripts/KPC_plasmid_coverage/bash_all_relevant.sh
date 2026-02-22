#!/bin/sh
cd /scratch/esnitkin_root/esnitkin1/kgontjes/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance/KPC_plasmid_coverage/CoverageStats

for i in $(ls *CoverageStats_*.sbat | grep -v 'model' | grep -v 'pretrimmed')
do
    echo $i
    sbatch $i
done
