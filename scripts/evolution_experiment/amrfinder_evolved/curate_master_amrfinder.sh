#!/bin/sh
cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/evolution_experiment/amrfinder_evolved

# Penn KPC 
awk '(NR==1)' PCMP_H19_amrfinder.out > evolved_amrfinder_results.txt

for isolate in $(ls PCMP_H* )
do 
        echo $isolate
        awk '(NR>1)' $isolate  >> evolved_amrfinder_results.txt
done