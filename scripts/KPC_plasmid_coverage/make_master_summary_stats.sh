#!/bin/sh
cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/KPC_plasmid_coverage/results
touch master_KPC_plasmid_coverage_summary.txt
head -n1 PCMP_H1_CoverageStats_summary.txt >> master_KPC_plasmid_coverage_summary.txt
for i in $(ls PCMP_H*CoverageStats_summary.txt | sed 's/_CoverageStats_summary.txt//' | sort)
do
    echo $i
    tail -n +2 ${i}_CoverageStats_summary.txt >> master_KPC_plasmid_coverage_summary.txt
done