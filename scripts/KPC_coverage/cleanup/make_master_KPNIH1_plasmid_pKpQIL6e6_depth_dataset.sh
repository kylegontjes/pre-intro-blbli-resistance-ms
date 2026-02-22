#!/bin/sh
cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/KPC_coverage/KPNIH1_plasmid_pKpQIL6e6
cat PCMP_H200_summary  | cut -f1 -d, > master_KPNIH1_plasmid_pKpQIL6e6_depth_matrix.csv
isolates=$(ls PCMP_H*_summary | grep -v ".sample" | sed 's/_summary//g' | sort | uniq)

for i in $isolates
do
    echo $i
    isolate_name=$(echo $i  )
    isolate_entry=$(cat ${isolate_name}_summary | cut -f2 -d,| tail -n +2)
    printf "%s\n" $isolate_name $isolate_entry > temp.csv 
    paste -d',' <(sed '/^$/d' master_KPNIH1_plasmid_pKpQIL6e6_depth_matrix.csv) <(sed '/^$/d' temp.csv) > merged.csv
    rm temp.csv
    mv merged.csv master_KPNIH1_plasmid_pKpQIL6e6_depth_matrix.csv
done