#!/bin/sh
cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/mobsuite_ONT

awk '(NR==1)' ./PCMP_H148/mobtyper_results.txt > mobtyper_results_ONT.txt

for isolate in $(ls -d */ | cut -d/ -f1)
do 
        echo $isolate
        awk '(NR>1)' ./$isolate/mobtyper_results.txt  >> mobtyper_results_ONT.txt

done

awk '(NR==1)' ./PCMP_H148/contig_report.txt > mobsuite_contig_report_ONT.txt

for isolate in $(ls -d */ | cut -d/ -f1)
do
        echo $isolate
        awk '(NR>1)' ./$isolate/contig_report.txt >> mobsuite_contig_report_ONT.txt
done