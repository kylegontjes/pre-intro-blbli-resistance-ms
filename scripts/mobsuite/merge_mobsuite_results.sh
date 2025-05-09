#!/bin/sh
cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/mobsuite

awk '(NR==1)' ./PCMP_H1/mobtyper_results.txt > CRKP_mobtyper_results.txt

for isolate in $(ls -d */ | cut -d/ -f1)
do 
        echo $isolate
        awk '(NR>1)' ./$isolate/mobtyper_results.txt  >> CRKP_mobtyper_results.txt

done

awk '(NR==1)' ./PCMP_H1/contig_report.txt > CRKP_contig_report.txt

for isolate in $(ls -d */ | cut -d/ -f1)
do
        echo $isolate
        awk '(NR>1)' ./$isolate/contig_report.txt >> CRKP_contig_report.txt
done