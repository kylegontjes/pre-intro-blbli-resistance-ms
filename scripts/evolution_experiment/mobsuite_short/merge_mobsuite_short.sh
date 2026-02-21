#!/bin/sh
cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/evolution_experiment/mobsuite_short

awk '(NR==1)' ./PCMP_H116/mobtyper_results.txt > mobtyper_results_evolved_short.txt

for isolate in $(ls -d */ | cut -d/ -f1)
do 
        echo $isolate
        awk '(NR>1)' ./$isolate/mobtyper_results.txt  >> mobtyper_results_evolved_short.txt

done

awk '(NR==1)' ./PCMP_H116/contig_report.txt > mobsuite_contig_report_evolved_short.txt

for isolate in $(ls -d */ | cut -d/ -f1)
do
        echo $isolate
        awk '(NR>1)' ./$isolate/contig_report.txt >> mobsuite_contig_report_evolved_short.txt
done