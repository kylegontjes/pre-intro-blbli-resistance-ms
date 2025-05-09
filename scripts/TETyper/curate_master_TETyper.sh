#!/bin/bash
cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/TETyper
 
for specimen in `ls -d */ | cut -d/ -f1`
do 
echo $specimen
cd $specimen
awk -v isolate="$specimen" 'BEGIN{ FS = OFS = "\t" } { print $0 , (NR==1? "isolate_no" : isolate)}' $specimen\_summary.txt  > $specimen\_results.txt
cd ..
done

awk '(NR==1)' ./PCMP_H3/PCMP_H3_results.txt > TETyper_results.txt

for isolate in $(ls | grep "PCMP_H" )
do
        echo $isolate
         awk '(NR>1)' ./$isolate/$isolate\_results.txt >> TETyper_results.txt
done