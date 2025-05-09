# Curate ISfinder results 
cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/panISa

## 0. Create file
touch ISFinder_master.txt
# 1. Add first line of first isolate
first_isolate=$(ls *_ISFinder.txt | head -n1 | sed 's/_ISFinder.txt//')

cat $first_isolate\_ISFinder.txt | head -n1 >> ISFinder_master.txt

# 2. Add all isolates
for isolate in `ls *_ISFinder.txt | sed 's/_ISFinder.txt//'`
do 
echo $isolate
tail -n +2 $isolate\_ISFinder.txt >> ISFinder_master.txt
done
