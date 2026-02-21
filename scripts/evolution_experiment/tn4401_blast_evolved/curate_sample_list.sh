echo "sample_id" > /scratch/esnitkin_root/esnitkin1/kgontjes/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance/tn4401_blast_evolved/GeneScreener/config/sample_tn4401_evolved.tsv

for isolate in $(cat /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/scripts/evolution_experiment/evolution_experiment_strains.txt)
do 
name="${isolate}_flye_medaka_polypolish"
echo $name >>  /scratch/esnitkin_root/esnitkin1/kgontjes/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance/tn4401_blast_evolved/GeneScreener/config/sample_tn4401_evolved.tsv
done