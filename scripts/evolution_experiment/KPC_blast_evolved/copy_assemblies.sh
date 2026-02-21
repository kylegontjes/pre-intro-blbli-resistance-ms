for isolate in $(cat /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/scripts/evolution_experiment/evolution_experiment_strains.txt)
do 
echo $isolate
cp /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Sequence_data/assembly/hybrid/polypolish/$isolate/$isolate/$isolate\_flye_medaka_polypolish.fasta /scratch/esnitkin_root/esnitkin1/kgontjes/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance/blast_KPC_evolved/GeneScreener/assemblies
done