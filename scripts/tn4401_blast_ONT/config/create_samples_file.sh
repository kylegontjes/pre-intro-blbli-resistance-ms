#!/bin/sh
# Script to create samples file
path="/scratch/esnitkin_root/esnitkin1/kgontjes/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance/tn4401_blast_ONT/GeneScreener/assemblies/"
sample_id="sample_id"
sample_names=$(ls -1 $path | grep fasta | cut -d. -f1 | sort | uniq)
echo -e\n $sample_id $sample_names | tr ' ' '\n' > /scratch/esnitkin_root/esnitkin1/kgontjes/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance/tn4401_blast_ONT/GeneScreener/config/sample.tsv