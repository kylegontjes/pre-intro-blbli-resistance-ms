#!/bin/sh

# Create directories
cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/
## Download CARD Dataset
mkdir CARD
cd CARD
wget https://card.mcmaster.ca/download/0/broadstreet-v4.0.0.tar.bz2
tar -xvf broadstreet-v4.0.0.tar.bz2

# Create KPC Specific dataset 
mkdir KPC_database
cd KPC_database
## Accessions
grep 'KPC beta-lactamase' ../aro_index.tsv | cut -f1 | cut -d: -f2 > kpc_aro_accessions
## Fasta files
touch kpc_nuc_db.fasta
for accession in `cat kpc_aro_accessions`;
do
  sed -n "/$accession/,/>/{p;}" ../nucleotide_fasta_protein_homolog_model.fasta | sed '$d' >> kpc_nuc_db.fasta
done