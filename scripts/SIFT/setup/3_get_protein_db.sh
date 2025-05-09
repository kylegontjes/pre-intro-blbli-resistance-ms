#!/bin/sh
#1. Add protein db
## Note: I used the suggested default protein database: uniprot's sprot database (https://www.uniprot.org/)

cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/SIFT
mkdir proteindb
cd proteindb
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz 
