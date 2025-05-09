#!/bin/sh
#1. Make config directory
cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/scripts/SIFT
mkdir config
cd config 
cp /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/scripts/SIFT/sift_annotation/scripts_to_build_SIFT_db/test_files/candidatus_carsonella_ruddii_pv_config.txt kpnih1_config.txt

# Edit the config directory as follows: a. Change genetic code to 11 and mito code to 0; b. Change parent directory and organism name; c. specify path to SIFT4g file and protein database 

#PARENT_DIR=/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/SIFT/kpnih1
#ORG=klebsiella_pneumoniae
#ORG_VERSION=CP008827.1

#Running SIFT 4G
#SIFT4G_PATH=/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/scripts/SIFT/sift_annotation/sift4g/bin/sift4g
#PROTEIN_DB=/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/SIFT/proteindb/uniprot_sprot.fasta

