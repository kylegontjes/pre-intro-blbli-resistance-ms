#!/bin/sh
cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/references/KPNIH1

# BWA INDEX
module load Bioinformatics
module load bwa
bwa index KPNIH1.fasta

# Samtools faidx
module purge 
module load Bioinformatics
module load samtools
samtools faidx KPNIH1.fasta

# GATK Dictionary
module purge
module load Bioinformatics
module load gatk
gatk CreateSequenceDictionary -R  KPNIH1.fasta