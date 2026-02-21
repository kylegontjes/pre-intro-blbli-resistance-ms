#!/bin/bash
cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/evolution_experiment/coverage_statistics/database/
module load Bioinformatics samtools
samtools faidx PCMP_H261.fasta
samtools faidx PCMP_H318.fasta
samtools faidx PCMP_H116.fasta
samtools faidx PCMP_H19.fasta
samtools faidx PCMP_H226.fasta
samtools faidx PCMP_H323.fasta
samtools faidx PCMP_H414.fasta
samtools faidx PCMP_H464.fasta
module purge
module load Bioinformatics bwa
bwa index PCMP_H261.fasta
bwa index PCMP_H318.fasta
bwa index PCMP_H116.fasta
bwa index PCMP_H19.fasta
bwa index PCMP_H226.fasta
bwa index PCMP_H323.fasta
bwa index PCMP_H414.fasta
bwa index PCMP_H464.fasta
module purge
module load Bioinformatics gatk
gatk CreateSequenceDictionary -R PCMP_H261.fasta
gatk CreateSequenceDictionary -R PCMP_H318.fasta
gatk CreateSequenceDictionary -R PCMP_H116.fasta
gatk CreateSequenceDictionary -R PCMP_H19.fasta
gatk CreateSequenceDictionary -R PCMP_H226.fasta
gatk CreateSequenceDictionary -R PCMP_H323.fasta
gatk CreateSequenceDictionary -R PCMP_H414.fasta
gatk CreateSequenceDictionary -R PCMP_H464.fasta
module purge
