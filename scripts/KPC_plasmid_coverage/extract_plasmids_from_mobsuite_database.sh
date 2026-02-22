#!/bin/bash
cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/KPC_plasmid_coverage/database/data
module load Bioinformatics samtools
samtools faidx ncbi_plasmid_full_seqs.fas
samtools faidx ncbi_plasmid_full_seqs.fas CP006919 -o CP006919.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP007729 -o CP007729.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP008833 -o CP008833.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP009773 -o CP009773.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP009776 -o CP009776.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP010362 -o CP010362.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP011977 -o CP011977.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP011978 -o CP011978.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP018430 -o CP018430.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP020075 -o CP020075.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP021548 -o CP021548.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP021752 -o CP021752.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP021778 -o CP021778.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP021853 -o CP021853.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP021900 -o CP021900.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP025009 -o CP025009.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP025010 -o CP025010.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP027152 -o CP027152.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP028181 -o CP028181.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP029136 -o CP029136.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP032209 -o CP032209.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP036441 -o CP036441.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP036449 -o CP036449.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP039971 -o CP039971.fasta
samtools faidx ncbi_plasmid_full_seqs.fas CP039975 -o CP039975.fasta
samtools faidx ncbi_plasmid_full_seqs.fas FJ223607 -o FJ223607.fasta
samtools faidx ncbi_plasmid_full_seqs.fas HG969997 -o HG969997.fasta
samtools faidx ncbi_plasmid_full_seqs.fas HQ589350 -o HQ589350.fasta
samtools faidx ncbi_plasmid_full_seqs.fas JX430448 -o JX430448.fasta
samtools faidx ncbi_plasmid_full_seqs.fas JX442974 -o JX442974.fasta
samtools faidx ncbi_plasmid_full_seqs.fas KF954759 -o KF954759.fasta
samtools faidx ncbi_plasmid_full_seqs.fas KJ721789 -o KJ721789.fasta
samtools faidx ncbi_plasmid_full_seqs.fas KU665642 -o KU665642.fasta
samtools faidx ncbi_plasmid_full_seqs.fas KU934011 -o KU934011.fasta
samtools faidx ncbi_plasmid_full_seqs.fas KY930324 -o KY930324.fasta
samtools faidx ncbi_plasmid_full_seqs.fas KY940546 -o KY940546.fasta
samtools faidx ncbi_plasmid_full_seqs.fas NC_019151 -o NC_019151.fasta
samtools faidx CP006919.fasta
samtools faidx CP007729.fasta
samtools faidx CP008833.fasta
samtools faidx CP009773.fasta
samtools faidx CP009776.fasta
samtools faidx CP010362.fasta
samtools faidx CP011977.fasta
samtools faidx CP011978.fasta
samtools faidx CP018430.fasta
samtools faidx CP020075.fasta
samtools faidx CP021548.fasta
samtools faidx CP021752.fasta
samtools faidx CP021778.fasta
samtools faidx CP021853.fasta
samtools faidx CP021900.fasta
samtools faidx CP025009.fasta
samtools faidx CP025010.fasta
samtools faidx CP027152.fasta
samtools faidx CP028181.fasta
samtools faidx CP029136.fasta
samtools faidx CP032209.fasta
samtools faidx CP036441.fasta
samtools faidx CP036449.fasta
samtools faidx CP039971.fasta
samtools faidx CP039975.fasta
samtools faidx FJ223607.fasta
samtools faidx HG969997.fasta
samtools faidx HQ589350.fasta
samtools faidx JX430448.fasta
samtools faidx JX442974.fasta
samtools faidx KF954759.fasta
samtools faidx KJ721789.fasta
samtools faidx KU665642.fasta
samtools faidx KU934011.fasta
samtools faidx KY930324.fasta
samtools faidx KY940546.fasta
samtools faidx NC_019151.fasta
module purge
module load Bioinformatics bwa
bwa index CP006919.fasta
bwa index CP007729.fasta
bwa index CP008833.fasta
bwa index CP009773.fasta
bwa index CP009776.fasta
bwa index CP010362.fasta
bwa index CP011977.fasta
bwa index CP011978.fasta
bwa index CP018430.fasta
bwa index CP020075.fasta
bwa index CP021548.fasta
bwa index CP021752.fasta
bwa index CP021778.fasta
bwa index CP021853.fasta
bwa index CP021900.fasta
bwa index CP025009.fasta
bwa index CP025010.fasta
bwa index CP027152.fasta
bwa index CP028181.fasta
bwa index CP029136.fasta
bwa index CP032209.fasta
bwa index CP036441.fasta
bwa index CP036449.fasta
bwa index CP039971.fasta
bwa index CP039975.fasta
bwa index FJ223607.fasta
bwa index HG969997.fasta
bwa index HQ589350.fasta
bwa index JX430448.fasta
bwa index JX442974.fasta
bwa index KF954759.fasta
bwa index KJ721789.fasta
bwa index KU665642.fasta
bwa index KU934011.fasta
bwa index KY930324.fasta
bwa index KY940546.fasta
bwa index NC_019151.fasta
module purge
module load Bioinformatics gatk
gatk CreateSequenceDictionary -R CP006919.fasta
gatk CreateSequenceDictionary -R CP007729.fasta
gatk CreateSequenceDictionary -R CP008833.fasta
gatk CreateSequenceDictionary -R CP009773.fasta
gatk CreateSequenceDictionary -R CP009776.fasta
gatk CreateSequenceDictionary -R CP010362.fasta
gatk CreateSequenceDictionary -R CP011977.fasta
gatk CreateSequenceDictionary -R CP011978.fasta
gatk CreateSequenceDictionary -R CP018430.fasta
gatk CreateSequenceDictionary -R CP020075.fasta
gatk CreateSequenceDictionary -R CP021548.fasta
gatk CreateSequenceDictionary -R CP021752.fasta
gatk CreateSequenceDictionary -R CP021778.fasta
gatk CreateSequenceDictionary -R CP021853.fasta
gatk CreateSequenceDictionary -R CP021900.fasta
gatk CreateSequenceDictionary -R CP025009.fasta
gatk CreateSequenceDictionary -R CP025010.fasta
gatk CreateSequenceDictionary -R CP027152.fasta
gatk CreateSequenceDictionary -R CP028181.fasta
gatk CreateSequenceDictionary -R CP029136.fasta
gatk CreateSequenceDictionary -R CP032209.fasta
gatk CreateSequenceDictionary -R CP036441.fasta
gatk CreateSequenceDictionary -R CP036449.fasta
gatk CreateSequenceDictionary -R CP039971.fasta
gatk CreateSequenceDictionary -R CP039975.fasta
gatk CreateSequenceDictionary -R FJ223607.fasta
gatk CreateSequenceDictionary -R HG969997.fasta
gatk CreateSequenceDictionary -R HQ589350.fasta
gatk CreateSequenceDictionary -R JX430448.fasta
gatk CreateSequenceDictionary -R JX442974.fasta
gatk CreateSequenceDictionary -R KF954759.fasta
gatk CreateSequenceDictionary -R KJ721789.fasta
gatk CreateSequenceDictionary -R KU665642.fasta
gatk CreateSequenceDictionary -R KU934011.fasta
gatk CreateSequenceDictionary -R KY930324.fasta
gatk CreateSequenceDictionary -R KY940546.fasta
gatk CreateSequenceDictionary -R NC_019151.fasta
module purge
