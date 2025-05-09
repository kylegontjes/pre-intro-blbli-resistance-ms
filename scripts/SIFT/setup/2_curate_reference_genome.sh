#!/bin/sh
# Part 2: Get Reference Genome from NCBI
#1. Get Reference Genome Data from NCBI

cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/SIFT
mkdir genomic_data
cd genomic_data

#1a. Get fasta file of reference genome
## Note: I found that NCBI's GenBank assembly worked. The RefSeq assembly did not run properly
mkdir fasta
cd fasta 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/281/535/GCF_000281535.2_ASM28153v2/GCF_000281535.2_ASM28153v2_genomic.fna.gz
gunzip GCF_000281535.2_ASM28153v2_genomic.fna.gz
sed -i "1s/.*/>CP008827.1/" GCF_000281535.2_ASM28153v2_genomic.fna
cp GCF_000281535.2_ASM28153v2_genomic.fna kpnih1.fa 
gzip kpnih1.fa

#1b. Download GTF file from NCBI
cd ..
mkdir gtf
cd gtf 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/281/535/GCF_000281535.2_ASM28153v2/GCF_000281535.2_ASM28153v2_genomic.gtf.gz
gunzip GCF_000281535.2_ASM28153v2_genomic.gtf.gz

#2.For SIFT4g database to be generated, the gtf file must include an exon and transcript entry for each locus_tag. I had to manually update the gtf file using the script below. 
#2a. Update kpnih1_gtf.R script
module load R/4.2.0
module load Rtidyverse/4.2.0
module list
####  Commands your job should run follow this line
echo "Running R"

R CMD BATCH --no-restore --no-save --quiet /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/scripts/SIFT/setup/r_script_add_exon-transcript.R

#2b. gzip
gzip kpnih1.gtf

#3.  Make parent directory & add fasta and gtf files per github page instructions 
## Note: Ensure that the gtf and fasta files are gunzip-ed 

cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/SIFT
mkdir kpnih1
cd kpnih1
mkdir gene-annotation-src
cp ../genomic_data/gtf/kpnih1.gtf.gz ./gene-annotation-src
mkdir chr-src
cp ../genomic_data/fasta/kpnih1.fa.gz ./chr-src
mkdir dbSNP