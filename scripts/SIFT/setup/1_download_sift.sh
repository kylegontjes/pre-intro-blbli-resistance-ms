#!/bin/sh
# Name: Download Sift

# Date: 02-24-25
# Notes: Here is the code that I used to generate the SIFT4g 
# 1. SIFT4g Dabase Generation Github Page: https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB
# 2. Annotator Github Page: https://github.com/pauline-ng/SIFT4G_Annotator

cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/scripts/SIFT
mkdir sift_annotation
cd sift_annotation

#1.  Download the SIFT Annotator jar file and the instructions to your directory
wget https://github.com/pauline-ng/SIFT4G_Annotator/raw/master/SIFT4G_Annotator.jar
wget https://sift.bii.a-star.edu.sg/sift4g/AnnotateVariants.html 

#2.  Download Sift4g algorithm from Github
module load Bioinformatics
module load gcc/11.2.0
module load bioperl
module list

git clone --recursive https://github.com/rvaser/sift4g.git sift4g
cd sift4g/
make

cd ..

#3. Get Scripts to Build SIFT databse for Reference Genome of Choice
git clone https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB.git scripts_to_build_SIFT_db