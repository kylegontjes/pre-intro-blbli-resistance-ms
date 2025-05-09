
# Hogwash Run
# Start Date: 04/21/25
# Activate
library(hogwash)
library(tidyverse)
library(ape)
Sys.info()
sessionInfo()

setwd('/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/')

set.seed(45)

#Load Dataset & Tree
tr <- read.tree(file=paste0('/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/tree/tree.treefile'))
df <-  readRDS(file=paste0('/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/dataset/df.RDS')) %>% .[match(tr$tip.label,rownames(.)),] 
rownames(df) <- df$isolate_no

#Phenotype
phenotype <- df %>% select(MVB_log_2) %>% as.matrix

#Genotypes
geno <- readRDS(file=paste0('./data/GWAS/matrices/','model3.2_burden','.RDS')) %>% subset(rownames(.) %in% rownames(df)) %>% .[match(tr$tip.label,rownames(.)),] %>% as.data.frame %>% select_if(.,colSums(.) > 1 & colSums(.) < nrow(.) - 2) %>% as.matrix

print('running hogwash for model: analysis')
hogwash(
  pheno = phenotype,
  geno = geno,
  tree = tr,
  fdr = 0.10,
  file_name = paste0('MVB_log_2','_','model3.2_burden'),
  dir = '/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/GWAS/results') 
