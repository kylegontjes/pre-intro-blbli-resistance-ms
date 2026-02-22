library(tidyverse)
library(data.table)
library(vroom)
library(parallel)
library(pbmcapply)

setwd('/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/')
df <- readRDS("./data/dataset/df.RDS")
KPNIH1_CDS <- readRDS("./data/references/KPNIH1/KPNIH1_features.RDS")

# Load in data
setwd("/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/KPC_coverage/KPNIH1/")
no <- df$isolate_no  
filenames <-  paste0("./",no,"_summary")
kpnih1_depth_datafiles <- lapply(filenames,FUN=function(x){
    print(x)
    fread(x)
})

kpnih1_depth_datafiles <- lapply(kpnih1_depth_datafiles,as.data.frame) 

# Calculate Depth by Feature
depth_by_feature <- function(featureset,depth_data){
  features <- featureset$locus_tag
  dataset <- lapply(features,FUN=function(x,depth_data){
    locus_tag_info <- subset(featureset,locus_tag==x)
    depth_data_subset <- depth_data[locus_tag_info$start:locus_tag_info$end,"Total_Depth"] 
    #mean_depth <-  sum(depth_data_subset)/length(depth_data_subset)
    median_depth <- median(depth_data_subset)
    results <- cbind.data.frame(locus_tag=x,median_depth)
  },depth_data=depth_data)  %>% do.call(rbind, .) %>% as.data.frame() 
  colnames(dataset) <- c("locus_tag","median_depth")
  return(dataset)
} 

median_KPNIH1_gene_depth_matrix <- lapply(kpnih1_depth_datafiles,FUN=function(depth,featureset){
  dataset <- depth_by_feature(featureset,depth)
  test <- dataset %>% t() %>% as.data.frame %>%`colnames<-`(.[1,]) %>% .[2,]  %>% do.call(cbind.data.frame, .)
  return(test)
},featureset = KPNIH1_CDS) %>% data.table::rbindlist(., use.names=TRUE, fill=TRUE) %>% as.data.frame()  %>% mutate(isolate_no = no)

saveRDS(median_KPNIH1_gene_depth_matrix,"/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/KPC_coverage/KPNIH1//KPNIH1_depth_by_feature_median.RDS")


# Calculate Depth by Feature
mean_depth_by_feature <- function(featureset,depth_data){
  features <- featureset$locus_tag
  dataset <- lapply(features,FUN=function(x,depth_data){
    locus_tag_info <- subset(featureset,locus_tag==x)
    depth_data_subset <- depth_data[locus_tag_info$start:locus_tag_info$end,"Total_Depth"] 
    mean_depth <-  sum(depth_data_subset)/length(depth_data_subset)
    #median_depth <- median(depth_data_subset)
    results <- cbind.data.frame(locus_tag=x,mean_depth)
  },depth_data=depth_data)  %>% do.call(rbind, .) %>% as.data.frame() 
  colnames(dataset) <- c("locus_tag","mean_depth")
  return(dataset)
} 

mean_KPNIH1_gene_depth_matrix <- lapply(kpnih1_depth_datafiles,FUN=function(depth,featureset){
  dataset <- mean_depth_by_feature(featureset,depth)
  test <- dataset %>% t() %>% as.data.frame %>%`colnames<-`(.[1,]) %>% .[2,]  %>% do.call(cbind.data.frame, .)
  return(test)
},featureset = KPNIH1_CDS)  %>% do.call(rbind, .) %>% as.data.frame()  %>% mutate(isolate_no = no)

saveRDS(mean_KPNIH1_gene_depth_matrix,"/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/KPC_coverage/KPNIH1/KPNIH1_depth_by_feature_mean.RDS")