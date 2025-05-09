# Generate KPC plasmid depth
library(tidyverse)
library(data.table)

setwd('/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/')
df <- readRDS("./data/dataset/df.RDS")

KPNIH1_CDS <- readRDS("./data/references/KPNIH1/KPNIH1_features.RDS")
KPNIH1_KPC_Plasmid_CDS <- readRDS("./data/references/KPNIH1_plasmid_pKpQIL6e6/KPNIH1_KPC_Plasmid_features.RDS")

# Load in data
setwd("/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/KPC_coverage/KPNIH1_plasmid_pKpQIL6e6")

no <- df$isolate_no  
filenames <-  paste0("./",no,"_summary")
plasmid_depth_datafiles <- lapply(filenames,fread) %>% lapply(.,as.data.frame)

calculate_overall_depth <- function(x){
  mean <-  sum(x$Total_Depth) / nrow(x)
  median <- median(x$Total_Depth)  
  non_zero <- x$Total_Depth %>% subset(.,.>0)
  mean_non_zero <-  sum(non_zero) / length(non_zero)
  median_non_zero <- median(non_zero)  
  return(cbind.data.frame(mean,median,mean_non_zero,median_non_zero))
}

calculated_plasmid_depth <- lapply(plasmid_depth_datafiles,calculate_overall_depth)  %>% data.table::rbindlist(., use.names=TRUE, fill=TRUE) %>% as.data.frame() %>% `colnames<-`(c("KPNIH1_KPC_plasmid_mean_depth","KPNIH1_KPC_plasmid_median_depth","KPNIH1_KPC_plasmid_mean_depth_no_zero","KPNIH1_KPC_plasmid__median_depth_no_zero")) %>% mutate(isolate_no = no)
saveRDS(calculated_plasmid_depth,"/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/KPC_coverage/KPNIH1_plasmid_pKpQIL6e6/KPNIH1_KPC_plasmid_depth_summary.RDS")

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


median_plasmid_gene_depth_matrix <- lapply(plasmid_depth_datafiles,FUN=function(depth,featureset){
  dataset <- depth_by_feature(featureset,depth)
  test <- dataset %>% t() %>% as.data.frame %>%`colnames<-`(.[1,]) %>% .[2,]  %>% do.call(cbind.data.frame, .)
  return(test)
},featureset = KPNIH1_KPC_Plasmid_CDS) %>% data.table::rbindlist(., use.names=TRUE, fill=TRUE) %>% as.data.frame()  %>% mutate(isolate_no = no)

saveRDS(median_plasmid_gene_depth_matrix,"/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/KPC_coverage/KPNIH1_plasmid_pKpQIL6e6/KPNIH1_KPC_plasmid_depth_by_feature_median.RDS")


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


mean_plasmid_gene_depth_matrix <- lapply(plasmid_depth_datafiles,FUN=function(depth,featureset){
  dataset <- mean_depth_by_feature(featureset,depth)
  test <- dataset %>% t() %>% as.data.frame %>%`colnames<-`(.[1,]) %>% .[2,]  %>% do.call(cbind.data.frame, .)
  return(test)
},featureset = KPNIH1_KPC_Plasmid_CDS)  %>% do.call(rbind, .) %>% as.data.frame()  %>% mutate(isolate_no = no)

saveRDS(mean_plasmid_gene_depth_matrix,"/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/KPC_coverage/KPNIH1_plasmid_pKpQIL6e6/KPNIH1_KPC_plasmid_depth_by_feature_mean.RDS")