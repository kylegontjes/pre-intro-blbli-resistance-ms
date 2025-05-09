library(tidyverse)
library(data.table)
library(vroom)
library(parallel)
library(pbmcapply)

setwd('/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/')
df <- readRDS("./data/dataset/df.RDS")

# Load in data
setwd("/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/KPC_coverage/KPNIH1/")

no <- df$isolate_no  
filenames <-  paste0("./",no,"_summary")
kpnih1_depth_datafiles <- lapply(filenames,vroom) %>% lapply(.,as.data.frame)

calculate_overall_depth <- function(x){
  mean <-  sum(x$Total_Depth) / nrow(x)
  median <- median(x$Total_Depth)  
  non_zero <- x$Total_Depth %>% subset(.,.>0)
  mean_non_zero <-  sum(non_zero) / length(non_zero)
  median_non_zero <- median(non_zero)  
  return(cbind.data.frame(mean,median,mean_non_zero,median_non_zero))
}

calculated_kpnih1_depth <- pbmclapply(kpnih1_depth_datafiles,calculate_overall_depth,mc.cores=10) %>% data.table::rbindlist(., use.names=TRUE, fill=TRUE) %>% as.data.frame() %>% `colnames<-`(c("KPNIH1_mean_depth","KPNIH1_median_depth","KPNIH1_mean_depth_non_zero","KPNIH1_median_depth_non_zero")) %>% mutate(isolate_no = no)  
saveRDS(calculated_kpnih1_depth,"/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/KPC_coverage/KPNIH1/KPNIH1_depth_summary.RDS")
