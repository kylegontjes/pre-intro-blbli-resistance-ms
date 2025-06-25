numerical_data <- function(variable,df){
  isolates_w_variant <- df %>% subset(get(variable) ==1)
  representation_in_clusters <-  table(isolates_w_variant %>% .$blbli_asr_cluster_renamed_simple)
  true_cluster_size <- table(df %>% .$blbli_asr_cluster_renamed_simple)
  
  res_prop_df <-lapply(names(true_cluster_size),FUN=function(x){
    if(x %in% names(representation_in_clusters) == F){
      prop =0
    }
    if(x %in% names(representation_in_clusters) == T){
      prop = representation_in_clusters[[x]] / true_cluster_size[[x]]
    }
    return(prop)
  }) %>% do.call(cbind.data.frame,.) %>% `colnames<-`(paste0(names(true_cluster_size),"_prop"))
  
  res_ct_df <-lapply(names(true_cluster_size),FUN=function(x){
    if(x %in% names(representation_in_clusters) == F){
      ct =0
    }
    if(x %in% names(representation_in_clusters) == T){
      ct = representation_in_clusters[[x]]  
    }
    return(ct)
  }) %>% do.call(cbind.data.frame,.) %>% `colnames<-`(paste0(names(true_cluster_size),"_ct"))
  
  res_prop_df$Cluster_prop <- nrow(isolates_w_variant %>% subset(!blbli_asr_cluster_renamed_simple %in% c("Singleton","Susceptible"))) / nrow(df %>% subset(!blbli_asr_cluster_renamed_simple %in% c("Singleton","Susceptible")))   
  res_ct_df$Cluster_ct <- nrow(isolates_w_variant %>% subset(!blbli_asr_cluster_renamed_simple %in% c("Singleton","Susceptible")))  
  results <- cbind(res_prop_df,res_ct_df)
  return(results)
}