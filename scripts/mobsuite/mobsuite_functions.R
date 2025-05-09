# Mobsuite functions
read_mobtyper_results <- function(path,df){
  report <- read.delim(file=path) %>% mutate(isolate_no = gsub("_contigs_l1000","",word(sample_id,1,1,sep=":")),sequencing = "short") %>% subset(.,isolate_no %in% df$isolate_no) %>% `colnames<-`(gsub("\\.","",colnames(.)))
  return(report)
}

get_mobtyper_cluster_matrix <- function (df, mobtyper){
  plasmid_cluster <- mobtyper$primary_cluster_id %>% unique
  results <- lapply(plasmid_cluster, FUN = function(x, df, mobtyper) {
    isolates_with_plasmid <- subset(mobtyper, primary_cluster_id == 
                                      x) %>% select(isolate_no) %>% unlist
    row <- assign(paste0(x), ifelse(df$isolate_no %in% isolates_with_plasmid, 
                                    1, 0))
    return(row)
  }, df, mobtyper) %>% do.call(bind_cols, .) %>% as.data.frame %>% 
    `colnames<-`(plasmid_cluster) %>% `rownames<-`(df$isolate_no)
  return(results)
}

get_mobtyper_mash_nn_matrix <- function(df, mobtyper){
  nn_plasmid <- mobtyper$mash_nearest_neighbor %>% unique
  results <- lapply(nn_plasmid, FUN = function(x, df, mobtyper) {
    isolates_with_plasmid <- subset(mobtyper, mash_nearest_neighbor == 
                                      x) %>% select(isolate_no) %>% unlist
    row <- assign(paste0(x), ifelse(df$isolate_no %in% isolates_with_plasmid, 
                                    1, 0))
    return(row)
  }, df, mobtyper) %>% do.call(bind_cols, .) %>% as.data.frame %>% 
    `colnames<-`(nn_plasmid) %>% `rownames<-`(df$isolate_no)
  return(results)
}

recode_matrix_by_kpc_content <- function(isolate,type,KPC_matrix,KPC_location_data){
  kpc_containing_plasmid <- KPC_location_data %>% subset(isolate_no == isolate) %>% select(paste0(type)) %>% unlist
  isolate_KPC_matrix <- KPC_matrix %>% subset(rownames(.) %in% isolate) 
  if(length(kpc_containing_plasmid)>0){
    for(i in 1:ncol(isolate_KPC_matrix)){ 
      isolate_KPC_matrix[,i] <- ifelse(isolate_KPC_matrix[,i] == 1 & colnames(isolate_KPC_matrix[i]) == kpc_containing_plasmid, "KPC Plasmid",      ifelse(isolate_KPC_matrix[,i]==1,"Non-KPC plasmid","Not Present"))
    }
  }
  
  if(length(kpc_containing_plasmid)==0){
    for(i in 1:ncol(isolate_KPC_matrix)){ 
      isolate_KPC_matrix[,i] <- ifelse(isolate_KPC_matrix[,i]==1,"Non-KPC plasmid","Not Present")
    }
  }
  
  return(isolate_KPC_matrix)
} 
