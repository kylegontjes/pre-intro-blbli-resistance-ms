# GWAS analysis scripts
# Hogwash output analysis script functions

# Setup
## 1. Load matrix
left_join_matrix <- function(matrix,df,number){
  df_fin <- df %>% left_join(.,matrix %>% `colnames<-`(paste0(colnames(.),number)) %>% mutate(isolate_no = rownames(.)))
  return(df_fin)
}

## 2. Load hogwash files
load_file <- function(name){
  load(file = paste0("./",name))%>% get(.)
}

# General data analyses
## 1. Create summary dataset
hogwash_sum_data <- function(object,object_list,test){
  name <- object
  object <- object_list[[paste0(name)]]
  tot_variants <- nrow(object$hit_pvals) + length(object$dropped_genotypes) + length(object$no_convergence_genotypes)
  num_dropped_genotypes <- length(object$dropped_genotypes)
  num_no_convergent_genotypes <- length(object$no_convergence_genotypes)
  num_convergent_genotypes <-  nrow(object$hit_pvals)

  if(test=="binary"){
    num_sig_hits <- nrow(object$sig_pvals)
    sig_hits <- ifelse(nrow(object$sig_pvals) >0 ,rownames(object$sig_pvals) %>% paste0(.,collapse="",sep=",") %>% gsub(",$", "", .),NA)
  }
  if(test=="continuous"){
    num_sig_hits <- nrow(object$sig_hits)
    sig_hits <- ifelse(nrow(object$sig_hits) >0 ,rownames(object$sig_hits) %>% paste0(.,collapse="",sep=",") %>% gsub(",$", "", .),NA)
  }
  phylogenetic_signal <- round(object$phylogenetic_signal,3)
  data_row <- cbind(name,tot_variants,num_dropped_genotypes,num_no_convergent_genotypes,num_convergent_genotypes,num_sig_hits,sig_hits,phylogenetic_signal)
  return(data_row)
}

## 2. Get metadata data from file name
get_metadata_from_filename <- function(x){
  name <- x
  test <- str_split(x,"_")  %>% lapply(.,FUN =function(x){x[1]}) %>% unlist
  phenotype <-  str_split(x,"_")  %>% lapply(.,FUN =function(x){x[2]}) %>% unlist
  model <-str_split(x,"_") %>% lapply(.,FUN =function(x){paste0(x[6:length(x)-1],collapse="_")})  %>% unlist
  results <- cbind(name,phenotype,test,model) %>% as.data.frame
  return(results)
}

## 3. Update core genome hits with eggnog data
get_core_gene_hit_annotations <- function(name,eggnog_gff_df){
  locus_tag = name
  if(grepl("-",locus_tag) ==F ){
    eggnog_info <- eggnog_gff_df %>% subset(locus_tag == name) %>% select(locus_tag,gene,start,end,strand,product)
  }

  if(grepl("-",locus_tag) ==T ){
    locus_tags <- str_split(pattern = "-",name,simplify=T)
    locus_tag_1 <- locus_tags[,1]
    locus_tag_2 <- locus_tags[,2]

    eggnog_info_1 <- eggnog_gff_df %>% subset(locus_tag == locus_tag_1) %>% select(gene,start,end,strand,product)
    eggnog_info_2 <- eggnog_gff_df %>% subset(locus_tag == locus_tag_2) %>% select(gene,start,end,strand,product)
    eggnog_data <- rbind(eggnog_info_1,eggnog_info_2) %>% apply(.,MARGIN = 2,FUN=paste0,collapse=";") %>% t %>% as.data.frame()
    eggnog_info <- cbind(locus_tag,eggnog_data)
  }
  return(eggnog_info)
}

## 4. Get phenotype data
get_pheno <- function(locus_tag,sum_dataset){
  locus_tag_df <- sum_dataset[grepl(locus_tag,sum_dataset$sig_hits,fixed = T),]
  test_name <- apply(locus_tag_df,1,FUN=function(x){paste(x["phenotype"],x["test"],x["model"],sep = "_")}) %>% paste0(.,collapse=";")
  data <- cbind(locus_tag,test_name) %>% as.data.frame
  return(data)
}

## 5. Get phenotype summary data
get_pheno_sum <- function(locus_tag,sum_dataset){
  locus_tag_df <- sum_dataset[grepl(locus_tag,sum_dataset[,"sig_hits"],fixed = TRUE),]
  ## blbli
  blbli_all <- subset(locus_tag_df,grepl("blbli",phenotype)) %>% as.data.frame
  blbli_all <- ifelse(nrow(blbli_all)>0,apply(blbli_all,1,FUN=function(x){paste(x["phenotype"],x["test"],x["model"],x["population"],sep = "_")}) %>% paste0(.,collapse=";"),"")
  blbli_sync <- str_split(blbli_all,";") %>% unlist %>% .[grepl("synchronous",.)] %>% gsub("synchronous_","",.) %>% paste0(.,collapse=";")
  blbli_phyc <- str_split(blbli_all,";") %>% unlist %>% .[grepl("phyc",.)] %>% gsub("phyc_","",.)  %>% paste0(.,collapse=";")
  blbli_cont <- str_split(blbli_all,";") %>% unlist %>% .[grepl("continuous",.)] %>% gsub("continuous_","",.)  %>% paste0(.,collapse=";")
  ##MVB
  MVB_all <- subset(locus_tag_df,grepl("MVB",phenotype)) %>% as.data.frame
  MVB_all <- ifelse(nrow(MVB_all)>0,apply(MVB_all,1,FUN=function(x){paste(x["phenotype"],x["test"],x["model"],x["population"],sep = "_")}) %>% paste0(.,collapse=";"),"")
  MVB_sync <- str_split(MVB_all,";") %>% unlist %>% .[grepl("synchronous",.)] %>% gsub("synchronous_","",.) %>% paste0(.,collapse=";")
  MVB_phyc <- str_split(MVB_all,";") %>% unlist %>% .[grepl("phyc",.)] %>% gsub("phyc_","",.)  %>% paste0(.,collapse=";")
  MVB_cont <- str_split(MVB_all,";") %>% unlist %>% .[grepl("continuous",.)] %>% gsub("continuous_","",.)  %>% paste0(.,collapse=";")
  ##IR
  IR_all <- subset(locus_tag_df,grepl("IR",phenotype)) %>% as.data.frame
  IR_all <- ifelse(nrow(IR_all)>0,apply(IR_all,1,FUN=function(x){paste(x["phenotype"],x["test"],x["model"],x["population"],sep = "_")}) %>% paste0(.,collapse=";"),"")
  IR_sync <- str_split(IR_all,";") %>% unlist %>% .[grepl("synchronous",.)] %>% gsub("synchronous_","",.) %>% paste0(.,collapse=";")
  IR_phyc <- str_split(IR_all,";") %>% unlist %>% .[grepl("phyc",.)] %>% gsub("phyc_","",.)  %>% paste0(.,collapse=";")
  IR_cont <- str_split(IR_all,";") %>% unlist %>% .[grepl("continuous",.)] %>% gsub("continuous_","",.)  %>% paste0(.,collapse=";")
  ##Any Variable
  blbli_any_sig <-  ifelse(nchar(blbli_all)>0,1,0) %>% as.numeric
  MVB_any_sig <- ifelse(nchar(MVB_all)>0,1,0) %>% as.numeric
  IR_any_sig <- ifelse(nchar(IR_all)>0,1,0) %>% as.numeric
  both_any_sig <- ifelse(MVB_any_sig ==1 & IR_any_sig ==1,1,0)
  ##Any Synchronous
  blbli_sync_sig <- ifelse(nchar(blbli_sync)>1,1,0)
  MVB_sync_sig <- ifelse(nchar(MVB_sync)>1,1,0)
  IR_sync_sig <- ifelse(nchar(IR_sync)>1,1,0)
  Both_sync_sig <- ifelse(MVB_sync_sig ==1 & IR_sync_sig ==1,1,0)
  any_sync <-ifelse(blbli_sync_sig==1 | MVB_sync_sig ==1 | IR_sync_sig ==1,1,0)
  ##Any Phyc
  blbli_phyc_sig <- ifelse(nchar(blbli_phyc)>1,1,0)
  MVB_phyc_sig <- ifelse(nchar(MVB_phyc)>1,1,0)
  IR_phyc_sig <- ifelse(nchar(IR_phyc)>1,1,0)
  Both_phyc_sig <- ifelse(MVB_phyc_sig ==1 & IR_phyc_sig ==1,1,0)
  any_phyc <-ifelse(blbli_phyc_sig==1 | MVB_phyc_sig ==1 | IR_phyc_sig ==1,1,0)
  ##Any Continuous
  blbli_cont_sig <- ifelse(nchar(blbli_cont)>1,1,0)
  MVB_cont_sig <- ifelse(nchar(MVB_cont)>1,1,0)
  IR_cont_sig <- ifelse(nchar(IR_cont)>1,1,0)
  Both_cont_sig <- ifelse(MVB_cont_sig ==1 & IR_cont_sig ==1,1,0)
  any_cont <-ifelse(MVB_cont_sig ==1 | IR_cont_sig ==1,1,0)
  ##Cont & Binary
  cont_only <- ifelse(sum(MVB_cont_sig,IR_cont_sig,blbli_cont_sig)>0 & sum(blbli_phyc_sig,MVB_phyc_sig,IR_phyc_sig,blbli_sync_sig,MVB_sync_sig,IR_sync_sig)==0,1,0)
  bin_only <- ifelse(sum(MVB_cont_sig,IR_cont_sig,blbli_cont_sig)==0 & sum(blbli_phyc_sig,MVB_phyc_sig,IR_phyc_sig,blbli_sync_sig,MVB_sync_sig,IR_sync_sig)>0,1,0)
  cont_bin_both <-  ifelse(sum(MVB_cont_sig,IR_cont_sig,blbli_cont_sig)>0 & sum(blbli_phyc_sig,MVB_phyc_sig,IR_phyc_sig,blbli_sync_sig,MVB_sync_sig,IR_sync_sig)>0,1,0)
  ##Model specific
  model1_sig <- ifelse(sum(grepl("model1",locus_tag_df$model))>0,1,0)
  model1.2_sig <- ifelse(sum(grepl("model1.2",locus_tag_df$model))>0,1,0)
  model2_sig <- ifelse(sum(grepl("model2",locus_tag_df$model))>0,1,0)
  model2.2_sig <- ifelse(sum(grepl("model2.2",locus_tag_df$model))>0,1,0)
  model3_sig <- ifelse(sum(grepl("model3",locus_tag_df$model))>0,1,0)
  model3.2_sig <- ifelse(sum(grepl("model3.2",locus_tag_df$model))>0,1,0)
  across_all_models <- ifelse(sum(model1_sig,model1.2_sig)>0 & sum(model2_sig,model2.2_sig)>0 & sum(model3_sig,model3.2_sig)>0,1,0)
  core_sig <- ifelse(sum(model1_sig,model1.2_sig,model2_sig,model2.2_sig,model3_sig,model3.2_sig)>0,1,0)
  kleborate_sig <- ifelse(sum(grepl("kleborate",locus_tag_df$model))>0,1,0)
  kpc_plasmid_sig <- ifelse(sum(grepl("KPC",locus_tag_df$model))>0,1,0)

  summary_df <- cbind(locus_tag,blbli_any_sig,MVB_any_sig,IR_any_sig,both_any_sig,cont_only,bin_only,cont_bin_both,blbli_sync_sig,MVB_sync_sig,IR_sync_sig,Both_sync_sig,any_sync,blbli_phyc_sig,MVB_phyc_sig,IR_phyc_sig,Both_phyc_sig,any_phyc,blbli_cont_sig,MVB_cont_sig,IR_cont_sig,Both_cont_sig,any_cont,model1_sig,model1.2_sig,model2_sig,model2.2_sig,model3_sig,model3.2_sig,across_all_models,core_sig,kleborate_sig,kpc_plasmid_sig) %>% as.data.frame

  return(summary_df)
}

## 6. Get sig tests' max model and sig models
get_sig_tests <- function(data_row){
  models_info <- data_row %>% str_split(.,";",simplify = T) %>% str_split(.,"_",simplify = T) %>% as.data.frame
  models <- models_info  %>% .[,ncol(.)] %>% gsub("model","",.) %>% sort %>% unique()
  max_model <- models  %>% max()
  sig_models <- models %>% paste0(.,collapse=";")
  data_row <- cbind(sig_models,max_model) %>% as.data.frame
  return(data_row)
}

# Generate additional data
## 1. Get marginal stats (i.e., sens, spec, npv, ppv)# Marginal Statistics
get_marginal_stats <- function(feature,outcome){
  if(class(outcome) != "numeric"){
    stop("Outcome is not numeric")
  }
  if(class(feature) != "numeric"){
    stop("Feature is not numeric")
  }

  n <- length(outcome)
  phenotype_count = sum(outcome)
  phenotype_freq = phenotype_count / n
  feature_count = sum(feature)
  feature_freq = feature_count / n

  # Cells
  tp = sum(outcome==1 & feature==1)
  fp = sum(outcome==0 & feature==1)
  tn = sum(outcome==0 & feature==0)
  fn = sum(outcome==1 & feature==0)
  # Summary stats
  sens = round((tp / c(tp + fn)),3)
  spec = round((tn / c(tn + fp)),3)
  ppv = round((tp / c(tp + fp)),3)
  npv = round((tn / c(tn + fn)),3)
  accuracy =  round(((tp+tn)/n),3)

  data <- as.data.frame(cbind(phenotype = colnames(outcome),phenotype_count,phenotype_freq,feature_count,feature_freq,tp,fp,fn,tn,sens,spec,ppv,npv,accuracy))
  rownames(data) <- colnames(feature)
  data$name <-colnames(feature)
  return(data)
}

## Get sensitivity and ppv for
### a. core
get_sensitivity_ppv_core <- function(locus_tag,df){
  models <- c("_1","_1.2","_2","_2.2","_3","_3.2")
  models <- paste0(locus_tag,models)
  dataframe <- df %>% select(any_of(c(models,"blbli_dich_num","MVB_dich_num","IR_dich_num"))) %>% `rownames<-`(.$isolate_no)
  models_in_df <- colnames(dataframe) %>% subset(!. %in% c("blbli_dich_num","MVB_dich_num","IR_dich_num"))
  ## blbli SENS & PPV
  blbli_sens <- lapply(models_in_df,FUN=function(x){get_marginal_stats(dataframe[,x],dataframe[,"blbli_dich_num"])}) %>% do.call(rbind,.)%>% select(sens)%>% t %>% `colnames<-`(paste0("sens_blbli_model",gsub(locus_tag,"",models_in_df))) %>% as.data.frame %>%  mutate(locus_tag = locus_tag)
  blbli_ppv <- lapply(models_in_df,FUN=function(x){get_marginal_stats(dataframe[,x],dataframe[,"blbli_dich_num"])}) %>% do.call(rbind,.)%>% select(ppv) %>% t %>% `colnames<-`(paste0("ppv_blbli_model",gsub(locus_tag,"",models_in_df)))%>% as.data.frame%>% mutate(locus_tag = locus_tag)
  ## MVB SENS & PPV
  MVB_sens <- lapply(models_in_df,FUN=function(x){get_marginal_stats(dataframe[,x],dataframe[,"MVB_dich_num"])}) %>% do.call(rbind,.)%>% select(sens) %>% t %>% `colnames<-`(paste0("sens_MVB_model",gsub(locus_tag,"",models_in_df)))%>% as.data.frame %>% mutate(locus_tag = locus_tag)
  MVB_ppv <- lapply(models_in_df,FUN=function(x){get_marginal_stats(dataframe[,x],dataframe[,"MVB_dich_num"])}) %>% do.call(rbind,.)%>% select(ppv) %>% t %>% `colnames<-`(paste0("ppv_MVB_model",gsub(locus_tag,"",models_in_df)))%>% as.data.frame %>% mutate(locus_tag = locus_tag) ## IR SENS & PPV
  IR_sens <- lapply(models_in_df,FUN=function(x){get_marginal_stats(dataframe[,x],dataframe[,"IR_dich_num"])}) %>% do.call(rbind,.)%>% select(sens) %>% t %>% `colnames<-`(paste0("sens_IR_model",gsub(locus_tag,"",models_in_df)))%>% as.data.frame %>% mutate(locus_tag = locus_tag)
  IR_ppv <- lapply(models_in_df,FUN=function(x){get_marginal_stats(dataframe[,x],dataframe[,"IR_dich_num"])}) %>% do.call(rbind,.)%>% select(ppv)%>% t %>% `colnames<-`(paste0("ppv_IR_model",gsub(locus_tag,"",models_in_df)))%>% as.data.frame %>% mutate(locus_tag = locus_tag) ## IR SENS & PPV
  result <-left_join(blbli_sens,blbli_ppv) %>% left_join(.,MVB_sens) %>% left_join(.,MVB_ppv) %>% left_join(.,IR_sens) %>% left_join(.,IR_ppv)
  return(result)
}

get_ppv_sens_datasets <- function(sig_hits_df){
  ppv_blbli_data <- sig_hits_df[,grepl("ppv_blbli",colnames(sig_hits_df))]
  ppv_MVB_data <- sig_hits_df[,grepl("ppv_MVB",colnames(sig_hits_df))]
  ppv_IR_data <- sig_hits_df[,grepl("ppv_IR",colnames(sig_hits_df))]
  sens_blbli_data <- sig_hits_df[,grepl("sens_blbli",colnames(sig_hits_df))]
  sens_MVB_data <- sig_hits_df[,grepl("sens_MVB",colnames(sig_hits_df))]
  sens_IR_data <- sig_hits_df[,grepl("sens_IR",colnames(sig_hits_df))]
  models <- str_split(colnames(ppv_MVB_data),"_") %>% lapply(.,FUN=function(x){x %>% unlist %>% .[length(.)]}) %>% as.vector
  ppv_blbli_df <- cbind(ppv_blbli_data %>% t %>% as.vector,models)%>% as.data.frame  %>% `colnames<-`(c("blbli_PPV","model"))  %>% mutate_if(is.list,unlist)
  ppv_MVB_df <- cbind(ppv_MVB_data %>% t %>% as.vector,models)%>% as.data.frame  %>% `colnames<-`(c("MVB_PPV","model"))  %>% mutate_if(is.list,unlist)
  ppv_IR_df <- cbind(ppv_MVB_data %>% t %>% as.vector,models)%>% as.data.frame %>% `colnames<-`(c("IR_PPV","model")) %>% mutate_if(is.list,unlist)
  sens_blbli_df <- cbind(sens_blbli_data %>% t %>% as.vector,models)%>% as.data.frame  %>% `colnames<-`(c("blbli_sens","model"))  %>% mutate_if(is.list,unlist)
  sens_MVB_df <- cbind(sens_MVB_data %>% t %>% as.vector,models)%>% as.data.frame  %>% `colnames<-`(c("MVB_sens","model"))  %>% mutate_if(is.list,unlist)
  sens_IR_df <- cbind(sens_IR_data %>% t %>% as.vector,models)%>% as.data.frame %>% `colnames<-`(c("IR_sens","model")) %>% mutate_if(is.list,unlist)
  ppv_sens_dataset <- left_join(ppv_blbli_df,ppv_MVB_df)  %>% left_join(.,ppv_IR_df) %>% left_join(.,sens_MVB_df) %>% left_join(.,sens_IR_df)%>% left_join(.,sens_blbli_df)
  return(ppv_sens_dataset)
}

get_max_model <- function(dataset){
  dataset <- dataset %>% subset(MVB_PPV != "NA")
  locus_tag  <-  dataset[["locus_tag"]]  %>% unique
  # blbli
  blbli_max_sens <- dataset[max.col(t(dataset$blbli_sens),ties.method = "last"),"model"]
  blbli_max_sens_ties <- ifelse(nrow(subset(dataset,blbli_sens == max(dataset$blbli_sens)))>1,"yes","no")
  blbli_max_ppv <- dataset[max.col(t(dataset$blbli_PPV),ties.method = "last"),"model"]
  blbli_max_ppv_ties <- ifelse(nrow(subset(dataset,blbli_PPV == max(dataset$blbli_PPV)))>1,"yes","no")
  blbli_max_sens_max_ppv_same <- ifelse(blbli_max_sens == blbli_max_ppv,1,0)
  # MVB
  MVB_max_sens <- dataset[max.col(t(dataset$MVB_sens),ties.method = "last"),"model"]
  MVB_max_sens_ties <- ifelse(nrow(subset(dataset,MVB_sens == max(dataset$MVB_sens)))>1,"yes","no")
  MVB_max_ppv <- dataset[max.col(t(dataset$MVB_PPV),ties.method = "last"),"model"]
  MVB_max_ppv_ties <- ifelse(nrow(subset(dataset,MVB_PPV == max(dataset$MVB_PPV)))>1,"yes","no")
  MVB_max_sens_max_ppv_same <- ifelse(MVB_max_sens == MVB_max_ppv,1,0)
  # IR
  IR_max_sens <- dataset[max.col(t(dataset$IR_sens),ties.method = "last"),"model"]
  IR_max_sens_ties <- ifelse(nrow(subset(dataset,IR_sens == max(dataset$IR_sens)))>1,"yes","no")
  IR_max_ppv <- dataset[max.col(t(dataset$IR_PPV),ties.method = "last"),"model"]
  IR_max_ppv_ties <- ifelse(nrow(subset(dataset,IR_PPV == max(dataset$IR_PPV)))>1,"yes","no")
  IR_max_sens_max_ppv_same <- ifelse(IR_max_sens == IR_max_ppv,1,0)
  # MVB IR
  MVB_IR_max_sens_same <- ifelse(MVB_max_sens == IR_max_sens,1,0)
  MVB_IR_max_ppv_same <- ifelse(MVB_max_ppv == IR_max_ppv,1,0)
  MVB_IR_max_sens_max_ppv_same <- ifelse(c(MVB_max_sens,MVB_max_ppv,IR_max_ppv,IR_max_sens) %>% unlist %>% unique %>% length(.) ==1,1,0)
  IR_max_sens_val <- max(dataset$IR_sens)
  MVB_max_sens_val <- max(dataset$MVB_sens)

  df_row <- cbind(locus_tag = locus_tag,best_blbli_sens = blbli_max_sens,blbli_max_sens_ties, best_blbli_ppv=blbli_max_ppv,blbli_max_ppv_ties, blbli_max_sens_max_ppv_same,best_MVB_sens = MVB_max_sens,MVB_max_sens_ties, best_MVB_ppv=MVB_max_ppv,MVB_max_ppv_ties, MVB_max_sens_max_ppv_same,best_IR_sens =IR_max_sens, IR_max_sens_ties,best_IR_ppv = IR_max_ppv,IR_max_ppv_ties,IR_max_sens_max_ppv_same,MVB_IR_max_sens_same,MVB_IR_max_ppv_same,MVB_IR_max_sens_max_ppv_same,MVB_max_sens_val,IR_max_sens_val) %>% as.data.frame
  return(df_row)
}

### b. generic
get_sensitivity_ppv_generic <- function(locus_tag,df){
  dataframe <- df %>% select(any_of(c(locus_tag,"blbli_dich_num","MVB_dich_num","IR_dich_num"))) %>% `rownames<-`(.$isolate_no)
  ## blbli SENS & PPV
  blbli_sens <- get_marginal_stats(dataframe[,locus_tag],dataframe[,"blbli_dich_num"]) %>% select(sens) %>% t %>% `colnames<-`('blbli_sens') %>% as.data.frame %>% mutate_all(as.numeric) %>% mutate(locus_tag = locus_tag)  %>% as.data.frame
  blbli_ppv <- get_marginal_stats(dataframe[,locus_tag],dataframe[,"blbli_dich_num"]) %>% select(ppv) %>% t %>% `colnames<-`('blbli_ppv') %>% as.data.frame  %>% mutate_all(as.numeric) %>% mutate(locus_tag = locus_tag)  %>% as.data.frame
  ## MVB SENS & PPV
  MVB_sens <- get_marginal_stats(dataframe[,locus_tag],dataframe[,"MVB_dich_num"]) %>% select(sens) %>% t %>% `colnames<-`('MVB_sens') %>% as.data.frame %>% mutate_all(as.numeric) %>% mutate(locus_tag = locus_tag)  %>% as.data.frame
  MVB_ppv <- get_marginal_stats(dataframe[,locus_tag],dataframe[,"MVB_dich_num"]) %>% select(ppv) %>% t %>% `colnames<-`('MVB_ppv') %>% as.data.frame  %>% mutate_all(as.numeric) %>% mutate(locus_tag = locus_tag)  %>% as.data.frame
  ## IR SENS & PPV
  IR_sens <- get_marginal_stats(dataframe[,locus_tag],dataframe[,"IR_dich_num"]) %>% select(sens) %>% t %>% `colnames<-`('IR_sens') %>% as.data.frame  %>% mutate_all(as.numeric) %>% mutate(locus_tag = locus_tag)  %>% as.data.frame
  IR_ppv <- get_marginal_stats(dataframe[,locus_tag],dataframe[,"IR_dich_num"]) %>% select(ppv) %>% t %>% `colnames<-`('IR_ppv') %>% as.data.frame  %>% mutate_all(as.numeric) %>% mutate(locus_tag = locus_tag)  %>% as.data.frame
  result <-left_join(blbli_sens,blbli_ppv) %>% left_join(.,MVB_sens) %>% left_join(.,MVB_ppv) %>% left_join(.,IR_sens) %>% left_join(.,IR_ppv)
  return(result)
}

## Get proportion in clustering
get_prop_in_res_isolates <- function(variable,df){
  isolates_w_variant <- df %>% subset(get(variable) ==1)
  representation_in_clusters <-  table(isolates_w_variant %>% .$blbli_asr_cluster_renamed)
  true_cluster_size <- table(df %>% .$blbli_asr_cluster_renamed)

  res_prop_df <-lapply(names(true_cluster_size),FUN=function(x){
    if(x %in% names(representation_in_clusters) == F){
      prop =0
    }
    if(x %in% names(representation_in_clusters) == T){
      prop = representation_in_clusters[[x]] / true_cluster_size[[x]]
    }
    return(prop)
  }) %>% do.call(cbind.data.frame,.) %>% `colnames<-`(paste0(names(true_cluster_size),"_prop"))

  res_prop_df$Cluster_prop <- nrow(isolates_w_variant %>% subset(!blbli_asr_cluster_renamed %in% c("Singleton","No Feature"))) / nrow(df %>% subset(!blbli_asr_cluster_renamed %in% c("Singleton","No Feature")))
  return(res_prop_df)
}

## Get nearest neighbor best model
generate_nn_data_best_model <- function(variable,nn_data){
  NN_data <- nn_data[[variable]]$nn_comps_with_variant_sum %>% as.data.frame
  NN_data$locus_tag <-  gsub("_1|_2|_1.2|_2.2|_3|_3.2","",variable)
  return(NN_data)
}

# Cleaning the data
## Resistance category
get_resistance_category_data <- function(gwas_table){
  gwas_table$resistance_category  <- ifelse(gwas_table$cont_only ==1 | gwas_table$any_phyc==1 & gwas_table$any_sync ==0, "Modulator",
                                            ifelse(gwas_table$any_sync == 1,"Mediator",""))
  return(gwas_table)
}


## Get significant tests
get_significant_tests <- function(gwas_table){
  gwas_table$sig_tests <-
    cbind(ifelse(gwas_table$any_cont==1,"Continuous",""),
          ifelse(gwas_table$any_phyc==1,"PhyC",""),
          ifelse(gwas_table$any_sync==1,"Synchronous","")) %>% apply(.,MARGIN=1,paste0,collapse =" ") %>% trimws(.,which="both") %>% ifelse(.=="Continuous PhyC Synchronous","All tests",.) %>% ifelse(. == "Continuous  Synchronous","Continuous & Synchronous",.) %>% ifelse(.=="PhyC Synchronous","PhyC & Synchronous",.)
  return(gwas_table)
}

## Get significant model
get_significant_model <- function(gwas_table){
  gwas_table$sig_models_simple <-  cbind(ifelse(gwas_table$model1_sig ==1 | gwas_table$model1.2_sig ==1,"All variants",""),
                                         ifelse(gwas_table$model2_sig ==1 | gwas_table$model2.2_sig ==1,"Non-synonymous",""),
                                         ifelse(gwas_table$model3_sig ==1 | gwas_table$model3.2_sig ==1,"High-risk variants",""),
                                         ifelse(gwas_table$kleborate_sig==1,"Kleborate",""),ifelse(gwas_table$kpc_plasmid_sig==1,"KPC-associated plasmids","")) %>% apply(.,MARGIN=1,paste0,collapse =" ") %>% trimws(.,which="both") %>% ifelse(.=="All variants Non-synonymous High-risk variants","All variant models",.) %>% gsub("s N","s & N",.)
  gwas_table$sig_models_simple <- factor(gwas_table$sig_models_simple,levels = c("All variants","Non-synonymous",'All variants & Non-synonymous',"High-risk variants",'All variant models',"Kleborate","KPC-associated plasmids"))
  return(gwas_table)
}

## Sig Phenotypes
get_sig_phenotypes <- function(gwas_table){
  gwas_table$sig_phenotypes <- cbind(ifelse(gwas_table$blbli_any_sig==1,"BL/BLI",""),
                                    ifelse(gwas_table$IR_any_sig ==1,"IR",""),
                                    ifelse(gwas_table$MVB_any_sig==1,"MVB",""))  %>% apply(.,MARGIN=1,paste0,collapse=' ') %>% trimws(.,which="both")  %>% gsub("  "," ",.) %>% gsub(" "," & ",.)
  return(gwas_table)
}

## NN flag
get_nn_flag <- function(gwas_table){
  gwas_table$nn_qc <- ifelse(c(abs(gwas_table$MVB_num_log_2_diff_median) >1 | abs(gwas_table$IR_num_log_2_diff_median) >1),"Yes","No")
  return(gwas_table)
}

# Get variant frequency in our phenotypes
variants_freq <- function(locus_tags,df){
  MVB_sus<- subset(df,MVB_dich != "Non-Susceptible")
  MVB_res<- subset(df,MVB_dich == "Non-Susceptible")
  IR_sus<- subset(df,IR_dich != "Non-Susceptible")
  IR_res <- subset(df,IR_dich == "Non-Susceptible")
  blbli_sus<- subset(df,blbli_dich_num != 1)
  blbli_res <- subset(df,blbli_dich_num == 1)

  variant_freq <- function(locus_tag){
    overall_freq <- sum(df[,locus_tag])
    overall_prop <- overall_freq / nrow(df)
    blbli_sus_freq <- sum(blbli_sus[,locus_tag])
    blbli_sus_prop <- blbli_sus_freq / nrow(blbli_sus)
    blbli_res_freq <- sum(blbli_res[,locus_tag])
    blbli_res_prop <-    blbli_res_freq / nrow(blbli_res)
    MVB_sus_freq <- sum(MVB_sus[,locus_tag])
    MVB_sus_prop <- MVB_sus_freq / nrow(MVB_sus)
    MVB_res_freq <- sum(MVB_res[,locus_tag])
    MVB_res_prop <- MVB_res_freq / nrow(MVB_res)
    IR_sus_freq <- sum(IR_sus[,locus_tag])
    IR_sus_prop <- IR_sus_freq / nrow(IR_sus)
    IR_res_freq<- sum(IR_res[,locus_tag])
    IR_res_prop <- IR_res_freq / nrow(IR_res)
    var_freq_df <- cbind.data.frame(locus_tag,
                                    genotype=locus_tag,
                                    overall_freq,overall_prop,
                                    blbli_sus_freq,blbli_sus_prop,blbli_res_freq,blbli_res_prop,
                                    MVB_sus_freq,MVB_sus_prop,MVB_res_freq,MVB_res_prop,
                                    IR_sus_freq,IR_sus_prop,IR_res_freq,IR_res_prop)


    return(var_freq_df)
  }

  vars_freq_df <- lapply(locus_tags,FUN=variant_freq) %>% do.call(rbind,.)

  vars_freq_df$locus_tag <- vars_freq_df$locus_tag %>% gsub("_1|_2|_1.2|_2.2|_3|_3.2","",.)
  return(vars_freq_df)
}

freq_plot <- function(hits_data){
  melted <- hits_data %>%  select(locus_tag,genotype,blbli_sus_prop,blbli_res_prop,`Significant Hit`) %>% reshape2::melt()
  figure <- ggplot(data=melted,aes(fill=variable,y=`Significant Hit`,x = value)) + geom_bar(position="dodge", stat="identity") + xlim(0,1) + resistance_prop_scale + theme_bw()+ theme(legend.position="bottom",axis.text = element_text(size=56,colour = "black"),axis.title = element_text(size=56,colour = "black"),legend.title = element_text(size=56,colour = "black"),legend.text = element_text(size=56,colour = "black"),legend.title.align=0.5 ,legend.key.size = unit(1, "cm"),legend.key.width = unit(1, "cm")) + xlab("Proportion of Isolates with Genotype")  +ylab("")
  return(figure)
}

# NN Data
generate_nn_melt_data <- function(variable,nn_data){
  NN_data <- nn_data[[variable]][['nn_comps_with_variant']]   %>% mutate(locus_tag = rep(variable,nrow(.))) %>% select(locus_tag,MVB_num_log_2_diff,IR_num_log_2_diff)
  nn_data_melt <- NN_data %>% reshape2::melt()
  nn_data_melt$locus_tag <-  gsub("_1|_2|_1.2|_2.2|_3|_3.2","",nn_data_melt$locus_tag)
  return(nn_data_melt)
}

nn_plot <- function(nn_data_melt){
  nn_data_melt$variable <- factor(nn_data_melt$variable,levels = rev(unique(nn_data_melt$variable)))
  figure <- ggplot(data=nn_data_melt,aes(fill=variable,y=locus_tag,x = value)) + geom_boxplot(outlier.size=3,color="black",linewidth = 0.5) + xlab("Log-2 Fold Change in MIC Relative to Nearest Neighbor") + ylab("") + xlim(min(nn_data_melt$value)-1,max(nn_data_melt$value)+1) + MVB_IR_scale + theme_bw() + theme(legend.position="bottom",axis.text = element_text(size=56,colour = "black"),axis.title = element_text(size=56,colour = "black"),legend.title = element_text(size=56,colour = "black"),legend.text = element_text(size=56,colour = "black"),legend.title.align=0.5,legend.key.size = unit(2.5, "cm"),legend.key.width = unit(2.5, "cm")) + geom_vline(xintercept = 0,colour = "black")
  return(figure)
}

# Figures
gwas_figures <- function(df,tr,gwas_mat,sig_hits_name){
  rownames(df) <- df$isolate_no
  #Step #1: Clades
  p.1 <- gheatmap(ggtree(tr),df %>% select(clade_I) %>% `colnames<-`("ST258 Clade"), colnames_position = "top",colnames_angle=90, colnames_offset_y = 0.25, hjust = 0, color = NA, font.size = 20, width = 0.1)   + clade_colors_scale_v + consistent_theme_GWAS

  #Step #3: BL/BLI Cluster
  p.2 <- p.1 + new_scale_fill()
  p.3 <-  gheatmap(p.2,df %>% select(blbli_asr_cluster_renamed) %>% mutate_all(as.factor) %>% `colnames<-`("BL/BLI Clustering"), colnames_position = "top",colnames_angle=90, colnames_offset_y = 0.25, hjust = 0, color = NA, font.size = 20, width = 0.1,offset =.000025) + cluster_scale + consistent_theme_GWAS

  #Step #4: MVB Binary
  p.4 <- p.3 + new_scale_fill()
  p.5 <-  gheatmap(p.4,df %>% select(blbli_dich,MVB_dich,IR_dich) %>% `colnames<-`(c("BL/BLI Resistance","MVB Resistance","IR Resistance")), colnames_position = "top",colnames_angle=90, colnames_offset_y = 0.25, hjust = 0, color = NA, font.size = 20, width = 0.3,offset =.00005) + resistance_scale + consistent_theme_GWAS

  #Step #6: MVB MIC
  p.6 <- p.5 + new_scale_fill()
  p.7 <-  gheatmap(p.6,df %>% select(MVB_log_2,IR_log_2) %>%   mutate_all(as.factor)  %>%`colnames<-`(c("MVB MIC","IR MIC")), colnames_position = "top",colnames_angle=90, colnames_offset_y = 0.25, hjust = 0, color = NA, font.size = 20, width = 0.2,offset =.000125) + Log2_scale + consistent_theme_GWAS

  # Step #7: Add
  p.7.1 <- p.7 + new_scale_fill()
  p.8 <-  gheatmap(p.7.1,gwas_mat %>% mutate_all(as.factor) %>% `colnames<-`(sig_hits_name), colnames_position = "top",colnames_angle=90, colnames_offset_y = 0.25, hjust = 0, color = NA, font.size = 20, width = 0.1*ncol(gwas_mat),offset =.000175)  + potential_scale_general + consistent_theme_GWAS
  return(p.8)
}

get_gwas_hit_matrix <- function(hits,df){
  resistance_modulator <- hits %>% subset(resistance_category == "Modulator") %>% select(genotype) %>% unlist %>% `names<-`(NULL)
  resistance_mediator <-  hits %>% subset(resistance_category == "Mediator") %>% select(genotype) %>% unlist  %>% `names<-`(NULL)

  df_mat <- cbind(df %>% select(resistance_modulator) %>% {ifelse(.==1,"Modulator","None")},
                  df %>% select(resistance_mediator) %>% {ifelse(.==1,"Mediator","None")}) %>% as.data.frame %>% `rownames<-`(df$isolate_no)

  return(df_mat)
}

