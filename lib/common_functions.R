# Common functions
left_join_genomic_matrix <- function(matrix,df,number){
  df_fin <- df %>% left_join(.,matrix %>% `colnames<-`(paste0(colnames(.),number)) %>% mutate(isolate_no = rownames(.)))
  return(df_fin)
}

# Get presence absence data
get_presence_absence_matrix <- function(variable,df){
  presences <- str_split(unlist(c(df %>% dplyr::select(paste0(variable)))),";") %>% unlist(.)  %>% unique(.) %>% .[!. %in% "-"] %>% sort
  vector <- as.vector(tidyr::unite(df %>% dplyr::select(paste0(variable)),determinants,sep =";")) %>% unlist() %>% as.vector
  results <- matrix(nrow = nrow(df),ncol=length(presences))  %>% `colnames<-`(presences)
  for(i in presences){
    results[,paste0(i)] <- as.vector(ifelse(unlist(lapply(str_split(as.matrix(df[,paste0(variable)]),";"), function(x) i %in% x))==TRUE,1,0))
  }
  rownames(results) <- rownames(df)
  return(results)
}

# Functions for figure 1 & 2
## Create MIC histogram
create_MIC_histogram <- function(phenotype_cont,phenotype_cat = NULL,phenotype_name,title){ 
  data <- bind_cols(phenotype_cont,phenotype_cat,quiet = TRUE)
  colnames(data) <- c("phenotype_cont","phenotype_cat","variant_vector")
  histogram <-  ggplot(data) +geom_bar(aes(x=phenotype_cont,fill=phenotype_cat)) + ylab("Number of Isolates") + labs(fill="Resistance Profile") +ylab("Number of Isolates") +xlab(paste0(phenotype_name," MIC (µg/mL)"))  +  ggtitle(paste0(title))  +theme_bw_me 
  return(histogram)
} 

## Figure manipulations
# Get legend - comparable to cowplot's "get_legend"
get_legend <- function(plot){
  cowplot::get_plot_component(plot,"guide-box",return_all=T) %>% .[[3]]
}

## Create tableone output into a dataframe for export as a kable output
convert_tableone_into_df <- function(dataset,vars,strata=NULL,argsNormal=NULL,factorVars=NULL,outcome_names=NULL,exact=NULL){
  if(is.null(strata)==T){
    overall <- tableone::CreateTableOne(vars = vars, data=dataset,argsNormal = argsNormal,factorVars=factorVars)
    bound_table <-  capture.output(x <- print(overall, quote = FALSE, noSpaces = TRUE ))
    names <- x %>% as.matrix()  %>% rownames
    values <- x %>% as.data.frame %>% `rownames<-`(NULL)
    table <- cbind.data.frame(names,x %>% as.data.frame()) %>% `rownames<-`(NULL)
    colnames(table) <- c("Variable",paste0("Overall (n=",table[1,2],")"))
    rownames(table) <- NULL
    table <- table[-1,]
    return(table)
  }
  if(is.null(strata)==F){
    overall <- tableone::CreateTableOne(vars = vars, data=dataset,strata=strata,addOverall=T,argsNormal = argsNormal,factorVars=factorVars)
    bound_table <-  capture.output(x <- print(overall, quote = FALSE, noSpaces = TRUE,exact=exact))
    table <- x[,-ncol(x)]  %>% as.data.frame()
    names <- x %>% as.matrix()  %>% rownames
    table <- cbind.data.frame(names,table %>% as.data.frame()) %>% `rownames<-`(NULL) 
    colnames(table) <- c("Variable",paste0(c("Overall",outcome_names)," (n=",table[1,2:c(ncol(table)-1)],")"),"p-value")
    rownames(table) <- NULL
    table <- table[-1,]
    return(table)
  }
}

median_IQR <- function(pheno,df){
  median <- median(df[,pheno])
  IQR <- IQR(df[,pheno])
  median_IQR <- paste0(median," (",IQR,")") %>% `names<-`(paste0(pheno))
  return(median_IQR)
}

# Phylogeny
easy_resistance_figures <- function(tr,df){
  rownames(df) <- df$isolate_no
  #Step #1: Clades
  p.1 <- gheatmap(ggtree(tr),df %>% select(clade_I) %>% `colnames<-`("ST258 Clade"), colnames_position = "top",colnames_angle=90, colnames_offset_y = 0.25, hjust = 0, color = NA, font.size = 20, width = 0.1)   + clade_colors_scale_v
  
  #Step #3: BL/BLI Cluster
  p.2 <- p.1 + new_scale_fill()
  p.3 <-  gheatmap(p.2,df %>% select(blbli_asr_cluster_renamed) %>% mutate_all(as.factor) %>% `colnames<-`("BL/BLI Clustering"), colnames_position = "top",colnames_angle=90, colnames_offset_y = 0.25, hjust = 0, color = NA, font.size = 20, width = 0.1,offset =.000003) + cluster_scale + consistent_theme  
  
  #Step #4: MVB Binary
  p.4 <- p.3 + new_scale_fill()
  p.5 <-  gheatmap(p.4,df %>% select(blbli_dich,MVB_dich,IR_dich) %>% `colnames<-`(c("BL/BLI Non-Susceptible","MVB Non-Susceptible","IR Non-Susceptible")), colnames_position = "top",colnames_angle=90, colnames_offset_y = 0.25, hjust = 0, color = NA, font.size = 20, width = 0.3,offset =.000006) + NonSus_scale + consistent_theme  
  
  #Step #6: MVB MIC
  p.6 <- p.5 + new_scale_fill()
  p.7 <-  gheatmap(p.6,df %>% select(MVB_log_2,IR_log_2) %>%   mutate_all(as.factor)  %>%`colnames<-`(c("MVB MIC","IR MIC")), colnames_position = "top",colnames_angle=90, colnames_offset_y = 0.25, hjust = 0, color = NA, font.size = 20, width = 0.2,offset =.000015) + Log2_scale + consistent_theme  
  
  return(p.7)
}

# Marginal stats
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

# get pattern
get_pattern <- function(names,df){
  convert_to_name <- function(name,df){
    vector <- df %>% select(name)
    result <- ifelse(vector ==1,name,NA)
    return(result)
  }
  pattern <- lapply(names,convert_to_name,df=df) %>% do.call(cbind,.) %>%  apply(.,1,function(x) paste0(x,collapse = ';')) %>%  gsub("NA","",.) %>% trimws(.,which="both",whitespace = ";") 
  presence <-  ifelse(pattern=="",0,1)
  results <- cbind.data.frame(pattern,presence)
  return(results)
}