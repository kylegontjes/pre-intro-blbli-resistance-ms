
MIC_ompk36_prop_hist <- function(dataset){
  # n per resistance category (summed across OmpK36 groups)
  n_cat <- dataset %>%
    dplyr::group_by(phenotype_cat) %>%
    dplyr::summarise(n = sum(n), .groups = "drop")   # needs a per-row count column `n`
  
  # Build x labels: padded category name + n on a second line
  padding <- max(nchar(as.character(dataset[['phenotype_cat']])))
  lab_df <- n_cat %>%
    dplyr::mutate(lab = paste0(str_pad(as.character(phenotype_cat), width = padding, side = "both"),
                               "\n(n = ", n, ")"))
  lab_vec <- setNames(lab_df$lab, as.character(lab_df$phenotype_cat))
  
  histogram <- ggplot(data = dataset, aes(x = phenotype_cat, fill = OmpK36_mutations, y = freq)) +
    geom_bar(position = "stack", stat = "identity", width = 0.55) +
    scale_x_discrete(labels = lab_vec) +
    ylab("Proportion of Isolates") + xlab("Resistance Profile")
  return(histogram)
}

KPC_copy_number_ompk_hist <- function(dataset,y,x,color,title,colors){
  data <- dataset %>% select(y,x,color) %>% `colnames<-`(c("y","x","color")) 
  figure <- ggplot(data= data) + geom_point(aes(y=y, x=log2(x),color=color) ) + theme(legend.position = "bottom")  + theme_bw_me       + geom_smooth(method='lm',mapping = aes(x=log2(x),y=y,color=color)) + ggtitle(title)
  return(figure)
}

KPC_copy_number_ompk_hist_facet <- function(dataset,y,x,color,title,colors){
  data <- dataset %>% select(y,x,color) %>% `colnames<-`(c("y","x","color")) 
  figure <- ggplot(data= data) + geom_point(aes(y=y, x=log2(x),color=color) ) + theme(legend.position = "bottom")  + theme_bw_me       + geom_smooth(method='lm',mapping = aes(x=log2(x),y=y,color=color)) + ggtitle(title) + facet_wrap(~color)
  return(figure)
}