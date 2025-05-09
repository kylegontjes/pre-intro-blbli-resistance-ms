# Figure 2 functions
MIC_prop_hist <- function(dataset,title,name){ 
  # Histogram
  histogram <-  ggplot(data=dataset,aes(x=clade_I,fill=phenotype_cat,y=freq)) + geom_bar(position = position_fill(reverse = TRUE),stat="identity")  +ylab("Proportion of Isolates") +xlab("") + labs(fill=" Resistance Profile")     + ggtitle(paste0(title))  + theme_bw_me
  return(histogram)
}

MIC_prop_hist_count <- function(dataset,title,name){ 
  # Histogram
  histogram <-  ggplot(data=dataset,aes(x=clade_I,fill=phenotype_cat,y=n))+geom_col()  +ylab("Proportion of Isolates") +xlab("") + labs(fill=" Resistance Profile")     + ggtitle(paste0(title))   + theme_bw_me
  return(histogram)
}

