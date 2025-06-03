# Figure 2 functions
MIC_prop_hist <- function(dataset){ 
  # Histogram
  histogram <-  ggplot(data=dataset,aes(x=clade_I,fill=fct_rev(phenotype_cat),y=freq)) + geom_bar(position ='stack',stat="identity") +ylab("Proportion (%)") +xlab("") + labs(fill=" Resistance Profile") + theme_bw_me
  return(histogram)
}

MIC_prop_hist_count <- function(dataset){ 
  # Histogram
  histogram <-  ggplot(data=dataset,aes(x=clade_I,fill=fct_rev(phenotype_cat),y=n)) + geom_bar(position = 'stack',stat="identity") +ylab("No. isolates") +xlab("") + labs(fill=" Resistance Profile") + theme_bw_me
  return(histogram)
}

