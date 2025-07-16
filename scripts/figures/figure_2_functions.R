# Figure 2 functions
MIC_prop_hist <- function(dataset){ 
  # Histogram
  histogram <-  ggplot(data=dataset,aes(x=clade_I,fill=fct_rev(phenotype_cat),y=freq)) + geom_bar(position ='stack',stat="identity") +ylab("Percentage") +xlab("") + labs(fill="Resistance profile") + theme_bw_me
  return(histogram)
}

MIC_prop_hist_count <- function(dataset){ 
  # Histogram
  histogram <-  ggplot(data=dataset,aes(x=clade_I,fill=fct_rev(phenotype_cat),y=n)) + geom_bar(position = 'stack',stat="identity") +ylab("Isolates") +xlab("") + labs(fill="Resistance profile") + theme_bw_me
  return(histogram)
}

