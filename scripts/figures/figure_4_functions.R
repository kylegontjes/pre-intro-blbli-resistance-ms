
MIC_ompk36_prop_hist <- function(dataset){
  ## Color Scale 
  histogram <-  ggplot(data=dataset,aes(x=phenotype_cat,fill=OmpK36_mutations,y=freq))+geom_bar(position="stack",stat="identity", width = 0.75)  +ylab("Proportion of Isolates") +xlab("Resistance Profile") 
  return(histogram)
}

KPC_copy_number_ompk_hist <- function(dataset,y,x,color,title,colors){
  data <- dataset %>% select(y,x,color) %>% `colnames<-`(c("y","x","color")) 
  figure <- ggplot(data= data) + geom_point(aes(y=y, x=log2(x),color=color) ) + theme(legend.position = "bottom")  + theme_bw_me       + geom_smooth(method='lm',mapping = aes(x=log2(x),y=y,color=color)) + ggtitle(title)
  return(figure)
}