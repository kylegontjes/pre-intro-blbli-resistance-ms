
MIC_ompk36_prop_hist <- function(dataset,title,colors){
  ## Color Scale
  colors_scale <- scale_fill_manual(values=colors, name="ompK36 Status", guide = guide_legend(nrow=1, title.position = "top", label.position = "right")) 
  histogram <-  ggplot(data=dataset,aes(x=phenotype_cat,fill=OmpK36_mutations,y=freq))+geom_bar(position="stack",stat="identity")  +ylab("Proportion of Isolates") +xlab("Resistance Profile") + ggtitle(paste0(title)) + theme(legend.position = "bottom",
                                                                                                                                                                                                                                                                                                      axis.text =   element_text(size=20,color="black"),
                                                                                                                                                                                                                                                                                                      axis.title = element_text(size = 24,color="black"),
                                                                                                                                                                                                                                                                                                      legend.text =   element_text(size=22,color="black"),
                                                                                                                                                                                                                                                                                                      legend.title = element_text(size = 24,color="black"),
                                                                                                                                                                                                                                                                                                      plot.title = element_text(size = 26,color="black"),
  ) + colors_scale
  return(histogram)
}

KPC_copy_number_ompk_hist <- function(dataset,y,x,color,title,colors){
  data <- dataset %>% select(y,x,color) %>% `colnames<-`(c("y","x","color")) 
  figure <- ggplot(data= data) + geom_point(aes(y=y, x=log2(x),color=color) ) + theme(legend.position = "bottom")  + theme_bw_me  + theme(legend.position = "bottom",
                                                                                                                                          axis.text =   element_text(size=20,color="black"),
                                                                                                                                          axis.title = element_text(size = 24,color="black"),
                                                                                                                                          legend.text =   element_text(size=22,color="black"),
                                                                                                                                          legend.title = element_text(size = 24,color="black"),
                                                                                                                                          plot.title = element_text(size = 26,color="black"),
  )   + geom_smooth(method='lm',mapping = aes(x=log2(x),y=y,color=color)) + ggtitle(title)
  return(figure)
}