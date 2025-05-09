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