# Figure 2 functions
MIC_prop_hist <- function(dataset){
  # n per clade for axis labels
  n_lab <- dataset %>%
    dplyr::group_by(clade_I) %>%
    dplyr::summarise(n = sum(n), .groups = "drop")
  clade_labs <- setNames(
    paste0(n_lab$clade_I, "\n(n = ", n_lab$n, ")"),
    as.character(n_lab$clade_I)
  )
  # Histogram
  histogram <- ggplot(data=dataset, aes(x=clade_I, fill=fct_rev(phenotype_cat), y=freq)) +
    geom_bar(position='stack', stat="identity") +
    ylab("Percentage") + xlab("") +
    labs(fill="Resistance profile") +
    scale_x_discrete(labels = clade_labs) +
    theme_bw_me
  return(histogram)
}

MIC_prop_hist_count <- function(dataset){
  # n per clade for axis labels
  n_lab <- dataset %>%
    dplyr::group_by(clade_I) %>%
    dplyr::summarise(n = sum(n), .groups = "drop")
  clade_labs <- setNames(
    paste0(n_lab$clade_I, "\n(n = ", n_lab$n, ")"),
    as.character(n_lab$clade_I)
  )
  # Histogram
  histogram <- ggplot(data=dataset, aes(x=clade_I, fill=fct_rev(phenotype_cat), y=n)) +
    geom_bar(position='stack', stat="identity") +
    ylab("Isolates") + xlab("") +
    labs(fill="Resistance profile") +
    scale_x_discrete(labels = clade_labs) +
    theme_bw_me
  return(histogram)
}