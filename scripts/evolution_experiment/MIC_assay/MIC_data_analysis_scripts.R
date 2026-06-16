# MIC assay functions
## Process plate 
process_plate <- function(start, end,plate) {
  # Plate name
  plate_name <- plate[[1]][start]
  # Extract block (rows A–H + maybe header row)
  block <- plate[(start + 1):end, ]
  # If your first row is column headers (1–12), keep this:
  colnames(block) <- c("row", as.character(1:12))
  # Remove any accidental header row 
  block <- block %>%  filter(!is.na(row))
  # Reshape
  block_long <- block %>%
    pivot_longer(
      cols = -row,
      names_to = "col",
      values_to = "value"
    ) %>%
    mutate(
      plate = plate_name %>% gsub('plate_',"",.) %>% as.numeric,
      col = as.integer(col),
      well = paste0(row, col)
    )
  return(block_long)
}

## Process data frame spreadsheet
process_dataframe <- function(plate){
  # Get starts and ends
  plate_starts <- c(which(grepl("^plate_", MIC_assay_MVB[[1]])))
  plate_ends <- c(plate_starts[-1] - 1, nrow(MIC_assay_MVB))
  
  # Get plate data
  plate_data <- map2_dfr(
    plate_starts,
    plate_ends,
    process_plate,
    plate = plate
  )
  return(plate_data)
}

# Get MIC value
get_MIC_value <- function(entry){ 
  isolate_df <- MIC_MVB_data %>% subset(isolate_no == entry)
  
  # Get first no growth
  first_no_growth <- isolate_df %>% group_by(isolate_no,isolate_type,bio_rep) %>% arrange(M_concentration) %>% subset(growth == F) %>% filter(row_number() == 1)
  # Ret replicate data
  replicate_summary <-  isolate_df %>% group_by(isolate_no,isolate_type,bio_rep) %>% summarise(wells = length(growth), growth_ct = sum(growth) , no_growth_ct = sum(growth==F)) 
  # Get check if no growth or all cells had growth 
  ## All growth (>256)
  replicate_summary$all_growth <- replicate_summary$wells == replicate_summary$growth_ct
  ## No growth (<=0.25)
  replicate_summary$all_no_growth <- replicate_summary$wells == replicate_summary$no_growth_ct 
  replicate_datas <- left_join(replicate_summary,first_no_growth %>% select(-isolate))
  replicate_datas$MIC <- ifelse(replicate_datas$all_growth, ">256",
                                ifelse(replicate_datas$all_no_growth, "<=0.25",
                                       replicate_datas$M_concentration)) %>% as.character()
  # Numeric MIC
  replicate_datas$MIC_num <-  gsub("<=|>", "", replicate_datas$MIC) %>%  as.numeric()
  # Log2 MIC
  replicate_datas$MIC_log_2 <-  log2(replicate_datas$MIC_num)
  
  return(replicate_datas)                            
}

## Summarize MIC data
MIC_data_summarize <- function(MIC_data){
  results <- MIC_data %>% group_by(isolate_no,isolate_type) %>% summarise(
  # Median values
  median_log_2_MIC = median(MIC_log_2),
  median_MIC = median(MIC_num),
  # Max values
  max_log_2_MIC = max(MIC_log_2),
  max_MIC = max(MIC_num),
  # Min values
  min_log_2_MIC = min(MIC_log_2),
  min_MIC = min(MIC_num),
  # Other clarifies
  unique_values = n_distinct(MIC_num),
  all_growth = sum(all_growth) ,
  all_no_growth = sum(all_no_growth))  
  return(results)
}

# MIC fold change analysis
MIC_fold_change_analysis <- function(isolate,MIC_data_summary){
  # Parent
  isolate_MIC_data_parent <- subset(MIC_data_summary,isolate_type=="parent" & isolate_no == isolate) %>% ungroup %>% select(-isolate_type)
  colnames(isolate_MIC_data_parent) <- c('isolate_no',paste0("parent_",colnames(isolate_MIC_data_parent %>% select(-isolate_no,))))
  
  # Resistant
  isolate_MIC_data_resistant <- subset(MIC_data_summary,isolate_type!="parent"& isolate_no ==isolate) %>% ungroup %>% select(-isolate_type)
  colnames(isolate_MIC_data_resistant) <- c('isolate_no',paste0("resistant_",colnames(isolate_MIC_data_resistant %>% select(-isolate_no,))))  
  
  # MIC data
  MIC_data <- left_join(isolate_MIC_data_parent,isolate_MIC_data_resistant) %>% mutate(
    # Median FC
    median_MIC_fc = resistant_median_MIC / parent_median_MIC,
    median_MIC_log_2_fc = log2(median_MIC_fc),
    change_median_log_2_MIC = resistant_median_log_2_MIC - parent_median_log_2_MIC,
    # Max FC
    max_MIC_fc = resistant_max_MIC / parent_max_MIC,
    max_MIC_log_2_fc = log2(max_MIC_fc),
    change_max_log_2_MIC = resistant_max_log_2_MIC - parent_max_log_2_MIC,
    # Min FC
    min_MIC_fc =  resistant_min_MIC / parent_min_MIC,
    min_MIC_log_2_fc = log2(min_MIC_fc),
    change_min_log_2_MIC = resistant_min_log_2_MIC - parent_min_log_2_MIC) 
  return(MIC_data)
}