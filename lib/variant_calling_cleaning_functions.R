# Variant calling cleaning functions
return_grouped_matrix <- function(mat_list,annots_list,tr){
  merged_matrix <- bind_cols(lapply(mat_list,get)) %>% t %>% as.matrix()
  merged_annots <-  bind_rows(lapply(annots_list,get))
  merged_gene_list <- select(merged_annots,locus_tag) %>% as.vector() %>% unlist %>%  as.character()
  merged_grouped_matrix <- prewas::collapse_snps_into_genes(merged_matrix,merged_gene_list) %>% t %>% as.data.frame
  return(merged_grouped_matrix)
}

clean_snpkit_datasets <- function(variant_matrix){
  #Add Rownames Because Column Row is The Variant Data
  variant_matrix <- column_to_rownames(variant_matrix, var = colnames(variant_matrix[,1]))
  #Update column names
  colnames(variant_matrix) <- gsub('_R1.fastq.gz','',colnames(variant_matrix))
  return(variant_matrix)
}

remove_rows_with_no_variants_in_subset <- function(parsed){
  bin_mat <- parsed$bin$mat %>% subset(.,rowSums(.)>0) %>% t %>% as.data.frame
  annots <- parsed$bin$annots %>% subset(.,raw_rownames %in% colnames(bin_mat)) %>% as.data.frame
  list(bin = bin_mat,annots = annots)
}
remove_rows_with_no_variants_in_subset_unparsed <- function(matrix,annots){
  bin_mat <- matrix %>% subset(.,rowSums(.)>0) %>% t %>% as.data.frame
  annots <- annots %>% subset(.,raw_rownames %in% colnames(bin_mat)) %>% as.data.frame
  list(bin = bin_mat,annots = annots)
}

clean_snpkit_output_SNP <- function(code,allele,tr){
  ## Change colnames so that subset it to the isolates in ST258
  SNP_code_clean <- clean_snpkit_datasets(code) %>% .[,match(tr$tip.label,colnames(.))]
  SNP_allele_clean <- clean_snpkit_datasets(allele) %>% .[,match(tr$tip.label,colnames(.))]
  ## Parse annotations and keep variants with code 0,1, and -1 in code matrix
  SNP_parsed <- parse_snps(SNP_code_clean,SNP_allele_clean,tree=NULL, og=NULL, remove_multi_annots = TRUE,ref_to_anc=FALSE,keep_conf_only = T)
  ## Remove variants that are not present in an isolate, subset annotation dataset, and format binary matrix for easy manipulation as heatmap or dataset
  SNP_parsed_post_processing <- remove_rows_with_no_variants_in_subset(SNP_parsed)
  return(SNP_parsed_post_processing)
}

clean_snpkit_output_INDEL <- function(code,allele,tr){
  # Change colnames so that subset it to the isolates in ST258
  Indel_code <- clean_snpkit_datasets(code) %>% .[,match(tr$tip.label,colnames(.))]
  Indel_allele <- clean_snpkit_datasets(allele) %>% .[,match(tr$tip.label,colnames(.))]
  # Parse annotations and keep variants with code 0,1, and -1 in code matrix
  Indel_parsed <- parse_indels(Indel_code,Indel_allele,tree=NULL, og=NULL, remove_multi_annots = TRUE, ref_to_anc=FALSE, keep_conf_only = T)
  # Remove variants that are not present in an isolate, subset annotation dataset, and format binary matrix for easy manipulation as heatmap or dataset
  Indel_parsed_post_processing <- remove_rows_with_no_variants_in_subset(Indel_parsed)
  return(Indel_parsed_post_processing)
}

# Functions to parse panISa output
# panISa on GitHub: https://github.com/bvalot/panISa

# Load libraries
library(genbankr)
library(IRanges)
library(dplyr)

# Function to parse is results into a matrix. Takes the output of the panisa isfinder script as input.
# Input: ismat: panISa isfinder script output text file
# Output: matrix where rows are IS element named by start_end position and columns are sample. 1 indicates presence, 0 indicates absence
make_ismat = function(isfinder){
  if(is.character(isfinder)){
    isfinder = read.delim(isfinder)
  }
  isfinder = isfinder[!grepl('Sample',isfinder$Sample),]
  uniq_samp = unique(isfinder$Sample)
  pos = paste0(isfinder$Start_Position,'_',isfinder$Stop_Position)
  uniq_pos = unique(pos)
  ismat = matrix(0,nrow=length(uniq_pos),ncol=length(uniq_samp))
  rownames(ismat) = uniq_pos
  colnames(ismat) = uniq_samp
  for(i in 1:nrow(isfinder)){
    ismat[pos[i],isfinder$Sample[i]] = 1
  }
  colnames(ismat) = gsub('_panISa','',colnames(ismat))
  return(ismat)
}

# Function to annotate row names of the IS matrix (made with parse_is) to include locus_tag information
# Input: ismat: insertion matrix from parse_is
#        gb: genbank file of a reference genome
# Output: matrix where row names contain locus_tag information.
#         & indicates the IS element spans more than 1 gene
#         | indicates that the IS element is intergenic
annotate_is = function(ismat,gb,accession=F){
  if(is.character(gb)){
    if(accession){
      gb = GBAccession(gb)
    }
    print('Reading in genbank file')
    gb = readGenBank(gb)
  }
  print('Getting gene info')
  gene_info = genes(gb)
  locus_tag = gene_info$locus_tag
  ranges = ranges(gene_info)
  gene_strand = as.character(strand(gene_info))
  is_locus_tags = sapply(rownames(ismat), function(x){
    is_pos = strsplit(x,'_')[[1]]
    is_start = as.numeric(is_pos[1])
    is_end = as.numeric(is_pos[2])
    is_lt_start = locus_tag[is_start >= start(ranges) & is_start <= end(ranges)]
    is_lt_end = locus_tag[is_end >= start(ranges) &  is_end <= end(ranges)]
    ig = sum(length(is_lt_start),length(is_lt_end)) == 0
    if(ig){
      ind = which(is_start >= start(ranges) & is_start <= lead(end(ranges)) | is_end >= start(ranges) & is_end <= lead(end(ranges)))
      is_lt_igs = locus_tag[c(ind,ind+1)]
      is_lt_ig = paste(is_lt_igs,collapse = '-')
    }else{
      is_lt_igs = NULL
      is_lt_ig = NULL
    }
    strand = gene_strand[locus_tag %in% c(is_lt_start,is_lt_igs,is_lt_end)]
    if(ig){
      paste0(is_lt_ig,'|',paste(gene_strand[c(ind,ind+1)],collapse='.')) # intergenic
    }else{
      paste0(paste(unique(c(is_lt_start,is_lt_end)),collapse='&'),'|',paste(strand,collapse='&')) # & means in more than 1 gene
    }
  })
  rownames(ismat) = paste0(rownames(ismat),'|',is_locus_tags)
  return(ismat)
}

parse_is = function(ismat,gb,accession=F){

  # MAKE ISFNDER OUTPUT INTO MATRIX
  print('Making IS mat')
  ismat = make_ismat(ismat)
  # ANNOTATE WITH START AND END
  print('Annotating IS mat')
  ismat = annotate_is(ismat,gb,accession=accession)

  # SEPARATE ANNOTATIONS
  print('Getting annotations')
  all_info = strsplit(rownames(ismat),'\\|')

  # GET START AND END
  start_end = sapply(all_info, function(x) x[[1]])
  start = gsub('_.*','',start_end)
  end = gsub('.*_','',start_end)

  # GET LOCUS TAGS
  locus_tag = sapply(all_info, function(x) x[[2]])
  # INTERGENIC
  locus_tag_ig_gene1 = sapply(strsplit(locus_tag,'-|&'),function(x) x[1])
  locus_tag_ig_gene2 = sapply(strsplit(locus_tag,'-'),function(x) x[2])
  # SPANS TWO GENES
  locus_tag_gene1 = sapply(strsplit(locus_tag,'-|&'),function(x) x[1])
  locus_tag_gene2 = sapply(strsplit(locus_tag,'&'),function(x) x[2])

  # INTERGENIC BOOLEAN
  intergenic = !is.na(locus_tag_ig_gene2)
  # SPANS TWO GENES BOOLEAN
  two_genes = !is.na(locus_tag_gene2)

  # GET STRAND
  strand = sapply(all_info, function(x) x[[3]])
  # INTERGENIC
  strand_ig_gene1 = sapply(strsplit(strand,'\\.|&'),function(x) x[1])
  strand_ig_gene2 = sapply(strsplit(strand,'\\.'),function(x) x[2])
  # SPANS TWO GENES
  strand_gene1 = sapply(strsplit(strand,'\\.|&'),function(x) x[1])
  strand_gene2 = sapply(strsplit(strand,'&'),function(x) x[2])

  # INTERGENIC UPSTREAM
  upstream_ig_gene1 = intergenic & strand_ig_gene1 == '-'
  upstream_ig_gene2 = intergenic & strand_ig_gene2 == '+'

  # high impact if intragenic or upstream of a gene
  high_impact = !intergenic | upstream_ig_gene1 | upstream_ig_gene2

  parsed = list(mat = ismat,
                start = start,
                end = end,
                locus_tag = locus_tag,
                locus_tag_ig_gene1 = locus_tag_ig_gene1,
                locus_tag_ig_gene2 = locus_tag_ig_gene2,
                locus_tag_gene1 = locus_tag_gene1,
                locus_tag_gene2 = locus_tag_gene2,
                intergenic = intergenic,
                two_genes = two_genes,
                strand = strand,
                strand_ig_gene1 = strand_ig_gene1,
                strand_ig_gene2 = strand_ig_gene2,
                strand_gene1 = strand_gene1,
                strand_gene2 = strand_gene2,
                upstream_ig_gene1 = upstream_ig_gene1,
                usptream_ig_gene2 = upstream_ig_gene2,
                high_impact = high_impact
  )

  return(parsed)
}

return_grouped_matrix <- function(mat_list,annots_list,tr){
  merged_matrix <- bind_cols(lapply(mat_list,get)) %>% t %>% as.matrix()
  merged_annots <-  bind_rows(lapply(annots_list,get))
  merged_gene_list <- select(merged_annots,locus_tag) %>% as.vector() %>% unlist %>%  as.character()
  merged_grouped_matrix <- prewas::collapse_snps_into_genes(merged_matrix,merged_gene_list) %>% t %>% as.data.frame
  return(merged_grouped_matrix)
}
