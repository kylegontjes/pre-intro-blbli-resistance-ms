setwd("/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/data/SIFT/genomic_data/gtf")
library(tidyverse)

Sys.info()
sessionInfo()

#read.gtf function
read.gtf <- function(gtf.file, trackline = FALSE, sep = "\t"){
    # read gtf file content
  gtf.input <- readr::read_delim(gtf.file,delim = sep, col_names = FALSE, comment = "#" )
    if (!trackline){
    gtf.input <- readr::read_delim(gtf.file,delim = sep, col_names = FALSE, comment = "#" )
  } else {
    gtf.input <- readr::read_delim(gtf.file,delim = sep, col_names = FALSE, comment = "#", skip = 1)
  }
    # name standardized columns
  gffNames <- c("seqname","source","feature",
                "start","end","score","strand",
                "frame","attribute")
  names(gtf.input)[1:ncol(gtf.input)] <- gffNames[1:ncol(gtf.input)]
    if (ncol(gtf.input) > 9)
    stop ("The gff file format can not store more than 9 columns!")
    return(gtf.input)
}

kpnih1_gtf <- read.gtf('./GCF_000281535.2_ASM28153v2_genomic.gtf') %>% subset(seqname == "NZ_CP008827.1")
kpnih1_gtf$seqname <- gsub("NZ_","",kpnih1_gtf$seqname)
#add transcript and exon
kpnih1_gtf_transcript <- kpnih1_gtf %>% subset(feature == "CDS") %>% mutate(feature = "transcript") 
kpnih1_gtf_exon <- kpnih1_gtf %>% subset(feature == "CDS") %>% mutate(feature = "exon")
#merge and sort
kpnih1_gtf_new_exon_transcript <- rbind(kpnih1_gtf,kpnih1_gtf_transcript,kpnih1_gtf_exon) %>% arrange(start)

#write.gtf function
write.gtf <- function(gtf.data, file.name = "file.gtf"){
  write.table(gtf.data,file.name,col.names = FALSE, row.names = FALSE, quote = FALSE,sep = "\t")
}

write.gtf(kpnih1_gtf_new_exon_transcript,"./kpnih1.gtf") 