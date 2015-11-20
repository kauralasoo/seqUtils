#' Import EBI GWAS catalog into R.
#'
#' Imports tsv version of the GWAS catalog into R. Parses total sample size from the 
#' INITIAL SAMPLE DESCRIPTION field and also flags each study that contains European samples.
#' 
#' @param path Path to the EBI GWAS catalog export file (must contain MAPPED_TRAIT field as well).
#' @return data frame containing the gwas catalog.
#' @author Kaur Alasoo
#' @export 
importGwasCatalog <- function(path){
  
  #Load gwas catlog from disk
  gwas_catalog = readr::read_delim(path, delim = "\t")
  
  #Construct a data frame of studies
  studies = gwas_catalog[,c("PUBMEDID","INITIAL SAMPLE DESCRIPTION")] %>% unique()
  colnames(studies) = c("pubmed_id", "description")
  
  #Extract sample size from the gwas catalog
  sample_sizes = lapply(as.list(studies$description), function(x){
    stringr::str_extract_all(x, "(\\d+,\\d+)|(\\d+)") %>%
      unlist() %>%
      stringr::str_replace(",","") %>%
      as.numeric() %>%
      sum()
  })
  
  #Add sample size and European status to study df
  studies_df = dplyr::mutate(studies, sample_size = unlist(sample_sizes)) %>%
    dplyr::mutate(is_european = grepl("european",description,ignore.case = TRUE)) %>%
    dplyr::select(pubmed_id, sample_size, is_european)
  
  #Extract snps form catalog
  gwas_columns = gwas_catalog[,c("PUBMEDID","CHR_ID","CHR_POS","SNPS","DISEASE/TRAIT","MAPPED_TRAIT")]
  colnames(gwas_columns) = c("pubmed_id", "chr", "pos", "snp_id", "trait", "mapped_trait")
  
  #Join study data with SNPS and remove NAs
  gwas_table = dplyr::left_join(studies_df, gwas_columns, by = "pubmed_id") %>%
    dplyr::filter(!is.na(chr), !is.na(pos))
  
  return(gwas_table)
}
