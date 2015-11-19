importGwasCatalog <- function(path){
  
  #Load gwas catlog from disk
  gwas_catalog = readr::read_delim(path, delim = "\t")
  
  #Construct a data frame of studies
  studies = gwas_catalog[,c("PUBMEDID","INITIAL SAMPLE DESCRIPTION")] %>% unique()
  colnames(studies) = c("pubmed_id", "description")
  
  #Extract sample size from the gwas catalog
  sample_sizes = lapply(as.list(studies$description), function(x){
    str_extract_all(x, "(\\d+,\\d+)|(\\d+)") %>%
      unlist() %>%
      str_replace(",","") %>%
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
  
  #Filter by ethinicty and 
  gwas_table = dplyr::left_join(studies_df, gwas_columns, by = "pubmed_id") %>%
    dplyr::filter(!is.na(chr), !is.na(pos))
  
  return(gwas_table)
}
