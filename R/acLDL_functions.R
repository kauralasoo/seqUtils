constructDesignMatrix_acLDL <- function(sample_ids){
  #Construct a design matrix from a list of sample ids
  
  design = data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
  design = design %>% 
    tidyr::separate(sample_id, into = c("donor", "time_point", "condition"), remove = FALSE) %>% 
    dplyr::mutate(donor = tolower(donor))
  rownames(design)= design$sample_id
  return(design)
}