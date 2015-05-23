constructDesignMatrix_SL1344 <- function(sample_ids){
  #Construct a design matrix from a list of sample ids
  
  design = data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
  condition_names = data.frame(condition = c("A","B","C","D"), 
                               condition_name = c("naive", "IFNg", "SL1344", "IFNg+SL1344"), 
                               stringsAsFactors = FALSE)
  design = design %>% 
    tidyr::separate(sample_id, into = c("donor", "condition", "replicate"), extra = "drop", remove = FALSE) %>%
    dplyr::mutate(replicate = ifelse(is.na(replicate), 1, replicate)) %>% #if no replicate in file name then set it to 1
    dplyr::mutate(replicate = as.numeric(replicate)) %>%
    dplyr::left_join(condition_names, by = "condition") %>%
    dplyr::mutate(SL1344 = ifelse(condition %in% c("C","D"), "infected", "control")) %>%
    dplyr::mutate(IFNg = ifelse(condition %in% c("B","D"), "primed", "naive"))
  rownames(design)= design$sample_id
  return(design)
}