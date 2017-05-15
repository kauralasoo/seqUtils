calculateMeanProportions <- function(proportions_df, se){
  
  #Extract metadata from SE
  gene_meta = SummarizedExperiment::rowData(se) %>% tbl_df2() %>%
    dplyr::select(transcript_id, gene_id)
  sample_meta = SummarizedExperiment::colData(se) %>% tbl_df2() %>%
    dplyr::select(sample_id, condition_name) %>% as.data.frame()
  
  #Remove gene_id and feature_id columns
  rownames(proportions_df) = proportions_df$feature_id
  proportions_df = proportions_df[,-c(1,2)]
  
  #Calculate mean proportions
  mean_proportions = calculateMean(proportions_df, sample_meta, factor = "condition_name")
  mean_prop_df = dplyr::mutate(mean_proportions, transcript_id = rownames(mean_proportions)) %>% tbl_df() %>%
    dplyr::left_join(gene_meta, by = "transcript_id") %>%
    dplyr::select(gene_id, transcript_id, everything())
  return(mean_prop_df)
}