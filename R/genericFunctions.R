extractSubset <- function(sample_meta, data, old_column_names = "sample_id", new_column_names = "donor"){
  #Extract subset of columns from a data.frame and rename column names
  
  sample_meta = as.data.frame(sample_meta)               #Ensure that its a data.frame
  subset_data = data[,sample_meta[,old_column_names]]    #Keep only sample that are in the metadata column
  colnames(subset_data) = sample_meta[,new_column_names] #Rename the columns with donor id
  
  return(subset_data)
}

idVectorToList <- function(id_vector){
  #Convert a vector of IDs into a list of ids where the elements have the same name
  id_list = as.list(id_vector)
  names(id_list) = id_vector
  return(id_list)
}

#' Convert named vector into a tidy data_frame.
#' 
#' @param named_vector Named vector.
#' @return data_frame with two columns: value, sample_id.
#' @author Kaur Alasoo
#' @export 
tidyVector <- function(named_vector, value_id = "value", sample_id = "sample_id"){
  res = dplyr::data_frame(sample = names(named_vector), value = named_vector)
  colnames(res) = c(sample_id, value_id)
  return(res)
}

listUnion <- function(granges_list){
  #Calculated the union of a GRangesList object
  union_obj = granges_list[[1]]
  if(length(granges_list) > 1){
    for(i in 2:length(granges_list)){
      union_obj = GenomicRanges::union(union_obj, granges_list[[i]]) 
    } 
  }
  return(union_obj)
}

#' Split vector of n elements into batches
#' 
#' Given a total number of elements n and batch_size, contruct a vector of length n where
# each element occurs at most batch_size times.
#'
#' @param n Length of the vector
#' @param batch_size Size of each batch
#'
#' @return vector of integers where each integer occurs at most batch_size times.
#' @export
#'
#' @examples
#' splitIntoBatches(11,3)
splitIntoBatches <- function(n, batch_size){
  n_batches = ceiling(n/batch_size)
  batch_ids = rep(seq(1:n_batches), each = batch_size)[1:n]
  return(batch_ids)
}

#Split ids into batches
constructIdBatches <- function(batch_string, id_vector){
  
  #Extract batch number in n_batches from batch_string
  b_vector = as.numeric(unlist(strsplit(batch_string, " ")))
  batch_id = b_vector[1]
  n_batches = b_vector[2]
  
  #Construct batches
  batch_size = ceiling(length(id_vector)/n_batches)
  batches = splitIntoBatches(length(id_vector), batch_size)
  
  #Extrat ids belonging to batch
  selection = batches == batch_id
  selected_ids = id_vector[selection]
  
  return(selected_ids)
}

matrixExtractPairs <- function(row_name, col_name, matrix){
  return(matrix[row_name, col_name])
}

#' Colour plaette of four colours to represent experimental conditions
conditionPalette <- function(){
  c("#67a9cf","#2166ac","#ef8a62","#b2182b")
}

friendlyNames <- function(){
  friendly_names = data.frame(condition_name = c("naive","IFNg","SL1344", "IFNg_SL1344"), 
                              friendly_name = c("Naive","IFNg","Salmonella","Both")) %>%
    dplyr::mutate(friendly_name = factor(friendly_name, levels = friendly_name))
}

addMafFromVariantInfo <- function(variants_df, variant_info){
  filtered_variants = dplyr::filter(variant_info, snp_id %in% variants_df$snp_id) %>% 
    dplyr::select(snp_id, MAF)
  result = dplyr::left_join(variants_df, filtered_variants, by = "snp_id")
  return(result)
}

#' Convert a data frame with an id column into a matrix with row names
tibbleToNamedMatrix <- function(tibble, row_names = "transcript_id"){
  assertthat::assert_that(assertthat::has_name(tibble, row_names))
  
  rows = as.vector(unlist(tibble[,row_names]))
  selected_columns = which(!colnames(tibble) == row_names)
  matrix = dplyr::select(tibble, selected_columns) %>% as.matrix()
  rownames(matrix) = rows
  return(matrix)
}

#' Convert bioconductor DataFrame object into tibble
tbl_df2 <- function(dataframe){
  return(tibble::as_tibble(BiocGenerics::as.data.frame(dataframe)))
}

#' Construct matching between minor allele count and genotype text
constructGenotypeText <- function(snp_id, variant_information){
  
  #Extraxt SNP entry
  var_info = dplyr::filter(variant_information, snp_id == "rs7594476")
  
  #Make an new tibble
  df = dplyr::data_frame(genotype_value = c("0","1","2"), 
                         genotype_text = c(paste0(var_info$ref, var_info$ref), 
                                           paste0(var_info$ref, var_info$alt), 
                                           paste0(var_info$alt, var_info$alt))) %>%
    dplyr::mutate(genotype_text = factor(genotype_text, levels = genotype_text))
  return(df)
}
