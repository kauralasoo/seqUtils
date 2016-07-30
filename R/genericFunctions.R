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

#' Force a vector of values into standard normal distribution
#'
#' @param x numeric vector with arbitrary distribution
#'
#' @return Vector with a standard normal distribution
#' @export
quantileNormaliseVector = function(x){
  qnorm(rank(x,ties.method = "random")/(length(x)+1))
}


quantileNormaliseMatrix <- function(matrix){
  quantile_matrix = matrix(0, nrow(matrix), ncol(matrix))
  for (i in seq_along(matrix[1,])){
    quantile_matrix[,i] = quantileNormaliseVector(matrix[,i])
  }
  #Add names
  rownames(quantile_matrix) = rownames(matrix)
  colnames(quantile_matrix) = colnames(matrix)
  return(quantile_matrix)
}

