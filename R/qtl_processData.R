#' Save list of expression matrices in a format suitable for runPEER.py
#'
#' @param expression_list list of expression matrices
#' @param outdir output directory
#' @param file_suffix suffix of the output file
#'
#' @export
#'
savePEERData <- function(expression_list, outdir, file_suffix = "exprs"){
  #Save eQTL gene expression data to be used with PEER
  sample_names = names(expression_list)
  for (sample in sample_names){
    matrix = expression_list[[sample]]
    path = file.path(outdir, paste(sample, file_suffix ,"txt", sep = "."))
    print(path)
    write.table(t(matrix), path, row.names = FALSE, col.names = FALSE, sep = ",")
  }
}

importPEERFactors <- function(file_path, design_matrix, remove_mean = TRUE){
  #Import output from runPEER.py script back into R
  peer_factors = read.table(file_path, sep =",")
  if(remove_mean){ #Remove the first column that contains the mean
    peer_factors = peer_factors[,2:ncol(peer_factors)]
  }
  colnames(peer_factors) = paste("PEER_factor_", c(1:ncol(peer_factors)), sep = "")
  peer_factors = dplyr::mutate(peer_factors, sample_id = design_matrix$sample_id) %>%
    dplyr::select(sample_id, everything())
  return(peer_factors)
}
