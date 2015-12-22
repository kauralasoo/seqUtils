#' Save a list of matrices into a suitable format for FastQTL
#'
#' Works with expression and covariates matrices.
#' 
#' @param data_list list of matrices
#' @param output_dir relative path to the output dir
#' @param file_suffix suffix added to each file after their name in the list.
#' @return None
#' @author Kaur Alasoo
#' @export 
saveFastqtlMatrices <- function(data_list, output_dir, file_suffix = "expression", col_names = TRUE){
  #Save data for FastQTL to disk

  #Save each matrix as a separate  txt file
  for (sn in names(data_list)){
    file_path = file.path(output_dir, paste(sn,file_suffix, "txt", sep = "."))
    print(file_path)
    write.table(data_list[[sn]], file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = col_names)
  }
}

#' Import fastQTL output table into R.
#'
#' Detect if the table is from nominal run or permutation run and add proper column names.
#' 
#' @param file_path Path to the fastQTL output file
#' @return data_frame containing gene_ids, snp ids and p-values.
#' @author Kaur Alasoo
#' @export 
importFastQTLTable <- function(file_path){
  table = read.table(file_path, stringsAsFactors = FALSE)
  if(ncol(table) == 11){
    colnames(table) = c("gene_id", "n_cis_snps", "beta1", "beta2", "dummy", "snp_id", "distance","p_nominal", "slope","p_perm","p_beta")
  }
  return(table)
}

#' Construct gene position data.frame suitable for FastQTL
#' 
#' @param gene_metadata Gene metadata data frame from expression list (columns: chr, start, end, gene_id)
#' @return data_frame suitable for FastQTL
#' @author Kaur Alasoo
#' @export 
constructFastQTLGenePos <- function(gene_metadata){
	genepos = dplyr::select(gene_metadata, chr, start, end, gene_id) %>% 
    dplyr::rename_("left" = "start", "right" = "end", "geneid" = "gene_id", "#chr" = "chr")
  return(genepos)
}

#' Add gene position to the FastQTL expression matrix
#' 
#' @param matrix expression matrix
#' @param genepos data frame form constructFastQTLGenePos
#' @return data_frame suitable for FastQTL
#' @author Kaur Alasoo
#' @export 
prepareFastqtlMatrix <- function(matrix, genepos){
  res = dplyr::mutate(as.data.frame(matrix), geneid = rownames(matrix)) %>%
    dplyr::select(geneid, everything()) %>%
    dplyr::left_join(genepos, ., by = "geneid")
}
