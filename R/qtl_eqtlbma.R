#' Save a list of matrices into a suitable format for eqtlbma
#'
#' Works with expression, genotypes and covariates. Also creates a text file with absolute paths to each individual file.
#' 
#' @param data_list list of matrices
#' @param output_dir relative path to the output dir
#' @param file_suffix suffix added to each file after their name in the list.
#' @param project_root Absolute path to the prohect root directory.
#' @return None
#' @author Kaur Alasoo
#' @export 
saveEqtlbmaMatrices <- function(data_list, output_dir, file_suffix = "exression", project_root){
  #Save data for eqtlbma to disk
  
  #Construct file list
  file_list = data_frame(sample_name = names(data_list), 
                         file_path = file.path(project_root, output_dir, 
                                               paste(names(data_list),file_suffix, "txt.gz", sep = ".")))
  file_list_path = file.path(output_dir, 
                             paste("list",file_suffix, "txt", sep = "."))
  write.table(file_list, file_list_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  #Save each matrix as a separate gzipped txt file
  print(file_list)
  for (sn in file_list$sample_name){
    file_path = file.path(output_dir, paste(sn,file_suffix, "txt.gz", sep = "."))
    print(file_path)
    file = gzfile(file_path, "w")
    write.table(data_list[[sn]], file, sep = "\t", quote = FALSE)
    close(file)
  }
}

#' Save eqtl_data_list into files suitable for eqtlbma
#'
#' Helper function to export processed eQTL data into multiple text files required by eqtlbma.
#' 
#' @param eqtl_data_list List of all data required for eQTL mapping with eqtlbma.
#' @param output_dir relative path to the output dir
#' @param project_root Absolute path to the prohect root directory.
#' @return None
#' @author Kaur Alasoo
#' @export 
saveEqtlbmaData <- function(eqtl_data_list, output_dir, project_root){
  
  #Save gene coords into a bed file
  genepos_path = file.path(output_dir, "gene_coords.bed")
  write.table(eqtl_data_list$genepos, genepos_path, sep ="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #Save SNP coords to disk
  snpspos_path = file.path(output_dir, "snp_coords.bed")
  snpspos_eqtlbma = dplyr::transmute(eqtl_data_list$snpspos, chr, start  = pos-1, end = pos, snpid)
  write.table(snpspos_eqtlbma, snpspos_path, sep ="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #Save expression data
  saveEqtlbmaMatrices(eqtl_data_list$exprs_cqn_list, output_dir, "expression", project_root)
  
  #Save covariates
  saveEqtlbmaMatrices(eqtl_data_list$covariates_list, output_dir, "covariates", project_root)
  
  #Save genotypes
  genotypes_path = file.path(output_dir, "genotypes.txt")
  write.table(eqtl_data_list$genotypes, genotypes_path, sep ="\t", quote = FALSE)
  
  #Save genotypes list
  genotypes_list_path = file.path(output_dir, "list.genotypes.txt")
  genotypes_list = data_frame(sample_name = names(eqtl_data_list$exprs_cqn_list), 
                         file_path = file.path(project_root, output_dir, "genotypes.txt"))
  write.table(genotypes_list, genotypes_list_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}


