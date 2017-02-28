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
    table = table %>% tbl_df() %>% 
      dplyr::filter(!is.na(p_beta)) %>%
      dplyr::mutate(p_fdr = p.adjust(p_beta, method = "fdr")) %>%
      dplyr::mutate(qvalue = qvalue::qvalue(p_beta)$qvalues) %>%
      dplyr::arrange(p_fdr)
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
	  dplyr::arrange(chr, start) %>%
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
    dplyr::left_join(genepos, ., by = "geneid") %>%
    dplyr::arrange()
}

fastqtlMetadataToCovariates <- function(metadata){
  genotype_ids = metadata$genotype_id
  genotype_position = which(colnames(metadata) == "genotype_id")
  covariate_ids = colnames(metadata)[-genotype_position]
  cov_matrix = t(metadata[-genotype_position]) %>% as.data.frame()
  colnames(cov_matrix) = genotype_ids
  cov_matrix = dplyr::mutate(cov_matrix, id = covariate_ids) %>% dplyr::select(id, everything())
  return(cov_matrix)
}

fastQTLCorrectEigenMT <- function(fastqtl_results, n_tests){
  res = dplyr::left_join(fastqtl_results, n_tests, by = "gene_id") %>% 
    dplyr::mutate(n_tests = ifelse(is.na(n_tests), n_cis_snps, n_tests)) %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::mutate(p_eigen = p.adjust(p_nominal, "bonferroni", n = n_tests)) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(p_eigen_fdr = p.adjust(p_eigen, method = "fdr"))
  return(res)
}


#' Fetch particular genes from tabix indexed FastQTL output file.
#'
#' @param phenotype_ranges GRanges object with coordinates of the cis regions around genes.
#' @param tabix_file Tabix-indexed fastqtl output file.
#'
#' @return List of data frames containing Rasqual results for each gene.
#' @export
fastqtlTabixFetchGenes <- function(phenotype_ranges, tabix_file){
  
  #Assertions about input
  assertthat::assert_that(class(phenotype_ranges) == "GRanges")
  assertthat::assert_that(assertthat::has_name(GenomicRanges::elementMetadata(phenotype_ranges), "phenotype_id"))
  
  #Set column names for rasqual
  fastqtl_columns = c("phenotype_id","chr","pos","snp_id","distance","p_nominal","beta")
  fastqtl_coltypes = "ccicidd"
  
  result = list()
  for (i in seq_along(phenotype_ranges)){
    selected_phenotype_id = phenotype_ranges[i]$phenotype_id
    print(i)
    tabix_table = scanTabixDataFrame(tabix_file, phenotype_ranges[i], col_names = fastqtl_columns, col_types = fastqtl_coltypes)[[1]] %>%
      dplyr::filter(phenotype_id == selected_phenotype_id)
    
    #Add additional columns
    result[[selected_phenotype_id]] = tabix_table
  }
  return(result)
}

fastqtlTabixFetchGenesQuick <- function(gene_ids, tabix_file, gene_metadata, cis_window = 5e5){
  gene_df = dplyr::data_frame(gene_id = gene_ids)
  
  #If gene_metadata is already a GRanges object then just filter it based on gene_ids
  if(class(gene_metadata) == "GRanges"){
    gene_ranges = gene_metadata[gene_metadata$gene_id %in% gene_ids]
  } else { #Otherwise construct a gene ranges object
    gene_ranges = rasqualTools::constructGeneRanges(gene_df, gene_metadata, cis_window = cis_window)
  }
  tabix_data = fastqtlTabixFetchGenes(gene_ranges, tabix_file)
  return(tabix_data)
}
