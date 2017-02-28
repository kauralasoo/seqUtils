#' Construct gene position data.frame suitable for QTLtools
#' 
#' @param gene_metadata Gene metadata data frame from expression list (columns: chr, start, end, gene_id, transcript_id, strand)
#' @return data_frame suitable for FastQTL
#' @author Kaur Alasoo
#' @export 
constructQTLtoolsGenePos <- function(gene_metadata){
  genepos = dplyr::select(gene_metadata, chr, start, end, transcript_id, gene_id, strand) %>% 
    dplyr::arrange(chr, start) %>%
    dplyr::rename_("left" = "start", "right" = "end", "#chr" = "chr") %>%
    dplyr::mutate(strand = ifelse(strand == 1, "+", "-"))
  return(genepos)
}


#' Add gene position to the QTLtools expression matrix
#' 
#' @param matrix expression matrix
#' @param genepos data frame form constructQTLtoolsGenePos
#' @return data_frame suitable for QTLtools
#' @author Kaur Alasoo
#' @export 
prepareQTLtoolsMatrix <- function(matrix, genepos){
  res = dplyr::mutate(as.data.frame(matrix), transcript_id = rownames(matrix)) %>%
    dplyr::select(transcript_id, everything()) %>%
    dplyr::left_join(genepos, ., by = "transcript_id") %>%
    dplyr::arrange()
}

#' Import QTLtools output table from permutation run into R.
#'
#' 
#' @param file_path Path to the QTLtools output file.
#' @return data_frame containing gene_ids, snp ids and p-values.
#' @author Kaur Alasoo
#' @export 
importQTLtoolsTable <- function(file_path){
  col_names = c("group_id", "pheno_chr", "pheno_start", "pheno_end", "strand", "phenotype_id", "group_size", "n_cis_snps", 
                "distance", "snp_id", "snp_chr", "snp_start", "snp_end", "df", "dummy", "beta1", 
                "beta2", "p_nominal","slope","p_perm","p_beta")
  col_types = "cciicciiicciiiddddddd"
  table = readr::read_delim(file_path, col_names = col_names, delim = " ", col_types = col_types) %>%
    dplyr::filter(!is.na(p_beta)) %>%
    dplyr::mutate(p_bonferroni = p_nominal*group_size*n_cis_snps) %>%
    dplyr::mutate(p_bonferroni = pmin(p_bonferroni,1)) %>%
    dplyr::mutate(p_fdr = p.adjust(p_beta, method = "fdr")) %>%
    dplyr::mutate(qvalue = qvalue::qvalue(p_beta)$qvalues) %>%
    dplyr::arrange(p_fdr)
  return(table)
}

#' Fetch particular genes from tabix indexed QTLtools output file.
#'
#' @param gene_ranges GRanges object with coordinates of the cis regions around genes.
#' @param tabix_file Tabix-indexed fastqtl output file.
#'
#' @return List of data frames containing QTLtools results for each gene.
#' @export
qtltoolsTabixFetchPhenotypes <- function(phenotype_ranges, tabix_file){
  
  #Assertions about input
  assertthat::assert_that(class(phenotype_ranges) == "GRanges")
  assertthat::assert_that(assertthat::has_name(GenomicRanges::elementMetadata(phenotype_ranges), "phenotype_id"))
  
  #Set column names for rasqual
  fastqtl_columns = c("phenotype_id","pheno_chr","pheno_start", "pheno_end",
                      "strand","n_snps", "distance", "snp_id", "snp_chr",
                      "snp_start", "snp_end", "p_nominal","beta", "is_lead")
  fastqtl_coltypes = "cciiciicciiddi"
  
  result = list()
  for (i in seq_along(phenotype_ranges)){
    selected_phenotype_id = phenotype_ranges[i]$phenotype_id
    print(i)
    tabix_table = scanTabixDataFrame(tabix_file, phenotype_ranges[i], 
                                     col_names = fastqtl_columns, col_types = fastqtl_coltypes)[[1]] %>%
      dplyr::filter(phenotype_id == selected_phenotype_id)
    
    #Add additional columns
    result[[selected_phenotype_id]] = tabix_table
  }
  return(result)
}


