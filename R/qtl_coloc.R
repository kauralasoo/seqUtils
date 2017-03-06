#' Wrapper around coloc.abf to directly work with tidy data frames
#'
#' Group QTL database by gene, sort by p_nominal, keep SNP with smallest p-value,
#' correct that using bonferroni correction and then apply FDR correction across genes.
#'
#' @param df1 Summary statistics from the first dataset
#' (required columns: snp_id, p_nominal, beta, MAF).
#' @param df1 Summary statistics from the second dataset
#' (required columns: snp_id, p_nominal, beta, MAF).
#' @param n1 Sample size of the first dataset.
#' @param n2 Sample size of the second dataset.
#' @param p1 prior probability a SNP is associated with trait 1, default 1e-4
#' @param p2 prior probability a SNP is associated with trait 2, default 1e-4
#' @param p12 prior probability a SNP is associated with both traits, default 1e-5
#'
#' @return List of colocalisation results (Same as coloc.abf).
#' @export
testColoc <- function(df1, df2, n1, n2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5){
  
  #Test for colocalisation between two sets of p-values
  result = tryCatch({
    result = coloc::coloc.abf(dataset1 = list(pvalues = df1$p_nominal, N = n1, 
                                            beta = df1$beta, MAF = df1$MAF, 
                                            snp = df1$snp_id, type = "quant"), 
                            dataset2 = list(pvalues = df2$p_nominal, N = n2, 
                                            beta = df2$beta, MAF = df2$MAF,
                                            snp = df2$snp_id, type = "quant"),
                            p1 = p1, p2 = p2, p12 = p12)
  }, error = function(err) {
    print(paste("ERROR:",err))
    result = NULL
  }
  )
  return(result)
}

#' Calculate posterior probabilities from association summary stats
#'
#' @param dataset Association summary statistics
#' (required columns: snp_id, p_nominal, effect_size, MAF).
#' @param n Sample size of the dataset
#'
#' @return Original data frame with one additional column containing posterior probabilities.
#' @export
addAssociationPosterior <- function(dataset, n){
  
  #Use the coloc package to calculate ABFs
  coloc_res = testColoc(dataset, dataset, n, n)
  post_df = dplyr::transmute(coloc_res$results, snp_id = snp, lABF = lABF.df1) %>% 
    dplyr::mutate(posterior = exp(lABF)/sum(exp(lABF))) %>% 
    dplyr::select(snp_id, lABF, posterior) %>%
    dplyr::mutate(snp_id = as.character(snp_id))
  result = dplyr::left_join(dataset, post_df, by = "snp_id")
  return(result)
}


#' Construct the credible set of associated SNPs for each gene.
#'
#' @param dataset Output from addAssociationPosterior
#' @param threshold Cumulative posterior porbability of most associated SNPs.
#'
#' @return Dataset that only contains the credible set of associated SNPs.
#' @export
constructCredibleSet <- function(dataset, threshold = 0.95){
  result = dplyr::arrange(dataset, -posterior) %>% 
    dplyr::mutate(cum_posterior = cumsum(posterior))
  #Identify the number of top associated SNPs whose cumulative posterior is above the threshold
  n_top_rows = nrow(result) - (length(which(result$cum_posterior >= threshold))-1)
  result = dplyr::filter(result, row_number() <= n_top_rows)
  return(result)
}


#' Construct a GRanges obejct corresponding to a cis region around one variant.
#'
#' @param variant_df data frame with at least snp_id column
#' @param variant_information data.frame from importVariantInformation() function
#' @param cis_dist Number of basepairs upstream and downstream of the variant.
#'
#' @return GRanges 
#' @export
constructVariantRanges <- function(variant_df, variant_information, cis_dist){
  
  #Make key assertions
  assertthat::assert_that(assertthat::has_name(variant_df, "snp_id"))
  assertthat::assert_that(assertthat::has_name(variant_information, "snp_id"))
  assertthat::assert_that(assertthat::has_name(variant_information, "chr"))
  assertthat::assert_that(assertthat::has_name(variant_information, "pos"))
  
  #Filter variant information to contain only required snps
  var_info = dplyr::filter(variant_information, snp_id %in% variant_df$snp_id) %>%
    dplyr::select(snp_id, chr, pos, MAF)
  
  #Add variant info to variant df
  var_df = dplyr::left_join(variant_df, var_info, by = "snp_id")
  
  #Make a ranges object
  var_ranges = var_df %>%
    dplyr::rename(seqnames = chr) %>%
    dplyr::mutate(start = pos - cis_dist, end = pos + cis_dist, strand = "*") %>%
    dataFrameToGRanges()
  
  return(var_ranges)
}

summaryReplaceCoordinates <- function(summary_df, variant_information){
  
  #Make key assertions
  assertthat::assert_that(assertthat::has_name(summary_df, "snp_id"))
  assertthat::assert_that(assertthat::has_name(summary_df, "pos"))
  assertthat::assert_that(assertthat::has_name(summary_df, "chr"))
  
  #Filter variant information to contain only required snps
  var_info = dplyr::filter(variant_information, snp_id %in% summary_df$snp_id) %>%
    dplyr::select(snp_id, chr, pos, MAF)
  
  #Remove MAF if it is present
  if(assertthat::has_name(summary_df, "MAF")){
    summary_df = dplyr::select(summary_df, -MAF)
  }
  
  #Add new coordinates
  new_coords = dplyr::select(summary_df, -chr, -pos) %>%
    dplyr::left_join(var_info, by = "snp_id") %>%
    dplyr::filter(!is.na(pos)) %>%
    dplyr::arrange(pos)
  
  return(new_coords)
}

summaryReplaceSnpId <- function(summary_df, variant_information){
  
  #Make key assertions
  assertthat::assert_that(assertthat::has_name(summary_df, "snp_id"))
  assertthat::assert_that(assertthat::has_name(summary_df, "pos"))
  assertthat::assert_that(assertthat::has_name(summary_df, "chr"))
  
  #Filter variant information to contain only required snps
  var_info = dplyr::filter(variant_information, pos %in% summary_df$pos) %>%
    dplyr::select(snp_id, chr, pos, MAF)
  
  #Remove MAF if it is present
  if(assertthat::has_name(summary_df, "MAF")){
    summary_df = dplyr::select(summary_df, -MAF)
  }
  
  #Add new coordinates
  new_coords = dplyr::select(summary_df, -snp_id) %>%
    dplyr::left_join(var_info, by = c("chr","pos")) %>%
    dplyr::filter(!is.na(snp_id)) %>%
    dplyr::arrange(pos)
  
  return(new_coords)
}

#' Test colocalisation between molecular QTL and GWAS summary stats
#'
#' @param qtl QTL summary stats (p_nominal, MAF, beta, snp_id)
#' @param gwas GWAS summary stats(beta, se, MAF, log_OR)
#' @param N_qtl Sample size of the QTL mapping study
#'
#' @return coloc.abf result object
#' @export
colocQtlGWAS <- function(qtl, gwas, N_qtl){
  
  #Check that QTL df has all correct names
  assertthat::assert_that(assertthat::has_name(qtl, "snp_id"))
  assertthat::assert_that(assertthat::has_name(qtl, "beta"))
  assertthat::assert_that(assertthat::has_name(qtl, "MAF"))
  assertthat::assert_that(assertthat::has_name(qtl, "p_nominal"))
  
  #Check that GWAS df has all correct names
  assertthat::assert_that(assertthat::has_name(gwas, "beta"))
  assertthat::assert_that(assertthat::has_name(gwas, "se"))
  assertthat::assert_that(assertthat::has_name(gwas, "snp_id"))
  assertthat::assert_that(assertthat::has_name(gwas, "log_OR"))
  assertthat::assert_that(assertthat::has_name(gwas, "MAF"))
  
  #Count NAs for log_OR and beta
  log_OR_NA_count = length(which(is.na(gwas$log_OR)))
  beta_NA_count = length(which(is.na(gwas$beta)))
  
  #Remove GWAS SNPs with NA std error
  gwas = dplyr::filter(gwas, !is.na(se))
  
  #If beta is not specified then use log_OR
  if(beta_NA_count <= log_OR_NA_count){
    coloc_res = coloc::coloc.abf(dataset1 = list(pvalues = qtl$p_nominal, 
                                                 N = N_qtl, 
                                                 MAF = qtl$MAF, 
                                                 type = "quant", 
                                                 beta = qtl$beta,
                                                 snp = qtl$snp_id),
                                 dataset2 = list(beta = gwas$beta, 
                                                 varbeta = gwas$se^2, 
                                                 type = "cc", 
                                                 snp = gwas$snp_id, 
                                                 MAF = gwas$MAF))
  } else{
    coloc_res = coloc::coloc.abf(dataset1 = list(pvalues = qtl$p_nominal, 
                                                 N = N_qtl, 
                                                 MAF = qtl$MAF, 
                                                 type = "quant", 
                                                 beta = qtl$beta,
                                                 snp = qtl$snp_id),
                                 dataset2 = list(beta = gwas$log_OR, 
                                                 varbeta = gwas$se^2, 
                                                 type = "cc", 
                                                 snp = gwas$snp_id, 
                                                 MAF = gwas$MAF))
  }
  
  return(coloc_res)
}

colocMolecularQTLs <- function(qtl_df, qtl_summary_path, gwas_summary_path, 
                               GRCh37_variants, GRCh38_variants,
                               N_qtl = 84, cis_dist = 1e5, QTLTools = TRUE){
  
  #Assertions
  assertthat::assert_that(assertthat::has_name(qtl_df, "phenotype_id"))
  assertthat::assert_that(assertthat::has_name(qtl_df, "snp_id"))
  assertthat::assert_that(nrow(qtl_df) == 1)
  
  assertthat::assert_that(is.numeric(cis_dist))
  assertthat::assert_that(is.numeric(N_qtl))
  
  
  result = tryCatch({
    #Make GRanges object to fetch data
    qtl_ranges = constructVariantRanges(qtl_df, GRCh38_variants, cis_dist = cis_dist)
    gwas_ranges = constructVariantRanges(qtl_df, GRCh37_variants, cis_dist = cis_dist)

    #Fetch QTL summary stats
    if(QTLTools){
      qtl_summaries = qtltoolsTabixFetchPhenotypes(qtl_ranges, qtl_summary_path)[[1]] %>%
        dplyr::transmute(snp_id, chr = snp_chr, pos = snp_start, p_nominal, beta)
    } else{
      qtl_summaries = fastqtlTabixFetchGenes(qtl_ranges, qtl_summary_path)[[1]]
    }

    #Fetch GWAS summary stats
    gwas_summaries = tabixFetchGWASSummary(gwas_ranges, gwas_summary_path)[[1]]
    
    #Substitute coordinate for the eqtl summary stats and add MAF
    qtl = summaryReplaceCoordinates(qtl_summaries, GRCh37_variants)
    
    #Substitute snp_id for the GWAS summary stats and add MAF
    gwas = summaryReplaceSnpId(gwas_summaries, GRCh37_variants)
    
    #Extract minimal p-values for both traits
    qtl_min = dplyr::arrange(qtl, p_nominal) %>% dplyr::filter(row_number() == 1)
    gwas_min = dplyr::arrange(gwas, p_nominal) %>% dplyr::filter(row_number() == 1)
    
    #Perform coloc analysis
    coloc_res = colocQtlGWAS(qtl, gwas, N_qtl = N_qtl)
    coloc_summary = dplyr::tbl_df(t(data.frame(coloc_res$summary))) %>%
      dplyr::mutate(qtl_pval = qtl_min$p_nominal, gwas_pval = gwas_min$p_nominal,
                    qtl_lead = qtl_min$snp_id, gwas_lead = gwas_min$snp_id) #Add minimal pvalues
    
    #Summary list
    data_list = list(qtl = qtl, gwas = gwas)
    
    result = list(summary = coloc_summary, data = data_list)
  }, error = function(err) {
    print(paste("ERROR:",err))
    result = list(summary = NULL, data = NULL)
  }
  )
  return(result)
}

colocMolecularQTLsByRow <- function(qtl_df, ...){
  result = purrr::by_row(qtl_df, ~colocMolecularQTLs(.,...)$summary, .collate = "rows")
}

#Make coloc plot
makeColocPlot <- function(data_list){
  #Join data together
  trait_df = purrr::map_df(data_list, ~dplyr::select(.,chr, pos, p_nominal), .id = "trait") %>%
    dplyr::mutate(log10p = -log(p_nominal, 10))
  
  #Keep only matching varaints
  trait_df = dplyr::mutate(dplyr::group_by(trait_df, pos), pos_count = length(pos)) %>% 
    dplyr::filter(pos_count == 2) %>% 
    dplyr::ungroup()
  
  #Make a plot
  plot = ggplot(trait_df, aes(x = pos, y = log10p)) + 
    geom_point() + facet_wrap(~trait, ncol = 1, scales = "free_y")
  return(plot)
}


#' Perform a quick pre-filtering between QTLs and GWAS hits to reduce the number of coloc tests
#'
#' @param qtl_min_pvalues List of data frames with QTL lead pvalues. Each data frame must contain
#' gene_id, snp_id and p_fdr and should not contain other columns.
#' @param gwas_prefix Prefix of the GWAS summarystats file
#' @param GRCh37_variants QTL variant information in GWAS coordinates.
#' @param fdr_thresh Minimal QTL FDR threshold
#' @param overlap_dist Max distance between GWAS and QTL variants.
#' @param gwas_thresh Minimal GWAS p-value threshold
#'
#' @return List of data.frames with phenotype_ids and snp_ids to be tested with coloc.
#' @export
prefilterColocCandidates <- function(qtl_min_pvalues, gwas_prefix, GRCh37_variants, 
                                     fdr_thresh = 0.1, overlap_dist = 1e5, gwas_thresh = 1e-5){
  
  #Make sure that the qtl_df has neccessary columns
  assertthat::assert_that(assertthat::has_name(qtl_min_pvalues[[1]], "phenotype_id"))
  assertthat::assert_that(assertthat::has_name(qtl_min_pvalues[[1]], "snp_id"))
  assertthat::assert_that(assertthat::has_name(qtl_min_pvalues[[1]], "p_fdr"))
  assertthat::assert_that(ncol(qtl_min_pvalues[[1]]) == 3)
  
  
  #Import top GWAS p-values
  gwas_pvals = importGWASSummary(paste0(gwas_prefix,".top_hits.txt.gz")) %>%
    dplyr::filter(p_nominal < gwas_thresh) %>%
    dplyr::transmute(chr = chr, gwas_pos = pos)

  #Filter lead variants
  qtl_hits = purrr::map(qtl_min_pvalues, ~dplyr::filter(., p_fdr < fdr_thresh))
  lead_variants = purrr::map_df(qtl_hits, identity) %>% unique()
  selected_variants = dplyr::filter(GRCh37_variants, snp_id %in% lead_variants$snp_id) %>% 
    dplyr::select(chr, pos, snp_id)
  
  #Add GRCh37 coordinates
  qtl_pos = purrr::map(qtl_hits, ~dplyr::left_join(., selected_variants, by = "snp_id") %>%
                         dplyr::filter(!is.na(pos)))

  #Identify genes that have associated variants nearby (ignoring LD)
  qtl_df_list = purrr::map(qtl_pos, ~dplyr::left_join(., gwas_pvals, by = "chr") %>%
                             dplyr::mutate(distance = abs(gwas_pos - pos)) %>%
                             dplyr::filter(distance < overlap_dist) %>%
                             dplyr::select(phenotype_id, snp_id) %>% unique())
  
}

colocGeneAgainstPeaks <- function(gene_df, peaks_df, eqtl_summaries, caqtl_summaries, 
                                  n_eqtl, n_caqtl, variant_information, cis_dist = 2e5){
  print(gene_df$phenotype_id)
  
  #Construct gene ranges
  eqtl_ranges = constructVariantRanges(gene_df, variant_information, cis_dist = cis_dist)
  
  #Fetch summary stats
  eqtl_summaries = fastqtlTabixFetchGenes(eqtl_ranges, eqtl_summaries)[[1]] %>%
    summaryReplaceSnpId(variant_information)
  
  #Construct peak ranges
  peak_ranges =  purrr::by_row(peaks_df, 
                               function(x, eqtl_ranges){eqtl_ranges$phenotype_id = x$peak_id; return(eqtl_ranges)}, eqtl_ranges)
  peak_ranges_list = setNames(peak_ranges$.out, peak_ranges$peak_id)
  
  #Fetch peak summaries
  caqtl_summaries_list = purrr::map(peak_ranges_list, ~fastqtlTabixFetchGenes(., caqtl_summaries)[[1]] %>%
                                      summaryReplaceSnpId(variant_information)) 

  #Perform coloc
  coloc_res = purrr::map(caqtl_summaries_list, 
                            ~testColoc(eqtl_summaries, ., n1 = n_eqtl, n2 = n_caqtl)$summary)
  coloc_res_filtered = coloc_res[!purrr::map_lgl(coloc_res, is.null)]#Remove NULL entries from list
  coloc_res_df = coloc_res_filtered %>%
    purrr::map_df(., ~dplyr::tbl_df(t(.)), .id = "peak_id") %>%
    dplyr::mutate(gene_id = gene_df$phenotype_id) %>%
    dplyr::select(gene_id, everything())
  
  return(coloc_res_df)
}
