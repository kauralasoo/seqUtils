
#' Add ATAC peak annotations to the eQTL credible set
#'
#' @param eqtl_credible_set 
#' @param atac_peak_metadata 
#' @param atac_tabix_path 
#'
#' @return Credible set with ATAC peak annotations
#' @export
annotateCredibleSet <- function(eqtl_credible_set, atac_peak_metadata, atac_tabix_path){
  
  #Construct Granges for the credible set SNPs
  variant_ranges = dplyr::transmute(eqtl_credible_set, seqnames = chr, start = pos, end = pos, strand = "+", snp_id) %>%
    dataFrameToGRanges()
  
  #Add overlapping peaks
  peak_finemapping = addOverlappingPeaks(eqtl_credible_set, atac_peak_metadata)
  
  #Find most associated peaks for each SNP
  associated_peaks = tabixFetchSNPs(variant_ranges, atac_tabix_path)
  if (!is.null(associated_peaks)){
    associated_peaks = associated_peaks %>% 
      dplyr::group_by(snp_id) %>% 
      dplyr::arrange(p_nominal) %>% 
      dplyr::filter(row_number() == 1) %>% 
      dplyr::ungroup() %>% 
      dplyr::transmute(assoc_peak_id = gene_id, snp_id, p_nominal_peak = p_nominal)
    
    #Add associated peaks
    peak_finemapping = dplyr::left_join(peak_finemapping, associated_peaks, by = "snp_id")
  }else{
    peak_finemapping = dplyr::mutate(peak_finemapping, assoc_peak_id = NA)
  }
  
  return(peak_finemapping)
}

#' Add overlapping ATAC peaks to a credible set of QTL variants
#'
#' @param qtl_credible_set Data frame with credible set of QTL variants 
#' (columns: snp_id, chr, pos)
#' @param atac_peak_metadata Data frame of start end coordinates of ATAC peaks 
#' (columns: gene_id, chr, start, end, strand)
#' @param extend Number of bp added to peak start and end coordinates to increase the chance of finding an overlap.
#'
#' @return Modifed qtl_credible_set data frame with added overlap_peak_id column.
#' @export
addOverlappingPeaks <- function(qtl_credible_set, atac_peak_metadata, extend = 0) {
  
  #Some basic assertions
  #Check that atac_peak_metadata has required columns
  assertthat::assert_that(assertthat::has_name(atac_peak_metadata, "gene_id"))
  assertthat::assert_that(assertthat::has_name(atac_peak_metadata, "chr"))
  assertthat::assert_that(assertthat::has_name(atac_peak_metadata, "start"))
  assertthat::assert_that(assertthat::has_name(atac_peak_metadata, "end"))
  assertthat::assert_that(assertthat::has_name(atac_peak_metadata, "strand"))
  
  #Check that qtl_credible_set has required columns
  assertthat::assert_that(assertthat::has_name(qtl_credible_set, "snp_id"))
  assertthat::assert_that(assertthat::has_name(qtl_credible_set, "chr"))
  assertthat::assert_that(assertthat::has_name(qtl_credible_set, "pos"))
  
  assertthat::assert_that(is.numeric(extend))
  
  #Extend peaks by some buffer sequence
  if(extend > 0){
    atac_peak_metadata = dplyr::mutate(atac_peak_metadata, start = start - extend, end = end + extend)
  }
  
  #construct Granges object with atac peak data
  peak_ranges = dataFrameToGRanges(dplyr::rename(atac_peak_metadata, seqnames = chr))
  
  #Construct Granges for the credible set SNPs
  variant_ranges = dplyr::transmute(qtl_credible_set, seqnames = chr, start = pos, end = pos, strand = "+", snp_id) %>%
    dataFrameToGRanges()
  
  #Find overlapping peaks for each variant
  olaps = GenomicRanges::findOverlaps(variant_ranges, peak_ranges, ignore.strand = TRUE)
  variants_overlaps = variant_ranges[GenomicRanges::queryHits(olaps)]
  peak_overlaps = peak_ranges[GenomicRanges::subjectHits(olaps),]
  overlapping_peaks = dplyr::data_frame(snp_id = variants_overlaps$snp_id, overlap_peak_id = peak_overlaps$gene_id)
  
  #Add overlapping peaks
  peak_finemapping = dplyr::left_join(qtl_credible_set, overlapping_peaks, by = "snp_id")
  return(peak_finemapping)
}


#Calculate summary stats on the credible sets
summariseCredibleSets <- function(credible_sets_df){
  
  #Calculate credible set size
  credible_set_size = dplyr::group_by(credible_sets_df, gene_id) %>% 
    dplyr::arrange(p_nominal) %>% 
    dplyr::summarise(credible_set_size = length(snp_id))
  
  #Remove SNPs that do not overlap peaks
  peaks_filtered = dplyr::filter(credible_sets_df, !is.na(overlap_peak_id))
  
  #Find peaks where at least one variant lies within the peak
  overlap_snp_count = dplyr::filter(peaks_filtered, gene_id == overlap_peak_id) %>% 
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(overlap_snp_count = length(snp_id)) %>%
    dplyr::ungroup()
  
  #Count the number of peaks overlapping SNPs in the credible set
  overlap_peak_count = peaks_filtered %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::summarise(overlap_peak_count = length(overlap_peak_id)) %>% 
    ungroup()
  
  #Find lead SNPs
  lead_snps = dplyr::group_by(peak_cs_df, gene_id) %>% dplyr::arrange(p_nominal) %>% 
    dplyr::filter(row_number() == 1) %>% dplyr::ungroup() %>%
    dplyr::left_join(credible_set_size, by = "gene_id") %>%
    dplyr::left_join(overlap_snp_count, by = "gene_id") %>%
    dplyr::mutate(overlap_snp_count = ifelse(is.na(overlap_snp_count), 0, overlap_snp_count)) %>%
    dplyr::left_join(overlap_peak_count, by = "gene_id") %>%
    dplyr::mutate(overlap_peak_count = ifelse(is.na(overlap_peak_count), 0, overlap_peak_count))
  
  return(lead_snps)
}

constructClustersFromGenePairs <- function(gene_pairs, cluster_name_prefix = "gene_cluster_"){
  assertthat::assert_that(ncol(gene_pairs) == 2)
  
  #Construct igraph object and extract clusters
  graph = igraph::graph_from_data_frame(gene_pairs, directed = FALSE)
  clusters = igraph::clusters(graph)
  
  #Make df
  clusters_df = data_frame(gene_id = names(clusters$membership), cluster_number = clusters$membership) %>% 
    arrange(cluster_number) %>% 
    dplyr::mutate(cluster_id = paste0(cluster_name_prefix, cluster_number)) %>% 
    dplyr::select(gene_id, cluster_id)
  
  #Count genes in clusters
  gene_counts = dplyr::group_by(clusters_df, cluster_id) %>% 
    dplyr::summarise(gene_count = length(gene_id))
  clusters_df = dplyr::left_join(clusters_df, gene_counts, by = "cluster_id")
  return(clusters_df)
  
}

credibleSetJaccard <- function(row_df, master_credible_sets, dependent_credible_sets){
  
  #Extract credible sets
  master_cs = master_credible_sets[[row_df$master_condition]][[row_df$master_id]] %>% dplyr::arrange(p_nominal)
  dependent_cs = dependent_credible_sets[[row_df$dependent_condition]][[row_df$dependent_id]] %>% dplyr::arrange(p_nominal)

  #Calculate two types of jaccard indexes
  jaccard = length(intersect(master_cs$snp_id, dependent_cs$snp_id)) / length(union(master_cs$snp_id, dependent_cs$snp_id))
  min_jaccard = length(intersect(master_cs$snp_id, dependent_cs$snp_id)) / 
    min(length(master_cs$snp_id), length(dependent_cs$snp_id))
  res = data_frame(jaccard = jaccard, min_jaccard = min_jaccard, snp_id = master_cs$snp_id[1], p_nominal = master_cs$p_nominal[1])
  return(res)
}

#Convert credible sets into gigantic data frame
credibleSetsToDf <- function(credible_sets){
  result = purrr::map_df(credible_sets, ~purrr::map_df(., 
                  ~dplyr::mutate(.,chr = as.character(chr))) %>% 
                  dplyr::filter(chr != "X"), .id = "condition_name")
  return(result)
}

credibleSetsToGranges <- function(credible_sets_df){
  granges = credible_sets_df %>% 
    dplyr::transmute(condition_name, gene_id, snp_id, chr, seqnames = chr, start = pos, end = pos, strand = "+") %>% 
    dataFrameToGRanges()
  return(granges)
}

findCredibleSetOverlaps <- function(master_credible_set, dependent_credible_set){
  #Find overlaps
  olaps = findOverlaps(master_credible_set, dependent_credible_set)
  
  #Eextract query overlaps
  queries = master_credible_set[queryHits(olaps),] %>% 
    elementMetadata() %>% as.data.frame() %>% 
    dplyr::transmute(master_condition = condition_name, master_id = gene_id, chr1 = chr) %>% 
    tbl_df()
  #Extract subject overlaps
  subjects = dependent_credible_set[subjectHits(olaps),] %>% 
    elementMetadata() %>% as.data.frame() %>% 
    dplyr::transmute(dependent_condition = condition_name, dependent_id = gene_id, chr2 = chr) %>% 
    tbl_df()
  
  #Identify all genes that share at least one SNP in their credible set
  pairwise_shared = dplyr::bind_cols(queries, subjects) %>% 
    unique() %>% 
    dplyr::filter(chr1 == chr2) %>% 
    dplyr::filter(master_id != dependent_id) %>% 
    dplyr::select(master_id, master_condition, dependent_id, dependent_condition) %>% 
    unique()
  return(pairwise_shared)
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

tabixFetchGWASSummary <- function(granges, summary_path){
  gwas_col_names = c("snp_id", "chr", "pos", "effect_allele", "MAF", 
                     "p_nominal", "beta", "OR", "log_OR", "se", "z_score", "trait", "PMID", "used_file")
  gwas_col_types = c("ccicdddddddccc")
  gwas_pvalues = scanTabixDataFrame(summary_path, granges, col_names = gwas_col_names, col_types = gwas_col_types)
  return(gwas_pvalues)
}

importGWASSummary <- function(summary_path){
  gwas_col_names = c("snp_id", "chr", "pos", "effect_allele", "MAF", 
                     "p_nominal", "beta", "OR", "log_OR", "se", "z_score", "trait", "PMID", "used_file")
  gwas_col_types = c("ccicdddddddccc")
  gwas_pvals = readr::read_tsv(summary_path,
                               col_names = gwas_col_names, col_types = gwas_col_types, skip = 1)
  return(gwas_pvals)
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
#' @param qtl QTL summary stats (p_nominal, MAF, neta, snp_id)
#' @param gwas GWAS summary stats(beta, se, MAF)
#' @param N_qtl Sample size of the QTL mapping study
#'
#' @return
#' @export
colocQtlGWAS <- function(qtl, gwas, N_qtl){
  
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
                               GRCh37_variants, GRCh38_variants, qtl_type = "rasqual",
                               N_qtl = 84, cis_dist = 1e5){
  
  result = tryCatch({
    #Make GRanges object to fetch data
    qtl_ranges = constructVariantRanges(qtl_df, GRCh38_variants, cis_dist = cis_dist)
    gwas_ranges = constructVariantRanges(qtl_df, GRCh37_variants, cis_dist = cis_dist)
    
    #Fetch QTL summary stats
    if(qtl_type == "rasqual"){
      qtl_summaries = rasqualTools::tabixFetchGenes(qtl_ranges, qtl_summary_path)[[1]]
    } else{
      qtl_summaries = fastqtlTabixFetchGenes(qtl_ranges, qtl_summary_path)[[1]]
    }
    #Fetch GWAS summary stats
    gwas_summaries = tabixFetchGWASSummary(gwas_ranges, gwas_summary_path)[[1]]
    
    #Substitute coordinate for the eqtl summary stats and add MAF
    qtl = summaryReplaceCoordinates(qtl_summaries, GRCh37_variants)

    #Substitute snp_id for the GWAS summary stats and add MAF
    gwas = summaryReplaceSnpId(gwas_summaries, GRCh37_variants)
    
    #Perform coloc analysis
    coloc_res = colocQtlGWAS(qtl, gwas, N_qtl = N_qtl)
    coloc_summary = dplyr::tbl_df(t(data.frame(coloc_res$summary))) %>%
      dplyr::mutate(qtl_pval = min(qtl$p_nominal), gwas_pval = min(gwas$p_nominal)) #Add minimal pvalues
    
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

prefilterColocCandidates <- function(qtl_min_pvalues, gwas_prefix, GRCh37_variants, fdr_thresh = 0.1, overlap_dist = 1e5, gwas_thresh = 1e-5){
  
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
                             dplyr::select(gene_id, snp_id) %>% unique())
  
}





