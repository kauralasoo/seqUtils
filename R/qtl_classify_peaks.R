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
  variants_overlaps = variant_ranges[S4Vectors::queryHits(olaps)]
  peak_overlaps = peak_ranges[S4Vectors::subjectHits(olaps),]
  overlapping_peaks = dplyr::data_frame(snp_id = variants_overlaps$snp_id, overlap_peak_id = peak_overlaps$gene_id)
  
  #Add overlapping peaks
  peak_finemapping = dplyr::left_join(qtl_credible_set, overlapping_peaks, by = "snp_id")
  return(peak_finemapping)
}

#' Import credible sets from disk and add overlapping peaks to each QTL
importCredibleSets <- function(credible_set_path, peak_metadata){
  
  #Import raw credible sets
  credible_sets = readRDS(credible_set_path)
  
  #Add overlapping peaks to each variant in the credible set
  credible_sets_df = purrr::map(credible_sets, ~purrr::map_df(., 
                                 ~dplyr::mutate(.,chr = as.character(chr))) %>%
                                  dplyr::filter(chr != "X") %>%
                                  addOverlappingPeaks(peak_metadata) %>%
                                  unique()) #For some weird reason the credible sets contain duplicate SNPs
  return(credible_sets_df)
}

#' Identify all peak-peak pairs linked by credible sets
credibleSetPeakOverlaps <- function(credible_sets_df){
  shared_credible_sets = purrr::map(credible_sets_df, ~dplyr::filter(.,!is.na(overlap_peak_id))) %>% 
    purrr::map_df(., identity) %>% 
    dplyr::select(gene_id, overlap_peak_id) %>% 
    unique() %>% 
    dplyr::rename(master_peak_id = gene_id)
  return(shared_credible_sets)
}

#' Identify master peaks that share credible sets with each other
identifyAmbiguousMasters <- function(shared_credible_sets, master_peak_list){
  
  shared_masters = dplyr::semi_join(shared_credible_sets, master_peak_list, by = "master_peak_id") %>% 
    dplyr::semi_join(master_peak_list, by = c("overlap_peak_id" = "master_peak_id")) %>% 
    dplyr::filter(master_peak_id != overlap_peak_id)
  
  shared_master_peaks = dplyr::data_frame(master_peak_id = unique(unlist(shared_masters)))
  return(shared_master_peaks)
}

identifyUniqueMasters <- function(all_master_peaks, ambiguous_master_peaks, credible_sets_df){
  
  #Identify unique master peaks
  credible_sets_by_condition = purrr::map_df(credible_sets_df, identity, .id = "condition_name")
  unique_master_peaks = dplyr::anti_join(all_master_peaks, ambiguous_master_peaks, by = "master_peak_id")
  unique_master_pairs = dplyr::semi_join(credible_sets_by_condition, unique_master_peaks, by = c("gene_id" = "master_peak_id")) %>%
    dplyr::select(gene_id, snp_id, overlap_peak_id) %>% unique()
  
  #Count SNPs in the credible set and SNPs overlapping peak
  snp_count = dplyr::group_by(unique_master_pairs, gene_id) %>% dplyr::summarise(snp_count = length(snp_id))
  overlap_snp_count = dplyr::filter(unique_master_pairs, gene_id == overlap_peak_id) %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::summarise(overlap_snp_count = length(snp_id))
  
  #Add counts to the unique masters
  unique_masters_counted = dplyr::left_join(unique_master_pairs, snp_count, by = "gene_id") %>%
    dplyr::left_join(overlap_snp_count, by = "gene_id")
  
  #Find lead SNPs for unique masters accross conditions
  snp_stats = dplyr::left_join(unique_masters_counted, credible_sets_by_condition, by = c("gene_id","snp_id", "overlap_peak_id"))
  unique_lead_snps = group_by(snp_stats, gene_id) %>% 
    dplyr::arrange(p_nominal) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::select(gene_id, snp_id, snp_count, overlap_snp_count, condition_name, R2, chr, pos) %>%
    dplyr::ungroup()
  
  return(unique_lead_snps)
}

identifyLinkedLeadSNPs <- function(unique_master_peaks, genotypes){
  #Find out some additional potential multi-peak QTLs
  r2_filtered = filterHitsR2(dplyr::mutate(unique_master_peaks, gene_id = "gene"), genotypes)
  overlapping_snps = dplyr::anti_join(unique_master_peaks, r2_filtered, by = "snp_id")
  
  #For each SNP, find the other overlapping SNPs
  genotype_matrix = t(genotypes[unique_master_peaks$snp_id,])
  r2 = cor(genotype_matrix, use = "pairwise.complete.obs")^2
  olap_snps = r2[overlapping_snps$snp_id,]
  
  results = list()
  for(snp_id in overlapping_snps$snp_id) {
    olap_ids = names(which(olap_snps[snp_id, ] > 0.8))
    df = data_frame(master_snp_id = snp_id, overlap_snp_id = olap_ids)
    results[[snp_id]] = df
  }
  
  #Identify linked SNP pairs
  shared_snps = purrr::map_df(results, identity) %>% 
    dplyr::filter(master_snp_id != overlap_snp_id)
  return(shared_snps)
}

clusterAmbiguousMasters <- function(all_shared_masters, credible_sets_df){
  #Convert shared peaks into clusters
  atac_clusters = constructClustersFromGenePairs(all_shared_masters, cluster_name_prefix = "ATAC_cluster_")
  
  #Count the number of SNPs in each peak and cluster
  potential_master_peaks = purrr::map_df(credible_sets_df, identity, .id = "condition_name") %>% 
    dplyr::filter(gene_id == overlap_peak_id)
  atac_clusters_counted = dplyr::left_join(atac_clusters, potential_master_peaks, by = "gene_id") %>% 
    dplyr::group_by(gene_id, cluster_id) %>% 
    dplyr::summarise(peak_snp_count = length(snp_id)) %>% 
    dplyr::group_by(cluster_id) %>% 
    dplyr::mutate(cluster_snp_count = sum(peak_snp_count), peak_count = length(gene_id)) %>% 
    dplyr::ungroup() %>%
    arrange(cluster_snp_count)
  
  #Find lead SNPs for each cluster of peaks
  all_credible_sets = dplyr::filter(purrr::map_df(credible_sets_df, identity, .id = "condition_name"), 
                                    gene_id %in% atac_clusters_counted$gene_id)
  cluster_lead_snps = dplyr::left_join(atac_clusters_counted, all_credible_sets, by = "gene_id") %>% 
    dplyr::group_by(cluster_id) %>% 
    dplyr::arrange(cluster_id, p_nominal) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::select(cluster_id, gene_id, snp_id, R2, peak_count, cluster_snp_count, chr, pos) %>%
    dplyr::ungroup()
  
  return(list(cluster_memberships = atac_clusters_counted, lead_snps = cluster_lead_snps))
}
