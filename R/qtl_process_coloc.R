#Functions
identifyColocHits <- function(coloc_df, PP_power_thresh = 0.8, PP_coloc_thresh = 0.9, nsnps_thresh = 10){
  coloc_hits = dplyr::filter(coloc_df, PP_power > PP_power_thresh) %>% 
    dplyr::group_by(trait, gwas_lead, phenotype_id) %>% 
    dplyr::arrange(trait, phenotype_id, gwas_lead, -PP.H4.abf) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::filter(PP_coloc > PP_coloc_thresh) %>%
    dplyr::filter(nsnps > nsnps_thresh) %>%
    dplyr::ungroup()
  return(coloc_hits)
}

countConditionSpecificOverlaps <- function(coloc_filtered, PP_power_thresh = 0.8, PP_coloc_thresh = 0.9){
  #Count the number of overlaps added by each additional condition
  coloc_counts = dplyr::mutate(coloc_filtered, is_hit = ifelse(PP_power > PP_power_thresh & PP_coloc > PP_coloc_thresh, 1, 0)) %>% 
    dplyr::select(summarised_trait, phenotype_id, snp_id, condition_name, is_hit) %>% 
    tidyr::spread(condition_name, is_hit) %>%
    dplyr::mutate(naive_IFNg = pmax(naive, IFNg)) %>%
    dplyr::mutate(naive_IFNg_SL1344 = pmax(naive_IFNg, SL1344)) %>%
    dplyr::mutate(all = pmax(naive_IFNg_SL1344, IFNg_SL1344)) %>%
    dplyr::mutate(IFNg_added = naive_IFNg - naive) %>%
    dplyr::mutate(SL1344_added = naive_IFNg_SL1344 - naive_IFNg) %>%
    dplyr::mutate(IFNg_SL1344_added = all - naive_IFNg_SL1344) %>%
    dplyr::select(summarised_trait, phenotype_id, snp_id, naive, IFNg_added, SL1344_added, IFNg_SL1344_added) %>%
    dplyr::rename(IFNg = IFNg_added, SL1344 = SL1344_added, IFNg_SL1344 = IFNg_SL1344_added) %>%
    tidyr::gather("condition_name", "is_hit", naive:IFNg_SL1344) %>% 
    dplyr::filter(is_hit == 1) %>%
    dplyr::left_join(figureNames(), by = "condition_name")
}

importAndFilterColocHits <- function(gwas_stats, coloc_suffix = ".eQTL.1e+05.coloc.txt", coloc_prefix = "results/SL1344/coloc/coloc_lists/",
                                     PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, gwas_pval_thresh = 1e-6, mhc_phenotypes){
  #Import coloc hits
  
  #Name all files
  file_names = as.list(paste0(coloc_prefix, gwas_stats_labeled$trait, coloc_suffix))
  names(file_names) = gwas_stats_labeled$trait
  
  #Import enrichments
  coloc_df = purrr::map_df(file_names, ~readr::read_delim(., delim = "\t"), .id = "trait") %>%
    dplyr::mutate(PP_power = (PP.H4.abf + PP.H3.abf), PP_coloc = PP.H4.abf/PP_power) %>%
    dplyr::mutate(summarised_trait = ifelse(trait %in% c("IBD","UC","CD"), "IBD", trait))
  
  #Identify one overlap GWAS lead varaint
  coloc_hits = identifyColocHits(coloc_df, PP_power_thresh, PP_coloc_thresh, nsnps_thresh) %>%
    dplyr::filter(gwas_pval < 1e-6) %>%
    dplyr::filter(!(phenotype_id %in% mhc_phenotypes$phenotype_id))
  
  #Merge IBD overlaps together, keep only stronges association per gene
  coloc_hits_merged = dplyr::group_by(coloc_hits, summarised_trait, phenotype_id) %>% 
    dplyr::arrange(summarised_trait, phenotype_id, -PP.H4.abf) %>% 
    dplyr::filter(row_number() == 1) %>%
    dplyr::ungroup()
  
  #Retreive posterior probabilitites in all conditions
  coloc_filtered = dplyr::semi_join(coloc_df, coloc_hits_merged, by = c("trait", "phenotype_id", "snp_id"))
  return(list(coloc_filtered = coloc_filtered, coloc_df = coloc_df))
}



importSummariesForPlotting <- function(qtl_df, gwas_stats_labeled, gwas_dir, 
                                       qtl_paths, GRCh37_variants, GRCh38_variants, cis_dist = 1e5, use_rasqual = FALSE){
  
  assertthat::assert_that(assertthat::has_name(qtl_df, "phenotype_id"))
  assertthat::assert_that(assertthat::has_name(qtl_df, "snp_id"))
  assertthat::assert_that(assertthat::has_name(qtl_df, "trait"))
  
  #Print gene_id
  print(qtl_df$phenotype_id)
  
  #Make GRanges objectsto fetch data
  qtl_ranges = constructVariantRanges(qtl_df, GRCh38_variants, cis_dist = cis_dist)
  gwas_ranges = constructVariantRanges(qtl_df, GRCh37_variants, cis_dist = cis_dist)
  
  #Set Identify the location of the gwas file
  gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == qtl_df$trait)$file_name
  gwas_prefix = file.path(gwas_dir, gwas_file_name)
  
  #Import GWAS pvalues
  trait_name = as.vector(qtl_df$trait)
  gwas_summaries = tabixFetchGWASSummary(gwas_ranges, summary_path = paste0(gwas_prefix, ".sorted.txt.gz"))[[1]] %>%
    summaryReplaceSnpId(GRCh37_variants) %>% 
    summaryReplaceCoordinates(GRCh38_variants) %>%
    dplyr::transmute(condition_name = trait_name, snp_id, chr, pos, p_nominal, log10p = -log(p_nominal, 10))
  
  #Fetch QTL summary stats
  if(use_rasqual){
    qtl_fetch_function = rasqualTools::tabixFetchGenes
  } else{
    qtl_fetch_function = fastqtlTabixFetchGenes
  }
  qtl_summaries = purrr::map_df(qtl_paths, ~qtl_fetch_function(qtl_ranges, .)[[1]], .id = "condition_name") %>%
    dplyr::transmute(condition_name, snp_id, chr, pos, p_nominal, log10p = -log(p_nominal, 10))
  
  #Merge data
  full_data = dplyr::bind_rows(gwas_summaries, qtl_summaries) %>% 
    dplyr::mutate(condition_name = factor(condition_name, levels = c(trait_name, names(qtl_paths))))
  
  return(full_data)
  
}


