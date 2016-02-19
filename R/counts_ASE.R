constructExonRanges <- function(gene_data){
  gene_df = data_frame(gene_id = gene_data$gene_id, 
                       seqnames = gene_data$chr, 
                       strand = gene_data$strand, 
                       start = as.numeric(unlist(strsplit(gene_data$exon_starts, ","))), 
                       end = as.numeric(unlist(strsplit(gene_data$exon_ends, ","))))
  granges = dataFrameToGRanges(gene_df)
  return(granges)
}

tabixFetchASE <- function(ranges, tabix_file){
  assertthat::assert_that(class(ranges)[1] == "GRanges")
  
  #Load tabix column names
  tabix_header = scan(tabix_file, nlines = 1, what = "character")
  #Specify column types
  column_types = paste("dd", paste(rep("c",length(tabix_header)-2), collapse = ""), sep = "")
  result = scanTabixDataFrame(tabix_file, param = ranges, col_names = tabix_header, col_types = column_types) %>%
    plyr::ldply(.id = "query")
  if(nrow(result) == 0){ #Output is an emtpy DF
    return(result)
  } else{
   result = dplyr::mutate(result,query = as.character(query)) %>%
     dplyr::rename(chr = contig, pos = position, feature_snp_id = variantID)
  }
  return(result)
}

rangesQueryMeta <- function(ranges){
  #Add metadata to the tabix import
  df = as.data.frame(elementMetadata(ranges))
  query_string = GenomicRanges::as.data.frame(ranges) %>% 
    tbl_df() %>% 
    dplyr::mutate(query = paste(paste(seqnames, start, sep = ":"), end, sep = "-")) %>% 
    dplyr::select(query)
  query_metadata = dplyr::bind_cols(query_string, df)
  return(query_metadata)
}

tidyASECounts <- function(ase_counts){
  result = tidyr::gather(ase_counts, key = "sample_id", value = "ASE", 7:ncol(a)) %>% 
    tidyr::separate(ASE,c("ref_count", "alt_count"), sep = ",") %>%
    dplyr::mutate(ref_count = as.numeric(ref_count), alt_count = as.numeric(alt_count)) %>%
    dplyr::mutate(total_count = ref_count + alt_count, max_count = pmax(ref_count, alt_count)) %>%
    dplyr::mutate(ratio = ifelse(ref_count > 0, ref_count/total_count, ifelse(total_count > 0, 1, NA)))
  return(result)
}

addVariantGenotype <- function(variant_df, genotypes, snp_col = "snp_id", value_id = "lead_snp"){
  variant_df = as.data.frame(variant_df)
  snp_id = variant_df[,snp_col][1]
  variant_gt = tidyVector(genotypes[snp_id,], value_id = value_id, sample_id = "genotype_id")
  result = dplyr::left_join(variant_df, variant_gt, by = "genotype_id")
  return(result)
}

#Fetch ASE data for a single gene
fetchGeneASEData <- function(gene_snp_pair, tabix_file, genotypes, sample_metadata, gene_metadata){
  assertthat::assert_that(nrow(gene_snp_pair) == 1) #Must contain only one gene and SNP
  
  gene_ranges = dplyr::filter(gene_metadata, gene_id == gene_snp_pair$gene_id) %>% constructExonRanges()
  ase_data = tabixFetchASE(gene_ranges, tabix_file)
  if(nrow(ase_data) > 0){
    ase_data = tidyASECounts(ase_data) %>%
      dplyr::left_join(rangesQueryMeta(gene_ranges), by = "query") %>% dplyr::select(-query) %>% #Add gene id
      dplyr::left_join(gene_snp_pair, by = "gene_id") %>% #Add lead SNP id
      dplyr::left_join(sample_metadata, ., by = "sample_id") %>%
      addVariantGenotype(genotypes, snp_col = "snp_id", value_id = "lead_snp")
    
    #Split by variant
    variant_list = dlply(ase_data, .variables = "feature_snp_id")
    variant_list = lapply(variant_list, addVariantGenotype, genotypes, snp_col = "feature_snp_id", value_id = "feature_snp")
    variant_df = ldply(variant_list, .id = NULL)
    return(variant_df)
  }else{
    return(ase_data)
  }
}

filterASEforPlotting <-function(ase_data){
  #Keep only hets
  hets = dplyr::filter(ase_data, lead_snp == 1, feature_snp == 1)
  
  #Remove feature SNPs with very low coverage
  expressed_snps = dplyr::group_by(hets, feature_snp_id) %>% dplyr::summarise(mean_total_count = mean(total_count)) %>% dplyr::filter(mean_total_count > 3)
  selected_hets = dplyr::semi_join(hets, expressed_snps, by = "feature_snp_id")
  return(selected_hets)
}
