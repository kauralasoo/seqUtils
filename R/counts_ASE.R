constructExonRanges <- function(selected_gene_id, lead_snp_id, gene_metadata){
  gene_data = dplyr::filter(gene_metadata, gene_id == selected_gene_id)
  gene_df = data_frame(gene_id = gene_data$gene_id,
                       lead_snp_id = lead_snp_id,
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
  df = as.data.frame(GenomicRanges::elementMetadata(ranges))
  query_string = GenomicRanges::as.data.frame(ranges) %>% 
    tbl_df() %>% 
    dplyr::mutate(query = paste(paste(seqnames, start, sep = ":"), end, sep = "-")) %>% 
    dplyr::select(query)
  query_metadata = dplyr::bind_cols(query_string, df)
  return(query_metadata)
}

tidyASECounts <- function(ase_counts){
  result = tidyr::gather(ase_counts, key = "sample_id", value = "ASE", 7:ncol(ase_counts)) %>% 
    tidyr::separate(ASE,c("ref_count", "alt_count"), sep = ",") %>%
    dplyr::mutate(ref_count = as.numeric(ref_count), alt_count = as.numeric(alt_count)) %>%
    dplyr::mutate(total_count = ref_count + alt_count, max_count = pmax(ref_count, alt_count)) %>%
    dplyr::mutate(ratio = ifelse(ref_count > 0, ref_count/total_count, ifelse(total_count > 0, 1, NA)))
  return(result)
}

#Fetch ASE data for a single gene
fetchGeneASEData <- function(gene_ranges, tabix_file, sample_metadata){
  assertthat::assert_that(class(gene_ranges)[1] == "GRanges")

  ase_data = tabixFetchASE(gene_ranges, tabix_file)
  if(nrow(ase_data) > 0){
    ase_data = tidyASECounts(ase_data) %>% #Tidyfy
      dplyr::left_join(rangesQueryMeta(gene_ranges), by = "query") %>% dplyr::select(-query) %>% #Add gene id and lead snp id
      dplyr::mutate(sample_id = as.character(sample_id)) %>%
      dplyr::left_join(sample_metadata, ., by = "sample_id") #Add sample metadata
    return(ase_data)
  }else{
    return(ase_data)
  }
}

filterASEforPlotting <-function(ase_data){
  #Keep only hets
  hets = dplyr::filter(ase_data, feature_snp_value == 1)
  
  #Remove feature SNPs with very low coverage
  expressed_snps = dplyr::group_by(hets, feature_snp_id) %>% dplyr::summarise(mean_total_count = mean(total_count)) %>% dplyr::filter(mean_total_count > 10)
  selected_hets = dplyr::semi_join(hets, expressed_snps, by = "feature_snp_id")
  return(selected_hets)
}

#' Convert genotype matrix from gdsToMatrix into a tidy data frame
#'
#' @param genotypes Matrix with variants in rows and individuals in columms.
#' @param selected_snp_id Name for the snp_id column (default: snp_id).
#' @param value_id Name for the allele dosage column (default: snp_value).
#'
#' @return Tidy data frame with one genotype per row.
#' @export
tidyGenotypeMatrix <- function(genotypes, selected_snp_id = "snp_id", value_id = "snp_value"){
  rows = rownames(genotypes)
  result = tbl_df(as.data.frame(genotypes)) %>% 
    dplyr::mutate(snp_id = rows) %>% 
    tidyr::gather(genotype_id, value_id, 1:ncol(genotypes))  %>%
    dplyr::mutate(genotype_id = as.character(genotype_id)) %>%
    dplyr::rename_(.dots = setNames("snp_id", selected_snp_id)) %>%
    dplyr::rename_(.dots = setNames("value_id", value_id))
  return(result)
}

#' Title Add genotype data to tidy ASE data frame.
#'
#' @param ase_data Tidy ASE data frame with at least the followign columns: feature_snp_id, lead_snp_id, genptype_id.
#' @param genotypes Standard genotype matrix (variants in rows, individuals in columns)
#'
#' @return ASE data frame that contains genotypes for each individual for each feature snp and lead snp.
#' @export
aseDataAddGenotypes <- function(ase_data, genotypes){
  assertthat::assert_that(nrow(ase_data) > 0)
  
  gt_matrix = genotypes[unique(ase_data$feature_snp_id),,drop = FALSE]
  tidy_gt_matrix = tidyGenotypeMatrix(gt_matrix, value_id = "feature_snp_value", selected_snp_id = "feature_snp_id")
  ase_data_gt = dplyr::left_join(ase_data, tidy_gt_matrix, by = c("feature_snp_id", "genotype_id")) %>%
    dplyr::left_join(tidyGenotypeMatrix(genotypes[unique(ase_data$lead_snp_id),,drop = FALSE], 
        value_id = "lead_snp_value", selected_snp_id = "lead_snp_id"), by = c("lead_snp_id", "genotype_id"))
  return(ase_data_gt)
}
