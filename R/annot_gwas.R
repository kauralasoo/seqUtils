#' Import EBI GWAS catalog into R.
#'
#' Imports tsv version of the GWAS catalog into R. Parses total sample size from the 
#' INITIAL SAMPLE DESCRIPTION field and also flags each study that contains European samples.
#' 
#' @param path Path to the EBI GWAS catalog export file (must contain MAPPED_TRAIT field as well).
#' @return data frame containing the gwas catalog.
#' @author Kaur Alasoo
#' @export 
importGwasCatalog <- function(path){
  
  #Load gwas catlog from disk
  gwas_catalog = readr::read_delim(path, delim = "\t")
  
  #Construct a data frame of studies
  studies = gwas_catalog[,c("PUBMEDID","INITIAL SAMPLE DESCRIPTION")] %>% unique()
  colnames(studies) = c("pubmed_id", "description")
  
  #Extract sample size from the gwas catalog
  sample_sizes = lapply(as.list(studies$description), function(x){
    stringr::str_extract_all(x, "(\\d+,\\d+)|(\\d+)") %>%
      unlist() %>%
      stringr::str_replace(",","") %>%
      as.numeric() %>%
      sum()
  })
  
  #Add sample size and European status to study df
  studies_df = dplyr::mutate(studies, sample_size = unlist(sample_sizes)) %>%
    dplyr::mutate(is_european = grepl("european",description,ignore.case = TRUE)) %>%
    dplyr::select(pubmed_id, sample_size, is_european)
  
  #Extract snps form catalog
  gwas_columns = gwas_catalog[,c("PUBMEDID","CHR_ID","CHR_POS","SNPS","DISEASE/TRAIT","MAPPED_TRAIT","P-VALUE")]
  colnames(gwas_columns) = c("pubmed_id", "chr", "pos", "snp_id", "trait", "mapped_trait","gwas_pvalue")
  
  #Join study data with SNPS and remove NAs
  gwas_table = dplyr::left_join(studies_df, gwas_columns, by = "pubmed_id") %>%
    dplyr::filter(!is.na(chr), !is.na(pos)) %>%
    dplyr::mutate(chr = as.character(chr)) %>%
    dplyr::mutate(chr = ifelse(chr == "23","X",chr))
  
  return(gwas_table)
}

#' Calculate the enrichment of GWAS traits among eQTLs
#'
#' Based on Supplementary Figure S16 in (Kumasaka et al 2015).  
#' 
#' @param gwas_catalog EBI GWAS catalog imported using the importGwasCatalog function.
#' @param gwas_snp_pvalues FastQTL eQTL pvalues for each SNP in the GWAS catalog, from importGwasPvalue.
#' @param eqtl_lead_snps Lead eQTL SNPs for each significant gene.
#' @param gene_id_name_map - data frame with two columns: gene_id, gene_name.
#' @param n_genes Total numner of genes used in the eQTL analysis.
#' @param p_diff_thresh Log10 p-value difference between lead eQTL SNP and lead GWAS SNP.
#' @param min_genes_per_trait Minimal number of unique genes per trait to be used for enrichment analysis.
#' @return list(enrichment = traits enriched among eQTLs, overlap = table of overlapping genes and traits).
#' @author Kaur Alasoo
#' @export 
calculateEnrichment <- function(gwas_catalog, gwas_snp_pvalues, eqtl_lead_snps, gene_id_name_map, 
                                n_genes = 14163, p_diff_thresh = 3, min_genes_per_trait = 20){
  
  #For each trait-SNP pair find the most associated eQTL gene
  genes_with_traits = dplyr::inner_join(gwas_snp_pvalues, gwas_catalog, by = c("chr", "pos")) %>%
    dplyr::select(gene_id, chr, pos, pvalue, trait, mapped_trait) %>% 
    dplyr::group_by(chr, pos, mapped_trait) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(gene_id, chr, pos, pvalue, mapped_trait) %>% 
    dplyr::group_by(gene_id, mapped_trait) %>% 
    dplyr::arrange(pvalue) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::rename(trait_chr = chr, trait_pos = pos, trait_pvalue = pvalue)
  
  #Count the number of idendentified genes per trait
  genes_per_trait = genes_with_traits %>% 
    group_by(mapped_trait) %>% 
    summarize(trait_gene_count = length(gene_id)) %>% 
    arrange(-trait_gene_count)
  
  #Join QTLs and GWAS hits
  trait_qtl_map = dplyr::inner_join(eqtl_lead_snps, genes_with_traits, by = "gene_id") %>%
    dplyr::mutate(pvalue_diff = -log10(p_nominal) +log10(trait_pvalue)) %>% 
    dplyr::filter(pvalue_diff < p_diff_thresh) %>% #Difference in eQTL pvalues at most 3 orders of magnitude
    dplyr::left_join(gene_id_name_map, by = "gene_id") %>% 
    dplyr::select(gene_id, gene_name, snp_id, p_nominal, trait_pvalue, pvalue_diff, mapped_trait, everything())
  
  #Count the number QTL genes that overlap a trait
  trait_qtl_overlap = trait_qtl_map %>% 
    group_by(mapped_trait) %>% 
    dplyr::summarise(trait_qtl_overlap = length(gene_id)) %>% 
    arrange(-trait_qtl_overlap)
  
  #Calculate odd-ratios of enrichment
  n_qtls = nrow(fastqtl_callset$naive)
  n_genes = nrow(eqtl_data_list$exprs_cqn)
  n_genes_with_traits = length(unique(genes_with_traits$gene_id))
  n_qtls_with_traits = length(unique(trait_qtl_map$gene_id))
  overlap_df = dplyr::left_join(genes_per_trait, trait_qtl_overlap, by = "mapped_trait") %>% 
    dplyr::mutate(trait_qtl_overlap = ifelse(is.na(trait_qtl_overlap), 0, trait_qtl_overlap)) %>%
    dplyr::mutate(OR_all_qtls = (trait_qtl_overlap/(n_qtls - trait_qtl_overlap)) / 
                    ((trait_gene_count-trait_qtl_overlap)/(n_genes-n_qtls))) %>%
    dplyr::mutate(OR_trait_qtls = (trait_qtl_overlap/(n_qtls_with_traits - trait_qtl_overlap)) / 
                    ((trait_gene_count-trait_qtl_overlap)/(n_genes_with_traits-n_qtls_with_traits))) %>%
    dplyr::arrange(-OR_trait_qtls) %>%
    dplyr::filter(trait_gene_count > min_genes_per_trait)
  
  return(list(enrichment = overlap_df, overlap = trait_qtl_map))
}

importGwasPvalue <- function(path, gwas_catalog){
  n_pval = readr::read_delim(path, delim = " ", col_names = FALSE)
  colnames(n_pval) = c("gene_id", "chr","pos", "snp_id", "distance", "pvalue", "effect_size")
  gwas_pvalues = dplyr::semi_join(n_pval, gwas_catalog, by = c("chr","pos"))
  
  return(gwas_pvalues)
}