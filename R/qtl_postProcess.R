
filterEQTLs <- function(data_frame, gene_id_name_map, fdr_cutoff = 0.1){
  dat = dplyr::filter(data_frame, FDR < fdr_cutoff) %>% 
    dplyr::rename(gene_id = gene, snp_id = snps) %>%
    dplyr::group_by(gene_id) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::filter(row_number() == 1) %>%
    dplyr::ungroup() %>% 
    dplyr::arrange(pvalue) %>%
    dplyr::left_join(gene_id_name_map, by = "gene_id")
  return(dat)
}

testInteraction <- function(gene_id, genotype_id, expression_dataset, genotype_dataset, line_metadata){
  genotype_ids = dplyr::select(line_metadata, donor, replicate, genotype_id)
  geno = genotype_dataset$genotypes[genotype_id,]
  geno_df = data_frame(genotype_id = names(geno), genotype = geno)
    
  exprs = expression_dataset$exprs_cqn[gene_id,]
  exprs_df = data_frame(sample_id = names(exprs), expression = exprs)
  
  data = dplyr::left_join(expression_dataset$design, exprs_df, by = "sample_id") %>% 
    dplyr::left_join(genotype_ids, by = c("donor","replicate")) %>%
    dplyr::left_join(geno_df, by = "genotype_id")
  
  no_interaction = lm(expression~genotype + condition_name, as.data.frame(data))
  interaction = lm(expression~genotype + condition_name + condition_name:genotype, as.data.frame(data))
  res = anova(no_interaction, interaction)
  return(res)
}

testMultipleInteractions <- function(snps_df, expression_dataset, genotype_dataset, line_metadata){
  #Plot eQTL results for a list of gene and SNP pairs.
  result = list()
  for(i in 1:nrow(snps_df)){
    gene_id = snps_df[i,]$gene_id
    snp_id = snps_df[i,]$snp_id
    print(gene_id)
    test = testInteraction(gene_id, snp_id, expression_dataset, genotype_dataset, line_metadata)
    result[[gene_id]] = test
  }
  return(result)
}
