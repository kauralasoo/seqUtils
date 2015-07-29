vcfToMatrix <- function(file, genome){
  #Read vcf file into R and convert it onto matrix of SNP positions and matrix of genotypes
  genotypes_vcf = VariantAnnotation::readVcf(file, genome)
  
  # Extract SNP positions from the VCF file
  variant_granges = GenomicRanges::rowData(genotypes_vcf)
  elementMetadata(variant_granges) = c()
  print(head(variant_granges))
  snp_positions = GenomicRanges::as.data.frame(variant_granges)
  snpspos = dplyr::mutate(snp_positions, snpid = rownames(snp_positions)) %>% 
    dplyr::select(snpid, seqnames, start) %>%
    dplyr::rename(chr = seqnames, pos = start)
  
  #Extract genotype matrix
  genotypes = VariantAnnotation::geno(genotypes_vcf)$GT
  genotypes[genotypes == "1/1"] = 2
  genotypes[genotypes == "0/1"] = 1
  genotypes[genotypes == "1/0"] = 1
  genotypes[genotypes == "0/0"] = 0
  genotypes[genotypes == "."] = "NA"
  mode(genotypes) = "numeric"
  
  return(list(snpspos = snpspos, genotypes = genotypes))
}

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

plotEQTL <- function(selected_gene_id, genotype_id, expression_dataset, genotype_dataset, line_metadata){
  
  #Extraxt gene_name
  gene_name = dplyr::filter(expression_dataset$gene_metadata, gene_id == selected_gene_id)$gene_name
  print(gene_name)
  
  #extract genotypes
  geno_vector = genotype_dataset$genotypes[genotype_id,]
  genotype_df = data.frame(genotype_id = names(geno_vector), genotype_value = as.character(geno_vector), stringsAsFactors = FALSE, row.names = NULL)
  
  #expression
  expression_vector = expression_dataset$exprs_cqn[selected_gene_id,]
  exprs_df = data.frame(sample_id = names(expression_vector), norm_exp = expression_vector, stringsAsFactors = FALSE, row.names = NULL)
  
  #Map genotype ids to donor names
  donor_genotype_map = dplyr::select(line_metadata, donor, genotype_id) %>% unique()
  
  plot_df = dplyr::left_join(expression_dataset$design, donor_genotype_map, by ="donor") %>%
    left_join(exprs_df, by = "sample_id") %>% 
    left_join(genotype_df, by = "genotype_id") %>%
    dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg", "SL1344", "IFNg_SL1344")))
  
  plot = ggplot(plot_df, aes(x = genotype_value, y = norm_exp)) + 
    facet_wrap(~ condition_name) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(position = position_jitter(width = .1)) + 
    ylab("Normalized expression") +
    xlab(genotype_id) + 
    labs(title = gene_name)
  
  return(plot)
}

plotEQTL_acLDL <- function(selected_gene_id, genotype_id, expression_dataset, genotype_dataset, line_metadata){
  
  #Extraxt gene_name
  gene_name = dplyr::filter(expression_dataset$gene_metadata, gene_id == selected_gene_id)$gene_name
  print(gene_name)
  
  #extract genotypes
  geno_vector = genotype_dataset$genotypes[genotype_id,]
  genotype_df = data.frame(genotype_id = names(geno_vector), genotype_value = as.character(geno_vector), stringsAsFactors = FALSE, row.names = NULL)
  
  #expression
  expression_vector = expression_dataset$exprs_cqn[selected_gene_id,]
  exprs_df = data.frame(sample_id = names(expression_vector), norm_exp = expression_vector, stringsAsFactors = FALSE, row.names = NULL)
  
  #Map genotype ids to donor names
  donor_genotype_map = dplyr::select(line_metadata, donor, genotype_id) %>% unique()
  
  plot_df = dplyr::left_join(expression_dataset$design, donor_genotype_map, by ="donor") %>%
    left_join(exprs_df, by = "sample_id") %>% 
    left_join(genotype_df, by = "genotype_id") %>%
    dplyr::mutate(condition_name = factor(condition_name, levels = c("Ctrl","AcLDL")))
  
  plot = ggplot(plot_df, aes(x = genotype_value, y = norm_exp)) + 
    facet_wrap(~ condition_name) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(position = position_jitter(width = .1)) + 
    ylab("Normalized expression") +
    xlab(genotype_id) + 
    labs(title = gene_name)
  
  return(plot)
}
makeMultiplePlots <- function(snps_df, expression_dataset, genotype_dataset, line_metadata){
  #Plot eQTL results for a list of gene and SNP pairs.
  result = list()
  for(i in 1:nrow(snps_df)){
    gene_id = snps_df[i,]$gene_id
    snp_id = snps_df[i,]$snp_id
    print(gene_id)
    plot = plotEQTL(gene_id, snp_id, expression_dataset, genotype_dataset, line_metadata)
    result[[gene_id]] = plot 
  }
  return(result)
}

makeMultiplePlots_acLDL <- function(snps_df, expression_dataset, genotype_dataset, line_metadata){
  #Plot eQTL results for a list of gene and SNP pairs.
  result = list()
  for(i in 1:nrow(snps_df)){
    gene_id = snps_df[i,]$gene_id
    snp_id = snps_df[i,]$snp_id
    print(gene_id)
    plot = plotEQTL_acLDL(gene_id, snp_id, expression_dataset, genotype_dataset, line_metadata)
    result[[gene_id]] = plot 
  }
  return(result)
}

savePlots <- function(plot_list, path, width, height){
  #Save a list of plots into the folder specified by path
  for (plot in plot_list){
    gene_name = plot$labels$title
    file_name = file.path(path, paste(gene_name, ".pdf", sep = ""))
    ggsave(file_name, plot, width = width, height = height)
  }
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