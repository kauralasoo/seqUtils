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