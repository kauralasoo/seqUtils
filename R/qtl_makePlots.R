#' Make a boxplot of expression QTL in multiple conditions.
#'
#' Plots are faceted by condition_name column in sample_metadata.
#' 
#' @param selected_gene_id ID if the gene of interest.
#' @param gentype_id ID of the variant of interest.
#' @param expression_matrix Matrix of normalised gene expression values (genes in rows, samples in columns). 
#' Column names correspond to sample_id in sample_metadata. 
#' @param genotype_matrix Matrix of genotypes (variants in rows, individuals in columns). 
#' Column names correspond to genotype_id in sample_metadata.
#' @param sample_metadata data frame linking samples to genotypes (Required columns: sample_id, genotype_id, condition_name)
#' @param gene_metadata data frame for linking gene ids to gene names. (Required columns: gene_id, gene_name) 
#' @return ggplot2 object
#' @author Kaur Alasoo
#' @export 
plotEQTL <- function(selected_gene_id, genotype_id, expression_matrix, genotype_matrix, sample_metadata, 
                     gene_metadata, return_df = FALSE){
  
  #Extraxt gene_name
  gene_name = dplyr::filter(gene_metadata, gene_id == selected_gene_id)$gene_name
  print(gene_name)
  
  #extract genotypes
  geno_vector = genotype_matrix[genotype_id,]
  genotype_df = data_frame(genotype_id = names(geno_vector), genotype_value = as.character(geno_vector))
  
  #expression
  expression_vector = expression_matrix[selected_gene_id,]
  exprs_df = dplyr::data_frame(sample_id = names(expression_vector), norm_exp = expression_vector)
  
  #Join all of the data together
  g_id = genotype_id
  plot_df = dplyr::left_join(sample_metadata, exprs_df, by ="sample_id") %>%
    dplyr::left_join(genotype_df, by = "genotype_id") %>%
    dplyr::mutate(snp_id = g_id, gene_name = gene_name)
  
  #Either return raw data frame or make plot and return ggplot2 object
  if (return_df){
    return(plot_df)
  } else{
    plot = ggplot2::ggplot(plot_df, ggplot2::aes(x = genotype_value, y = norm_exp)) + 
      ggplot2::facet_wrap(~ condition_name) + 
      ggplot2::geom_boxplot(outlier.shape = NA) + 
      ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .1)) + 
      ggplot2::ylab("Normalized expression") +
      ggplot2::xlab(genotype_id) + 
      ggplot2::labs(title = gene_name)
    
    return(plot)
  }
}

plotQtlRow <- function(qtl_df){
  plot = ggplot2::ggplot(qtl_df, ggplot2::aes(x = genotype_value, y = norm_exp)) + 
    ggplot2::facet_grid(~condition_name) + 
    ggplot2::geom_boxplot(outlier.shape = NA) + 
    ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .2), size = 0.5) + 
    ggplot2::ylab("Normalized expression") +
    ggplot2::xlab(qtl_df$snp_id[1]) + 
    ggplot2::labs(title = qtl_df$gene_name[1])
  return(plot)
}

#' Make a list of plotEQTL plots.
#'
#' Plots are faceted by condition_name column in sample_metadata.
#' 
#' @param snps_df Data frame with at least two columns (gene_id, snp_id) corresponding to snps and genes to be plotted.
#' @param expression_matrix Matrix of normalised gene expression values (genes in rows, samples in columns). 
#' Column names correspond to sample_id in sample_metadata. 
#' @param genotype_matrix Matrix of genotypes (variants in rows, individuals in columns). 
#' Column names correspond to genotype_id in sample_metadata.
#' @param sample_metadata data frame linking samples to genotypes (Required columns: sample_id, genotype_id, condition_name)
#' @param gene_metadata data frame for linking gene ids to gene names. (Required columns: gene_id, gene_name) 
#' @return List of ggplot2 objects.
#' @author Kaur Alasoo
#' @export 
makeMultiplePlots <- function(snps_df, expression_matrix, genotype_matrix, sample_metadata, gene_metadata){
  #Plot eQTL results for a list of gene and SNP pairs.
  result = list()
  for(i in 1:nrow(snps_df)){
    gene_id = snps_df[i,]$gene_id
    snp_id = snps_df[i,]$snp_id
    print(gene_id)
    plot = plotEQTL(gene_id, snp_id, expression_matrix, genotype_matrix, sample_metadata, gene_metadata)
    plot_name = paste(gene_id, "-",snp_id, sep = "")
    result[[plot_name]] = plot 
  }
  return(result)
}

#' Save a list of ggplot2 plots into a folder
#'
#' File name is extracted automatically from the plot title.
#' 
#' @param plot_list List of ggplot2 objects.
#' @param path Path to the output directory, created automatically if does not exist.
#' @param width Width of the plot.
#' @param height Height of the plot.
#' @return None
#' @author Kaur Alasoo
#' @export 
savePlots <- function(plot_list, path, width, height){
  #Save a list of plots into the folder specified by path
  for (plot in plot_list){
    gene_name = plot$labels$title
    snp_id = plot$labels$x
    if (!dir.exists(path)){ dir.create(path) } #Create dir if not there
    file_name = file.path(path, paste(gene_name, "-",snp_id, ".pdf", sep = ""))
    ggplot2::ggsave(file_name, plot, width = width, height = height)
  }
}

savePlotList <- function(plot_list, output_folder, suffix = ".pdf", ...){
  for (name in names(plot_list)){
    path = file.path(output_folder, paste0(name, suffix))
    ggsave(path, plot_list[[name]], ...)
  }
}