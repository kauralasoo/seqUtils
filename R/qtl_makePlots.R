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
                     gene_metadata){
  
  #Construct plot data frame
  plot_df = constructQtlPlotDataFrame(selected_gene_id, genotype_id, expression_matrix, genotype_matrix, sample_metadata, 
                                      gene_metadata)
  
  #Make plot and return ggplot2 object
  plot = ggplot2::ggplot(plot_df, ggplot2::aes(x = genotype_value, y = norm_exp)) + 
    ggplot2::facet_wrap(~ condition_name) + 
    ggplot2::geom_boxplot(outlier.shape = NA) + 
    ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .1)) + 
    ggplot2::ylab("Normalized expression") +
    ggplot2::xlab(genotype_id) + 
    ggplot2::labs(title = plot_df$gene_name[[1]])
    
  return(plot)
}

constructQtlPlotDataFrame <- function(selected_gene_id, genotype_id, expression_matrix, genotype_matrix, sample_metadata, 
                                      gene_metadata){
  
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
  
  return(plot_df)
}

plotQtlRow <- function(qtl_df, ylabel = "Normalized expression"){
  plot = ggplot2::ggplot(qtl_df, ggplot2::aes(x = genotype_value, y = norm_exp, color = condition_name)) + 
    ggplot2::facet_grid(~condition_name) + 
    ggplot2::geom_boxplot(outlier.shape = NA) + 
    ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .2), size = 0.5) + 
    ggplot2::ylab(ylabel) +
    ggplot2::xlab(qtl_df$snp_id[1]) + 
    ggplot2::labs(title = qtl_df$gene_name[1]) + 
    ggplot2::theme_light() + 
    ggplot2::scale_color_manual(values = conditionPalette(), guide=FALSE)
  return(plot)
}

plotQtlCol <- function(qtl_df){
  plot = ggplot2::ggplot(qtl_df, ggplot2::aes(x = genotype_value, y = norm_exp, color = condition_name)) + 
    ggplot2::facet_wrap(~condition_name, ncol = 1) + 
    ggplot2::geom_boxplot(outlier.shape = NA) + 
    ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .2), size = 0.5) + 
    ggplot2::ylab("Normalized expression") +
    ggplot2::xlab(qtl_df$snp_id[1]) + 
    ggplot2::labs(title = qtl_df$gene_name[1]) + 
    ggplot2::theme_light() + 
    ggplot2::scale_color_manual(values = conditionPalette(), guide=FALSE)
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


#' Construct metadata data frame for wiggleplotr 
#'
#' @param count_matrix Feature count matrix used to for calculating library size
#' @param sample_metadata Sample metadata matrix (required columns: sample_id, genotype_id, condition_name)
#' @param bigWig_dir Path to the directory with bigWig files
#' @param condition_name_levels Levels of the condition name factor
#'
#' @return Wiggleplotr metadata data frame
#' @export
wiggleplotrConstructMetadata <- function(count_matrix, sample_metadata, bigWig_dir, bigWig_suffix = ".bw", 
                                         condition_name_levels = c("naive","IFNg", "SL1344", "IFNg_SL1344") ){
  assertthat::assert_that(assertthat::has_name(sample_metadata, "sample_id"))
  assertthat::assert_that(assertthat::has_name(sample_metadata, "genotype_id"))
  assertthat::assert_that(assertthat::has_name(sample_metadata, "condition_name"))
  
  
  #Calculate library sizes
  library_sizes = data_frame(sample_id = colnames(count_matrix), scaling_factor = colSums(count_matrix)/1e6)
  
  #Make a df with metadata
  plotting_meta = sample_metadata %>%
    dplyr::select(sample_id, genotype_id, condition_name) %>%
    dplyr::mutate(bigWig = file.path(bigWig_dir, paste(sample_id, bigWig_suffix, sep = ""))) %>%
    dplyr::mutate(track_id = factor(condition_name, levels = condition_name_levels)) %>%
    dplyr::left_join(library_sizes, by = "sample_id")
  return(plotting_meta)
}

wiggleplotrGenotypeColourGroup <- function(metadata, variant_id, genotype_matrix, beta_sign){
  
  #Set the correct sign for the colour group levels
  if(beta_sign >= 0){
    colour_group_levels = c(2,1,0)
  } else{
    colour_group_levels = c(0,1,2)
  }
  
  #Extract genotypes from the genotype matrix
  genotype_df = vcf_file$genotypes[variant_id,] %>% 
    tidyVector(sample_id = "genotype_id", value_id = "colour_group") %>%
    dplyr::mutate(colour_group = factor(colour_group, levels = colour_group_levels))
  
  #Add colour group to the metadata df
  new_meta = dplyr::left_join(metadata, genotype_df, by = "genotype_id")
  return(new_meta)
}

wiggpleplotrConstructPeakAnnotations <- function(selected_peaks){
  peak_list = list(ATAC = dplyr::transmute(selected_peaks, seqnames = chr, start, end, strand) %>% dataFrameToGRanges())
  peak_annot = data_frame(transcript_id = "ATAC", gene_id = "ATAC", gene_name = "ATAC-seq", strand = "+")
  
  return(list(peak_list = peak_list, peak_annot = peak_annot))
}



