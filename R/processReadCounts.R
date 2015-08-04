calculateTPM <- function(counts_matrix, lengths, selected_genes = NULL, fragment_length = 250){
  #Normalize gene expression using TPM method
  
  if (is.null(selected_genes)){
    selected_genes = rownames(counts_matrix)
  }
  
  #Add rownames to lengths df 
  rownames(lengths) = lengths$gene_id
  
  #Calculated scaling factor
  selected_counts = counts_matrix[intersect(rownames(counts_matrix), selected_genes),]
  selected_lengths = lengths[rownames(selected_counts),]$length
  scaling_factor = colSums((selected_counts*fragment_length)/selected_lengths)
  
  #Calculate TPM for each gene
  length_vector = lengths[rownames(counts_matrix),]$length
  tpm = t(t(counts_matrix * fragment_length * 1e6)/scaling_factor)
  tpm = tpm/length_vector
  
  return(tpm)
}

calculateCQN <- function(counts_matrix, gene_metadata){
  #Normalize read counts using the CQN method.
  expression_cqn = cqn(counts = counts_matrix[gene_metadata$gene_id,], x = gene_metadata$percentage_gc_content, 
                       lengths = gene_metadata$length, verbose = TRUE)
  expression_norm = expression_cqn$y + expression_cqn$offset
  return(expression_norm)
}

loadCounts <- function(sample_dir, sample_names, counts_suffix = ".counts.txt", sub_dir = TRUE){
  #Load featureCounts output into R
  matrix = c()
  for (i in c(1:length(sample_names))){
    if (sub_dir == TRUE){
      path = file.path(sample_dir, sample_names[i], paste(sample_names[i], counts_suffix, sep = ""))
    } else {
      path = file.path(sample_dir, paste(sample_names[i], counts_suffix, sep = ""))      
    }
    print(sample_names[i])
    table = read.table(path, header = TRUE)
    print(head(table))
    if (i == 1){
      matrix = table[,c(1,6,7)]
    }
    else{
      matrix = cbind(matrix, table[,7])
    }
  }
  colnames(matrix) = c("gene_id", "length", sample_names)
  return(matrix)
}

zScoreNormalize <- function(matrix){
  #Normalize expression matrix by z-score
  matrix = matrix - rowMeans(matrix)
  matrix = matrix / apply(matrix, 1, sd)
  return(matrix)
}

performPCA <- function(matrix, design, ...){
  #Perform PCA of gene expression matrix add experimental design metadata to the results
  pca = prcomp(t(matrix), ...)
  pca_matrix = as.data.frame(pca$x) %>% 
    dplyr::mutate(sample_id = rownames(pca$x)) %>%
    dplyr::left_join(design, by = "sample_id")
  return(list(pca_matrix = pca_matrix, pca_object = pca))
}

calculateMean <- function(matrix, design, factor, sample_id_col = "sample_id"){
  #Calculate the mean value in matrix over all possible factor values.
  
  #If the factor is not a factor then make it a factor.
  if(!is.factor(design[,factor])){
    design[,factor] = factor(design[,factor])
  }
  
  #Set sample_id column as rownames
  rownames(design) = design[,sample_id_col]
  factor = design[,factor]
  levs = levels(factor)
  result = c()
  for (lev in levs){
    filter = factor == lev
    samples = rownames(design[filter,])
    mat = matrix[,samples]
    mat = rowMeans(mat)
    result = cbind(result, mat)
  }
  colnames(result) = levs
  return(data.frame(result))
}

plotGene <- function(gene_id, matrix, design, gene_metadata, colors = c("#d95f02","#1b9e77","#7570b3")){
  #Plot the expression values of each gene in different conditions
  matrix = matrix[,match(rownames(design), colnames(matrix))]
  gene_expression = matrix[gene_id,]
  print(gene_expression)
  design$expression = as.numeric(gene_expression)
  gene_name = as.vector(gene_metadata[gene_metadata$gene_id == gene_id,]$gene_name)
  
  #Plot results
  plot = ggplot(design, aes(x = condition_name, y = expression)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(position = position_jitter(width = .1)) + 
    ylab(expression(paste(Log[2], " expression"))) + 
    ggtitle(gene_name) +
    #scale_fill_manual(values = colors) + 
    theme(legend.position="none", text = element_text(size=20), axis.text.x = element_text(angle = 20), axis.title.x = element_blank())
  return(plot)
}

filterDESeqResults <- function(results,gene_metadata, min_padj = 0.01, min_fc = 1, biotype_filter = NULL){
  #Add gene name to the DESeq result and filter out up and downregulated genes.
  
  #Construct a results table
  result_table = results %>% 
    data.frame() %>% 
    dplyr::mutate(gene_id = rownames(results)) %>% 
    tbl_df() %>% 
    dplyr::left_join(gene_metadata, by = "gene_id") %>% 
    dplyr::arrange(padj)
  
  #Find up and down-regulated genes
  up_genes = dplyr::filter(result_table, padj < min_padj, log2FoldChange > min_fc) %>% 
    arrange(-log2FoldChange)
  down_genes = dplyr::filter(result_table, padj < min_padj, log2FoldChange < -min_fc) %>% 
    arrange(log2FoldChange)
  
  #Filter up and down-regulated genes by biotype
  if(!is.null(biotype_filter)){
    up_genes = dplyr::filter(up_genes, gene_biotype == biotype_filter)
    down_genes = dplyr::filter(down_genes, gene_biotype == biotype_filter)
  }
  return(list(up_genes = up_genes, down_genes = down_genes, results_table = result_table))
}

calculateNormFactors <- function(counts_matrix, method = "RLE", output = "rasqual"){
  #Calculate norm factors for a counts matrix
  dge = edgeR::DGEList(counts = counts_matrix)
  dge = edgeR::calcNormFactors(dge, method = method)
  sample_info = dge$samples[,-1]
  if (output == "rasqual"){
    size_matrix = matrix(rep(sample_info$norm.factors, nrow(counts_matrix)), nrow = nrow(counts_matrix), byrow = TRUE)
    rownames(size_matrix) = rownames(counts_matrix)
    return(size_matrix)
  }else{
    return(sample_info)
  }
}

loadIntronEventCounts <- function(sample_dir, sample_names, counts_suffix = ".intron_events.txt", sub_dir = TRUE){
  #Load featureCounts output into R
  matrix = c()
  for (i in c(1:length(sample_names))){
    if (sub_dir == TRUE){
      path = file.path(sample_dir, sample_names[i], paste(sample_names[i], counts_suffix, sep = ""))
    } else {
      path = file.path(sample_dir, paste(sample_names[i], counts_suffix, sep = ""))      
    }
    print(sample_names[i])
    table = readr::read_tsv(path, col_names = FALSE, col_types = "ccc")
    if (i == 1){
      matrix = table
    }
    else{
      matrix = cbind(matrix, table[,3])
    }
  }
  colnames(matrix) = c("gene_id", "event_id", sample_names)
  return(matrix)
}