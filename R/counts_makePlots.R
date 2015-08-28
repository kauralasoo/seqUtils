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