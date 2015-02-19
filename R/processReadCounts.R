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

loadCounts <- function(sample_dir, sample_names, counts_suffix = ".counts.txt" ){
  #Load featureCounts output into R
  matrix = c()
  for (i in c(1:length(sample_names))){
    path = file.path(sample_dir, sample_names[i], paste(sample_names[i], counts_suffix, sep = ""))
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
