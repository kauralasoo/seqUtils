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

loadVerifyBamID <- function(sample_names, sample_dir, suffix = ".verifyBamID.bestSM"){
  matrix = c()
  for (i in c(1:length(sample_names))){
    path = file.path(sample_dir, sample_names[i], paste(sample_names[i], suffix, sep = ""))
    table = read.table(path, comment.char = "", sep ="\t", header = TRUE) %>%
      dplyr::select(CHIP_ID, FREEMIX) %>%
      dplyr::mutate(sample_id = sample_names[i]) %>%
      dplyr::rename(genotype_id = CHIP_ID, freemix = FREEMIX)
    matrix = rbind(matrix, table)
  }
  return(matrix)
}
