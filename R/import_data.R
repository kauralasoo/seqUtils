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

loadChrCounts <- function(sample_dir, sample_names, counts_suffix = ".chr_counts"){
  #Import read counts per chromosome for each sample (used for estimating MT fraction in ATAC-Seq)
  matrix = c()
  for (i in c(1:length(sample_names))){
    path = file.path(sample_dir, sample_names[i], paste(sample_names[i], counts_suffix, sep = ""))
    print(path)
    table = read.table(path, header = FALSE, stringsAsFactors = FALSE)
    colnames(table) = c(sample_names[i], "chr_name")
    if (i == 1){
      matrix = table[,c(2,1)]
      print(head(matrix))
    }
    else{
      matrix = dplyr::full_join(matrix, table, by = "chr_name")
      print(head(matrix))
    }
  }
  return(matrix)
}

loadFeaturCountsSummary <- function(sample_dir, sample_names, counts_suffix = ".counts.txt.summary"){
  #Import featureCounts summary files
  matrix = c()
  for (i in c(1:length(sample_names))){
    path = file.path(sample_dir, sample_names[i], paste(sample_names[i], counts_suffix, sep = ""))
    print(path)
    table = read.table(path, header = TRUE, stringsAsFactors = FALSE)
    colnames(table) = c("status", sample_names[i])
    if (i == 1){
      matrix = table[,c(1,2)]
      print(head(matrix))
    }
    else{
      matrix = dplyr::full_join(matrix, table, by = "status")
      print(head(matrix))
    }
  }
  rownames(matrix) = matrix$status
  matrix = matrix[,-1]
  matrix = as.data.frame(t(matrix))
  matrix = dplyr::mutate(matrix, sample_id = rownames(matrix)) %>% dplyr::select(sample_id, everything())
  return(matrix)
}

loadMarkDuplicates <- function(sample_dir, sample_names, counts_suffix = ".MarkDuplicates.txt"){
  #Import MarkDuplicates summary files
  matrix = c()
  for (i in c(1:length(sample_names))){
    path = file.path(sample_dir, sample_names[i], paste(sample_names[i], counts_suffix, sep = ""))
    print(path)
    data = read.table(path, skip = 6, nrows = 2, stringsAsFactors = FALSE)
    table = as.data.frame(t(data),stringsAsFactors = FALSE)
    colnames(table) = c("statistic", sample_names[i])
    if (i == 1){
      matrix = table[,c(1,2)]
      print(head(matrix))
    }
    else{
      matrix = dplyr::full_join(matrix, table, by = "statistic")
      print(head(matrix))
    }
  }
  rownames(matrix) = matrix$statistic
  matrix = matrix[,-1]
  matrix = as.data.frame(t(matrix), stringsAsFators = FALSE)
  matrix = dplyr::mutate(matrix, sample_id = rownames(matrix)) %>% dplyr::select(sample_id, everything())
  return(matrix)
}

loadNarrowPeaks <- function(sample_dir, sample_names, peaks_suffix = "_peaks.narrowPeak"){
  #Import narrowPeak files into a list
  result = list()
  extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
  for (i in c(1:length(sample_names))){
    path = file.path(sample_dir, sample_names[i], paste(sample_names[i], peaks_suffix, sep = ""))
    print(path)
    peaks = rtracklayer::import(path, format = "BED", extraCols = extraCols_narrowPeak)
    result[[sample_names[i]]] = peaks
  }
  return(result)
}

loadBroadPeaks <- function(sample_dir, sample_names, peaks_suffix = "_peaks.broadPeak"){
  #Import broakPeak files into a list
  result = list()
  extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric")
  for (i in c(1:length(sample_names))){
    path = file.path(sample_dir, sample_names[i], paste(sample_names[i], peaks_suffix, sep = ""))
    print(path)
    peaks = rtracklayer::import(path, format = "BED", extraCols = extraCols_narrowPeak)
    result[[sample_names[i]]] = peaks
  }
  return(result)
}

#' Import genotypes in GDS file to a matrix
#'
#' Open GDS file (created with SNPRelate::snpgdsVCF2GDS) into R and convert it into a matrix of 
#' variant positions and matrix of allele dosage (0,1,2)
#' 
#' @param file Path to the GDS file.
#' @return List containing snpspos and genotypes matrix.
#' @author Kaur Alasoo
#' @export 
gdsToMatrix <- function(gds_file){
  
  #Extract genotypes
  gds <- GWASTools::GdsGenotypeReader(gds_file)
  genotypes = GWASTools::getGenotype(gds)
  sample_ids = GWASTools::getVariable(gds, "sample.id")
  snp_ids = GWASTools::getVariable(gds, "snp.rs.id")
  colnames(genotypes) = sample_ids
  rownames(genotypes) = snp_ids
  
  #Extract SNP coordinates
  snpspos = data_frame(snpid = GWASTools::getVariable(gds, "snp.rs.id"), 
             chr = GWASTools::getVariable(gds, "snp.chromosome"), 
             pos = GWASTools::getVariable(gds, "snp.position"))
  GWASTools::close(gds)
  
  return(list(snpspos = snpspos, genotypes = genotypes))
}
