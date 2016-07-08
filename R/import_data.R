#' Import featureCounts output table into R.
#'
#' Skips the first comment line.
#' 
#' @param sample_dir Path to the directory containing the count files.
#' @param sample_names Vector of sample names.
#' @param counts_suffix Suffix of the counts file.
#' @param sub_dir If TRUE, count files are nested in subfolders named after sample names, otherwise counts files
#' are directly in sample_dir.
#' @return Counts matrix where the first two columns are gene_id and feature length.
#' @author Kaur Alasoo
#' @export 
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
    table = readr::read_tsv(path, skip = 1, col_types = "cccccii")
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

loadExonCounts <- function(sample_dir, sample_names, counts_suffix = ".exon_counts.txt", sub_dir = TRUE){
  #Load featureCounts output into R
  matrix = c()
  for (i in c(1:length(sample_names))){
    if (sub_dir == TRUE){
      path = file.path(sample_dir, sample_names[i], paste(sample_names[i], counts_suffix, sep = ""))
    } else {
      path = file.path(sample_dir, paste(sample_names[i], counts_suffix, sep = ""))      
    }
    print(sample_names[i])
    table = readr::read_tsv(path, skip = 1, col_types = "cccccii")
    print(head(table))
    if (i == 1){
      matrix = table[,c(1,3,4,6,7)]
    }
    else{
      matrix = cbind(matrix, table[,7])
    }
  }
  colnames(matrix) = c("gene_id", "start", "end", "length", sample_names)
  matrix = dplyr::tbl_df(matrix)
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

#' Import narrowPeak files into a list of GRanges objects.
#'
#' All narrowPeak files are assumed to be nested in subdirectories named adter sample names.
#' Based on rtracklayer::import.bed, but adds additional columns present only in narrowPeak files.
#' 
#' @param sample_dir Path to the directory containing the narroPeak files.
#' @param sample_names Vector of sample names.
#' @param peaks_suffix Suffix of the narrowPeak files.
#' @return List of GRanges objects corresponding to peak calls from each sample.
#' @param sub_dir If TRUE, count files are nested in subfolders named after sample names, otherwise counts files
#' are directly in sample_dir.
#' @author Kaur Alasoo
#' @export 
loadNarrowPeaks <- function(sample_dir, sample_names, peaks_suffix = "_peaks.narrowPeak", sub_dir = TRUE){
  #Import narrowPeak files into a list
  result = list()
  extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
  for (i in c(1:length(sample_names))){
    if (sub_dir == TRUE){
      path = file.path(sample_dir, sample_names[i], paste(sample_names[i], peaks_suffix, sep = ""))
    } else {
      path = file.path(sample_dir, paste(sample_names[i], peaks_suffix, sep = ""))
    }
    print(path)
    peaks = rtracklayer::import(path, format = "BED", extraCols = extraCols_narrowPeak)
    result[[sample_names[i]]] = peaks
  }
  return(result)
}

#' Import broadPeak files into a list of GRanges objects.
#'
#' All broadPeak files are assumed to be nested in subdirectories named adter sample names.
#' Based on rtracklayer::import.bed, but adds additional columns present only in broadPeak files.
#' 
#' @param sample_dir Path to the directory containing the broadPeak files.
#' @param sample_names Vector of sample names.
#' @param peaks_suffix Suffix of the broadPeak files.
#' @param sub_dir If TRUE, count files are nested in subfolders named after sample names, otherwise counts files
#' are directly in sample_dir.
#' @return List of GRanges objects corresponding to peak calls from each sample.
#' @author Kaur Alasoo
#' @export 
loadBroadPeaks <- function(sample_dir, sample_names, peaks_suffix = "_peaks.broadPeak", sub_dir = TRUE){
  #Import broakPeak files into a list
  result = list()
  extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric")
  for (i in c(1:length(sample_names))){
    if (sub_dir == TRUE){
      path = file.path(sample_dir, sample_names[i], paste(sample_names[i], peaks_suffix, sep = ""))
    } else {
      path = file.path(sample_dir, paste(sample_names[i], peaks_suffix, sep = ""))
    }
    print(path)
    peaks = rtracklayer::import(path, format = "BED", extraCols = extraCols_broadPeak)
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
  snp_rs_ids = GWASTools::getVariable(gds, "snp.rs.id")
  snp_ids = GWASTools::getVariable(gds, "snp.id")
  
  #Invent id for snps that do not have an rs id
  new_snp_ids = paste("snp",snp_ids[snp_rs_ids == ""], sep = "")
  snp_rs_ids[snp_rs_ids == ""] = new_snp_ids
  colnames(genotypes) = sample_ids
  rownames(genotypes) = snp_rs_ids
  
  #Extract SNP coordinates
  snpspos = dplyr::data_frame(snpid = snp_rs_ids, 
             chr = GWASTools::getVariable(gds, "snp.chromosome"), 
             pos = GWASTools::getVariable(gds, "snp.position"))
  GWASTools::close(gds)
  return(list(snpspos = snpspos, genotypes = genotypes))
}

#' Import output files from bedCountFragmentLengths.py script into R.
#'
#' All fragment length files are assumed to be nested in subdirectories named after sample names.
#' 
#' @param sample_dir Path to the directory containing the fragment length files.
#' @param sample_names Vector of sample names.
#' @param file_suffix Suffix of the fragment lengths files.
#' @return List of GRanges objects corresponding to peak calls from each sample.
#' @author Kaur Alasoo
#' @export 
loadFragmentLengths <- function(sample_dir, sample_names, file_suffix = ".fragment_lengths.txt"){
  result = c()
  for (i in c(1:length(sample_names))){
    path = file.path(sample_dir, sample_names[i], paste(sample_names[i], file_suffix, sep = ""))
    table = read.table(path) %>% dplyr::tbl_df()
    colnames(table) = c("count", "fragment_length")
    sample_table = dplyr::mutate(table, sample_id = sample_names[i])
    if (i == 1){
      result = sample_table
    }
    else{
      result = rbind(result, sample_table)
    }
  }
  return(result)
}

#' A general function to quickly import tabix indexed tab-separated files into data_frame
#'
#' @param tabix_file Path to tabix-indexed text file
#' @param param A instance of GRanges, RangedData, or RangesList 
#' provide the sequence names and regions to be parsed. Passed onto Rsamtools::scanTabix()
#' @param ... Additional parameters to be passed on to readr::read_delim()
#'
#' @return List of data_frames, one for each entry in the param GRanges object.
#' @export
scanTabixDataFrame <- function(tabix_file, param, ...){
  tabix_list = Rsamtools::scanTabix(tabix_file, param = param)
  df_list = lapply(tabix_list, function(x,...){
    if(length(x) > 0){
      if(length(x) == 1){
        #Hack to make sure that it also works for data frames with only one row
        #Adds an empty row and then removes it
        result = paste(paste(x, collapse = "\n"),"\n",sep = "")
        result = readr::read_delim(result, delim = "\t", ...)[1,]
      }else{
        result = paste(x, collapse = "\n")
        result = readr::read_delim(result, delim = "\t", ...)
      }
    } else{
      #Return NULL if the nothing is returned from tabix file
      result = NULL
    }
    return(result)
  }, ...)
  return(df_list)
}

#' Import variant information extracted from VCF file into R
importVariantInformation <- function(path){
  info_col_names = c("chr","pos","snp_id","ref","alt","type","AC","AN")
  into_col_types = "ciccccii"
  snp_info = readr::read_delim(path, delim = "\t", col_types = into_col_types, col_names = info_col_names)
  snp_info = dplyr::mutate(snp_info, indel_length = pmax(nchar(alt), nchar(ref))) %>%
    dplyr::mutate(is_indel = ifelse(indel_length > 1, TRUE, FALSE)) %>%
    dplyr::mutate(MAF = pmin(AC/AN, 1-(AC/AN)))
  return(snp_info)
}

importGWASSummaryStats <- function(region, path){
  gwas_colnames = c("chr","pos","pos2","snp_id","A1","A2", "INFO","OR","SE","p_nominal")
  gwas_coltypes = "ciicccdddd"
  gwas_pvalues = scanTabixDataFrame("databases/GWAS/IGAP_summary_statistics/IGAP_stage_1.GRCh38.sorted.bed.gz", region, 
                                   col_names = gwas_colnames, col_types = gwas_coltypes)
  gwas_pvalues = purrr::map(gwas_pvalues, ~dplyr::select(.,-pos2))
  return(gwas_pvalues)
}

importMotifDisruptions <- function(path){
  motif_colnames = c("gene_id","snp_id","snp_count",".row","start","strand","motif_id","ref_abs_score","ref_rel_score",
                     "ref_match","alt_abs_score","alt_rel_score","alt_match","rel_diff","max_rel_score")
  motif_disruptions = readr::read_delim(path, delim = "\t", col_names = motif_colnames)
  return(motif_disruptions)
}

