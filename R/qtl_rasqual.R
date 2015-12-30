#' Save a list of matrices into a suitable format for RASQUAL
#'
#' Works with expression and covariates matrices.
#' 
#' @param data_list list of matrices
#' @param output_dir relative path to the output dir
#' @param file_suffix suffix added to each file after their name in the list.
#' @return None
#' @author Kaur Alasoo
#' @export 
saveRasqualMatrices <- function(data_list, output_dir, file_suffix = "expression"){
  #Save data for FastQTL to disk
  
  #Save each matrix as a separate  txt file
  for (sn in names(data_list)){
    file_path = file.path(output_dir, paste(sn,file_suffix, "txt", sep = "."))
    file_path_bin = file.path(output_dir, paste(sn,file_suffix, "bin", sep = "."))
    print(file_path)
    write.table(data_list[[sn]], file_path, quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")
    writeBin(as.double(c(t(data_list[[sn]]))), file_path_bin)
  }
}


rasqualSizeFactorsMatrix <- function(expression_list, factor_name){
  factor_vector = expression_list$norm_factor[,factor_name]
  size_matrix = matrix(rep(factor_vector, nrow(expression_list$counts)), nrow = nrow(expression_list$counts), byrow = TRUE)
  rownames(size_matrix) = rownames(expression_list$counts)
  return(size_matrix)
}

#' Count the number of feature SNPs and cis SNPs overlapping a set of peak calls
#'
#' @param peak_metadata Data frame with peak metadata (required columns: gene_id, chr, start, end)
#' @param snp_coords Data frame with SNP coordinates from a VCF file (required columns: chr, pos, snp_id)
#' @param cis_window Size of the cis window from both sides of the peak.
#'
#' @return Data frame with peak coordinates, cis region coordiantes as well as number of cis and feature snps.
#' @export
countSnpsOverlapingPeaks <- function(peak_metadata, snp_coords, cis_window = 500000){
  
  #Construct peak coords data frame
  peak_coords = dplyr::transmute(peak_metadata, gene_id, chromosome_name = chr, strand = 1, 
                                 exon_starts = start, exon_ends = end) %>%
    dplyr::mutate(range_start = pmax(0, exon_starts - cis_window), range_end = exon_ends + cis_window)
  
  #Construct GRanges objects
  peak_granges = GenomicRanges::GRanges(seqnames = peak_coords$chromosome_name, 
                         ranges = IRanges::IRanges(start = peak_coords$exon_starts, end = peak_coords$exon_ends))
  region_granges = GenomicRanges::GRanges(seqnames = peak_coords$chromosome_name, 
                           ranges = IRanges::IRanges(start = peak_coords$range_start, end = peak_coords$range_end))
  snp_granges = GenomicRanges::GRanges(seqnames = snp_coords$chr, ranges = IRanges::IRanges(start = snp_coords$pos, end = snp_coords$pos))
  
  #Count overlaps
  feature_snps = GenomicRanges::countOverlaps(peak_granges, snp_granges)
  cis_snps = GenomicRanges::countOverlaps(region_granges, snp_granges)
  
  new_peak_coords = dplyr::mutate(peak_coords, feature_snp_count = feature_snps, cis_snp_count = cis_snps)
  return(new_peak_coords)
}

#' Import rasqual output table into R
#'
#' Skipped gene-SNP pairs are automatically removed and 
#' chisq statistic is converted into p-value (p_nominal).
#'
#' @param path Bath to rasqual output file.
#'
#' @return data_frame
#' @export
importRasqualTable <- function(path){
  rasqual_results = readr::read_delim(path, delim = "\t", col_types = "cccidddiii",col_names = FALSE)
  colnames(rasqual_results) = c("gene_id", "snp_id", "chr", "pos", "allele_freq", "chisq", "pi", "n_feature_snps", "n_cis_snps", "converged")
  
  rasqual_pvalues = dplyr::filter(rasqual_results, snp_id != "SKIPPED") %>%
    dplyr::mutate(p_nominal = pchisq(chisq, df = 1, lower = FALSE))
  
  return(rasqual_pvalues)
}

#' Helper function for rasqualGcCorrection
#' 
#' @author Natsuhiko Kumasaka
Quantile <- function(x,k=20){
  x=rank(x,ties="random")
  z=rep(0,length(x))
  for(i in 1:k){
    z = z+as.numeric(x<=quantile(x,i/k,na.rm=T))
  }
  k-z
}

#' Correct rasqual size factors matrix for gc content
#'
#' @param Y Matrix of read counts
#' @param gcvec Vector of GC content percentages
#' @param PLOT 
#'
#' @return Matrix of GC-corrected library sizes.
#' @export
#' @author Natsuhiko Kumasaka
rasqualGcCorrection <- function(Y,gcvec,PLOT=F){
  bin=Quantile(gcvec,200);
  x=sort(unlist(lapply(split(gcvec,bin),mean)))
  S=apply(Y,2,function(y){unlist(lapply(split(y,bin),sum))[as.character(0:199)]});
  Fs=log(t(t(S)/apply(S,2,sum))/apply(S,1,sum)*sum(S));
  Gs=apply(Fs,2,function(y){smooth.spline(x,y,spar=1)$y}); 
  if(PLOT){
    par(mfcol=c(5,5),mar=c(2,2,2,2)); 
    for(i in 1:ncol(Y)){
      plot(Fs[,i])
      lines(Gs[,i],col=2)
    }
    matplot(x,Gs,type="l",col=2,lty=1)
  }
  exp(Gs[bin+1,])
}

#' Split gene ids into batches for runRasqual.py script
#'
#' @param gene_metadata Data frame with gene metadata (gene_id column required)
#' @param batch_size Number of genes in a batch
#'
#' @return Data frame of gene batches (columns: batch_id, gene_ids)
#' @export
rasqualConstructGeneBatches <- function(gene_metadata, batch_size){
  batch_df = dplyr::select(gene_metadata, gene_id) %>%
    dplyr::mutate(batch_number = splitIntoBatches(length(gene_id), batch_size)) %>%
    dplyr::mutate(batch_id = paste("batch", batch_number, sep = "_")) %>% 
    dplyr::group_by(batch_id) %>%
    dplyr::summarize(gene_ids = paste(gene_id, collapse = ","))
  return(batch_df)
}

#' Find SNP with minimal p-value per gene.
#'
#' Group QTL matrix by gene, sort by p_nominal, keep SNP with smalles p-value,
#' Correct that using bonferroni correction and then apply FDR correction across genes.
#'
#' @param qtl_df Data frame with QTL mapping results 
#' (required columns: gene_id, p_nominal, n_cis_snps)
#'
#' @return Only gene-SNP pairs with minimal p-values,
#'  added columns: p_bonferroni, p_fdr.
#' @export
findMinimalSnpPvalues <- function(qtl_df){
  result = dplyr::group_by(qtl_df, gene_id) %>%
    dplyr::arrange(p_nominal) %>% 
    dplyr::filter(row_number() == 1) %>%
    dplyr::mutate(p_bonferroni = p.adjust(p_nominal, "bonferroni", n_cis_snps)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(p_nominal) %>%
    dplyr::mutate(p_fdr = p.adjust(p_bonferroni, "fdr"))
  return(result)
}
