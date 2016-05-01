#' Filter peak calls based on how many other samples they are found in.
#'
#' For each set of peak calls, find the ones that are present in at least minOverlapCount samples.
#' 
#' @param peak_list list of GRanges objects corresponding to peak calls.
#' @param minOverlapCount Minimal number of samples the peak has to be present to be retained.
#' @return List of GRanges objects corresponding to filtered peak calls from each sample, only peaks 
#' that were present in minOverlapCount samples are retained.
#' @author Kaur Alasoo
#' @export 
filterOverlaps <- function(peak_list, minOverlapCount = 3){
  #For each set of peaks, find the ones that are present in at least minOverlapCount samples
  result = list()
  for (i in 1:length(peak_list)){
    query_peaks = peak_list[[i]]
    query_hits = c()
    for (j in 1:length(peak_list)){
      target_peaks = peak_list[[j]]
      overlaps = GenomicRanges::findOverlaps(query_peaks, target_peaks)
      unique_hits = unique(GenomicRanges::queryHits(overlaps))
      query_hits = c(query_hits, unique_hits)
    }
    query_hits = table(query_hits)
    query_hits = query_hits[query_hits >= minOverlapCount]
    query_filtered = query_peaks[as.numeric(rownames(query_hits))]
    result[i] = query_filtered
  }
  return(result)
}

#' Construct the union of list of peak calls and add gene_id and type metadata columns for read counting.
makeUnionPeaks <- function(peak_list, seqlevels, id_prefix){
  peaks = lapply(peak_list, function(x) keepSeqlevels(x, seqlevels))
  union_peaks = listUnion(peaks)
  metadata = dplyr::data_frame(type = "exon", gene_id = paste(id_prefix, c(1:length(union_peaks)), sep = ""))
  elementMetadata(union_peaks) = metadata
  return(union_peaks)
}