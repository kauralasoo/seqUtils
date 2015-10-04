filterOverlaps <- function(peak_list, minOverlapCount = 3){
  #For each set of peaks, find the ones that are present in at least minOverlapCount samples
  result = list()
  for (i in 1:length(peak_list)){
    query_peaks = peak_list[[i]]
    query_hits = c()
    for (j in 1:length(peak_list)){
      target_peaks = peak_list[[j]]
      overlaps = findOverlaps(query_peaks, target_peaks)
      unique_hits = unique(queryHits(overlaps))
      query_hits = c(query_hits, unique_hits)
    }
    query_hits = table(query_hits)
    query_hits = query_hits[query_hits >= minOverlapCount]
    query_filtered = query_peaks[as.numeric(rownames(query_hits))]
    result[i] = query_filtered
  }
  return(result)
}