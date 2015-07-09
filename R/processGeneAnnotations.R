listUnion <- function(granges_list){
  #Calculated the union of a GRangesList object
  union_obj = granges_list[[1]]
  if(length(granges_list) > 1){
    for(i in 2:length(granges_list)){
      union_obj = GenomicRanges::union(union_obj, granges_list[[i]]) 
    } 
  }
  return(union_obj)
}


reduceDataFrameGRanges <- function(df){
  #Convert data frame into a a Granges object and reduce overlapping exons, then convert back to df
  
  gr = GRanges(seqnames = df$seqnames, ranges = IRanges(start = df$start, end = df$end), strand = df$strand)
  #Reduce exons
  gr = IRanges::reduce(gr)
  #Convert back to data.frame
  gr = GenomicRanges::as.data.frame(gr)
  return(gr)
}

GRangesListToDataFrame <- function(gr_list, joint_id_df){
  #Convert GRangesList into a data frame of exon coordinates
  # gr_list - GRangesList with tr_ids as names
  # joint_id_df - data.frame(with rownames) containing three columns: gene_id, transcript_id and joint_id
  
  names(gr_list) = joint_id_df[names(gr_list),]$joint_id
  exons_df = plyr::ldply(gr_list, GenomicRanges::as.data.frame, .id = "joint_id") %>%
    dplyr::mutate(joint_id = as.vector(joint_id)) %>%
    tidyr::separate(joint_id, into = c("gene_id", "transcript_id"), sep = ":")
  return(exons_df)
}

extractExonsStartEnd <- function(gr_list, gene_tx_map){
  #Extract Exon start and end coorinates per gene from a GrangesList object of exons by transcript id.
  
  #Count the number of transcripts per gene
  transcript_counts = dplyr::group_by(gene_tx_map, gene_id) %>% 
    dplyr::summarize(transcript_count = length(transcript_id))
  
  #Extract genes with only one transcipt
  single_transcript_genes = dplyr::filter(transcript_counts, transcript_count == 1)
  print(dim(single_transcript_genes))
  single_tx_annot = dplyr::semi_join(gene_tx_map, single_transcript_genes, by = "gene_id")
  
  #Extract genes with more than one transcript
  multi_transcript_genes = dplyr::filter(transcript_counts, transcript_count > 1)
  print(dim(multi_transcript_genes))
  multi_tx_annot = dplyr::semi_join(gene_tx_map, multi_transcript_genes, by = "gene_id")
  
  #Construct joint id for multi tx genes
  joint_id_df = dplyr::mutate(gene_tx_map, joint_id = paste(gene_id, transcript_id, sep = ":")) %>% as.data.frame()
  rownames(joint_id_df) = joint_id_df$transcript_id
  
  #Extract exon coordinates for single tx genes
  single_tx_exons = exons[single_tx_annot$transcript_id]
  single_tx_exons_df = GRangesListToDataFrame(single_tx_exons, joint_id_df)
  single_tx_exons_coords = dplyr::group_by(single_tx_exons_df, gene_id) %>% 
    arrange(start) %>% 
    summarize(starts = paste(start, collapse = ","), ends = paste(end, collapse = ","))
  
  #Extract joint exon coordinates for multi tx genes
  multi_tx_exons = exons[multi_tx_annot$transcript_id]
  multi_tx_exons_df = GRangesListToDataFrame(multi_tx_exons, joint_id_df)
  multi_tx_exons_reduced_df = dlply(multi_tx_exons_df, .(gene_id), seqUtils::reduceDataFrameGRanges) %>%
    ldply(.id = "gene_id")
  multi_tx_exons_coords = dplyr::group_by(multi_tx_exons_reduced_df, gene_id) %>% 
    arrange(start) %>% 
    summarize(starts = paste(start, collapse = ","), ends = paste(end, collapse = ","))
  
  #bind results together
  result = rbind(single_tx_exons_coords, multi_tx_exons_coords)
  return(result)
  
}
