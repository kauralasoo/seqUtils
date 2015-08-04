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

constructIntronExonDf <- function(gene_id, gene_data, intron_gap = 10, type = "intron"){
  #Return data frames with intron and exon coordinates and names
  # exon_starts - comma-separated string of exon start coords.
  # exon_ends - comma-separated string of exon end coords.
  
  #extract gene data
  gene_record = gene_data[gene_data$gene_id == gene_id,]
  exon_starts = gene_record$exon_starts
  exon_ends = gene_record$exon_ends
  seqname = gene_record$chr
  strand = gene_record$strand
  
  #Extract cooridnates from string
  exon_start_coords = as.numeric(unlist(strsplit(exon_starts, ",")))
  exon_end_coords = as.numeric(unlist(strsplit(exon_ends, ",")))

  #Construct intron data frame
  if (type == "intron"){
    if(length(exon_start_coords) > 1){ #Only genes with more than one exon can have introns
      intron_df = data_frame(start = exon_end_coords[1:length(exon_end_coords)-1] + intron_gap +1,
                             end = exon_start_coords[2:length(exon_start_coords)] - intron_gap -1,
                             strand = strand,
                             seqnames = seqname,
                             Parent = gene_id,
                             type = "intron",
                             gene_id = gene_id)
      intron_df$ID = paste(gene_id, "intron", seq(1,nrow(intron_df)), sep = "_")
      return(intron_df)
    } else{
      return(NULL)
    } 
  } else if (type == "exon"){   #Construct exon data.frame
    exon_df = data_frame(start = exon_start_coords, 
                         end = exon_end_coords,
                         strand = strand,
                         seqnames = seqname,
                         Parent = gene_id,
                         type = "exon",
                         gene_id = gene_id)
    exon_df$ID = paste(gene_id, "exon", seq(1,nrow(exon_df)), sep = "_")
    return(exon_df)
  }
}

dataFrameToGRanges <- function(df){
  #Convert a data.frame into a GRanges object
  
  gr = GRanges(seqnames = df$seqnames, 
               ranges = IRanges(start = df$start, end = df$end),
               strand = df$strand)
  
  #Add metadata
  meta = dplyr::select(df, -start, -end, -strand, -seqnames)
  elementMetadata(gr) = meta
  
  return(gr)
}