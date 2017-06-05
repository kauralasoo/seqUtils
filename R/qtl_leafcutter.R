#Construct metadata df from leafcutter event ids
leafcutterConstructMeta <- function(leafcutter_ids){
  intron_metadata = dplyr::data_frame(transcript_id = leafcutter_ids) %>%
    tidyr::separate(transcript_id, c("chr","transcript_start","transcript_end","gene_id"), sep = ":", remove = FALSE) %>%
    dplyr::mutate(transcript_start = as.integer(transcript_start), transcript_end = as.integer(transcript_end)) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(start = min(transcript_start), end = max(transcript_end)) %>%
    dplyr::mutate(strand = "+") %>%
    dplyr::ungroup()
  return(intron_metadata)
}


#Convert leafcutter metadata into a Granges list of alternative exons
leafcutterMetaToGrangesList <- function(leafutter_metadata, exon_length = 75){
  
  assertthat::assert_that(assertthat::has_name(leafutter_metadata, "transcript_id"))
  assertthat::assert_that(assertthat::has_name(leafutter_metadata, "chr"))
  assertthat::assert_that(assertthat::has_name(leafutter_metadata, "transcript_start"))
  assertthat::assert_that(assertthat::has_name(leafutter_metadata, "transcript_end"))
  assertthat::assert_that(assertthat::has_name(leafutter_metadata, "strand"))
  
  intron_df = dplyr::transmute(leafutter_metadata, transcript_id, seqnames = chr, transcript_start, 
                               transcript_end, strand) %>% 
    tidyr::gather("type", "pos", transcript_start:transcript_end)
  intron_start = dplyr::filter(intron_df, type == "transcript_start") %>% 
    dplyr::mutate(start = pos - exon_length, end = pos) %>% 
    dplyr::select(-pos)
  intron_end = dplyr::filter(intron_df, type == "transcript_end") %>% 
    dplyr::mutate(start = pos, end = pos + exon_length) %>% 
    dplyr::select(-pos)
  
  #Constrcut GRanges list
  inton_both = bind_rows(intron_start, intron_end) %>%
    dplyr::select(-type) %>%
    dplyr::group_by(transcript_id) %>%
    purrrlyr::by_slice(~dataFrameToGRanges(.))
  intron_list = setNames(inton_both$.out, inton_both$transcript_id)
  return(intron_list)
}

