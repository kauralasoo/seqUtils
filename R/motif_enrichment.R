cisBPImportMotif <- function(motif_id, cisbp_dir){
  motif_path = file.path(cisbp_dir, paste(motif_id, ".txt", sep = ""))
  motif = read.table(motif_path, header = TRUE)
  motif = t(motif[,-1])
  motif = round(motif*1000)
  storage.mode(motif) = "integer"
  return(motif)
}

cisBPMotifToPFMatrix <- function(motif_record, cisbp_dir){
  pfm = cisBPImportMotif(motif_record$Motif_ID, cisbp_dir)
  if(ncol(pfm) == 0){
    warning(paste("Skipping ", motif_record$Motif_ID, " because the PFM is empty.", sep = ""))
    return(NULL)
  }
  pfmatrix = TFBSTools::PFMatrix(ID = motif_record$Motif_ID, name = motif_record$TF_Name,
                      tags = list(gene_id = motif_record$gene_id),
                      profileMatrix = pfm)
  return(pfmatrix)
}

asPWMatrixList <- function(pwm_list){
  result = TFBSTools::PWMatrixList()
  for (i in seq_along(pwm_list)){
    result[[i]] = pwm_list[[i]]
  }
  names(result) = names(pwm_list)
  return(result)
}

cisBPImportRecords <- function(motifs_df, cisbp_dir, toPWM = TRUE){
  
  #Convert motifs df to list
  motif_list = plyr::dlply(motifs_df, .(Motif_ID))

  #Import PWMs
  pwm_list = lapply(motif_list, function(motif_record, cisbp_dir) {
    pfm = cisBPMotifToPFMatrix(motif_record, cisbp_dir)
    if (!is.null(pfm)){
      if(toPWM == TRUE){
        pwm = TFBSTools::toPWM(pfm)
      } else{
        pwm = pfm
      }
    } else{
      pwm = NULL
    }
    return(pwm)
  }, cisbp_dir)
  
  #Remove empty matrices
  non_empty = which(!unlist(lapply(pwm_list, is.null)))
  pwm_list = pwm_list[non_empty]
  
  #Convert to PWMatrixtList object
  pwmatrix_list = asPWMatrixList(pwm_list)
  return(pwmatrix_list)
}

subsetDNAStringSet <- function(dna_string_set, batch_number, n_batches){
  n_sequences = length(dna_string_set)
  batch_size = ceiling(n_sequences/n_batches)
  batches = splitIntoBatches(n_sequences, batch_size)
  subset = dna_string_set[batches == batch_number]
  return(subset)
}

#' Calculate relative motif enrichment in foreground set of peaks compared to a specifed background.
#'
#' @param foreground Data frame of foreground peaks, should be a subset of background peaks.
#' @param background Data frame of background peaks (deafult: NULL). If NULL then uses all peaks in peak_metadata.
#' @param fimo_hits Data frame containg all FIMO motif matches in all peaks.
#' @param peak_metadata Data frame with peak start and end coordinates.
#'
#' @return Data frame with relative enrichment of each motif inf the foreground set compared to the background set.
#' @export
fimoRelativeEnrichment <- function(foreground, background = NULL, fimo_hits, peak_metadata){
  
  #Some assertions
  assertthat::assert_that(nrow(foreground) > 0)
  assertthat::assert_that(assertthat::has_name(foreground, "gene_id"))
  assertthat::has_name(peak_metadata, "gene_id")
  assertthat::has_name(peak_metadata, "start")
  assertthat::has_name(peak_metadata, "end")
  assertthat::has_name(fimo_hits, "gene_id")
  assertthat::has_name(fimo_hits, "motif_id")
  
  #Calculate peak widths
  width_matrix = dplyr::mutate(peak_metadata, width = end - start)

  #Use all peaks as background
  if (is.null(background)){
    background = data_frame(gene_id = peak_metadata$gene_id)
  }
  
  #Count background hits
  bg_hits = fimo_hits[fimo_hits$gene_id %in% background$gene_id,]
  bg_widths = width_matrix[width_matrix$gene_id %in% background$gene_id,]
  
  #Calculate total number of matches for each motif
  bg_matches = dplyr::group_by(bg_hits, motif_id) %>%
    dplyr::summarise(n_bg_matches = length(motif_id)) %>%
    dplyr::mutate(bg_length = sum(bg_widths$width)) %>%
    dplyr::ungroup()
  
  #Count foreground hits
  fg_lengths = sum(dplyr::filter(width_matrix, gene_id %in% foreground$gene_id)$width)
  fg_matches = dplyr::filter(fimo_hits, gene_id %in% foreground$gene_id) %>% 
    dplyr::group_by(motif_id) %>% 
    dplyr::summarise(n_fg_matches = length(motif_id)) %>%
    dplyr::mutate(fg_length = fg_lengths) %>%
    dplyr::ungroup()
  
  #Motif_enrichment
  enrichments = dplyr::left_join(bg_matches, fg_matches, by = "motif_id") %>% 
    dplyr::mutate(enrichment = (n_fg_matches/fg_length)/(n_bg_matches/bg_length)) %>% 
    arrange(-enrichment) %>% 
    dplyr::group_by(motif_id) %>%
    dplyr::mutate(p_hyper = phyper(n_fg_matches - 1, n_bg_matches, bg_length - n_bg_matches, fg_length, lower.tail = FALSE)) %>%
    dplyr::mutate(p_fdr = p.adjust(p_hyper, "fdr")) %>%
    dplyr::ungroup()
  
  return(enrichments)
}
