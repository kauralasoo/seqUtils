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
  assertthat::assert_that(assertthat::has_name(peak_metadata, "gene_id"))
  assertthat::assert_that(assertthat::has_name(peak_metadata, "start"))
  assertthat::assert_that(assertthat::has_name(peak_metadata, "end"))
  assertthat::assert_that(assertthat::has_name(fimo_hits, "gene_id"))
  assertthat::assert_that(assertthat::has_name(fimo_hits, "motif_id"))
  
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
    dplyr::filter(!is.na(n_fg_matches)) %>% #Not sure why this is necessary
    purrr::by_row(., ~subsampleFisherTest(.$n_fg_matches, .$fg_length, .$n_bg_matches, .$bg_length), .collate = "rows")
  return(enrichments)
}


#' Perform Fisher's Exact test when sample is a subsample of a larger population
#' 
#' Return a tidy data frame with results (p-value, OR, confidence intervals and their log2 coversions)
#'
#' @param sample_hits number of hits in the sample
#' @param sample_size sample size
#' @param total_hits number of hits in the total populaton
#' @param total_size population size
#'
#' @return data frame with p-value, OR and confidence intervals
#' @export
subsampleFisherTest <- function(sample_hits, sample_size, total_hits, total_size){
  assertthat::assert_that(assertthat::is.number(sample_hits))
  assertthat::assert_that(assertthat::is.number(sample_size))
  assertthat::assert_that(assertthat::is.number(total_hits))
  assertthat::assert_that(assertthat::is.number(total_size))
  
  #Set up data for the test
  bg_hits = total_hits - sample_hits
  sample_fails = sample_size - sample_hits
  bg_fails = total_size - total_hits - sample_size + sample_hits
  
  #Perform Fisher's Exact Test
  fet_matrix = matrix(c(sample_hits, bg_hits, sample_fails, bg_fails), 2,2)
  test = fisher.test(fet_matrix)
  
  #Extract results
  results_df = data_frame(fisher_pvalue = test$p.value, OR = test$estimate, 
                          ci_lower = test$conf.int[1], ci_higher = test$conf.int[2]) %>%
    dplyr::mutate(OR_log2 = log(OR,2), ci_lower_log2 = log(ci_lower,2), 
                  ci_higher_log2 = log(ci_higher,2))
  return(results_df)
  
  #hyper = phyper(sample_hits-1, motif_row$baseline_disruption, 
  #       motif_row$baseline_peak_count-motif_row$baseline_disruption, 
  #       motif_row$cluster_size, lower.tail = FALSE)
  #print(hyper)
}


#' Censor log2 odds ratios and their confidence intervals for easier plotting
#'
#' @param or_df data frame of ORs
#' @param min_OR_log2 
#' @param min_CI_lower_log2 
#'
#' @return Modified or_df with censored values
#' @export
censorLog2OR <- function(or_df, min_OR_log2 = -3, min_CI_lower_log2 = -4){
  result = dplyr::mutate(or_df, OR_log2 = ifelse(OR_log2 < min_OR_log2, min_OR_log2, OR_log2)) %>%
    dplyr::mutate(ci_lower_log2 = ifelse(ci_lower_log2 < min_CI_lower_log2, min_CI_lower_log2, ci_lower_log2)) %>%
    dplyr::mutate(ci_higher_log2 = ifelse(ci_higher_log2 < min_CI_lower_log2, min_CI_lower_log2, ci_higher_log2))
  return(result)
}

fimoFisherTest <- function(bg_motif_hits, fg_motif_hits, bg_seq_length, fg_seq_length){
  
  #Count the number of mathces in the background sequences
  bg_matches = dplyr::group_by(fimo_promoter_hits, motif_id) %>% 
    dplyr::summarise(n_bg_matches = length(seq_name)) %>% 
    dplyr::mutate(bg_length = bg_seq_length) %>%
    dplyr::ungroup()
  
  #Count the number of matches in the foreground sequences
  fg_matches = dplyr::group_by(fimo_hits, motif_id) %>%
    dplyr::summarise(n_fg_matches = length(motif_id)) %>%
    dplyr::mutate(fg_length = fg_seq_length) %>%
    dplyr::ungroup()
  
  #Put the two lists together
  joint_matches = dplyr::left_join(bg_matches, fg_matches, by = "motif_id") %>%
    dplyr::mutate(fold_enrichment = (n_fg_matches/fg_length)/(n_bg_matches/bg_length)) %>% 
    dplyr::arrange(-fold_enrichment)
  
  #Peform Fisher's Exact test on each motif
  fisher_tests = purrr::by_row(joint_matches, ~fisher.test(matrix(c(.$n_fg_matches, .$fg_length-.$n_fg_matches,
                              .$n_bg_matches, .$bg_length-.$n_bg_matches), 2, 2) ), .to = "fisher_test")
  
  #Extract p-values and ORs
  enrichment = purrr::by_row(fisher_tests, function(df_row){
    test = df_row$fisher_test[[1]]
    result = data_frame(fisher_pvalue = test$p.value,
                        OR = test$estimate, 
                        ci_lower = test$conf.int[1], 
                        ci_higher = test$conf.int[2])
  }, .collate = "rows") %>%
    dplyr::mutate(OR_log2 = log(OR,2), ci_lower_log2 = log(ci_lower,2), 
                  ci_higher_log2 = log(ci_higher,2)) %>%
    dplyr::select(-fisher_test, -.row)
  
  return(enrichment)
}

modifyDNAString <- function(dna_string, pos, ref_value, alt_value){
  
  #Convert ref and alt to DNAString
  ref = as(ref_value, "DNAString")
  alt = as(alt_value, "DNAString")
  
  #Extraxt upstream and downstream sequences from the dna_string
  if(pos <= 1){ #If variant is at first position
    upstream = as("","DNAString")
  } else{
    upstream = dna_string[1:pos-1]
  }
  
  #If variant is at last position
  if(pos == length(dna_string)){
    downstream = as("","DNAString")
  } else{
    downstream_start_coord = pos + length(ref)
    #Handle cases where indel is wider than the peak itself.
    if(downstream_start_coord <= length(dna_string)){
      downstream = dna_string[(pos + length(ref)):length(dna_string)]
    } else{
      warning("Reference indel length exceeds peak boundaries.")
      downstream = as("","DNAString")
    }
  }
  
  #Make alternate sequence
  alt_sequence = c(upstream, alt, downstream)
  return(alt_sequence)
}

#' Estimate how much PWM binding score is altered by the SNP
#'
#' @param pwm Position weight matrix in TFBSTools format
#' @param peak_id ID of the ATAC peak
#' @param snp_id ID of the variant
#' @param peak_metadata Peak metadata (required columns: gene_id, start)
#' @param peak_sequences DNAStringSet object containing the peak sequences.
#' @param snp_metadata SNP metadata (required columns: snp_id, chr, pos, ref, alt)
#' @param window_size Window size around the SNP
#'
#' @return Data frame of binding socres for the reference and alternate versions of the sequence.
#' @export
quantifyMotifDisruption <- function(pwm, peak_id, snp_id, peak_metadata, peak_sequences, snp_metadata, window_size = 25){
  
  #Some assertions
  assertthat::assert_that(assertthat::has_name(peak_metadata, "gene_id"))
  assertthat::assert_that(assertthat::has_name(peak_metadata, "start"))
  assertthat::assert_that(assertthat::has_name(snp_metadata, "snp_id"))
  assertthat::assert_that(assertthat::has_name(snp_metadata, "chr"))
  assertthat::assert_that(assertthat::has_name(snp_metadata, "pos"))
  assertthat::assert_that(assertthat::has_name(snp_metadata, "ref"))
  assertthat::assert_that(assertthat::has_name(snp_metadata, "alt"))
  assertthat::assert_that(class(sequences) == "DNAStringSet")
  
  #Extract data
  snp = snp_id
  ref_seq = sequences[[peak_id]]
  peak_info = dplyr::filter(peak_metadata, gene_id == peak_id)
  snp_info = dplyr::filter(snp_metadata, snp_id == snp)
  
  #Construct alternate sequence
  variant_pos = snp_info$pos - peak_info$start + 1
  alt_seq = modifyDNAString(ref_seq, variant_pos, snp_info$ref, snp_info$alt)
  length_diff = length(alt_seq) - length(ref_seq)

  #Keep only sequence around the SNP
  region_start = max(variant_pos - window_size, 0)
  min_length = min(length(ref_seq), length(alt_seq))
  region_end = min(variant_pos + window_size, min_length)
  
  #Construct regions around the variant
  ref_region = ref_seq[region_start:region_end]
  #Make sure that end of of the alternative sequence is not before the variant position (problem with long indels)
  alt_end = max(variant_pos, (region_end + length_diff))
  alt_region = alt_seq[region_start:alt_end]
  
  #Search for motif matches
  ref_matches = TFBSTools::searchSeq(pwm, ref_region, min.score="0%", strand="*") %>% TFBSTools::as.data.frame()
  alt_matches = TFBSTools::searchSeq(pwm, alt_region, min.score="0%", strand="*") %>% TFBSTools::as.data.frame()

  #Rename columns
  ref_scores = dplyr::transmute(ref_matches, start, strand, motif_id = ID, ref_abs_score = absScore, ref_rel_score = relScore, ref_match = siteSeqs)
  alt_scores = dplyr::transmute(alt_matches, start, strand, alt_abs_score = absScore, alt_rel_score = relScore, alt_match = siteSeqs)
  
  #Join both of the matches together
  if (length_diff == 0){ #For SNPs perform position by position comparison
   combined_scores = dplyr::left_join(ref_scores, alt_scores , by = c("start","strand")) %>% 
      dplyr::mutate(rel_diff = alt_rel_score - ref_rel_score, max_rel_score = pmax(ref_rel_score, alt_rel_score))
  } else {
    #For indels positions get shifted, better to look if the max score of the motif changes because of the indel
    ref_max = ref_scores %>% dplyr::arrange(-ref_rel_score) %>% dplyr::filter(row_number() == 1)
    alt_max = alt_scores %>% dplyr::arrange(-alt_rel_score) %>% dplyr::filter(row_number() == 1) %>% 
      dplyr::select(-start, -strand)
    combined_scores = cbind(ref_max, alt_max) %>%
      dplyr::mutate(rel_diff = alt_rel_score - ref_rel_score, max_rel_score = pmax(ref_rel_score, alt_rel_score))
  }

  return(combined_scores)
}

#' Applies quantifyMotifDisruption to a list of motifs
#' 
#' Converts the output into a data frame and performs some basic filtering
#'
#' @param max_score_thresh Maximal relative score across variants must be greater than this (default = 0.8)
#' @param rel_diff_thresh Difference in relative binding score must be greater than this (default = 0)
quantifyMultipleMotifs <- function(peak_id, snp_id, pwm_list, peak_metadata, peak_sequences, snp_metadata, window_size = 25, 
                                   max_score_thresh = 0.7, rel_diff_thresh = 0){
  motif_disruptions = purrr::map(as.list(pwm_list), ~quantifyMotifDisruption(., peak_id, snp_id, peak_metadata, 
                                                                             peak_sequences,snp_metadata, window_size) %>%
                                   dplyr::filter(max_rel_score >= max_score_thresh, abs(rel_diff) > rel_diff_thresh))
  result_df = purrr::map_df(motif_disruptions, ~dplyr::mutate(.,strand = as.character(strand), motif_id = as.character(motif_id)))
  return(result_df)
}

#' Applies quantifyMultipleMotifs to a data frame of peak and snp ids
#' 
#' @param peak_snp_list data frame of peak and snp ids (required columns: gene_id, snp_id) 
quantifyMultipleVariants <- function(peak_snp_list, motif_list, peak_metadata, peak_sequences, snp_metadata, window_size = 25){
  
  motif_disruptions = purrrlyr::by_row(peak_snp_list, function(x, ...){
    print(x$gene_id)
    quantifyMultipleMotifs(x$gene_id, x$snp_id, ...)
  },motif_list, peak_metadata, peak_sequences, snp_metadata, window_size, .collate = "rows")
  
  return(motif_disruptions)
}

saveMotifListLogos <- function(motif_list, outdir, width, height){
  motif_names = names(motif_list)
  for(motif_name in motif_names){
    motif = motif_list[[motif_name]]
    outfile = file.path(outdir, paste0(motif_name, ".pdf"))
    pdf(outfile, width, height)
    seqLogo::seqLogo(motif)
    dev.off()
  }
}

PWMSimilarityMatrix <- function(pwm_list, method){
  sim_matrix = matrix(0, length(pwm_list), length(pwm_list))
  colnames(sim_matrix) = names(pwm_list)
  rownames(sim_matrix) = names(pwm_list)
  for (i in seq_along(pwm_list)){
    for (j in seq_along(pwm_list)) {
      sim_matrix[i,j] = TFBSTools::PWMSimilarity(pwm_list[i], pwm_list[j], method)
    }
  }
  return(sim_matrix)
}