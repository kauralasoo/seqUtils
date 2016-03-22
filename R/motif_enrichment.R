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
