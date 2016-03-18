
cisBPImportMotif <- function(motif_id, cisbp_dir){
  motif_path = file.path(cisbp_dir, paste(motif_id, ".txt", sep = ""))
  motif = read.table(motif_path, header = TRUE)
  motif = t(motif[,-1])
  motif = round(motif*1000)
  return(motif)
}
