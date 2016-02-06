#Scripts to convert data to/from format suitable to for eigenMT.py script

#' Export genotype data in format suitable for eigenMT
#'
#'  Creates two text files in the eigenMT_dir directory: [prefix].snp_positions.txt and 
#'  [prefix].genotypes.txt.
#'  
#' @param vcf_file Genotype data imported with gdsToMatrix
#' @param eigenMT_dir Directory of the eigenMT input files
#' @param prefix Prefix of the output files,
#'
#' @export
eigenMTExportGenotypes <- function(vcf_file, eigenMT_dir, prefix){
  #Save SNP positions
  snp_pos_df = vcf_file$snpspos %>% dplyr::rename(snp = snpid, chr_snp = chr)
  write.table(snp_pos_df, file.path(eigenMT_dir, paste(prefix, ".snp_positions.txt", sep ="")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  #Save genotypes
  genotypes = dplyr::mutate(as.data.frame(vcf_file$genotypes), ID = rownames(vcf_file$genotypes)) %>% 
    dplyr::select(ID, everything())
  write.table(genotypes, file.path(eigenMT_dir, paste(prefix,".genotypes.txt",sep = "")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
}

eigenMTExportGenotypesByChr <- function(chromosomes, gds_dir, eigenMT_dir, gds_prefix = "chr_"){
  for (chr in chromosomes){
    gds_file = file.path(gds_dir, paste(gds_prefix, chr, ".gds", sep = ""))
    print(gds_file)
    vcf_file = gdsToMatrix(gds_file)
    eigenMT_prefix = paste(gds_prefix, chr, sep = "")
    eigenMTExportGenotypes(vcf_file, eigenMT_dir, eigenMT_prefix)
  }
}

eigenMTExportGeneMetadata <- function(gene_metadata, eigenMT_dir){
  #Save gene positions
  gene_data = dplyr::transmute(gene_metadata, gene_id, chrom_probe = chr, s1 = start, s2 = end) %>%
    dplyr::arrange(chrom_probe, s1)
  write.table(gene_data, file.path(eigenMT_dir, "gene_positions.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

eigenMTImportResults <- function(file_path){
  result = readr::read_delim(file_path, delim = "\t", col_types = "ccdddddd", 
                             col_names = c("snp_id", "gene_id","chisq","p_nominal", "FDR", "pi", "p_eigen", "n_tests")) %>%
    dplyr::select(-FDR, -pi) %>%
    dplyr::mutate(p_fdr = p.adjust(p_eigen, method = "fdr")) %>% 
    dplyr::arrange(p_eigen)
}

