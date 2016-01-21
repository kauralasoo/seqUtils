#Scripts to convert data to/from format suitable to for eigenMT.py script

eigenMTExportGenotypes <- function(vcf_file, eigenMT_dir){
  #Save SNP positions
  snp_pos_df = vcf_file$snpspos %>% dplyr::rename(snp = snpid, chr_snp = chr)
  write.table(snp_pos_df, file.path(eigenMT_dir, "snp_positions.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  #Save genotypes
  genotypes = dplyr::mutate(as.data.frame(vcf_file$genotypes), ID = rownames(vcf_file$genotypes)) %>% dplyr::select(ID, everything())
  write.table(genotypes, file.path(eigenMT_dir, "genotypes.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

eigenMTExportGeneMetadata <- function(gene_metadata, eigenMT_dir){
  #Save gene positions
  gene_data = dplyr::transmute(gene_metadata, gene_id, chrom_probe = chromosome_name, s1 = start_position, s2 = end_position) %>%
    dplyr::arrange(chrom_probe, s1)
  write.table(gene_data, file.path(eigenMT_dir, "gene_positions.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

eigenMTImportResults <- function(file_path){
  result = readr::read_delim(file_path, delim = "\t", col_types = "ccdddddd") %>%
    dplyr::transmute(gene_id = gene, snp_id = snps, chisq = statistic, p_nominal = pvalue, 
                     p_eigen = BF, n_tests = TESTS, p_fdr = p.adjust(p_eigen, method = "fdr")) %>% 
    dplyr::arrange(p_eigen)
}

