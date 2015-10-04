vcfToMatrix <- function(file, genome){
  #Read vcf file into R and convert it onto matrix of SNP positions and matrix of genotypes
  genotypes_vcf = VariantAnnotation::readVcf(file, genome)
  
  # Extract SNP positions from the VCF file
  variant_granges = GenomicRanges::rowData(genotypes_vcf)
  GenomicRanges::elementMetadata(variant_granges) = c()
  snp_positions = GenomicRanges::as.data.frame(variant_granges)
  snpspos = dplyr::mutate(snp_positions, snpid = rownames(snp_positions)) %>% 
    dplyr::select(snpid, seqnames, start) %>%
    dplyr::rename(chr = seqnames, pos = start)
  
  #Extract genotype matrix
  genotypes = VariantAnnotation::geno(genotypes_vcf)$GT
  genotypes[genotypes == "1/1"] = 2
  genotypes[genotypes == "0/1"] = 1
  genotypes[genotypes == "1/0"] = 1
  genotypes[genotypes == "0/0"] = 0
  genotypes[genotypes == "."] = "NA"
  mode(genotypes) = "numeric"
  
  return(list(snpspos = snpspos, genotypes = genotypes))
}