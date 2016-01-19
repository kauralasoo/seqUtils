#' Save a list of matrices into a suitable format for RASQUAL
#'
#' Works with expression and covariates matrices.
#' 
#' @param data_list list of matrices
#' @param output_dir relative path to the output dir
#' @param file_suffix suffix added to each file after their name in the list.
#' @return None
#' @author Kaur Alasoo
#' @export 
saveRasqualMatrices <- function(data_list, output_dir, file_suffix = "expression"){
  #Save data for FastQTL to disk
  
  #Save each matrix as a separate  txt file
  for (sn in names(data_list)){
    file_path = file.path(output_dir, paste(sn,file_suffix, "txt", sep = "."))
    file_path_bin = file.path(output_dir, paste(sn,file_suffix, "bin", sep = "."))
    print(file_path)
    write.table(data_list[[sn]], file_path, quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")
    writeBin(as.double(c(t(data_list[[sn]]))), file_path_bin)
  }
}

#' Calculate sample-specic offsets for RASQUAL
#'
#'  Calculate offset values for each gene in each condition that
#'  correct for library size and GC content bias.  
#'
#' @param counts Matrix of read counts (genes in rows, samples in columns)
#' @param gene_metadata Matrix of gene metadata (required columns: gene_id, precentage_gc_content). 
#' Used for GC content bias correction. 
#' @param gc_correct If true then correct for GC content bias in addition to library size.
#'
#' @return Matrix of gene-specifc offsets.
#' @export
rasqualCalculateSampleOffsets <- function(counts, gene_metadata, gc_correct = TRUE){
  
  #Calculate library sizes
  library_size = colSums(counts)
  size_factors = library_size/mean(library_size) #Standardise
  size_matrix = matrix(rep(size_factors, nrow(counts)), nrow = nrow(counts), byrow = TRUE)
  rownames(size_matrix) = rownames(counts)
  colnames(size_matrix) = colnames(counts)
  
  #Apply GC correction
  if(gc_correct == TRUE){
    gc_factor = rasqualGcCorrection(counts, gene_metadata)
    size_matrix = size_matrix * gc_factor
  }
  return(size_matrix)
}

#' Count the number of feature SNPs and cis SNPs overlapping a set of peak calls
#'
#' @param peak_metadata Data frame with peak metadata (required columns: gene_id, chr, start, end)
#' @param snp_coords Data frame with SNP coordinates from a VCF file (required columns: chr, pos, snp_id)
#' @param cis_window Size of the cis window from both sides of the peak.
#'
#' @return Data frame with peak coordinates, cis region coordiantes as well as number of cis and feature snps.
#' @export
countSnpsOverlapingPeaks <- function(peak_metadata, snp_coords, cis_window = 500000){
  
  #Construct peak coords data frame
  peak_coords = dplyr::transmute(peak_metadata, gene_id, chromosome_name = chr, strand = 1, 
                                 exon_starts = start, exon_ends = end) %>%
    dplyr::mutate(range_start = pmax(0, exon_starts - cis_window), range_end = exon_ends + cis_window)
  
  #Construct GRanges objects
  peak_granges = GenomicRanges::GRanges(seqnames = peak_coords$chromosome_name, 
                         ranges = IRanges::IRanges(start = peak_coords$exon_starts, end = peak_coords$exon_ends))
  region_granges = GenomicRanges::GRanges(seqnames = peak_coords$chromosome_name, 
                           ranges = IRanges::IRanges(start = peak_coords$range_start, end = peak_coords$range_end))
  snp_granges = GenomicRanges::GRanges(seqnames = snp_coords$chr, ranges = IRanges::IRanges(start = snp_coords$pos, end = snp_coords$pos))
  
  #Count overlaps
  feature_snps = GenomicRanges::countOverlaps(peak_granges, snp_granges, ignore.strand = TRUE)
  cis_snps = GenomicRanges::countOverlaps(region_granges, snp_granges, ignore.strand = TRUE)
  
  new_peak_coords = dplyr::mutate(peak_coords, feature_snp_count = feature_snps, cis_snp_count = cis_snps)
  return(new_peak_coords)
}

#' Count the number of feature SNPs and cis SNPs overlapping exons of all genes
#'
#' @param gene_metadata Data frame with gene metadata (required columns: gene_id, chr, strand, exon_starts, exon_ends)
#' -exon_stars: comma-separated list of exon start coordinates
#' -exon_ends: comma-separated list of exon end coordinates.
#' @param snp_coords Data frame with SNP coordinates from a VCF file (required columns: chr, pos, snp_id)
#' @param cis_window Size of the cis window from both sides of the gene.
#'
#' @return Data frame with exon coordinates, cis region coordiantes as well as the number of cis and feature snps.
#' @export
countSnpsOverlapingExons <- function(gene_metadata, snp_coords, cis_window = 5e5){
  
  #Split exon coordinates into separate rows
  gene_df_list = plyr::dlply(gene_metadata, .(gene_id), function(x){
                          data.frame(gene_id = x$gene_id, 
                                     seqnames = x$chr,
                                     strand = x$strand,
                                     start = as.numeric(unlist(strsplit(x$exon_starts,","))),
                                     end = as.numeric(unlist(strsplit(x$exon_ends,","))) )
                        })
  exon_df = plyr::ldply(gene_df_list, .id = NULL) %>% tbl_df()
  
  #Counts the number of feature SNPs
  exon_granges = dataFrameToGRanges(exon_df)
  snp_granges = GenomicRanges::GRanges(seqnames = snp_coords$chr, 
                                       ranges = IRanges::IRanges(start = snp_coords$pos, end = snp_coords$pos))
  n_feature_snps = GenomicRanges::countOverlaps(exon_granges, snp_granges, ignore.strand=TRUE)
  feature_snp_df = dplyr::mutate(exon_df, feature_snp_count = n_feature_snps) %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::summarise(seqnames = seqnames[1], strand = strand[1], start = min(start), end = max(end), 
                     feature_snp_count = sum(feature_snp_count))
  
  #Count the number of cis SNPs
  cis_df = dplyr::mutate(feature_snp_df, start = pmax(0, start - cis_window), end = end + cis_window)
  cis_granges = dataFrameToGRanges(cis_df)
  n_cis_snps = GenomicRanges::countOverlaps(cis_granges, snp_granges, ignore.strand=TRUE)
  result = dplyr::mutate(cis_df, cis_snp_count = n_cis_snps, gene_id = as.character(gene_id)) %>%
    dplyr::rename(range_start = start, range_end = end)
  
  #Add exon start and end coords and reorder columns
  start_end_df = dplyr::select(gene_metadata, gene_id, exon_starts, exon_ends)
  result = dplyr::left_join(result, start_end_df, by = "gene_id") %>%
    dplyr::transmute(gene_id, chromosome_name = seqnames, strand, exon_starts, exon_ends, 
                     range_start, range_end, feature_snp_count, cis_snp_count)
  return(result)
  }

#' Import rasqual output table into R
#'
#' Skipped gene-SNP pairs are automatically removed and 
#' chisq statistic is converted into p-value (p_nominal).
#'
#' @param path Bath to rasqual output file.
#'
#' @return data_frame
#' @export
importRasqualTable <- function(path){
  rasqual_results = readr::read_delim(path, delim = "\t", col_types = "cccidddiii",col_names = FALSE)
  colnames(rasqual_results) = c("gene_id", "snp_id", "chr", "pos", "allele_freq", "chisq", "pi", "n_feature_snps", "n_cis_snps", "converged")
  
  rasqual_pvalues = dplyr::filter(rasqual_results, snp_id != "SKIPPED") %>%
    dplyr::mutate(p_nominal = pchisq(chisq, df = 1, lower = FALSE))
  
  return(rasqual_pvalues)
}

#' Helper function for rasqualGcCorrection
#' 
#' @author Natsuhiko Kumasaka
Quantile <- function(x,k=20){
  x=rank(x,ties="random")
  z=rep(0,length(x))
  for(i in 1:k){
    z = z+as.numeric(x<=quantile(x,i/k,na.rm=T))
  }
  k-z
}

#' Estimate the effect of GC-bias for each feature in each sample.
#' 
#' This function does not correct for differences in library size between samples.
#'
#' @param Y Matrix of read counts
#' @param gene_metadata Data frame with gene metadata.
#' (required columnss: gene_id, percentage_gc_content)
#' @param PLOT 
#'
#' @return Matrix of GC-bias offsets for each gene in each condition.
#' @export
#' @author Natsuhiko Kumasaka
rasqualGcCorrection <- function(Y,gene_metadata,PLOT=F){
  
  #Extract GC vector from the metadata matrix
  gene_metadata = as.data.frame(gene_metadata)
  rownames(gene_metadata) = gene_metadata$gene_id
  gene_metadata = gene_metadata[rownames(Y),]
  gcvec = gene_metadata$percentage_gc_content
  
  #Perfrorm GC correction
  bin=Quantile(gcvec,200);
  x=sort(unlist(lapply(split(gcvec,bin),mean)))
  S=apply(Y,2,function(y){unlist(lapply(split(y,bin),sum))[as.character(0:199)]});
  Fs=log(t(t(S)/apply(S,2,sum))/apply(S,1,sum)*sum(as.numeric(S)));
  Gs=apply(Fs,2,function(y){smooth.spline(x,y,spar=1)$y}); 
  if(PLOT){
    par(mfcol=c(5,5),mar=c(2,2,2,2)); 
    for(i in 1:ncol(Y)){
      plot(Fs[,i])
      lines(Gs[,i],col=2)
    }
    matplot(x,Gs,type="l",col=2,lty=1)
  }
  result_matrix = exp(Gs[bin+1,])
  #Add gene ids to rows
  rownames(result_matrix) = rownames(Y)
  return(result_matrix)
}

#' Split gene ids into batches for runRasqual.py script
#'
#' @param gene_metadata Data frame with gene metadata (gene_id column required)
#' @param batch_size Number of genes in a batch
#'
#' @return Data frame of gene batches (columns: batch_id, gene_ids)
#' @export
rasqualConstructGeneBatches <- function(gene_metadata, batch_size){
  batch_df = dplyr::select(gene_metadata, gene_id) %>%
    dplyr::mutate(batch_number = splitIntoBatches(length(gene_id), batch_size)) %>%
    dplyr::mutate(batch_id = paste("batch", batch_number, sep = "_")) %>% 
    dplyr::group_by(batch_id) %>%
    dplyr::summarize(gene_ids = paste(gene_id, collapse = ","))
  return(batch_df)
}

#' Find SNP with minimal p-value per gene.
#'
#' Group QTL matrix by gene, sort by p_nominal, keep SNP with smalles p-value,
#' Correct that using bonferroni correction and then apply FDR correction across genes.
#'
#' @param qtl_df Data frame with QTL mapping results 
#' (required columns: gene_id, p_nominal, n_cis_snps)
#'
#' @return Only gene-SNP pairs with minimal p-values,
#'  added columns: p_bonferroni, p_fdr.
#' @export
findMinimalSnpPvalues <- function(qtl_df){
  result = dplyr::group_by(qtl_df, gene_id) %>%
    dplyr::arrange(p_nominal) %>% 
    dplyr::filter(row_number() == 1) %>%
    dplyr::mutate(p_bonferroni = p.adjust(p_nominal, "bonferroni", n_cis_snps)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(p_nominal) %>%
    dplyr::mutate(p_fdr = p.adjust(p_bonferroni, "fdr"))
  return(result)
}

#' Find SNP with minimal p-value per gene from SQLite database.
#'
#' Group QTL database by gene, sort by p_nominal, keep SNP with smallest p-value,
#' correct that using bonferroni correction and then apply FDR correction across genes.
#'
#' @param sqlite_path Path to the SQLlite database file with RASQUAL output. 
#' Table name rasqual (required columns: gene_id, p_nominal, n_cis_snps)
#'
#' @return Only gene-SNP pairs with minimal p-values,
#'  added columns: p_bonferroni, p_fdr.
#' @export
findMinimalSnpPvaluesSQLite <- function(sqlite_path){
  
  #Connect to database filter by maximal chisq value
  my_db <- dplyr::src_sqlite(sqlite_path)
  table = tbl(my_db, sql("SELECT *,MAX(chisq) AS chisq2 FROM rasqual GROUP BY gene_id")) %>% 
    collect() %>%
    dplyr::filter(snp_id != "SKIPPED") %>%
    dplyr::select(-chisq2) %>% #Remove the second copy of chisq column
    dplyr::mutate(p_nominal = pchisq(chisq, df = 1, lower = FALSE)) %>% #Calculate p-value
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(p_bonferroni = p.adjust(p_nominal, "bonferroni", n_cis_snps)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(p_nominal) %>%
    dplyr::mutate(p_fdr = p.adjust(p_bonferroni, "fdr"))
  return(table)
}

#' Extract rasqual results for particular gene_id or snp_id from the SQLite database
#'
#' Either selected_gene_id or selected_snp_id must be specified. 
#' 
#' @param db_table Connection to SQLite database table.
#' (required columns: gene_id, snp_id, p_nominal, allele_freq).
#' @param selected_gene_id Id of the selected gene. 
#' @param selected_snp_id Id of the selected SNP.
#'
#' @return data.frame with rasqual results for a particular gene or SNP.
#' @export
fetchSQLite <- function(db_table, selected_gene_id = NULL, selected_snp_id = NULL){
  
  if(is.null(selected_snp_id) & is.null(selected_gene_id)){
    stop("Need to specify at least selected_gene_id or selected_snp_id.")
  } else if (!is.null(selected_snp_id) & !is.null(selected_gene_id)){
    result = dplyr::filter(db_table, gene_id == selected_gene_id, snp_id == selected_snp_id) %>% dplyr::collect()
  } else if (!is.null(selected_gene_id)){
    result = dplyr::filter(db_table, gene_id == selected_gene_id) %>% dplyr::collect()
  } else {
    result = dplyr::filter(db_table, snp_id == selected_snp_id) %>% dplyr::collect()
  }
  
  #Add p-value, beta and MAF
  result = dplyr::mutate(result, p_nominal = pchisq(chisq, df = 1, lower = FALSE)) %>%
    dplyr::mutate(MAF = pmin(allele_freq, 1-allele_freq)) %>%
    dplyr::mutate(beta = -log(effect_size/(1-effect_size),2)) #Calculate beta from rasqual pi
  return(result)
}

#' Fetch multiple genes from SQLite database
#' NOTE: This implementation is extremely inefficient
fetchMultipleGenes <- function(gene_ids, db_table){
  
  gene_id_list = idVectorToList(gene_ids)
  result = lapply(gene_id_list, function(gene_id, db){ fetchSQLite(db, selected_gene_id = gene_id) }, db_table) %>%
    ldply(.id = NULL)
  return(result)
}

#' Write multiple files for RASQUAL onto disk.
#'
#' @param condition_list Named list of expression lists, each expression list needs to
#' contain at least the following elements: 'counts' matrix, 'sample_metadata' df,
#' @param rasqual_input_folder Path to the RASQUAL input folder.
#' @param max_batch_size Maximal number of feaures to be included in a single batch.
#'
#' @export
exportDataForRasqual <- function(condition_list, rasqual_input_folder, max_batch_size = 50){
  if(!is.list(condition_list)){
    stop("Input condiiton_list must be a list.")
  }
  if(length(condition_list) < 1){
    stop("The conditon_list must have at least one element.")
  }
  
  #Extract and save read count matrices
  counts_list = lapply(condition_list, function(x){x$counts})
  saveRasqualMatrices(counts_list, rasqual_input_folder, file_suffix = "expression")
  
  #Extract sample-genotype map for each condition
  sg_map = lapply(condition_list, function(x){ dplyr::select(x$sample_metadata, sample_id, genotype_id) })
  saveFastqtlMatrices(sg_map, rasqual_input_folder, file_suffix = "sg_map", col_names = FALSE)
  
  #Export library size
  library_size_list = lapply(counts_list, rasqualCalculateSampleOffsets, condition_list[[1]]$gene_metadata, gc_correct = FALSE)
  saveRasqualMatrices(library_size_list, rasqual_input_folder, file_suffix = "library_size")
  
  #Export GC-corrected library sizes
  gc_library_size_list = lapply(counts_list, rasqualCalculateSampleOffsets, condition_list[[1]]$gene_metadata)
  saveRasqualMatrices(gc_library_size_list, rasqual_input_folder, file_suffix = "gc_library_size")
  
  #Export offsets from the cqn pacakge
  
  #Save feature names to disk
  feature_names = rownames(counts_list[[1]])
  write.table(feature_names, file.path(rasqual_input_folder, "feature_names.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  #Construct a list of batches for rasqual
  feature_batches = rasqualConstructGeneBatches(condition_list[[1]]$gene_metadata, max_batch_size)
  write.table(feature_batches, file.path(rasqual_input_folder, "feature_batches.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}


