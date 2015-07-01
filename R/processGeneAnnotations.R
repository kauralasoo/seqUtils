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

findGeneExons <- function(gene_id, gene_tx_map, tx_exon_list){
  tx_vector = gene_tx_map[gene_tx_map$ensembl_gene_id == gene_id,]$ensembl_transcript_id
  selected_txs = tx_exon_list[tx_vector]
  gene_exons = listUnion(selected_txs)
  return(gene_exons)
}

