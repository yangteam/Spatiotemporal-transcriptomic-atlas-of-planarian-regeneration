BLAST_to_ID <- function(BLAST_cluster, type = "direct"){
  if (type == "direct") {
    return(gene_info_blast$Smed.ID[match(BLAST_cluster, gene_info_blast$BLAST)])
  }
  else if (type == "best e-value") {
    ID_cluster <- rep(0, length(BLAST_cluster))
    for (i in 1:length(BLAST_cluster)) {
      correspond_id <- gene_info_blast$Smed.ID[which(gene_info_blast$BLAST == BLAST_cluster[i])]
      correspond_evalue <- gene_info_blast$E.Value[match(correspond_id, gene_info_blast$Smed.ID)]
      correspond_id <- correspond_id[!is.na(correspond_evalue)]
      correspond_evalue <- correspond_evalue[!is.na(correspond_evalue)]
      correspond_id <- correspond_id[which(correspond_evalue == min(correspond_evalue))]
      if (length(correspond_id) == 0){
        ID_cluster[i] <- gene_info_blast$Smed.ID[match(BLAST_cluster[i],gene_info_blast$BLAST)]
      }
      else {
        ID_cluster[i] <- correspond_id
      }
    }
    return(ID_cluster)
  }
}

ID_to_BLAST <- function(ID_cluster){
  return(gene_info_blast$BLAST[match(ID_cluster, gene_info_blast$Smed.ID)])
}
