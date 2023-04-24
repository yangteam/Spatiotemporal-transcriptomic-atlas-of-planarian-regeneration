ST_sc_marker_enrich <- function(
        ST.level, 
        sc.level, 
        ST_marker.list, 
        sc_marker.list, 
        ST.logFC_thre, 
        sc.logFC_thre,
        gene_number=27664){
    
    enrich_matrix <- matrix(0, length(sc.level), length(ST.level))
    colnames(enrich_matrix) <- ST.level
    rownames(enrich_matrix) <- sc.level
    for (i in 1:length(sc.level)) {
        for (j in 1:length(ST.level)) {
            sc_marker <- sc_marker.list$gene[sc_marker.list$avg_log2FC > sc.logFC_thre & 
                                                 sc_marker.list$p_val_adj < 0.05 &
                                                 sc_marker.list$cluster == sc.level[i]]
            ST_marker <- ST_marker.list$gene[ST_marker.list$avg_log2FC > ST.logFC_thre & 
                                                 ST_marker.list$p_val_adj < 0.05 &
                                                 ST_marker.list$cluster == ST.level[j]]
            q <- length(intersect(sc_marker, ST_marker))-1
            m <- length(sc_marker)
            n <- gene_number-length(sc_marker)
            k <- length(ST_marker)
            enrich_matrix[i,j] <- round(-log10(phyper(q-1, m, n, k, lower.tail = F)), 1)
        }
    }
    return(enrich_matrix)
}

