library(Seurat)
library(dplyr)
library(Matrix)
library(monocle)


load("orderCells.Robj")


my_pseudotime_de <- differentialGeneTest(normal,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 10)

# write.table(my_pseudotime_de,file="my_pseudotime_de.txt")

c<-subset(my_pseudotime_de, qval < 0.01)

# write.table(c[order(c[,4]),],file="my_pseudotime_de_qval.txt")

sig_gene_names <- row.names(subset(my_pseudotime_de, qval < 0.01))


BEAM_res <- BEAM(normal, branch_point = 9, cores = 10)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
#save(BEAM_res,file="BEAM_res.Robj")

pdf(file="pseudotime_gene_branch.pdf")
plot_genes_branched_heatmap(normal[row.names(subset(BEAM_res, qval < 0.01)),], branch_point = 9, num_clusters = 3, cores = 10, use_gene_short_name = T, show_rownames = T)
dev.off()

#不同分支基因
cluster_gene_branch<-plot_pseudotime_heatmap(normal[sig_gene_names,], num_clusters = 3,cores = 10,use_gene_short_name = TRUE,show_rownames = TRUE,return_heatmap = TRUE)

gene_branch<-cutree(cluster_gene_branch$tree_row, k = 3)
write.table(gene,file="cutree_branch.txt",quote=F)

