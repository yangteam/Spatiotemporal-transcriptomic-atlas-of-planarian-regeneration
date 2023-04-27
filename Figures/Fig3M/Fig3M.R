library(Seurat)
library(dplyr)
library(Matrix)
library(monocle)


# data1 <- readRDS('adata_Neoblast.rds')
# sub <- subset(data1, subset=`leiden`  %in% c("0","7","15"))

# data <- as(as.matrix(sub@assays$RNA@data), 'sparseMatrix')
# fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
# save(data,file="data.Robj",version =2)
# save(fData,file="fData.Robj",version =2)
# save(sub,file="meta.Robj",version =2)

load("data.Robj")
load("fData.Robj")
load('meta.Robj')

pd <- new('AnnotatedDataFrame', data = sub@meta.data)
fd <- new('AnnotatedDataFrame', data = fData)


#Construct monocle cds
normal <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
                         
                        

normal <- estimateSizeFactors(normal)
normal <- estimateDispersions(normal)
expressed_genes <- row.names(subset(fData(normal)))


save(normal,file="sample_normal.Robj")

###monoclemarker

normal <- detectGenes(normal, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(normal), num_cells_expressed >= 10))

diff_test_res <- differentialGeneTest(normal[expressed_genes,], fullModelFormulaStr = "~leiden")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
normal <- setOrderingFilter(normal, ordering_genes = ordering_genes)


###max_components
normal <- reduceDimension(normal, max_components = 20, method = 'DDRTree')

##
normal <- orderCells(normal,reverse=F)
save(normal,file="orderCells.Robj")



pdf(file="trajectory.pdf")
plot_cell_trajectory(normal, color_by = "Pseudotime",show_branch_points=T,cell_size=1)+scale_color_gradient(low="#C64032",high="#FDF6BA")
plot_cell_trajectory(normal, color_by = "leiden", show_branch_points=T,cell_size=1.0) + scale_color_manual(breaks = c("0", "7", "15"), values=c("#85A8CD", "#C9E0A4", "#BE80AC")) + theme(legend.position = "right")
dev.off()

