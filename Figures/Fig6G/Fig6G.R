
gene<-read.table('gene.txt',header=F)
cd_genes <- as.character(gene$V1)
gene_name<-read.table('genename.txt',header=F)
cd_genes_name <- as.character(gene_name$V1)
cd_genes_name<-tolower(cd_genes_name)

pdf('dotplot_filter_red_1.pdf', heigh=10, width=5)
DotPlot(object = merged, features = cd_genes,group.by="type",cols=c("gray","#DC4B39"))+coord_flip()+scale_x_discrete(labels=cd_genes_name)+theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(face="italic"))
dev.off()
