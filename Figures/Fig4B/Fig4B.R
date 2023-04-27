
library(Seurat)
library(ggplot2)
library(patchwork)

load('cut12h.Robj')
DefaultAssay(cut12h)<-"Spatial"


# geneTable<-read.table('cluster.txt',sep=" ",header=1)
# rownames(geneTable)<-geneTable$Gene

#merge and rename similiar Modules
# geneTable[geneTable[,"Module"]==1,"Module_new"] <- 10
# geneTable[geneTable[,"Module"]==2,"Module_new"] <- 5
# geneTable[geneTable[,"Module"]==3,"Module_new"] <- 6
# geneTable[geneTable[,"Module"]==4,"Module_new"] <- 1
# geneTable[geneTable[,"Module"]==5,"Module_new"] <- 15
# geneTable[geneTable[,"Module"]==6,"Module_new"] <- 7
# geneTable[geneTable[,"Module"]==7,"Module_new"] <- 8
# geneTable[geneTable[,"Module"]==8,"Module_new"] <- 9
# geneTable[geneTable[,"Module"]==9,"Module_new"] <- 4
# geneTable[geneTable[,"Module"]==10,"Module_new"] <- 5
# geneTable[geneTable[,"Module"]==11,"Module_new"] <- 2
# geneTable[geneTable[,"Module"]==12,"Module_new"] <- 11
# geneTable[geneTable[,"Module"]==13,"Module_new"] <- 5
# geneTable[geneTable[,"Module"]==14,"Module_new"] <- 13
# geneTable[geneTable[,"Module"]==15,"Module_new"] <-16
# geneTable[geneTable[,"Module"]==16,"Module_new"] <-12
# geneTable[geneTable[,"Module"]==17,"Module_new"] <- 3
# geneTable[geneTable[,"Module"]==18,"Module_new"] <- 15
# geneTable[geneTable[,"Module"]==19,"Module_new"] <- 14

# write.table(geneTable,file='cluster_new.txt')

geneTable<-read.table('cluster_new.txt',sep=" ",header=1)
rownames(geneTable)<-geneTable$Gene

for (i in unique(geneTable$Module_new)){
    gene_plot<-subset(geneTable, subset=Module_new == i)
    modules <- AddModuleScore(object = cut12h, features = list(gene_plot$Gene), ctrl = 50,name = 'gene_plot')
    m=paste0("Moduled_",i,".pdf")
    pdf(m,width=10,heigh=10)
    p1=SpatialFeaturePlot(object = modules, features="gene_plot1",min.cutoff="q20",combine = FALSE,alpha=c(0.6,0.8))
    plot_list = list()
    for (i in 1:4) {
	    p=p1[[i]]+theme(legend.title=element_blank(),legend.position="right")+scale_fill_viridis()
	    plot_list[[i]] = p}
    p1=wrap_plots(plot_list , nrow=1)
    print(p1)
    dev.off()
}

