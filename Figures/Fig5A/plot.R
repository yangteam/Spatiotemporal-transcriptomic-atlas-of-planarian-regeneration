
library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)

load('./normalizeddata/cut12h.Robj')
DefaultAssay(cut12h)<- "Spatial"

#Load the Modules of cut12h data
geneTable<-read.table('./hotspot/cut12h/Module.txt',sep=" ",header=1)
rownames(geneTable)<-geneTable$Gene

# geneTable[geneTable[,"Module"]==1,"Module_new"] <- 10
# geneTable[geneTable[,"Module"]==2,"Module_new"] <- 5
# geneTable[geneTable[,"Module"]==3,"Module_new"] <- 6
# geneTable[geneTable[,"Module"]==4,"Module_new"] <- 1
# geneTable[geneTable[,"Module"]==5,"Module_new"] <- 14
# geneTable[geneTable[,"Module"]==6,"Module_new"] <- 7
# geneTable[geneTable[,"Module"]==7,"Module_new"] <- 8
# geneTable[geneTable[,"Module"]==8,"Module_new"] <- 9
# geneTable[geneTable[,"Module"]==9,"Module_new"] <- 4
# geneTable[geneTable[,"Module"]==10,"Module_new"] <- 5
# geneTable[geneTable[,"Module"]==11,"Module_new"] <- 2
# geneTable[geneTable[,"Module"]==12,"Module_new"] <- 11
# geneTable[geneTable[,"Module"]==13,"Module_new"] <- 5
# geneTable[geneTable[,"Module"]==14,"Module_new"] <- 13
# geneTable[geneTable[,"Module"]==15,"Module_new"] <-15
# geneTable[geneTable[,"Module"]==16,"Module_new"] <-12
# geneTable[geneTable[,"Module"]==17,"Module_new"] <- 3
# geneTable[geneTable[,"Module"]==18,"Module_new"] <- 14
# geneTable[geneTable[,"Module"]==19,"Module_new"] <- 13


#Calculate Score of Each Module
for (i in unique(geneTable$Module_new)){
    gene_plot<-subset(geneTable, subset=Module_new == i)
    cut12h <- AddModuleScore(object = cut12h, features = list(gene_plot$Gene), ctrl = 50,name = paste0("Module_",i))
}

cut12h_modules<-cut12h@meta.data[,c("Module_11","Module_21","Module_31","Module_41","Module_51","Module_61","Module_71","Module_81","Module_91","Module_101","Module_111","Module_121","Module_131","Module_141","Module_151")]

celltype_cell2location <- read.table("cell_location_section.txt",header=T)

rownames(celltype_cell2location)<-gsub("-1-cut12h_A1","-1_1",rownames(celltype_cell2location))
rownames(celltype_cell2location)<-gsub("-1-cut12h_B1","-1_2",rownames(celltype_cell2location))
rownames(celltype_cell2location)<-gsub("-1-cut12h_C1","-1_3",rownames(celltype_cell2location))
rownames(celltype_cell2location)<-gsub("-1-cut12h_D1","-1_4",rownames(celltype_cell2location))


cut12h_modules_celltype<-merge(celltype_cell2location,cut12h_modules,by="row.names")
rownames(cut12h_modules_celltype)<-cut12h_modules_celltype$Row.names
cut12h_modules_celltype<-cut12h_modules_celltype[,-1]

cut12h_modules_celltype_top8<-subset(cut12h_modules_celltype, subset=Section_id %in% c(1,2,3,4,5,6,7,8))
write.table(cut12h_modules_celltype_top8,file='cut12h_modules_celltype_top8.txt')

#Calculate the correlation between cell types and modules
library(corrplot)

cut12h_modules_celltype_top8<-read.table(file='cut12h_modules_celltype_top8.txt')

matrix <- cor(cut12h_modules_celltype_top8[,3:9],cut12h_modules_celltype_top8[,10:24],method="pearson",use="complete.obs")

my_colors <- brewer.pal(5, "Spectral")

my_colors <- colorRampPalette(my_colors)(100)[100:1]


library(pheatmap)
library(grid)

pdf('cor.pdf',width=10,heigh=5)
ph <- pheatmap(matrix, color=my_colors, scale="none", border="black")

ph$gtable$grobs[[1]]$gp <- gpar(lwd = 2)
ph$gtable$grobs[[2]]$gp <- gpar(lwd = 2)

grid.newpage()
grid.draw(ph$gtable)
dev.off()
