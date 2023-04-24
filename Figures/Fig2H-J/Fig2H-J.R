require(RNAMagnet)
require(ggplot2)
library(stringr)
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(circlize)
library(MASS)

setwd("G:/Work/planarian/github/Fig 2h-j/RNAmagnet/")

# There need Seurat object of scRNA-seq data of 8 time points

gene_info_blast <- read.csv("gene_info_blast.csv", row.names = 1)

source("./switch_function.R")

# Run RNAmagnet

keep_gene_blast <- unique(na.omit(gene_info_blast$BLAST))
keep_gene_ID <- BLAST_to_ID(keep_gene_blast, type = "best e-value")

keep_gene_blast <- keep_gene_blast[keep_gene_ID %in% rownames(sc_0d)]
keep_gene_ID <- keep_gene_ID[keep_gene_ID %in% rownames(sc_0d)]


ligrec <- getLigandsReceptors(version = "latest",
                              cellularCompartment = c("Membrane", "ECM", "Both", "Secreted"),
                              manualAnnotation = "Correct", 
                              ligandClass = c("Other", "Cytokine","Chemokine", "GrowthFactor", "Interleukin"))
head(ligrec)

# prepare seurat object, data process
# 0d
sc_0d_count.blast <- sc_0d@assays$RNA@counts[keep_gene_ID,]
rownames(sc_0d_count.blast) <- str_to_title(keep_gene_blast)
sc_0d_count.blast <- sc_0d_count.blast[rowSums(sc_0d_count.blast)>0,]
sc_0d.blast <- CreateSeuratObject(sc_0d_count.blast)
sc_0d.blast <- AddMetaData(sc_0d.blast, sc_0d@meta.data[,4:5])
Idents(sc_0d.blast) <- sc_0d.blast$cell_type
#Idents(sc_0d.blast) <- sc_0d.blast$cell_subtype
sc_0d.blast <- sc_0d.blast %>% 
    NormalizeData() %>% 
    FindVariableFeatures(nfeatures = 3000) %>% 
    ScaleData() %>%
    RunPCA() %>% 
    RunUMAP(dims = 1:30) %>%
    RunTSNE(dims = 1:30)
# 6h
sc_6h_count.blast <- sc_6h@assays$RNA@counts[keep_gene_ID,]
rownames(sc_6h_count.blast) <- str_to_title(keep_gene_blast)
sc_6h_count.blast <- sc_6h_count.blast[rowSums(sc_6h_count.blast)>0,]
sc_6h.blast <- CreateSeuratObject(sc_6h_count.blast)
sc_6h.blast <- AddMetaData(sc_6h.blast, sc_6h@meta.data[,4:5])
Idents(sc_6h.blast) <- sc_6h.blast$cell_type
#Idents(sc_6h.blast) <- sc_6h.blast$cell_subtype
sc_6h.blast <- sc_6h.blast %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData() %>%
  RunPCA() %>% 
  RunUMAP(dims = 1:30) %>%
  RunTSNE(dims = 1:30)
# 12h
sc_12h_count.blast <- sc_12h@assays$RNA@counts[keep_gene_ID,]
rownames(sc_12h_count.blast) <- str_to_title(keep_gene_blast)
sc_12h_count.blast <- sc_12h_count.blast[rowSums(sc_12h_count.blast)>0,]
sc_12h.blast <- CreateSeuratObject(sc_12h_count.blast)
sc_12h.blast <- AddMetaData(sc_12h.blast, sc_12h@meta.data[,4:5])
Idents(sc_12h.blast) <- sc_12h.blast$cell_type
#Idents(sc_12h.blast) <- sc_12h.blast$cell_subtype
sc_12h.blast <- sc_12h.blast %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData() %>%
  RunPCA() %>% 
  RunUMAP(dims = 1:30) %>%
  RunTSNE(dims = 1:30)
# 1d
sc_1d_count.blast <- sc_1d@assays$RNA@counts[keep_gene_ID,]
rownames(sc_1d_count.blast) <- str_to_title(keep_gene_blast)
sc_1d_count.blast <- sc_1d_count.blast[rowSums(sc_1d_count.blast)>0,]
sc_1d.blast <- CreateSeuratObject(sc_1d_count.blast)
sc_1d.blast <- AddMetaData(sc_1d.blast, sc_1d@meta.data[,4:5])
Idents(sc_1d.blast) <- sc_1d.blast$cell_type
#Idents(sc_1d.blast) <- sc_1d.blast$cell_subtype
sc_1d.blast <- sc_1d.blast %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData() %>%
  RunPCA() %>% 
  RunUMAP(dims = 1:30) %>%
  RunTSNE(dims = 1:30)
# 2d
sc_2d_count.blast <- sc_2d@assays$RNA@counts[keep_gene_ID,]
rownames(sc_2d_count.blast) <- str_to_title(keep_gene_blast)
sc_2d_count.blast <- sc_2d_count.blast[rowSums(sc_2d_count.blast)>0,]
sc_2d.blast <- CreateSeuratObject(sc_2d_count.blast)
sc_2d.blast <- AddMetaData(sc_2d.blast, sc_2d@meta.data[,4:5])
Idents(sc_2d.blast) <- sc_2d.blast$cell_type
#Idents(sc_2d.blast) <- sc_2d.blast$cell_subtype
sc_2d.blast <- sc_2d.blast %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData() %>%
  RunPCA() %>% 
  RunUMAP(dims = 1:30) %>%
  RunTSNE(dims = 1:30)
# 3d
sc_3d_count.blast <- sc_3d@assays$RNA@counts[keep_gene_ID,]
rownames(sc_3d_count.blast) <- str_to_title(keep_gene_blast)
sc_3d_count.blast <- sc_3d_count.blast[rowSums(sc_3d_count.blast)>0,]
sc_3d.blast <- CreateSeuratObject(sc_3d_count.blast)
sc_3d.blast <- AddMetaData(sc_3d.blast, sc_3d@meta.data[,4:5])
Idents(sc_3d.blast) <- sc_3d.blast$cell_type
#Idents(sc_3d.blast) <- sc_3d.blast$cell_subtype
sc_3d.blast <- sc_3d.blast %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData() %>%
  RunPCA() %>% 
  RunUMAP(dims = 1:30) %>%
  RunTSNE(dims = 1:30)
# 5d
sc_5d_count.blast <- sc_5d@assays$RNA@counts[keep_gene_ID,]
rownames(sc_5d_count.blast) <- str_to_title(keep_gene_blast)
sc_5d_count.blast <- sc_5d_count.blast[rowSums(sc_5d_count.blast)>0,]
sc_5d.blast <- CreateSeuratObject(sc_5d_count.blast)
sc_5d.blast <- AddMetaData(sc_5d.blast, sc_5d@meta.data[,4:5])
Idents(sc_5d.blast) <- sc_5d.blast$cell_type
#Idents(sc_5d.blast) <- sc_5d.blast$cell_subtype
sc_5d.blast <- sc_5d.blast %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData() %>%
  RunPCA() %>% 
  RunUMAP(dims = 1:30) %>%
  RunTSNE(dims = 1:30)
# 7d
sc_7d_count.blast <- sc_7d@assays$RNA@counts[keep_gene_ID,]
rownames(sc_7d_count.blast) <- str_to_title(keep_gene_blast)
sc_7d_count.blast <- sc_7d_count.blast[rowSums(sc_7d_count.blast)>0,]
sc_7d.blast <- CreateSeuratObject(sc_7d_count.blast)
sc_7d.blast <- AddMetaData(sc_7d.blast, sc_7d@meta.data[,4:5])
Idents(sc_7d.blast) <- sc_7d.blast$cell_type
#Idents(sc_7d.blast) <- sc_7d.blast$cell_subtype
sc_7d.blast <- sc_7d.blast %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 3000) %>% 
  ScaleData() %>%
  RunPCA() %>% 
  RunUMAP(dims = 1:30) %>%
  RunTSNE(dims = 1:30)


result_0d <- RNAMagnetSignaling(sc_0d.blast, .version = "latest", .minExpression = 10,
                                .cellularCompartment=c('Membrane', 'ECM', 'Secreted', 'Both'))
result_6h <- RNAMagnetSignaling(sc_6h.blast, .version = "latest", .minExpression = 10,
                                .cellularCompartment=c('Membrane', 'ECM', 'Secreted', 'Both'))
result_12h <- RNAMagnetSignaling(sc_12h.blast, .version = "latest", .minExpression = 10,
                                 .cellularCompartment=c('Membrane', 'ECM', 'Secreted', 'Both'))
result_1d <- RNAMagnetSignaling(sc_1d.blast, .version = "latest", .minExpression = 10,
                                .cellularCompartment=c('Membrane', 'ECM', 'Secreted', 'Both'))
result_2d <- RNAMagnetSignaling(sc_2d.blast, .version = "latest", .minExpression = 10,
                                .cellularCompartment=c('Membrane', 'ECM', 'Secreted', 'Both'))
result_3d <- RNAMagnetSignaling(sc_3d.blast, .version = "latest", .minExpression = 10,
                                .cellularCompartment=c('Membrane', 'ECM', 'Secreted', 'Both'))
result_5d <- RNAMagnetSignaling(sc_5d.blast, .version = "latest", .minExpression = 10,
                                .cellularCompartment=c('Membrane', 'ECM', 'Secreted', 'Both'))
result_7d <- RNAMagnetSignaling(sc_7d.blast, .version = "latest", .minExpression = 10,
                                .cellularCompartment=c('Membrane', 'ECM', 'Secreted', 'Both'))



dir.create("./LR_pair")

for (i in names(table(sc_0d.blast$cell_type))) {
  for (j in names(table(sc_0d.blast$cell_type))) {
    LR_pair <- data.frame(time=0, score=0, pair=0)
    LR_pair <- rbind(LR_pair, cbind(time="0d", getRNAMagnetGenes(result_0d, i, j, thresh = 0)))
    LR_pair <- rbind(LR_pair, cbind(time="6h", getRNAMagnetGenes(result_6h, i, j, thresh = 0)))
    LR_pair <- rbind(LR_pair, cbind(time="12h", getRNAMagnetGenes(result_12h, i, j, thresh = 0)))
    LR_pair <- rbind(LR_pair, cbind(time="1d", getRNAMagnetGenes(result_1d, i, j, thresh = 0)))
    LR_pair <- rbind(LR_pair, cbind(time="2d", getRNAMagnetGenes(result_2d, i, j, thresh = 0)))
    LR_pair <- rbind(LR_pair, cbind(time="3d", getRNAMagnetGenes(result_3d, i, j, thresh = 0)))
    LR_pair <- rbind(LR_pair, cbind(time="5d", getRNAMagnetGenes(result_5d, i, j, thresh = 0)))
    LR_pair <- rbind(LR_pair, cbind(time="7d", getRNAMagnetGenes(result_7d, i, j, thresh = 0)))
    LR_pair <- LR_pair[-1,]
    
    LR_pair <- LR_pair[LR_pair$score>=0.01,]
    if (dim(LR_pair)[1]>0){
      write.csv(LR_pair, file = paste0("./LR_pair/", i, "_to_", j, ".csv"), row.names = F)
    }
  }
}


colors <- c(Epidermal = "#aa40fc", Gut = "#279e68",
            Muscle = "#d62728", Neoblast = "#A9A9A9",
            Neuronal = "#ff7f0e", Parenchymal = "#e377c2",
            Secretory = "#1f77b4")


pops <- names(colors)


# plot chordDiagram for Supplementary Fig 4c
if (!dir.exists("./Fig S4c")) {dir.create("./Fig S4c")}

i=0
for (result in c(result_0d, result_6h, result_12h, result_1d, result_2d, result_3d, result_5d, result_7d)) {
  i=i+1
  mean_by_pop <- sapply(pops, function(id) {
    sapply(pops, function(id2) {
      mean(result@specificity[result@celltype == id2, id])
    })
  })
  
  mean_by_pop[mean_by_pop<=0.01]=0
  
  pdf(paste0("./Fig S4c/time points ", as.character(i), ".pdf"), width = 5, height = 5)
  
  diag(mean_by_pop)=0
  par(cex=2, mar=c(0,0,0,0))
  chordDiagram(mean_by_pop, symmetric = F, grid.col = colors, annotationTrackHeight = c(0.1),
               annotationTrack = c("grid"))
  dev.off()
}


# There need Seurat object of ST data for 6 time points

ST_0d_count.blast <- ST_0d@assays$Spatial@counts[keep_gene_ID,]
rownames(ST_0d_count.blast) <- str_to_title(keep_gene_blast)

ST_6h_count.blast <- ST_6h@assays$Spatial@counts[keep_gene_ID,]
rownames(ST_6h_count.blast) <- str_to_title(keep_gene_blast)

ST_12h_count.blast <- ST_12h@assays$Spatial@counts[keep_gene_ID,]
rownames(ST_12h_count.blast) <- str_to_title(keep_gene_blast)

ST_1d_count.blast <- ST_1d@assays$Spatial@counts[keep_gene_ID,]
rownames(ST_1d_count.blast) <- str_to_title(keep_gene_blast)

ST_3d_count.blast <- ST_3d@assays$Spatial@counts[keep_gene_ID,]
rownames(ST_3d_count.blast) <- str_to_title(keep_gene_blast)

ST_7d_count.blast <- ST_7d@assays$Spatial@counts[keep_gene_ID,]
rownames(ST_7d_count.blast) <- str_to_title(keep_gene_blast)

# Fig 2h
# spatial colocolization degree 
source("./spatial_colocalization_function.R")

l11 <- spatial_colocalization_degree("Ncam1", "Fgfr1", ST_0d_count.blast)
l12 <- spatial_colocalization_degree("Ncam1", "Fgfr1", ST_6h_count.blast)
l13 <- spatial_colocalization_degree("Ncam1", "Fgfr1", ST_12h_count.blast)
l14 <- spatial_colocalization_degree("Ncam1", "Fgfr1", ST_1d_count.blast)

l21 <- spatial_colocalization_degree("Bmp7", "Acvr1", ST_0d_count.blast)
l22 <- spatial_colocalization_degree("Bmp7", "Acvr1", ST_6h_count.blast)
l23 <- spatial_colocalization_degree("Bmp7", "Acvr1", ST_12h_count.blast)
l24 <- spatial_colocalization_degree("Bmp7", "Acvr1", ST_1d_count.blast)


df=data.frame(x=rep(c("0 hpa", "6 hpa", "12 hpa", "24 hpa"), 2), 
              y=c(l21, l22, l23, l24, l11, l12, l13, l14),
              LR_pair=c(rep("bmp7-acvr1",4), rep("ncam1-fgfr1",4)))

df$x=factor(df$x, levels = c("0 hpa", "6 hpa", "12 hpa", "24 hpa"))

if (!dir.exists("./Figures")) {dir.create("./Figures")}

gg <- ggplot(df, aes(x=x,y=y,group=LR_pair))+
    geom_point(aes(color=LR_pair), size=3)+
    geom_line(aes(color=LR_pair), size=1)+
    theme_test()+
    labs(x="",y="")+
    theme(#legend.position = c(0.75,0.7),
          legend.title = element_text(size = 8),
          legend.text  = element_text(size = 8),
          legend.key.size = unit(1.5, "lines"),
          axis.text.x = element_text(size = 12, color = "black", angle = 30, hjust = 1,vjust = 1),
          axis.text.y = element_text(size = 10, color = "black"),
          panel.border = element_rect(fill=NA, color="black", size=1.5, linetype="solid"))
ggsave(gg, filename = "./Figures/Fig 2h.pdf", width = 3.5, height = 3)

# Fig 2i
# dotplot
pdf("Figures/Fig 2i.pdf", width = 4, height = 3)
DotPlot(sc_0d, features = c("SMED30018873", "SMED30009204"), group.by = "cell_type")+
    scale_size(range = c(5, 5))+
    scale_color_gradientn(colors = rev(brewer.pal(11, 'RdYlBu')[2:10]))
dev.off()

# Fig 2j
# LR pair double spatial visualization

library(stringr)
library(ggplot2)
library(Seurat)

align_coord_0d=read.csv("./3D_align_result/0d.csv", row.names = 1)
align_coord_6h=read.csv("./3D_align_result/6h.csv", row.names = 1)
align_coord_12h=read.csv("./3D_align_result/12h.csv", row.names = 1)
align_coord_1d=read.csv("./3D_align_result/1d.csv", row.names = 1)
align_coord_3d=read.csv("./3D_align_result/3d.csv", row.names = 1)
align_coord_7d=read.csv("./3D_align_result/7d.csv", row.names = 1)

# the rownames of 3D align result must be same as spot name of ST seurat object

ST_0d <- AddMetaData(ST_0d, align_coord_0d)
ST_6h <- AddMetaData(ST_6h, align_coord_6h)
ST_12h <- AddMetaData(ST_12h, align_coord_12h)
ST_1d <- AddMetaData(ST_1d, align_coord_1d)
ST_3d <- AddMetaData(ST_3d, align_coord_3d)
ST_7d <- AddMetaData(ST_7d, align_coord_7d)


dimreduc_matrix <- as.matrix(ST_0d@meta.data[,c("align.y","align.x")])
colnames(dimreduc_matrix) <- c("align_1", "align_2")
rownames(dimreduc_matrix) <- rownames(ST_0d@meta.data)
ST_0d[["align"]] <- CreateDimReducObject(embeddings = dimreduc_matrix, key = "align_", assay = "Spatial")

dimreduc_matrix <- as.matrix(ST_6h@meta.data[,c("align.y","align.x")])
colnames(dimreduc_matrix) <- c("align_1", "align_2")
rownames(dimreduc_matrix) <- rownames(ST_6h@meta.data)
ST_6h[["align"]] <- CreateDimReducObject(embeddings = dimreduc_matrix, key = "align_", assay = "Spatial")

dimreduc_matrix <- as.matrix(ST_12h@meta.data[,c("align.y","align.x")])
colnames(dimreduc_matrix) <- c("align_1", "align_2")
rownames(dimreduc_matrix) <- rownames(ST_12h@meta.data)
ST_12h[["align"]] <- CreateDimReducObject(embeddings = dimreduc_matrix, key = "align_", assay = "Spatial")

dimreduc_matrix <- as.matrix(ST_1d@meta.data[,c("align.y","align.x")])
colnames(dimreduc_matrix) <- c("align_1", "align_2")
rownames(dimreduc_matrix) <- rownames(ST_1d@meta.data)
ST_1d[["align"]] <- CreateDimReducObject(embeddings = dimreduc_matrix, key = "align_", assay = "Spatial")




ligand="NCAM1"
receptor="FGFR1"
ligand.id=BLAST_to_ID(ligand,type = "best e-value")
receptor.id=BLAST_to_ID(receptor,type = "best e-value")
ligand=str_to_lower(ligand)
receptor=str_to_lower(receptor)


gg <- FeaturePlot(ST_0d[,which(ST_0d$align.z==1)], cols = c("grey", "red", "darkblue"), pt.size = 1.8,
                  features = c(ligand.id, receptor.id), reduction = "align", blend = T, blend.threshold = 0)
gg[[1]] <- gg[[1]]+coord_fixed()+theme_void()+NoLegend()+labs(title = ligand)
gg[[2]] <- gg[[2]]+coord_fixed()+theme_void()+NoLegend()+labs(title = receptor)
gg[[3]] <- gg[[3]]+coord_fixed()+theme_void()+NoLegend()+labs(title = paste0(ligand,"-",receptor))
gg[[4]] <- gg[[4]]+NoLegend()+coord_fixed(ratio = 3)+labs(title = "",x=ligand,y=receptor) +
    theme(axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25))
ggsave(paste0("./Figures/Fig 2j 0d sec 1.pdf"), gg, width = 9.7, height = 5.7)




gg <- FeaturePlot(ST_12h[,which(ST_12h$align.z==3)], cols = c("grey", "red", "darkblue"), pt.size = 2,
                  features = c(ligand.id, receptor.id), reduction = "align", blend = T,blend.threshold = 0)
gg[[1]] <- gg[[1]]+coord_fixed()+theme_void()+NoLegend()+labs(title = ligand)
gg[[2]] <- gg[[2]]+coord_fixed()+theme_void()+NoLegend()+labs(title = receptor)
gg[[3]] <- gg[[3]]+coord_fixed()+theme_void()+NoLegend()+labs(title = paste0(ligand,"-",receptor))
gg[[4]] <- gg[[4]]+NoLegend()+coord_fixed(ratio = 3)+labs(title = "",x=ligand,y=receptor) +
    theme(axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25))
ggsave(paste0("./Figures/Fig 2j 12h sec 3.pdf"), gg, width = 9.7, height = 5.7)



gg <- FeaturePlot(ST_1d[,which(ST_1d$align.z==2)], cols = c("grey", "red", "darkblue"), pt.size = 2.4,
                  features = c(ligand.id, receptor.id), reduction = "align", blend = T,blend.threshold = 0)
gg[[1]] <- gg[[1]]+coord_fixed()+theme_void()+NoLegend()+labs(title = ligand)
gg[[2]] <- gg[[2]]+coord_fixed()+theme_void()+NoLegend()+labs(title = receptor)
gg[[3]] <- gg[[3]]+coord_fixed()+theme_void()+NoLegend()+labs(title = paste0(ligand,"-",receptor))
gg[[4]] <- gg[[4]]+NoLegend()+coord_fixed(ratio = 3)+labs(title = "",x=ligand,y=receptor) +
    theme(axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25))
ggsave(paste0("./Figures/Fig 2j 1d sec 2.pdf"), gg, width = 9.7, height = 5.7)




