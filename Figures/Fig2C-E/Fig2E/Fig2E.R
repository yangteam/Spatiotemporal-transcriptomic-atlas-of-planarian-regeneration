library(Seurat)
library(ggplot2)
library(gridExtra)

file_path <- "../cell2location_result/"
file_list <- list.files(file_path)
tps_list <- substr(file_list, 1, nchar(file_list)-4)

deconv_score <- list()
for (i in 1:length(file_list)) {
  deconv_score[[tps_list[i]]] <- read.csv(paste0(file_path, file_list[[i]]), row.names = 1)
}

# There need Seurat object of ST data at 6 time points
# And spot name must be the same as the rownames of cell2location deconvolution results

ST_0d <- AddMetaData(ST_0d, deconv_score[["0d"]])
ST_6h <- AddMetaData(ST_6h, deconv_score[["6h"]])
ST_12h <- AddMetaData(ST_12h, deconv_score[["12h"]])
ST_1d <- AddMetaData(ST_1d, deconv_score[["1d"]])
ST_3d <- AddMetaData(ST_3d, deconv_score[["3d"]])
ST_7d <- AddMetaData(ST_7d, deconv_score[["7d"]])



mycol <- c(Neoblast = "#A9A9A9", Epidermal = "#aa40fc", Gut = "#279e68", Muscle = "#d62728", 
           Neuronal = "#ff7f0e", Parenchymal = "#e377c2", Secretory = "#1f77b4")

gg_list <- list()


gg_list[[1]] <- SpatialDimPlot(ST_0d, group.by = "hard_type", crop = T, cols = mycol)
gg_list[[2]] <- SpatialDimPlot(ST_6h, group.by = "hard_type", crop = T, cols = mycol)
gg_list[[3]] <- SpatialDimPlot(ST_12h, group.by = "hard_type", crop = T, cols = mycol)
gg_list[[4]] <- SpatialDimPlot(ST_1d, group.by = "hard_type", crop = T, cols = mycol)
gg_list[[5]] <- SpatialDimPlot(ST_3d, group.by = "hard_type", crop = T, cols = mycol)
gg_list[[6]] <- SpatialDimPlot(ST_7d, group.by = "hard_type", crop = T, cols = mycol)
for (i in 1:6) {
  for (j in 1:4) {
    gg_list[[i]][[j]] <- gg_list[[i]][[j]]+NoLegend()
  }
}


gg <- grid.arrange(gg_list[[1]][[1]], gg_list[[1]][[2]], gg_list[[1]][[3]], gg_list[[1]][[4]], 
                   gg_list[[2]][[1]], gg_list[[2]][[2]], gg_list[[2]][[3]], gg_list[[2]][[4]], 
                   gg_list[[3]][[1]], gg_list[[3]][[2]], gg_list[[3]][[3]], gg_list[[3]][[4]], 
                   gg_list[[4]][[1]], gg_list[[4]][[2]], gg_list[[4]][[3]], gg_list[[4]][[4]], 
                   gg_list[[5]][[1]], gg_list[[5]][[2]], gg_list[[5]][[3]], gg_list[[5]][[4]], 
                   gg_list[[6]][[1]], gg_list[[6]][[2]], gg_list[[6]][[3]], gg_list[[6]][[4]], 
                   ncol = 4)

if (!dir.exists("./Figures")) {dir.create("./Figures")}

ggsave("./Figures/hard_type.pdf", gg, width = 8, height = 16)

celltype.level <- c("Neoblast",
                    "Epidermal",
                    "Gut",
                    "Muscle",
                    "Neuronal",
                    "Parenchymal",
                    "Secretory")

gg_list <- list()
for (i in 1:length(celltype.level)) {
  
  gg_list[[1]] <- SpatialFeaturePlot(ST_0d, features = celltype.level[i])
  gg_list[[2]] <- SpatialFeaturePlot(ST_6h, features = celltype.level[i])
  gg_list[[3]] <- SpatialFeaturePlot(ST_12h, features = celltype.level[i])
  gg_list[[4]] <- SpatialFeaturePlot(ST_1d, features = celltype.level[i])
  gg_list[[5]] <- SpatialFeaturePlot(ST_3d, features = celltype.level[i])
  gg_list[[6]] <- SpatialFeaturePlot(ST_7d, features = celltype.level[i])
  
  gg <- grid.arrange(gg_list[[1]][[1]], gg_list[[1]][[2]], gg_list[[1]][[3]], gg_list[[1]][[4]], 
                     gg_list[[2]][[1]], gg_list[[2]][[2]], gg_list[[2]][[3]], gg_list[[2]][[4]], 
                     gg_list[[3]][[1]], gg_list[[3]][[2]], gg_list[[3]][[3]], gg_list[[3]][[4]], 
                     gg_list[[4]][[1]], gg_list[[4]][[2]], gg_list[[4]][[3]], gg_list[[4]][[4]], 
                     gg_list[[5]][[1]], gg_list[[5]][[2]], gg_list[[5]][[3]], gg_list[[5]][[4]], 
                     gg_list[[6]][[1]], gg_list[[6]][[2]], gg_list[[6]][[3]], gg_list[[6]][[4]], 
                     ncol = 4)
  ggsave(paste0("./Figures/", celltype.level[i], ".png"),gg, width = 10, height = 24)
}


