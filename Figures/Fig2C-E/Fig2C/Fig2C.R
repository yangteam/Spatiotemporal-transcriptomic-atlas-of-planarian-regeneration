library(pheatmap)
library(RColorBrewer)

source("./enrich_function.R")


file_path <- "../FindAllMarker_scRNAseq_results/"
file_list <- list.files(file_path)
tps_list <- substr(file_list, 1, nchar(file_list)-4)

all_marker.type <- list()
for (i in 1:length(file_list)) {
    all_marker.type[[tps_list[i]]] <- read.csv(paste0(file_path, file_list[[i]]), row.names = 1)
}


file_path <- "../FindAllMarker_ST_results/"
file_list <- list.files(file_path)
tps_list <- substr(file_list, 1, nchar(file_list)-4)

all_marker.ap <- list()
for (i in 1:length(file_list)) {
    all_marker.ap[[tps_list[i]]] <- read.csv(paste0(file_path, file_list[[i]]), row.names = 1)
}



interval.x.level <- unique(all_marker.ap[["0d"]]$cluster)
celltype.level <- unique(all_marker.type[["0d"]]$cluster)[-7]

type_ap_enrich.0d <- ST_sc_marker_enrich(rev(interval.x.level), celltype.level,
                                         all_marker.ap$`0d`, all_marker.type$`0d`, 0.5, 0.5)
type_ap_enrich.6h <- ST_sc_marker_enrich(rev(interval.x.level[2:8]), celltype.level, 
                                         all_marker.ap$`6h`, all_marker.type$`6h`, 0.5, 0.5)
type_ap_enrich.12h <- ST_sc_marker_enrich(rev(interval.x.level[2:8]), celltype.level, 
                                          all_marker.ap$`12h`, all_marker.type$`12h`, 0.5, 0.5)
type_ap_enrich.1d <- ST_sc_marker_enrich(rev(interval.x.level[2:8]), celltype.level, 
                                         all_marker.ap$`1d`, all_marker.type$`1d`, 0.5, 0.5)
type_ap_enrich.3d <- ST_sc_marker_enrich(rev(interval.x.level[2:8]), celltype.level, 
                                         all_marker.ap$`3d`, all_marker.type$`3d`, 0.5, 0.5)
type_ap_enrich.7d <- ST_sc_marker_enrich(rev(interval.x.level[2:8]), celltype.level, 
                                         all_marker.ap$`7d`, all_marker.type$`7d`, 0.5, 0.5)

display_type_ap_enrich.0d <- matrix("", nrow(type_ap_enrich.0d), ncol(type_ap_enrich.0d))
rownames(display_type_ap_enrich.0d) <- rownames(type_ap_enrich.0d)
colnames(display_type_ap_enrich.0d) <- colnames(type_ap_enrich.0d)
display_type_ap_enrich.0d[type_ap_enrich.0d>=-log10(0.05)] <- "*"

display_type_ap_enrich.6h <- matrix("", nrow(type_ap_enrich.6h), ncol(type_ap_enrich.6h))
rownames(display_type_ap_enrich.6h) <- rownames(type_ap_enrich.6h)
colnames(display_type_ap_enrich.6h) <- colnames(type_ap_enrich.6h)
display_type_ap_enrich.6h[type_ap_enrich.6h>=-log10(0.05)] <- "*"

display_type_ap_enrich.12h <- matrix("", nrow(type_ap_enrich.12h), ncol(type_ap_enrich.12h))
rownames(display_type_ap_enrich.12h) <- rownames(type_ap_enrich.12h)
colnames(display_type_ap_enrich.12h) <- colnames(type_ap_enrich.12h)
display_type_ap_enrich.12h[type_ap_enrich.12h>=-log10(0.05)] <- "*"

display_type_ap_enrich.1d <- matrix("", nrow(type_ap_enrich.1d), ncol(type_ap_enrich.1d))
rownames(display_type_ap_enrich.1d) <- rownames(type_ap_enrich.1d)
colnames(display_type_ap_enrich.1d) <- colnames(type_ap_enrich.1d)
display_type_ap_enrich.1d[type_ap_enrich.1d>=-log10(0.05)] <- "*"

display_type_ap_enrich.3d <- matrix("", nrow(type_ap_enrich.3d), ncol(type_ap_enrich.3d))
rownames(display_type_ap_enrich.3d) <- rownames(type_ap_enrich.3d)
colnames(display_type_ap_enrich.3d) <- colnames(type_ap_enrich.3d)
display_type_ap_enrich.3d[type_ap_enrich.3d>=-log10(0.05)] <- "*"

display_type_ap_enrich.7d <- matrix("", nrow(type_ap_enrich.7d), ncol(type_ap_enrich.7d))
rownames(display_type_ap_enrich.7d) <- rownames(type_ap_enrich.7d)
colnames(display_type_ap_enrich.7d) <- colnames(type_ap_enrich.7d)
display_type_ap_enrich.7d[type_ap_enrich.7d>=-log10(0.05)] <- "*"

# plot cell2location mean score
prop_ap.0d <- matrix(0, nrow(type_ap_enrich.0d), ncol(type_ap_enrich.0d))
rownames(prop_ap.0d) <- rownames(type_ap_enrich.0d)
colnames(prop_ap.0d) <- colnames(type_ap_enrich.0d)

prop_ap.6h <- matrix(0, nrow(type_ap_enrich.6h), ncol(type_ap_enrich.6h))
rownames(prop_ap.6h) <- rownames(type_ap_enrich.6h)
colnames(prop_ap.6h) <- colnames(type_ap_enrich.6h)

prop_ap.12h <- matrix(0, nrow(type_ap_enrich.12h), ncol(type_ap_enrich.12h))
rownames(prop_ap.12h) <- rownames(type_ap_enrich.12h)
colnames(prop_ap.12h) <- colnames(type_ap_enrich.12h)

prop_ap.1d <- matrix(0, nrow(type_ap_enrich.1d), ncol(type_ap_enrich.1d))
rownames(prop_ap.1d) <- rownames(type_ap_enrich.1d)
colnames(prop_ap.1d) <- colnames(type_ap_enrich.1d)

prop_ap.3d <- matrix(0, nrow(type_ap_enrich.3d), ncol(type_ap_enrich.3d))
rownames(prop_ap.3d) <- rownames(type_ap_enrich.3d)
colnames(prop_ap.3d) <- colnames(type_ap_enrich.3d)

prop_ap.7d <- matrix(0, nrow(type_ap_enrich.7d), ncol(type_ap_enrich.7d))
rownames(prop_ap.7d) <- rownames(type_ap_enrich.7d)
colnames(prop_ap.7d) <- colnames(type_ap_enrich.7d)



file_path <- "../cell2location_result/"
file_list <- list.files(file_path)
tps_list <- substr(file_list, 1, nchar(file_list)-4)

deconv_score <- list()
for (i in 1:length(file_list)) {
    deconv_score[[tps_list[i]]] <- read.csv(paste0(file_path, file_list[[i]]), row.names = 1)
}

for (it in colnames(prop_ap.0d)) {
    prop_ap.0d[,it] <- apply(deconv_score[["0d"]][(!is.na(deconv_score[["0d"]]$Neoblast))&(deconv_score[["0d"]]$interval.x==it), 
                                                  rownames(prop_ap.0d)], 2, mean)
}
for (it in colnames(prop_ap.6h)) {
    prop_ap.6h[,it] <- apply(deconv_score[["6h"]][which((!is.na(deconv_score[["6h"]]$Neoblast))&(deconv_score[["6h"]]$interval.x==it)), 
                                                  rownames(prop_ap.6h)], 2, mean)
    prop_ap.12h[,it] <- apply(deconv_score[["12h"]][which((!is.na(deconv_score[["12h"]]$Neoblast))&(deconv_score[["12h"]]$interval.x==it)), 
                                                    rownames(prop_ap.12h)], 2, mean)
    prop_ap.1d[,it] <- apply(deconv_score[["1d"]][which((!is.na(deconv_score[["1d"]]$Neoblast))&(deconv_score[["1d"]]$interval.x==it)), 
                                                  rownames(prop_ap.1d)], 2, mean)
    prop_ap.3d[,it] <- apply(deconv_score[["3d"]][which((!is.na(deconv_score[["3d"]]$Neoblast))&(deconv_score[["3d"]]$interval.x==it)), 
                                                  rownames(prop_ap.3d)], 2, mean)
    prop_ap.7d[,it] <- apply(deconv_score[["7d"]][which((!is.na(deconv_score[["7d"]]$Neoblast))&(deconv_score[["7d"]]$interval.x==it)), 
                                                  rownames(prop_ap.7d)], 2, mean)
}

if (!dir.exists("./Figures")) {dir.create("./Figures")}

# neoblast

display_type_ap_enrich.neoblast <- rbind(display_type_ap_enrich.0d["Neoblast",],
                                   c("", "", display_type_ap_enrich.6h["Neoblast",], ""),
                                   c("", "", display_type_ap_enrich.12h["Neoblast",], ""),
                                   c("", "", display_type_ap_enrich.1d["Neoblast",], ""),
                                   c("", "", display_type_ap_enrich.3d["Neoblast",], ""),
                                   c("", "", display_type_ap_enrich.7d["Neoblast",], ""))
rownames(display_type_ap_enrich.neoblast) <- c("0d", "6h", "12h", "1d", "3d", "7d")

prop_ap.neoblast <- rbind(prop_ap.0d["Neoblast",],
                          c(0, 0, prop_ap.6h["Neoblast",], 0),
                          c(0, 0, prop_ap.12h["Neoblast",], 0),
                          c(0, 0, prop_ap.1d["Neoblast",], 0),
                          c(0, 0, prop_ap.3d["Neoblast",], 0),
                          c(0, 0, prop_ap.7d["Neoblast",], 0))
rownames(prop_ap.neoblast) <- c("0d", "6h", "12h", "1d", "3d", "7d")


color_max <- max(prop_ap.neoblast)

pdf("Figures/Fig 2c Neoblast.pdf", width = 5, height = 3)
pheatmap(prop_ap.neoblast, cluster_cols = F, cluster_rows = F, angle_col = 45, border_color = NA, 
         display_numbers = display_type_ap_enrich.neoblast, fontsize_number = 20,
         legend_breaks = seq(0,color_max,1),legend_labels = seq(0,color_max,1),breaks = seq(0,color_max,0.1), 
         cellwidth = 20, cellheight = 12, number_format = "%.2f", fontsize = 8, number_color = "black", 
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(length(seq(0,color_max,0.1))),
         gaps_col = c(2,9), gaps_row = 1)
dev.off()



# neuronal

display_type_ap_enrich.neuronal <- rbind(display_type_ap_enrich.0d["Neuronal",],
                                         c("", "", display_type_ap_enrich.6h["Neuronal",], ""),
                                         c("", "", display_type_ap_enrich.12h["Neuronal",], ""),
                                         c("", "", display_type_ap_enrich.1d["Neuronal",], ""),
                                         c("", "", display_type_ap_enrich.3d["Neuronal",], ""),
                                         c("", "", display_type_ap_enrich.7d["Neuronal",], ""))
rownames(display_type_ap_enrich.neuronal) <- c("0d", "6h", "12h", "1d", "3d", "7d")

prop_ap.neuronal <- rbind(prop_ap.0d["Neuronal",],
                          c(0, 0, prop_ap.6h["Neuronal",], 0),
                          c(0, 0, prop_ap.12h["Neuronal",], 0),
                          c(0, 0, prop_ap.1d["Neuronal",], 0),
                          c(0, 0, prop_ap.3d["Neuronal",], 0),
                          c(0, 0, prop_ap.7d["Neuronal",], 0))
rownames(prop_ap.neuronal) <- c("0d", "6h", "12h", "1d", "3d", "7d")


color_max <- max(prop_ap.neuronal)

pdf("Figures/Fig 2c Neuronal.pdf", width = 5, height = 3)
pheatmap(prop_ap.neuronal, cluster_cols = F, cluster_rows = F, angle_col = 45, border_color = NA, 
         display_numbers = display_type_ap_enrich.neuronal, fontsize_number = 20,
         legend_breaks = seq(0,color_max,2),legend_labels = seq(0,color_max,2),breaks = seq(0,color_max,0.1), 
         cellwidth = 20, cellheight = 12, number_format = "%.2f", fontsize = 8, number_color = "black", 
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(length(seq(0,color_max,0.1))),
         gaps_col = c(2,9), gaps_row = 1)
dev.off()


# parenchymal

display_type_ap_enrich.parenchymal <- rbind(display_type_ap_enrich.0d["Parenchymal",],
                                         c("", "", display_type_ap_enrich.6h["Parenchymal",], ""),
                                         c("", "", display_type_ap_enrich.12h["Parenchymal",], ""),
                                         c("", "", display_type_ap_enrich.1d["Parenchymal",], ""),
                                         c("", "", display_type_ap_enrich.3d["Parenchymal",], ""),
                                         c("", "", display_type_ap_enrich.7d["Parenchymal",], ""))
rownames(display_type_ap_enrich.parenchymal) <- c("0d", "6h", "12h", "1d", "3d", "7d")

prop_ap.parenchymal <- rbind(prop_ap.0d["Parenchymal",],
                          c(0, 0, prop_ap.6h["Parenchymal",], 0),
                          c(0, 0, prop_ap.12h["Parenchymal",], 0),
                          c(0, 0, prop_ap.1d["Parenchymal",], 0),
                          c(0, 0, prop_ap.3d["Parenchymal",], 0),
                          c(0, 0, prop_ap.7d["Parenchymal",], 0))
rownames(prop_ap.parenchymal) <- c("0d", "6h", "12h", "1d", "3d", "7d")


color_max <- max(prop_ap.parenchymal)

pdf("Figures/Fig 2c Parenchymal.pdf", width = 5, height = 3)
pheatmap(prop_ap.parenchymal, cluster_cols = F, cluster_rows = F, angle_col = 45, border_color = NA, 
         display_numbers = display_type_ap_enrich.parenchymal, fontsize_number = 20, 
         legend_breaks = seq(0,color_max,1),legend_labels = seq(0,color_max,1),breaks = seq(0,color_max,0.1), 
         cellwidth = 20, cellheight = 12, number_format = "%.2f", fontsize = 8, number_color = "black", 
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(length(seq(0,color_max,0.1))),
         gaps_col = c(2,9), gaps_row = 1)
dev.off()


# muscle

display_type_ap_enrich.muscle <- rbind(display_type_ap_enrich.0d["Muscle",],
                                         c("", "", display_type_ap_enrich.6h["Muscle",], ""),
                                         c("", "", display_type_ap_enrich.12h["Muscle",], ""),
                                         c("", "", display_type_ap_enrich.1d["Muscle",], ""),
                                         c("", "", display_type_ap_enrich.3d["Muscle",], ""),
                                         c("", "", display_type_ap_enrich.7d["Muscle",], ""))
rownames(display_type_ap_enrich.muscle) <- c("0d", "6h", "12h", "1d", "3d", "7d")

prop_ap.muscle <- rbind(prop_ap.0d["Muscle",],
                          c(0, 0, prop_ap.6h["Muscle",], 0),
                          c(0, 0, prop_ap.12h["Muscle",], 0),
                          c(0, 0, prop_ap.1d["Muscle",], 0),
                          c(0, 0, prop_ap.3d["Muscle",], 0),
                          c(0, 0, prop_ap.7d["Muscle",], 0))
rownames(prop_ap.muscle) <- c("0d", "6h", "12h", "1d", "3d", "7d")


color_max <- max(prop_ap.muscle)

pdf("Figures/Fig 2c Muscle.pdf", width = 5, height = 3)
pheatmap(prop_ap.muscle, cluster_cols = F, cluster_rows = F, angle_col = 45, border_color = NA, 
         display_numbers = display_type_ap_enrich.muscle, fontsize_number = 20, 
         legend_breaks = seq(0,color_max,1),legend_labels = seq(0,color_max,1),breaks = seq(0,color_max,0.1), 
         cellwidth = 20, cellheight = 12, number_format = "%.2f", fontsize = 8, number_color = "black", 
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(length(seq(0,color_max,0.1))),
         gaps_col = c(2,9), gaps_row = 1)
dev.off()



