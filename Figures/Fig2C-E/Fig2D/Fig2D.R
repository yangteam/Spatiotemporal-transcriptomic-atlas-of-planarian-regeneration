library(ggplot2)
library(tidyverse)
library(ggrepel)
library(stringr)


source("./switch_function.R")

gene_info_blast=read.csv("./gene_info_blast.csv", row.names = 1)

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



d1 <- all_marker.type$`6h`[all_marker.type$`6h`$cluster %in% c("Muscle"),]
d1$cluster <- paste0(d1$cluster, "-6hpa")
d2 <- all_marker.type$`1d`[all_marker.type$`1d`$cluster %in% c("Muscle"),]
d2$cluster <- paste0(d2$cluster, "-1dpa")

d3 <- all_marker.type$`0d`[all_marker.type$`0d`$cluster %in% c("Neuronal"),]
d3$cluster <- paste0(d3$cluster, "-0hpa")
d4 <- all_marker.type$`7d`[all_marker.type$`7d`$cluster %in% c("Neuronal"),]
d4$cluster <- paste0(d4$cluster, "-7dpa")

d5 <- all_marker.type$`6h`[all_marker.type$`6h`$cluster %in% c("Neoblast","Parenchymal"),]
d5$cluster <- paste0(d5$cluster, "-6hpa")
d6 <- all_marker.type$`1d`[all_marker.type$`1d`$cluster %in% c("Neoblast","Parenchymal"),]
d6$cluster <- paste0(d6$cluster, "-1dpa")

df <- rbind(d1,d2,d3,d4,d5,d6)

a1 <- intersect(all_marker.type$`6h`$gene[all_marker.type$`6h`$cluster %in% c("Muscle", "Parenchymal")],
               all_marker.ap$`6h`$gene[all_marker.ap$`6h`$cluster == "AP-8"])
a2 <- intersect(all_marker.type$`1d`$gene[all_marker.type$`1d`$cluster %in% c("Muscle", "Parenchymal")],
               all_marker.ap$`1d`$gene[all_marker.ap$`1d`$cluster == "AP-8"])

a3 <- intersect(all_marker.type$`0d`$gene[all_marker.type$`0d`$cluster == "Neuronal"],
               all_marker.ap$`0d`$gene[all_marker.ap$`0d`$cluster == "AP-1"])
a4 <- intersect(all_marker.type$`7d`$gene[all_marker.type$`7d`$cluster == "Neuronal"],
               all_marker.ap$`7d`$gene[all_marker.ap$`7d`$cluster == "AP-2"])

a5 <- intersect(all_marker.type$`6h`$gene[all_marker.type$`6h`$cluster %in% c("Neoblast", "Parenchymal")],
               all_marker.ap$`6h`$gene[all_marker.ap$`6h`$cluster == "AP-2"])
a6 <- intersect(all_marker.type$`1d`$gene[all_marker.type$`1d`$cluster %in% c("Neoblast", "Parenchymal")],
               all_marker.ap$`1d`$gene[all_marker.ap$`1d`$cluster == "AP-2"])

select_gene <- unique(c(a1,a2,a3,a4,a5,a6))

df <- df[which(df$gene %in% select_gene),]


df$cluster[which(df$cluster=="Neuronal-0hpa")] <- "Neuronal-0hpa-head"
df$cluster[which(df$cluster=="Neuronal-7dpa")] <- "Neuronal-7dpa-head"
df$cluster[which(df$cluster=="Neoblast-6hpa")] <- "Neoblast-6hpa-wound-A"
df$cluster[which(df$cluster=="Neoblast-1dpa")] <- "Neoblast-1dpa-wound-A"
df$cluster[which(df$cluster=="Muscle-6hpa")] <- "Muscle-6hpa-wound-P"
df$cluster[which(df$cluster=="Muscle-1dpa")] <- "Muscle-1dpa-wound-P"
df$cluster[which(df$cluster=="Parenchymal-6hpa")] <- "Parenchymal-6hpa-wounds"
df$cluster[which(df$cluster=="Parenchymal-1dpa")] <- "Parenchymal-1dpa-wounds"

df=rbind(df,d1)

unique(df$cluster)

df <- df[abs(df$avg_log2FC)>0.5,]
df <- df[df$avg_log2FC>=(-4),]
df <- df[df$p_val_adj<0.05,]

df$name <- ID_to_BLAST(df$gene)
df$name <- str_to_lower(df$name)
df$cluster <- factor(df$cluster, levels = c("Neoblast-6hpa-wound-A", "Neoblast-1dpa-wound-A",
                                            "Neuronal-0hpa-head", "Neuronal-7dpa-head",
                                            "Muscle-6hpa", "Muscle-6hpa-wound-P", "Muscle-1dpa-wound-P", 
                                            "Parenchymal-6hpa-wounds", "Parenchymal-1dpa-wounds"))
df$name[df$name=="notc1"]="notch1"
df$name[df$name=="clat"]="chat"
df$name[df$name=="lozen"]="runt-1"
df$name[df$name=="myod1"]="myoD"
df$name[df$name=="piwi1"]="smedwi-1"

df <- df[which(!is.na(df$name)),]

df$label <- "ordinary"


key_marker_gene <- read.table("./key_gene_list.txt")
key_marker_gene$V3 <- paste0(key_marker_gene$V1, "-", key_marker_gene$V2)

df_key <- df[which(df$name %in% key_marker_gene$V1),]

df_key$cell_type <- NA
for (i in 1:dim(df_key)[1]) {
    df_key$cell_type[i] <- paste0(df_key$name[i], "-", str_split(df_key$cluster, "-")[[i]][1])
}

df_key <- df_key[which(df_key$cell_type %in% key_marker_gene$V3),]
df_key$label <- "key"

df_key <- df_key[df_key$avg_log2FC>0,]


df %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> df_top10_up

df_top10_up$cell_type <- NA
for (i in 1:dim(df_top10_up)[1]) {
    df_top10_up$cell_type[i] <- paste0(df_top10_up$name[i], "-", str_split(df_top10_up$cluster, "-")[[i]][1])
}

df_top10_up <- df_top10_up[!(df_top10_up$gene %in% df_key$gene),]
df_top10_up$label <- "ordinary"

df$label[df$name %in% df_key$name] <- "key"

df_label <- rbind(df_key, df_top10_up)
df_label <- df_label[sample(dim(df_label)[1],dim(df_label)[1]),]
df_label <- df_label[!duplicated(df_label$cell_type),]

dfbar<-data.frame(x=unique(df$cluster),y=rep(0, length(unique(df$cluster))))
dfbar1<-data.frame(x=unique(df$cluster),y=rep(0, length(unique(df$cluster))))

for (i in 1:length(unique(df$cluster))) {
    dfbar$y[i]=min(df$avg_log2FC[df$cluster==unique(df$cluster)[i]])
    dfbar1$y[i]=max(df$avg_log2FC[df$cluster==unique(df$cluster)[i]])
}


dfcol<-data.frame(x=levels(df$cluster),
                  y=0,
                  label=levels(df$cluster))

#mycol <- c(Neoblast = "#A9A9A9", Epidermal = "#aa40fc", Gut = "#279e68", Muscle = "#d62728", 
#           Neuronal = "#ff7f0e", Parenchymal = "#e377c2", Secretory = "#1f77b4")

mycol <- c(rep("#A9A9A9", 2), rep("#ff7f0e", 2), rep("#d62728",3), rep("#e377c2",2))
names(mycol) <- levels(df$cluster)

pdf("./Figures.pdf", width = 18, height = 24)
ggplot()+
    geom_col(data = dfbar,
             mapping = aes(x = x,y = y),
             fill = "#dcdcdc",alpha = 0.6)+
    geom_col(data = dfbar1,
             mapping = aes(x = x,y = y),
             fill = "#dcdcdc",alpha = 0.6)+
    geom_jitter(data = df, position = position_jitter(seed = 123), 
                aes(x = cluster, y = avg_log2FC, color = label, size=label),
    )+
    geom_jitter(data = df_label, position = position_jitter(seed = 123), 
                aes(x = cluster, y = avg_log2FC, color = label, size=label),
    ) + 
    geom_tile(data = dfcol,
              aes(x=x,y=y),
              height=0.6,
              color="black",
              fill = mycol,
              alpha = 1,
              show.legend = F)+
    geom_text_repel(
        data=df_label, position = position_jitter(seed = 123), 
        aes(x=cluster,y=avg_log2FC,label=name, color=label, fontface = "italic"), size=18,
        force = 2, force_pull = 2, 
        arrow = arrow(length = unit(0.008, "npc"),
                      type = "open", ends = "last"),
        max.overlaps = 1000) +
    #geom_text(
    #    data=df_key, position = position_jitter(seed = 1), 
    #   aes(x=cluster,y=avg_log2FC,label=name)) +
    scale_color_manual(name=NULL, values = c("red","grey50"))+
    scale_size_manual(values = c(3,1))+
    labs(x="",y="")+
    guides(color="none", size="none")+
    #geom_text(data=dfcol,
    #          aes(x=x,y=y,label=label),
    #          size = 3,
    #          color ="white")+
    theme_minimal()+
    theme(axis.title = element_text(size = 20,
                                    color = "black",
                                    face = "bold"),
          axis.line.y = element_line(color = "black",
                                     size = 1.2),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "top",
          legend.direction = "vertical",
          legend.justification = c(1,0),
          legend.text = element_text(size = 15),
          axis.text.y = element_text(size = 50))
dev.off()


