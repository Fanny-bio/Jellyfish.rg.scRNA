###Monocle 3 for Aurelia
library(monocle3)
library(Seurat)
library(patchwork)
library(dplyr)
library(ggraph)
library(ggplot2)
library(cowplot)
library(sf)

load(".../Moon.singlecell.mt5.cellranger_v7.1.0.dim18.vst.AfterFindClusters.RData")

##Create monocle objects 
data <- GetAssayData(regeneration.combined, 
                     assay = 'RNA', slot = 'counts')
cell_metadata <- regeneration.combined@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
##preprocess data
cds <- preprocess_cds(cds)
cds <- align_cds(cds, alignment_group = "stim")
cds <- reduce_dimension(cds)
#Import integrated umap coordinates from seurat
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(regeneration.combined, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

## Clustering
cds <- cluster_cells(cds)

#replace monocle clusters with seurat
cds@clusters$UMAP$clusters <- regeneration.combined@meta.data$seurat_clusters
names(cds@clusters$UMAP$clusters) <- rownames(regeneration.combined@meta.data)

#Learn graph
cds <- learn_graph(cds, use_partition=FALSE, close_loop=FALSE)
##Color cells by pseudotime
###for Rhopilema clusters(cds) == 11
cds <- order_cells(cds, root_cells = colnames(cds[,colData(cds)$stim=="0h" & clusters(cds) == 4]))

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
#plot
custom_order <- c("0","11","1","7","5","6","2","3","8","9","10","4")
custom_colors <- c("0" = "#FF99CC", #"Epidermal muscle",
                   "11" = "#FF99CC", #"Epidermal muscle",
                   "1" = "#FFCCFF",#"Gastrodermis",
                   "7" = "#D9D9D9", #"Potential gastrodermis",
                   "5" = "#61CBF4",#"Neural cell",
                   "6" = "#61CBF4",#"Neural cell",
                   "2" = "#D9D9D9",#"Potential neural cell",
                   "3" = "#8ED973", #"Cnidocyte",
                   "8" = "#8ED973", #"Cnidocyte",
                   "9" = "#8ED973", #"Cnidocyte",
                   "10" = "#FF9966", #"Gland cell",
                   "4" = "#C67BFF") #"Stem cell",

data.pseudo$seurat_clusters <- factor(data.pseudo$seurat_clusters, levels = custom_order)
p <- ggplot(data.pseudo, 
            aes(monocle3_pseudotime, 
                seurat_clusters, 
                fill = seurat_clusters))+
  geom_boxplot()+
  scale_fill_manual(values = custom_colors)+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        #panel.border = element_blank())
