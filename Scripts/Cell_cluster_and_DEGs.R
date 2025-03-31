###Cell cluster and DEGs
library(Seurat)
library(dplyr)
library(ggraph)
library(ggplot2)
library(cowplot)
library(clustree)
library(metap)

#Load Aurelia single-cell data.
AAU0 <- Read10X(data.dir = ".../count_AAU_0h/outs/filtered_feature_bc_matrix")
AAU6 <- Read10X(data.dir = ".../count_AAU_6h/outs/filtered_feature_bc_matrix")
AAU12 <- Read10X(data.dir = ".../count_AAU_12h/outs/filtered_feature_bc_matrix")
AAU18 <- Read10X(data.dir = ".../count_AAU_18h/outs/filtered_feature_bc_matrix")
AAU24 <- Read10X(data.dir = ".../count_AAU24hRNAlib/outs/filtered_feature_bc_matrix")
AAU0 <- CreateSeuratObject(counts = AAU0, project = "AAU0", min.cells = 3, min.features = 200)
AAU6 <- CreateSeuratObject(counts = AAU6, project = "AAU6", min.cells = 3, min.features = 200)
AAU12 <- CreateSeuratObject(counts = AAU12, project = "AAU12", min.cells = 3, min.features = 200)
AAU18 <- CreateSeuratObject(counts = AAU18, project = "AAU18", min.cells = 3, min.features = 200)
AAU24 <- CreateSeuratObject(counts = AAU24, project = "AAU24", min.cells = 3, min.features = 200)

#QC
AAU0 <- subset(AAU0, subset = nCount_RNA>100 & percent.mt < 5)
AAU6 <- subset(AAU6, subset = nCount_RNA>100 & percent.mt < 5)
AAU12 <- subset(AAU12, subset = nCount_RNA>100 & percent.mt < 5)
AAU18 <- subset(AAU18, subset = nCount_RNA>100 & percent.mt < 5)
AAU24 <- subset(AAU24, subset = nCount_RNA>100 & percent.mt < 5)

#Data process
AAU0.QC$stim <- "0h"
AAU6.QC$stim <- "6h"
AAU12.QC$stim <- "12h"
AAU18.QC$stim <- "18h"
AAU24.QC$stim <- "24h"
AAU0 <- NormalizeData(AAU0, normalization.method = "LogNormalize", scale.factor = 10000)
AAU0 <- FindVariableFeatures(AAU0, selection.method = "vst")
AAU6 <- NormalizeData(AAU6, normalization.method = "LogNormalize", scale.factor = 10000)
AAU6 <- FindVariableFeatures(AAU6, selection.method = "vst")
AAU12 <- NormalizeData(AAU12, normalization.method = "LogNormalize", scale.factor = 10000)
AAU12 <- FindVariableFeatures(AAU12, selection.method = "vst")
AAU18 <- NormalizeData(AAU18, normalization.method = "LogNormalize", scale.factor = 10000)
AAU18 <- FindVariableFeatures(AAU18, selection.method = "vst")
AAU24 <- NormalizeData(AAU24, normalization.method = "LogNormalize", scale.factor = 10000)
AAU24 <- FindVariableFeatures(AAU24, selection.method = "vst")

#Perform integration
features <- SelectIntegrationFeatures(object.list = list(AAU0,AAU6,AAU12,AAU18,AAU24))
regeneration.anchors <- FindIntegrationAnchors(object.list = list(AAU0, AAU6,AAU12,AAU18,AAU24), anchor.features = features)
regeneration.combined <- IntegrateData(anchorset = regeneration.anchors)

#scale data and run PCA
DefaultAssay(regeneration.combined) <- "integrated"
regeneration.combined
regeneration.combined <- ScaleData(regeneration.combined)
regeneration.combined <- RunPCA(regeneration.combined)

#t-SNE and Clustering
#find reduction dim ref to ElbowPlot.
##for Rhopilema dim=21
regeneration.combined <- RunUMAP(regeneration.combined, reduction = "pca", dims = 1:18)
regeneration.combined <- FindNeighbors(regeneration.combined, reduction = "pca", dims = 1:18)

#find resolution
regeneration.combined_FR <- FindClusters(object = regeneration.combined, resolution = c(seq(.0,1,.1)))
clustree(regeneration.combined_FR@meta.data, prefix = "integrated_snn_res.")
rm(regeneration.combined_FR)

#cell clustering
#for Rhopilema resolution = 0.2
regeneration.combined <- FindClusters(object = regeneration.combined, resolution = 0.1)
regeneration.combined$stim<-factor(regeneration.combined$stim, levels = c("0h","6h","12h","18h","24h"))

#plot
p1 <- DimPlot(regeneration.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(regeneration.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

#see proportion of each cell cluster
table(regeneration.combined@meta.data$seurat_clusters, regeneration.combined@meta.data$orig.ident)

#save
save(".../Moon.singlecell.mt5.cellranger_v7.1.0.dim18.vst.AfterFindClusters.RData")

#find marker genes of each cell cluster
#findallmarker
DefaultAssay(regeneration.combined) <- "RNA"

allmarker.regeneration.response <- FindAllMarkers(regeneration.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(allmarker.regeneration.response,
            file = ".../Moon.rg.findallmarker.txt", 
            sep = "\t", quote = FALSE,
            row.names = FALSE)

#plot
top20 <- allmarker.regeneration.response %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DefaultAssay(regeneration.combined) <- "integrated"
DoHeatmap(regeneration.combined, features = top20$gene, label = F)

#DEGs
regeneration.combined$celltype.stim <- paste(Idents(regeneration.combined), regeneration.combined$stim, sep = "_")
regeneration.combined$celltype <- Idents(regeneration.combined)
unique(regeneration.combined$celltype)
Idents(regeneration.combined) <- "celltype.stim"

#count cell numbers of each cluster, when run Findmarkers(), cell number should > 3.
unique(Idents(regeneration.combined))
cellidents <- as.data.frame(Idents(regeneration.combined))
head(cellidents)
xtabs(~.,cellidents)

#create matrix
regeneration.stim <- matrix(data = NA, ncol = 7, 
                                dimnames = list(c(1),c("Genes","p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cell.type")))
regeneration.stim

celltypenumber<-sort(unique(regeneration.combined$celltype))
treattime<-c("0h","6h","12h","18h","24h")

regeneration.stim.wide<-read.delim(".../Aurelia.gene_id", header = FALSE)
names(regeneration.stim.wide)[1]="Genes"
regeneration.stim.wide$Genes<-gsub("_","-",regeneration.stim.wide$Genes)
regeneration.stim.wide.blank<-regeneration.stim.wide
regeneration.stim.blank <- matrix(data = NA, ncol = 7, 
                                dimnames = list(c(1),c("Genes","p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cell.type")))
regeneration.stim.blank<-as.data.frame(regeneration.stim.blank)
regeneration.stim.long<-regeneration.stim.blank
regeneration.stim.long$treatment<-NA

for (j in 1:(length(treattime)-1)){
         for (k in (j+1):length(treattime)){
         treattime1<-treattime[j]
        treattime2<-treattime[k]
        regeneration.stim.sub <- regeneration.stim.blank
    
    names(regeneration.stim.sub)<-paste0(names(regeneration.stim.sub),"_",treattime2, "VS", treattime1)
    names(regeneration.stim.sub)[1]="Genes" 
    regeneration.stim.sub<-regeneration.stim.sub[,1:6]
    regeneration.stim.wide <- merge(regeneration.stim.wide, regeneration.stim.sub, by = "Genes", all.x = TRUE)
         }
    }

regeneration.stim.wide.all<-regeneration.stim.wide
regeneration.stim.wide.all<-regeneration.stim.wide.all[1,]
regeneration.stim.wide.all[1,1]<-NA
regeneration.stim.wide.all$cell.type<-NA

#run findmarker
for (i in 1:length(celltypenumber)){
            celltypesub<-celltypenumber[i]
    regeneration.stim <- regeneration.stim.blank
    regeneration.stim.wide<-regeneration.stim.wide.blank
    for (j in 1:(length(treattime)-1)){
         for (k in (j+1):length(treattime)){
         treattime1<-treattime[j]
        treattime2<-treattime[k]
        regeneration.stim.sub <- FindMarkers(regeneration.combined,
                                             ident.1 = paste(celltypesub, treattime2, sep = "_"),
                                             ident.2 = paste(celltypesub, treattime1, sep = "_"),
                                             verbose = FALSE)
    regeneration.stim.sub$Genes <- rownames(regeneration.stim.sub)
    rownames(regeneration.stim.sub) <- c()
    regeneration.stim.sub$cell.type <- paste("cell_cluster",celltypesub, sep = "_")
    
    #for table
    regeneration.stim<-regeneration.stim.sub
    regeneration.stim$treatment<-paste0(treattime2, "VS", treattime1)
    regeneration.stim.long<-rbind(regeneration.stim.long,regeneration.stim)
         }
    }
}

write.table(regeneration.stim.long, 
           file = ".../Moon.rg.findmarker.txt", 
           sep = "\t",
           row.names = FALSE)
