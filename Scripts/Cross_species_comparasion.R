###Cross species comparasion using KLD caculation
library(LaplacesDemon)
library(preprocessCore)
library(pheatmap)
library(dplyr)
library(tidyr)

setwd(".../Synteny_KLD")

##Seurat to Seurat Aceo(indirect comparasion)
species<-c("Hvul",
           "Spis",
           "Nvec", 
           "Xesp",
           "Chem",
           "Resc" 
           )

#py = Aceo cell type
py<-read.delim(paste0("all.marker.wide.", "Aceo.",".txt"), header = TRUE)
py$gene_id<-gsub("-", "_", py$gene_id)
py.bk<-py

for (n in 1:length(species)) {
  #ortholog pair
  ortholog<-read.delim(paste0("ortholog.Aurelia__v__", species[n], ".tsv.txt"), header = TRUE)
  ortholog<-separate(ortholog, gene1, c("gene","t"), sep = "-", remove = F)
  ortholog<-ortholog[,c("gene","gene2")]
  ortholog$gene<-gsub("-", "_", ortholog$gene)
  ortholog<-unique(ortholog)
  
  names(ortholog)<-c("gene_id", "gene2")

#px = ref cell type
  px<-read.delim(paste0("all.marker.wide.", species[n], ".txt"), header = TRUE)  
  names(px)[names(px)=="gene_id"]<-"gene"

if (species[n]=="Hvul") {
    ortholog<-separate(ortholog, gene2, c("sca","gene"), sep = "\\.g", remove = F)
    ortholog<-ortholog[,c("gene_id","gene")]
    ortholog$gene<-paste0("g",ortholog$gene)
    ortholog$gene<-gsub("\\.t","-",ortholog$gene)
    px$gene<-gsub("Hvul-", "", px$gene)
  }
  
  if (species[n]=="Spis") {
    names(ortholog)<-c("gene_id","gene")
    ortholog$gene<-gsub("\\.","_",ortholog$gene)
    px$gene<-gsub("Spis-", "", px$gene)
    px$gene<-gsub("-", "_", px$gene)
  }
  
  if (species[n]=="Nvec") {
    names(ortholog)<-c("gene_id","gene")
    #ortholog$gene<-gsub("\\.","_",ortholog$gene)
    #px$gene<-gsub("Nvec-", "", px$gene)
    #px$gene<-gsub("-", "_", px$gene)
  }
  
  if (species[n]=="Xesp") {
    ortholog<-separate(ortholog, gene2, c("gene"), sep = "-", remove = F)
    ortholog<-ortholog[,c("gene_id","gene")]
    ortholog<-unique(ortholog)
    px$gene<-gsub("Xesp-", "Xe_", px$gene)
  }
  
  if (species[n]=="Chem") {
    ortholog<-unique(ortholog)
    px$gene<-gsub("-", "_", px$gene)
    names(ortholog)<-c("gene_id", "gene")
  }
  
  if (species[n]=="Resc") {
    ortholog<-separate(ortholog, gene2, c("gene"), sep = "-", remove = F)
    ortholog<-ortholog[,c("gene_id","gene")]
    ortholog<-unique(ortholog)
    px$gene<-gsub("-", "_", px$gene)
  }
ortholog<-unique(ortholog)
  
  genecount<-as.data.frame(xtabs(~gene, ortholog))
  genecount<-genecount[genecount$Freq <=3,]
  ortholog<-ortholog[ortholog$gene %in% genecount$gene,]
  
  ortho.dup<-ortholog[duplicated(ortholog$gene)==T,]
  ortholog<-ortholog[-which(ortholog$gene %in% ortho.dup$gene),]
  
  py<-py.bk
  py<-merge(py, ortholog, by = "gene_id", all = F)
  py<-py[,-which(names(py)=="gene_id")]
  row.names(px)<-px$gene
  px<-px[,-which(names(px)=="gene")]
  names(px)<-paste0(species[n],".",names(px))

  row.names(py)<-py$gene
  py<-py[,-which(names(py)=="gene")]
  names(py)<-paste0("AuMv1.",names(py))
  
  #quantile normalization
  pxn<-normalize.quantiles(as.matrix(px))
  dimnames(pxn)<-list(row.names(px), names(px))
  pxn<-as.data.frame(pxn)
  
  pyn<-normalize.quantiles(as.matrix(py))
  dimnames(pyn)<-list(row.names(py), names(py))
  pyn<-as.data.frame(pyn)
  
  #merge 
  pn<-merge(pxn, pyn, by = "row.names", all = FALSE)
  row.names(pn)<-pn$Row.names
  pn<-pn[,-which(names(pn)=="Row.names")]
  
  kld.matrix<-as.data.frame(matrix(data = NA, nrow = ncol(pn), ncol = ncol(pn)))
  colnames(kld.matrix)=colnames(pn)
  rownames(kld.matrix)=colnames(pn)

#Calculation KLD for each pair
  for (i in 1:ncol(pn)) {
    for (j in 1:ncol(pn)) {
      kld=KLD(pn[,i], pn[,j])
      kld.matrix[i,j]=kld$sum.KLD.px.py
    }
  } 
  
  write.table(kld.matrix, 
              paste0("seurat-seurat.",species[n],"-AuMv1.txt"),
              row.names = T, quote = F, sep = "\t")

  kld.matrix.label<-kld.matrix
  #Label 5% top
  for (x in 1:ncol(kld.matrix.label)) {
    co1<-quantile(as.matrix(kld.matrix[1:ncol(pxn),x]), 0.05)
    co2<-quantile(as.matrix(kld.matrix[(ncol(pxn)+1):ncol(kld.matrix), x]), 0.05)
    for (y in 1:ncol(pxn)) {
      if (x==y) {
        kld.matrix.label[y,x]<-""
      } else {
        if (kld.matrix.label[y,x]<=co1) {
          kld.matrix.label[y,x] <- "*"
        } else {
          kld.matrix.label[y,x] <- ""
        }
      }
    }
    
    for (y in (ncol(pxn)+1):ncol(kld.matrix)) {
      if (x==y) {
        kld.matrix.label[y,x]<-""
      } else {
        if (kld.matrix.label[y,x]<=co2) {
          kld.matrix.label[y,x] <- "*"
        } else {
          kld.matrix.label[y,x] <- ""
        }
      }
    }
  }
  
  #Visualisation
  
  if (max(kld.matrix)>=0.1) {
    break.limit<-ceiling(7.5*max(kld.matrix))/10
  } else if (max(kld.matrix)>=0.01){
    break.limit<-ceiling(75*max(kld.matrix))/100
  } else {
    break.limit<-ceiling(750*max(kld.matrix))/1000
  }
  
  heatp<-pheatmap(kld.matrix,
                  color = c(colorRampPalette(colors = c('#E31A1C','#F2F2F0','#5A8FCA'))(100)),
                  legend = TRUE, 
                  breaks = seq(0, break.limit, break.limit/100),
                  scale = "none",
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  cutree_rows = NA, treeheight_row = 0,
                  cutree_cols = NA, treeheight_col = 0,
                  border_color = NA,
                  legend_labels = "Kullback-Leibler divergence",
                  display_numbers = kld.matrix.label,
                  number_color = "black",
                  fontsize_number = 10,
                  cellwidth = 12, cellheight = 12, fontsize = 10,
                  angle_col = c("90"),
                  width = ncol(pn)*0.25+1, height = ncol(pn)*0.25,
                  show_colnames = TRUE,
                  show_rownames = TRUE,
                  gaps_row = c(ncol(pxn)),
                  gaps_col = c(ncol(pxn)),
                  filename = paste0("seurat-seurat.",species[n],"-AuMv1.png")
                  )
}

##ref to Seurat Aceo(direct comparasion)
species<-c("Hvul", 
           "Spis",
           "Nvec",
           "Xesp",
           "Chem"
)

refcelltype<-c("","_broad")

py<-read.delim(paste0("all.marker.wide.", "Aceo.",".txt"), header = TRUE)
py$gene_id<-gsub("-", "_", py$gene_id)
py.bk<-py

for (n in 1:length(species)) {
  for (k in 1:length(refcelltype)){
    
    #ortholog pair
    ortholog<-read.delim(paste0("ortholog.Aurelia__v__", species[n],".tsv.txt"), header = TRUE)
    ortholog<-separate(ortholog, gene1, c("gene","t"), sep = "-", remove = F)
    ortholog<-ortholog[,c("gene","gene2")]
    ortholog$gene<-gsub("-", "_", ortholog$gene)
    ortholog<-unique(ortholog)
    names(ortholog)<-c("gene_id", "gene2")

#px = ref cell type
    px<-read.delim(paste0(species[n], refcelltype[k], "_cell_type_gene_FC"), header = TRUE)
    row.names(px)<-gsub(paste0(species[n],"_"), "", row.names(px))
    row.names(px)<-gsub("-", "_", row.names(px))
    
    if (species[n]=="Hvul") {
      ortholog<-separate(ortholog, gene2, c("sca","gene"), sep = "\\.g", remove = F)
      ortholog<-ortholog[,c("gene_id","gene")]
      ortholog$gene<-paste0("g",ortholog$gene)
      ortholog$gene<-gsub("\\.t","-",ortholog$gene)
      ortholog$gene<-gsub("-", "_", ortholog$gene)
      px$gene<-row.names(px)
    }
    
    if (species[n]=="Spis") {
      names(ortholog)<-c("gene_id","gene")
      ortholog$gene<-gsub("\\.","_",ortholog$gene)
      px$gene<-row.names(px)
      px$gene<-gsub("Spis-", "", px$gene)
      px$gene<-gsub("-", "_", px$gene)
      
    }
    if (species[n]=="Nvec") {
      names(ortholog)<-c("gene_id","gene")
      ortholog$gene<-gsub("\\.","_",ortholog$gene)
      px$gene<-row.names(px)
      px$gene<-gsub("Nvec-", "", px$gene)
      px$gene<-gsub("-", "_", px$gene)
      
    }
    
    if (species[n]=="Xesp") {
      ortholog<-separate(ortholog, gene2, c("gene"), sep = "-", remove = F)
      ortholog<-ortholog[,c("gene_id","gene")]
      ortholog<-unique(ortholog)
      px$gene<-row.names(px)
      px$gene<-paste0("Xe_", px$gene)
    }
    
    if (species[n]=="Chem") {
      ortholog<-unique(ortholog)
      px$gene<-row.names(px)
      px$gene<-gsub("-", "_", px$gene)
      names(ortholog)<-c("gene_id", "gene")
    }
    ortholog<-unique(ortholog)
    
    genecount<-as.data.frame(xtabs(~gene, ortholog))
    genecount<-genecount[genecount$Freq <=3,]
    ortholog<-ortholog[ortholog$gene %in% genecount$gene,]
    
    ortho.dup<-ortholog[duplicated(ortholog$gene)==T,]
    ortholog<-ortholog[-which(ortholog$gene %in% ortho.dup$gene),]
    
    py<-py.bk
    py<-merge(py, ortholog, by = "gene_id", all = F)
    py<-py[,-which(names(py)=="gene_id")]
    row.names(px)<-px$gene
    px<-px[,-which(names(px)=="gene")]
    names(px)<-paste0(species[n],".",names(px))
    row.names(py)<-py$gene
    py<-py[,-which(names(py)=="gene")]
    names(py)<-paste0("Aaur.",names(py))
    
    #quantile normalization
    pxn<-normalize.quantiles(as.matrix(log2(px)))
    dimnames(pxn)<-list(row.names(px), names(px))
    pxn<-as.data.frame(pxn)
    
    pyn<-normalize.quantiles(as.matrix(py))
    dimnames(pyn)<-list(row.names(py), names(py))
    pyn<-as.data.frame(pyn)
    
    #merge 
    pn<-merge(pxn, pyn, by = "row.names", all = FALSE)
    row.names(pn)<-pn$Row.names
    pn<-pn[,-which(names(pn)=="Row.names")]
    
    kld.matrix<-as.data.frame(matrix(data = NA, nrow = ncol(pn), ncol = ncol(pn)))
    colnames(kld.matrix)=colnames(pn)
    rownames(kld.matrix)=colnames(pn)
    
    #Calculation KLD for each pair
    for (i in 1:ncol(pn)) {
      for (j in 1:ncol(pn)) {
        kld=KLD(pn[,i], pn[,j])
        kld.matrix[i,j]=kld$sum.KLD.px.py
      }
    } 
    
    write.table(kld.matrix, 
                paste0("ref-seuratAMv1.",species[n],refcelltype[k],".txt"),
                row.names = T, quote = F, sep = "\t")

    kld.matrix.label<-kld.matrix
    #Label 5% top
    for (x in 1:ncol(kld.matrix.label)) {
      co1<-quantile(as.matrix(kld.matrix[1:ncol(pxn),x]), 0.05)
      co2<-quantile(as.matrix(kld.matrix[(ncol(pxn)+1):ncol(kld.matrix), x]), 0.05)
      for (y in 1:ncol(pxn)) {
        if (x==y) {
          kld.matrix.label[y,x]<-""
        } else {
          if (kld.matrix.label[y,x]<=co1) {
            kld.matrix.label[y,x] <- "*"
          } else {
            kld.matrix.label[y,x] <- ""
          }
        }
      }
      
      for (y in (ncol(pxn)+1):ncol(kld.matrix)) {
        if (x==y) {
          kld.matrix.label[y,x]<-""
        } else {
          if (kld.matrix.label[y,x]<=co2) {
            kld.matrix.label[y,x] <- "*"
          } else {
            kld.matrix.label[y,x] <- ""
          }
        }
      }
    }
    
    #Visualisation
    
    if (max(kld.matrix)>=0.1) {
      break.limit<-ceiling(7.5*max(kld.matrix))/10
    } else if (max(kld.matrix)>=0.01){
      break.limit<-ceiling(75*max(kld.matrix))/100
    } else {
      break.limit<-ceiling(750*max(kld.matrix))/1000
    }
    
    heatp<-pheatmap(kld.matrix,
                    color = c(colorRampPalette(colors = c('#E31A1C','#F2F2F0','#5A8FCA'))(100)),
                    legend = TRUE, 
                    breaks = seq(0, break.limit, break.limit/100),
                    scale = "none",
                    cluster_rows = FALSE, cluster_cols = FALSE,
                    cutree_rows = NA, treeheight_row = 0,
                    cutree_cols = NA, treeheight_col = 0,
                    border_color = NA,
                    legend_labels = "Kullback-Leibler divergence",
                    display_numbers = kld.matrix.label,
                    number_color = "black",
                    fontsize_number = 10,
                    cellwidth = 12, cellheight = 12, fontsize = 10,
                    angle_col = c("90"),
                    width = ncol(pn)*0.25+1, height = ncol(pn)*0.25+1,
                    show_colnames = TRUE,
                    show_rownames = TRUE,
                    gaps_row = c(ncol(pxn)),
                    gaps_col = c(ncol(pxn)),
                    filename = paste0("ref-seuratAuMv1.",species[n],refcelltype[k],".png"))
  }
}
