library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(pheatmap)
library("readxl")
set.seed(1) 

saveext = "~/Desktop/Data/Endometrial/InVitro/Matteo/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))
PrimedNaiveForm <- read_excel("/Users/christopherpenfold/Downloads/mmc3.xls")

PrimedNaive <- read_excel("/Users/christopherpenfold/Downloads/mmc2.xls")
PrimedNaiveForm <- read_excel("/Users/christopherpenfold/Downloads/mmc3.xls")


AmnionMarkers <- read_excel("/Users/christopherpenfold/Downloads/AmnionMarkers.xlsx")

D <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/FinalAno.rds")

uID <- as.character(Idents(D))
uID[which(Idents(D)=="EmDisc" & D$Genotype%in%c("Epithelial_G1","Epithelial_G3"))] <- "Unciliated epithelia"
uID[which(Idents(D)=="EmDisc" & D$Genotype%in%c("Stromal_G4"))] <- "Stromal fibroblasts"
uID[which(Idents(D)=="Am/EmDisc_d14" & D$Genotype%in%c("Epithelial_G3"))] <- "Unciliated epithelia"
uID[which(Idents(D)=="Am/EmDisc_d14" & D$Genotype%in%c("Stromal_G4"))] <- "Stromal fibroblasts"
uID[which(Idents(D)=="Am_d14" & D$Genotype%in%c("Epithelial_G3"))] <- "Unciliated epithelia"
uID[which(Idents(D)=="Am_d14" & D$Genotype%in%c("Stromal_G3","Stromal_G4"))] <- "Stromal fibroblasts"
uID[which(Idents(D)=="ExMes_d14" & D$Genotype%in%c("Epithelial_G3"))] <- "Unciliated epithelia"
uID[which(Idents(D)=="ExMes_d14" & D$Genotype%in%c("Stromal_G3","Stromal_G4"))] <- "Stromal fibroblasts"
uID[which(Idents(D)=="Hyp_d14" & D$Genotype%in%c("Stromal_G3","Stromal_G4"))] <- "Stromal fibroblasts"
uID[which(Idents(D)=="Hyp_d14" & D$Genotype%in%c("Stromal_G3","Stromal_G4"))] <- "Stromal fibroblasts"
uID[which(Idents(D)=="ExMes_d14" & D$Genotype%in%c("Epithelial_G3"))] <- "Unciliated epithelia"
uID[which(Idents(D)=="ExMes_d14" & D$Genotype%in%c("Stromal_G3","Stromal_G4"))] <- "Stromal fibroblasts"
uID[which(Idents(D)=="VE_d14" & D$Genotype%in%c("Stromal_G3","Stromal_G4"))] <- "Stromal fibroblasts"
Idents(D) <- uID

uID<- as.character(Idents(D))
uID[which(uID%in%c("Unciliated epithelia","Lumenal"))] <- "Lumenal"
uID[which(uID%in%c("VE_d14"))] <- "Hyp_d14"

Idents(D) <- uID
Idents(D) <- as.factor(paste(D$ID3,Idents(D),sep="_"))


#uID<- as.character(Idents(mammal.combined))
#uID[which(uID%in%c("Ciliated","Prolif","Unciliated epithelia","SOX9P","Lumenal","SOX9LRG5","Glandular","SOX9LRG5"))] <- "Epithelia"
#Idents(mammal.combined) <- uID\



exps <- c("C1_CTB_d14", 
          "C2_CTB_d14",
          "C6_CTB_d14",
          "C8_CTB_d14",
          "C1_STB_d14",
          "C2_STB_d14",
          "C3_STB_d14",
          "C4_STB_d14",
          "C6_STB_d14",
          "C7_STB_d14",
          "C8_STB_d14",
          "C2_EVT_d14",
          "C3_EVT_d14",
          "C4_EVT_d14",
          "C6_EVT_d14",
          "C7_EVT_d14",
          "C8_EVT_d14",
          "C1_ExMes_d14",
          "C2_ExMes_d14",
          "C3_ExMes_d14",
          "C6_ExMes_d14",
          "C7_ExMes_d14",
          "C8_ExMes_d14",
          "C1_Hyp_d14",
          "C7_Hyp_d14",
          "C2_Hyp_d14",
          "C1_Am_d14",
          "C2_Am_d14" ,
          "C6_Am_d14",
          "C7_Am_d14",
          "C1_Am/EmDisc_d14",
          "C2_Am/EmDisc_d14",
          "C4_Am/EmDisc_d14",
          "C6_Am/EmDisc_d14",
          "C7_Am/EmDisc_d14",
          "C1_EmDisc_d14",
          "C2_EmDisc_d14",
          "C4_EmDisc_d14",
          "C6_EmDisc_d14",
          "C7_EmDisc_d14",
          "C8_EmDisc_d14",
"C1_Lumenal",
"C2_Lumenal",
"C3_Lumenal",
"C4_Lumenal",
"C6_Lumenal",
"C7_Lumenal",
"C8_Lumenal",
"C1_Prolif",
"C2_Prolif",
"C3_Prolif",
"C4_Prolif",
"C6_Prolif",
"C7_Prolif",
"C8_Prolif",
"C1_SOX9P",
"C2_SOX9P",
"C3_SOX9P",
"C4_SOX9P",
"C6_SOX9P",
"C7_SOX9P",
"C8_SOX9P",
"C1_SOX9LRG5",
"C2_SOX9LRG5",
"C4_SOX9LRG5",
"C6_SOX9LRG5",
"C8_SOX9LRG5",
"C1_Glandular",
"C2_Glandular",
"C3_Glandular",
"C4_Glandular",
"C6_Glandular",
"C7_Glandular",
"C8_Glandular",
"C1_Ciliated",
"C2_Ciliated",
"C3_Ciliated",
"C4_Ciliated",
"C6_Ciliated",
"C7_Ciliated",
"C8_Ciliated",
"C1_Stromal fibroblasts",
"C2_Stromal fibroblasts",
"C3_Stromal fibroblasts",
"C4_Stromal fibroblasts",
"C6_Stromal fibroblasts",
"C7_Stromal fibroblasts",
"C8_Stromal fibroblasts")


DefaultAssay(D) <- "RNA"

Av <- AverageExpression(D)
Av <- Av$RNA


Markers1 <- c("JAM3",
              "GATA2",
              "GATA3",
              "TFAP2A",
              "TFAP2C",
              "CGB8",
              "CGB5",
              "CGA",
              "ERVW-1",
              "PRDM6",
              "TBX3",
              "DIO2",
              "NOTUM",
              "HLA-G",
              "ASCL2",
              "SNAI1",
              "HGF",
              "SNAI2",
              "HAND1",
              "HAND2",
              "PDGFRA",
              "BST2",
              "GATA6",
              "GATA4",
              "CER1",
              "NODAL",
              "LEFTY1",
              "LEFTY2",
              "APOA1",
              "WNT6",
              "VTCN1",
              "BAMBI",
              "AKAP12",
              "PODXL",
              "SFRP1",
              "SOX15",
              "PDGFA",
              "POU5F1")

Markers2 <- c("MKI67","MECOM","SPON1","MSLN","CLDN3","SOX9","LGR5","ECM1","FOXO1","PTGS1","IL6","FOXJ1", "PIFO", "TP73")

#Markers2 <- c("SOX9","PGR","ESR1","MMP7",
#"CPM","MKI67","HMGB2",
##"PLAU","IL32","TNF","WNT7A","KRT7","KRT5","LGR5","IHH","EMID1","PPARG","MUC12",
#"CDC20B","CCNO","HES6","FOXJ1","PIFO","TP73","HEY1","ABCG1","SCGB2A2","C2CD4A","SLC18A2","PAEP","CXCL14","SPP1","DPP4","PAX2","VTCN1","SLC26A7","MSLN")


#Markers3 <- c("ACTA2","C7","IGF1","PCOLCE","MMP11","ECM1","FOXO1","IL15","CFD","CEBPB","PDGFA")
Markers3 <- c("PCOLCE","MMP11","DPP6","NETO1")


Markers <- c(Markers1,Markers2,Markers3)

Dat <- Av[Markers,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_selectmarkers.pdf",sep=""),scale="row",width=20,height=20)

Markers2 <- c("SOX9","PGR","ESR1","MMP7","CPM","MKI67","HMGB2","PLAU","IL32","TNF","WNT7A","KRT7","KRT5","LGR5","IHH","EMID1","PPARG","MUC12","CDC20B","CCNO","HES6","FOXJ1","PIFO","TP73","HEY1","ABCG1","SCGB2A2","C2CD4A","SLC18A2","PAEP","CXCL14","SPP1","DPP4","PAX2","VTCN1","SLC26A7","MSLN")
Markers3 <- c("ACTA2","C7","IGF1","PCOLCE","MMP11","ECM1","FOXO1","IL15","CFD","CEBPB","PDGFA")

Markers <- c(Markers1,Markers2,Markers3)
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88),cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_selectmarkers_2.pdf",sep=""),scale="row",width=20,height=20)


List1 <- c("C1_CTB_d14", 
                "C1_STB_d14",
                "C1_ExMes_d14",
                "C1_Hyp_d14",
                "C1_Am/EmDisc_d14",
                "C1_EmDisc_d14",
                "C1_Lumenal",
                "C1_Prolif",
                "C1_SOX9P",
                "C1_SOX9LRG5",
                "C1_Glandular",
                "C1_Ciliated",
                "C1_Stromal fibroblasts")

List2 <- c(
          "C2_CTB_d14",
          "C2_STB_d14",
          "C2_EVT_d14",
          "C2_ExMes_d14",
          "C2_Hyp_d14",
          "C2_Am/EmDisc_d14",
          "C2_EmDisc_d14",
          "C2_Lumenal",
          "C2_Prolif",
          "C2_SOX9P",
          "C2_SOX9LRG5",
          "C2_Glandular",
          "C2_Ciliated",
          "C2_Stromal fibroblasts")

List3 <- c("C3_STB_d14",
          "C3_EVT_d14",
          "C3_ExMes_d14",
          "C3_Lumenal",
          "C3_Prolif",
          "C3_SOX9P",
          "C3_Glandular",
          "C3_Ciliated",
          "C3_Stromal fibroblasts")


List4 <- c("C4_STB_d14",
          "C4_EVT_d14",
          "C4_Am/EmDisc_d14",
          "C4_EmDisc_d14",
          "C4_Lumenal",
          "C4_Prolif",
          "C4_SOX9P",
          "C4_SOX9LRG5",
          "C4_Glandular",
          "C4_Ciliated",
          "C4_Stromal fibroblasts")


List6 <- c("C6_CTB_d14",
          "C6_STB_d14",
          "C6_EVT_d14",
          "C6_ExMes_d14",
          "C6_Am_d14",
          "C6_Am/EmDisc_d14",
          "C6_EmDisc_d14",
          "C6_Lumenal",
          "C6_Prolif",
          "C6_SOX9P",
          "C6_SOX9LRG5",
          "C6_Glandular",
          "C6_Ciliated",
          "C6_Stromal fibroblasts")




List7 <- c("C7_STB_d14",
          "C7_EVT_d14",
          "C7_ExMes_d14",
          "C7_Hyp_d14",
          "C7_Am_d14",
          "C7_Am/EmDisc_d14",
          "C7_EmDisc_d14",
          "C7_Lumenal",
          "C7_Prolif",
          "C7_SOX9P",
          "C7_Glandular",
          "C7_Ciliated",
          "C7_Stromal fibroblasts")



List8 <- c("C8_CTB_d14",
          "C8_STB_d14",
          "C8_EVT_d14",
          "C8_ExMes_d14",
          "C8_EmDisc_d14",
          "C8_Lumenal",
          "C8_Prolif",
          "C8_SOX9P",
          "C8_SOX9LRG5",
          "C8_Glandular",
          "C8_Ciliated",
          "C8_Stromal fibroblasts")


#First do stromal
DEListA1 <- FindMarkers(D, ident.1 = "C1_Stromal fibroblasts", ident.2 = setdiff(List1, "C1_Stromal fibroblasts"),test.use = "MAST", only.pos = TRUE)
DEListA2 <- FindMarkers(D, ident.1 = "C2_Stromal fibroblasts", ident.2 = setdiff(List2, "C2_Stromal fibroblasts"),test.use = "MAST", only.pos = TRUE)
DEListA3 <- FindMarkers(D, ident.1 = "C3_Stromal fibroblasts", ident.2 = setdiff(List3, "C3_Stromal fibroblasts"),test.use = "MAST", only.pos = TRUE)
DEListA4 <- FindMarkers(D, ident.1 = "C4_Stromal fibroblasts", ident.2 = setdiff(List4, "C4_Stromal fibroblasts"),test.use = "MAST", only.pos = TRUE)
DEListA6 <- FindMarkers(D, ident.1 = "C6_Stromal fibroblasts", ident.2 = setdiff(List6, "C6_Stromal fibroblasts"),test.use = "MAST", only.pos = TRUE)
DEListA7 <- FindMarkers(D, ident.1 = "C7_Stromal fibroblasts", ident.2 = setdiff(List7, "C7_Stromal fibroblasts"),test.use = "MAST", only.pos = TRUE)
DEListA8 <- FindMarkers(D, ident.1 = "C8_Stromal fibroblasts", ident.2 = setdiff(List8, "C8_Stromal fibroblasts"),test.use = "MAST", only.pos = TRUE)

DEListA1 <- DEListA1[which(DEListA1$p_val_adj<0.001),]
DEListA1 <- DEListA1[order(DEListA1$avg_log2FC,decreasing=TRUE),]
DEListA2 <- DEListA2[which(DEListA2$p_val_adj<0.001),]
DEListA2 <- DEListA2[order(DEListA2$avg_log2FC,decreasing=TRUE),]
DEListA3 <- DEListA3[which(DEListA3$p_val_adj<0.001),]
DEListA3 <- DEListA3[order(DEListA3$avg_log2FC,decreasing=TRUE),]
DEListA4 <- DEListA4[which(DEListA4$p_val_adj<0.001),]
DEListA4 <- DEListA4[order(DEListA4$avg_log2FC,decreasing=TRUE),]
DEListA6 <- DEListA6[which(DEListA6$p_val_adj<0.001),]
DEListA6 <- DEListA6[order(DEListA6$avg_log2FC,decreasing=TRUE),]
DEListA7 <- DEListA7[which(DEListA7$p_val_adj<0.001),]
DEListA7 <- DEListA7[order(DEListA7$avg_log2FC,decreasing=TRUE),]
DEListt8 <- DEListA8[which(DEListA8$p_val_adj<0.001),]
DEListA8 <- DEListA8[order(DEListA8$avg_log2FC,decreasing=TRUE),]

L1 <- Reduce(intersect, list(rownames(DEListA1),rownames(DEListA2),rownames(DEListA3),rownames(DEListA4),rownames(DEListA6),rownames(DEListA7),rownames(DEListA8)))

L1 <- unique( c(
  intersect(rownames(DEListA1),rownames(DEListA2)),
  intersect(rownames(DEListA1),rownames(DEListA3)),
  intersect(rownames(DEListA1),rownames(DEListA4)),
  intersect(rownames(DEListA1),rownames(DEListA6)),
  intersect(rownames(DEListA1),rownames(DEListA7)),
  intersect(rownames(DEListA1),rownames(DEListA8)),
  intersect(rownames(DEListA2),rownames(DEListA3)),
  intersect(rownames(DEListA2),rownames(DEListA4)),
  intersect(rownames(DEListA2),rownames(DEListA6)),
  intersect(rownames(DEListA2),rownames(DEListA7)),
  intersect(rownames(DEListA2),rownames(DEListA8)),
  intersect(rownames(DEListA3),rownames(DEListA4)),
  intersect(rownames(DEListA3),rownames(DEListA6)),
  intersect(rownames(DEListA3),rownames(DEListA7)),
  intersect(rownames(DEListA3),rownames(DEListA8)),
  intersect(rownames(DEListA4),rownames(DEListA6)),
  intersect(rownames(DEListA4),rownames(DEListA7)),
  intersect(rownames(DEListA4),rownames(DEListA8)),
  intersect(rownames(DEListA6),rownames(DEListA7)),
  intersect(rownames(DEListA6),rownames(DEListA8)),
  intersect(rownames(DEListA7),rownames(DEListA8))
) )
#205

#First do lumenal
DEListB1 <- FindMarkers(D, ident.1 = c("C1_Lumenal"), ident.2 = setdiff(List1, "C1_Lumenal"),test.use = "MAST", only.pos = TRUE)
DEListB2 <- FindMarkers(D, ident.1 = c("C2_Lumenal"), ident.2 = setdiff(List2, "C2_Lumenal"),test.use = "MAST", only.pos = TRUE)
DEListB3 <- FindMarkers(D, ident.1 = c("C3_Lumenal"), ident.2 = setdiff(List3, "C3_Lumenal"),test.use = "MAST", only.pos = TRUE)
DEListB4 <- FindMarkers(D, ident.1 = c("C4_Lumenal"), ident.2 = setdiff(List4, "C4_Lumenal"),test.use = "MAST", only.pos = TRUE)
DEListB6 <- FindMarkers(D, ident.1 = c("C6_Lumenal"), ident.2 = setdiff(List6, "C6_Lumenal"),test.use = "MAST", only.pos = TRUE)
DEListB7 <- FindMarkers(D, ident.1 = c("C7_Lumenal"), ident.2 = setdiff(List7, "C7_Lumenal"),test.use = "MAST", only.pos = TRUE)
DEListB8 <- FindMarkers(D, ident.1 = c("C8_Lumenal"), ident.2 = setdiff(List8, "C8_Lumenal"),test.use = "MAST", only.pos = TRUE)

DEListB1 <- DEListB1[which(DEListB1$p_val_adj<0.001),]
DEListB1 <- DEListB1[order(DEListB1$avg_log2FC,decreasing=TRUE),]
DEListB2 <- DEListB2[which(DEListB2$p_val_adj<0.001),]
DEListB2 <- DEListB2[order(DEListB2$avg_log2FC,decreasing=TRUE),]
DEListB3 <- DEListB3[which(DEListB3$p_val_adj<0.001),]
DEListB3 <- DEListB3[order(DEListB3$avg_log2FC,decreasing=TRUE),]
DEListB4 <- DEListB4[which(DEListB4$p_val_adj<0.001),]
DEListB4 <- DEListB4[order(DEListB4$avg_log2FC,decreasing=TRUE),]
DEListB6 <- DEListB6[which(DEListB6$p_val_adj<0.001),]
DEListB6 <- DEListB6[order(DEListB6$avg_log2FC,decreasing=TRUE),]
DEListB7 <- DEListB7[which(DEListB7$p_val_adj<0.001),]
DEListB7 <- DEListB7[order(DEListB7$avg_log2FC,decreasing=TRUE),]
DEListB8 <- DEListB8[which(DEListB8$p_val_adj<0.001),]
DEListB8 <- DEListB8[order(DEListB8$avg_log2FC,decreasing=TRUE),]

L2 <- Reduce(intersect, list(rownames(DEListB1),rownames(DEListB2),rownames(DEListB3),rownames(DEListB4),rownames(DEListB6),rownames(DEListB7),rownames(DEListB8)))

L2 <- unique( c(
  intersect(rownames(DEListB1),rownames(DEListB2)),
  intersect(rownames(DEListB1),rownames(DEListB3)),
  intersect(rownames(DEListB1),rownames(DEListB4)),
  intersect(rownames(DEListB1),rownames(DEListB6)),
  intersect(rownames(DEListB1),rownames(DEListB7)),
  intersect(rownames(DEListB1),rownames(DEListB8)),
  intersect(rownames(DEListB2),rownames(DEListB3)),
  intersect(rownames(DEListB2),rownames(DEListB4)),
  intersect(rownames(DEListB2),rownames(DEListB6)),
  intersect(rownames(DEListB2),rownames(DEListB7)),
  intersect(rownames(DEListB2),rownames(DEListB8)),
  intersect(rownames(DEListB3),rownames(DEListB4)),
  intersect(rownames(DEListB3),rownames(DEListB6)),
  intersect(rownames(DEListB3),rownames(DEListB7)),
  intersect(rownames(DEListB3),rownames(DEListB8)),
  intersect(rownames(DEListB4),rownames(DEListB6)),
  intersect(rownames(DEListB4),rownames(DEListB7)),
  intersect(rownames(DEListB4),rownames(DEListB8)),
  intersect(rownames(DEListB6),rownames(DEListB7)),
  intersect(rownames(DEListB6),rownames(DEListB8)),
  intersect(rownames(DEListB7),rownames(DEListB8))
) )
#only 23

DEListC1 <- FindMarkers(D, ident.1 = c("C1_Prolif"), ident.2 = setdiff(List1, "C1_Prolif"),test.use = "MAST", only.pos = TRUE)
DEListC2 <- FindMarkers(D, ident.1 = c("C2_Prolif"), ident.2 = setdiff(List2, "C2_Prolif"),test.use = "MAST", only.pos = TRUE)
DEListC3 <- FindMarkers(D, ident.1 = c("C3_Prolif"), ident.2 = setdiff(List3, "C3_Prolif"),test.use = "MAST", only.pos = TRUE)
DEListC4 <- FindMarkers(D, ident.1 = c("C4_Prolif"), ident.2 = setdiff(List4, "C4_Prolif"),test.use = "MAST", only.pos = TRUE)
DEListC6 <- FindMarkers(D, ident.1 = c("C6_Prolif"), ident.2 = setdiff(List6, "C6_Prolif"),test.use = "MAST", only.pos = TRUE)
DEListC7 <- FindMarkers(D, ident.1 = c("C7_Prolif"), ident.2 = setdiff(List7, "C7_Prolif"),test.use = "MAST", only.pos = TRUE)
DEListC8 <- FindMarkers(D, ident.1 = c("C8_Prolif"), ident.2 = setdiff(List8, "C8_Prolif"),test.use = "MAST", only.pos = TRUE)

DEListC1 <- DEListC1[which(DEListC1$p_val_adj<0.001),]
DEListC1 <- DEListC1[order(DEListC1$avg_log2FC,decreasing=TRUE),]
DEListC2 <- DEListC2[which(DEListC2$p_val_adj<0.001),]
DEListC2 <- DEListC2[order(DEListC2$avg_log2FC,decreasing=TRUE),]
DEListC3 <- DEListC3[which(DEListC3$p_val_adj<0.001),]
DEListC3 <- DEListC3[order(DEListC3$avg_log2FC,decreasing=TRUE),]
DEListC4 <- DEListC4[which(DEListC4$p_val_adj<0.001),]
DEListC4 <- DEListC4[order(DEListC4$avg_log2FC,decreasing=TRUE),]
DEListC6 <- DEListC6[which(DEListC6$p_val_adj<0.001),]
DEListC6 <- DEListC6[order(DEListC6$avg_log2FC,decreasing=TRUE),]
DEListC7 <- DEListC7[which(DEListC7$p_val_adj<0.001),]
DEListC7 <- DEListC7[order(DEListC7$avg_log2FC,decreasing=TRUE),]
DEListC8 <- DEListC8[which(DEListC8$p_val_adj<0.001),]
DEListC8 <- DEListC8[order(DEListC8$avg_log2FC,decreasing=TRUE),]

L3 <- Reduce(intersect, list(rownames(DEListC1),rownames(DEListC2),rownames(DEListC3),rownames(DEListC4),rownames(DEListC6),rownames(DEListC7),rownames(DEListC8)))
L3 <- unique( c(
  intersect(rownames(DEListC1),rownames(DEListC2)),
  intersect(rownames(DEListC1),rownames(DEListC3)),
  intersect(rownames(DEListC1),rownames(DEListC4)),
  intersect(rownames(DEListC1),rownames(DEListC6)),
  intersect(rownames(DEListC1),rownames(DEListC7)),
  intersect(rownames(DEListC1),rownames(DEListC8)),
  intersect(rownames(DEListC2),rownames(DEListC3)),
  intersect(rownames(DEListC2),rownames(DEListC4)),
  intersect(rownames(DEListC2),rownames(DEListC6)),
  intersect(rownames(DEListC2),rownames(DEListC7)),
  intersect(rownames(DEListC2),rownames(DEListC8)),
  intersect(rownames(DEListC3),rownames(DEListC4)),
  intersect(rownames(DEListC3),rownames(DEListC6)),
  intersect(rownames(DEListC3),rownames(DEListC7)),
  intersect(rownames(DEListC3),rownames(DEListC8)),
  intersect(rownames(DEListC4),rownames(DEListC6)),
  intersect(rownames(DEListC4),rownames(DEListC7)),
  intersect(rownames(DEListC4),rownames(DEListC8)),
  intersect(rownames(DEListC6),rownames(DEListC7)),
  intersect(rownames(DEListC6),rownames(DEListC8)),
  intersect(rownames(DEListC7),rownames(DEListC8))
) )
#Only 41


#First do Glandular
DEListD3 <- FindMarkers(D, ident.1 = c("C3_Glandular"), ident.2 = setdiff(List3, "C3_Glandular"),test.use = "MAST", only.pos = TRUE)
DEListD4 <- FindMarkers(D, ident.1 = c("C4_Glandular"), ident.2 = setdiff(List4, "C4_Glandular"),test.use = "MAST", only.pos = TRUE)
DEListD6 <- FindMarkers(D, ident.1 = c("C6_Glandular"), ident.2 = setdiff(List6, "C6_Glandular"),test.use = "MAST", only.pos = TRUE)
DEListD7 <- FindMarkers(D, ident.1 = c("C7_Glandular"), ident.2 = setdiff(List7, "C7_Glandular"),test.use = "MAST", only.pos = TRUE)
DEListD8 <- FindMarkers(D, ident.1 = c("C8_Glandular"), ident.2 = setdiff(List8, "C8_Glandular"),test.use = "MAST", only.pos = TRUE)

DEListD3 <- DEListD3[which(DEListD3$p_val_adj<0.01),]
DEListD3 <- DEListD3[order(DEListD3$avg_log2FC,decreasing=TRUE),]
DEListD4 <- DEListD4[which(DEListD4$p_val_adj<0.01),]
DEListD4 <- DEListD4[order(DEListD4$avg_log2FC,decreasing=TRUE),]
DEListD6 <- DEListD6[which(DEListD6$p_val_adj<0.01),]
DEListD6 <- DEListD6[order(DEListD6$avg_log2FC,decreasing=TRUE),]
DEListD7 <- DEListD7[which(DEListD7$p_val_adj<0.01),]
DEListD7 <- DEListD7[order(DEListD7$avg_log2FC,decreasing=TRUE),]
DEListD8 <- DEListD8[which(DEListD8$p_val_adj<0.01),]
DEListD8 <- DEListD8[order(DEListD8$avg_log2FC,decreasing=TRUE),]

L4 <- Reduce(intersect, list(rownames(DEListD3),rownames(DEListD4),rownames(DEListD6),rownames(DEListD7),rownames(DEListD8)))
L4 <- unique( c(
  intersect(rownames(DEListD3),rownames(DEListD4)),
  intersect(rownames(DEListD3),rownames(DEListD6)),
  intersect(rownames(DEListD3),rownames(DEListD7)),
  intersect(rownames(DEListD3),rownames(DEListD8)),
  intersect(rownames(DEListD4),rownames(DEListD6)),
  intersect(rownames(DEListD4),rownames(DEListD7)),
  intersect(rownames(DEListD4),rownames(DEListD8)),
  intersect(rownames(DEListD6),rownames(DEListD7)),
  intersect(rownames(DEListD6),rownames(DEListD8)),
  intersect(rownames(DEListD7),rownames(DEListD8))
) )
#No overlap
L4 <- rownames(DEListD4)

#First do Ciliated
DEListE1 <- FindMarkers(D, ident.1 = c("C1_Ciliated"), ident.2 = setdiff(List1, "C1_Ciliated"),test.use = "MAST", only.pos = TRUE)
DEListE2 <- FindMarkers(D, ident.1 = c("C2_Ciliated"), ident.2 = setdiff(List2, "C2_Ciliated"),test.use = "MAST", only.pos = TRUE)
DEListE3 <- FindMarkers(D, ident.1 = c("C3_Ciliated"), ident.2 = setdiff(List3, "C3_Ciliated"),test.use = "MAST", only.pos = TRUE)
DEListE4 <- FindMarkers(D, ident.1 = c("C4_Ciliated"), ident.2 = setdiff(List4, "C4_Ciliated"),test.use = "MAST", only.pos = TRUE)
DEListE6 <- FindMarkers(D, ident.1 = c("C6_Ciliated"), ident.2 = setdiff(List6, "C6_Ciliated"),test.use = "MAST", only.pos = TRUE)
#DEListE7 <- FindMarkers(D, ident.1 = c("C7_Ciliated"), ident.2 = setdiff(List7, "C7_Ciliated"),test.use = "MAST", only.pos = TRUE)
DEListE8 <- FindMarkers(D, ident.1 = c("C8_Ciliated"), ident.2 = setdiff(List8, "C8_Ciliated"),test.use = "MAST", only.pos = TRUE)

DEListE1 <- DEListE1[which(DEListE1$p_val_adj<0.001),]
DEListE1 <- DEListE1[order(DEListE1$avg_log2FC,decreasing=TRUE),]
DEListE2 <- DEListE2[which(DEListE2$p_val_adj<0.001),]
DEListE2 <- DEListE2[order(DEListE2$avg_log2FC,decreasing=TRUE),]
DEListE3 <- DEListE3[which(DEListE3$p_val_adj<0.001),]
DEListE3 <- DEListE3[order(DEListE3$avg_log2FC,decreasing=TRUE),]
DEListE4 <- DEListE4[which(DEListE4$p_val_adj<0.001),]
DEListE4 <- DEListE4[order(DEListE4$avg_log2FC,decreasing=TRUE),]
DEListE6 <- DEListE6[which(DEListE6$p_val_adj<0.001),]
DEListE6 <- DEListE6[order(DEListE6$avg_log2FC,decreasing=TRUE),]
DEListE8 <- DEListE8[which(DEListE8$p_val_adj<0.001),]
DEListE8 <- DEListE8[order(DEListE8$avg_log2FC,decreasing=TRUE),]
L5 <- Reduce(intersect, list(rownames(DEListE1),rownames(DEListE2),rownames(DEListE3),rownames(DEListE4),rownames(DEListE6),rownames(DEListE8)))
L5 <- unique( c(
  intersect(rownames(DEListE1),rownames(DEListE2)),
  intersect(rownames(DEListE1),rownames(DEListE3)),
  intersect(rownames(DEListE1),rownames(DEListE4)),
  intersect(rownames(DEListE1),rownames(DEListE6)),
  intersect(rownames(DEListE1),rownames(DEListE8)),
  intersect(rownames(DEListE2),rownames(DEListE3)),
  intersect(rownames(DEListE2),rownames(DEListE4)),
  intersect(rownames(DEListE2),rownames(DEListE6)),
  intersect(rownames(DEListE2),rownames(DEListE8)),
  intersect(rownames(DEListE3),rownames(DEListE4)),
  intersect(rownames(DEListE3),rownames(DEListE6)),
  intersect(rownames(DEListE3),rownames(DEListE8)),
  intersect(rownames(DEListE4),rownames(DEListE6)),
  intersect(rownames(DEListE4),rownames(DEListE8)),
  intersect(rownames(DEListE6),rownames(DEListE8))) )
#MANY



#First do SOX9P
DEListG1 <- FindMarkers(D, ident.1 = c("C1_SOX9P"), ident.2 = setdiff(List1, "C1_SOX9P"),test.use = "MAST", only.pos = TRUE)
DEListG2 <- FindMarkers(D, ident.1 = c("C2_SOX9P"), ident.2 = setdiff(List2, "C2_SOX9P"),test.use = "MAST", only.pos = TRUE)
DEListG3 <- FindMarkers(D, ident.1 = c("C3_SOX9P"), ident.2 = setdiff(List3, "C3_SOX9P"),test.use = "MAST", only.pos = TRUE)
DEListG4 <- FindMarkers(D, ident.1 = c("C4_SOX9P"), ident.2 = setdiff(List4, "C4_SOX9P"),test.use = "MAST", only.pos = TRUE)
DEListG6 <- FindMarkers(D, ident.1 = c("C6_SOX9P"), ident.2 = setdiff(List6, "C6_SOX9P"),test.use = "MAST", only.pos = TRUE)
DEListG7 <- FindMarkers(D, ident.1 = c("C7_SOX9P"), ident.2 = setdiff(List7, "C7_SOX9P"),test.use = "MAST", only.pos = TRUE)
DEListG8 <- FindMarkers(D, ident.1 = c("C8_SOX9P"), ident.2 = setdiff(List8, "C8_SOX9P"),test.use = "MAST", only.pos = TRUE)

DEListG1 <- DEListG1[which(DEListG1$p_val_adj<0.001),]
DEListG1 <- DEListG1[order(DEListG1$avg_log2FC,decreasing=TRUE),]
DEListG2 <- DEListG2[which(DEListG2$p_val_adj<0.001),]
DEListG2 <- DEListG2[order(DEListG2$avg_log2FC,decreasing=TRUE),]
DEListG3 <- DEListG3[which(DEListG3$p_val_adj<0.001),]
DEListG3 <- DEListG3[order(DEListG3$avg_log2FC,decreasing=TRUE),]
DEListG4 <- DEListG4[which(DEListG4$p_val_adj<0.001),]
DEListG4 <- DEListG4[order(DEListG4$avg_log2FC,decreasing=TRUE),]
DEListG6 <- DEListG6[which(DEListG6$p_val_adj<0.001),]
DEListG6 <- DEListG6[order(DEListG6$avg_log2FC,decreasing=TRUE),]
DEListG7 <- DEListG7[which(DEListG7$p_val_adj<0.001),]
DEListG7 <- DEListG7[order(DEListG7$avg_log2FC,decreasing=TRUE),]
DEListG8 <- DEListG8[which(DEListG8$p_val_adj<0.001),]
DEListG8 <- DEListG8[order(DEListG8$avg_log2FC,decreasing=TRUE),]
L7 <- Reduce(intersect, list(rownames(DEListG1),rownames(DEListG2),rownames(DEListG3),rownames(DEListG4),rownames(DEListG6),rownames(DEListG7),rownames(DEListG8)))
L7 <- unique( c(
  intersect(rownames(DEListG1),rownames(DEListG2)),
  intersect(rownames(DEListG1),rownames(DEListG3)),
  intersect(rownames(DEListG1),rownames(DEListG4)),
  intersect(rownames(DEListG1),rownames(DEListG6)),
  intersect(rownames(DEListG1),rownames(DEListG7)),
  intersect(rownames(DEListG1),rownames(DEListG8)),
  intersect(rownames(DEListG2),rownames(DEListG3)),
  intersect(rownames(DEListG2),rownames(DEListG4)),
  intersect(rownames(DEListG2),rownames(DEListG6)),
  intersect(rownames(DEListG2),rownames(DEListG7)),
  intersect(rownames(DEListG2),rownames(DEListG8)),
  intersect(rownames(DEListG3),rownames(DEListG4)),
  intersect(rownames(DEListG3),rownames(DEListG6)),
  intersect(rownames(DEListG3),rownames(DEListG7)),
  intersect(rownames(DEListG3),rownames(DEListG8)),
  intersect(rownames(DEListG4),rownames(DEListG6)),
  intersect(rownames(DEListG4),rownames(DEListG7)),
  intersect(rownames(DEListG4),rownames(DEListG8)),
  intersect(rownames(DEListG6),rownames(DEListG7)),
  intersect(rownames(DEListG6),rownames(DEListG8)),
  intersect(rownames(DEListG7),rownames(DEListG8))
) )
#one 

#First do SOX9LRG5
DEListH1 <- FindMarkers(D, ident.1 = c("C1_SOX9LRG5"), ident.2 = setdiff(List1, "C1_SOX9LRG5"),test.use = "MAST", only.pos = TRUE)
DEListH2 <- FindMarkers(D, ident.1 = c("C2_SOX9LRG5"), ident.2 = setdiff(List2, "C2_SOX9LRG5"),test.use = "MAST", only.pos = TRUE)
#DEListH3 <- FindMarkers(D, ident.1 = c("C3_SOX9LRG5"), ident.2 = setdiff(List3, "C3_SOX9LRG5"),test.use = "MAST", only.pos = TRUE)
#DEListH4 <- FindMarkers(D, ident.1 = c("C4_SOX9LRG5"), ident.2 = setdiff(List4, "C4_SOX9LRG5"),test.use = "MAST", only.pos = TRUE)
DEListH6 <- FindMarkers(D, ident.1 = c("C6_SOX9LRG5"), ident.2 = setdiff(List6, "C6_SOX9LRG5"),test.use = "MAST", only.pos = TRUE)
#DEListH7 <- FindMarkers(D, ident.1 = c("C7_SOX9LRG5"), ident.2 = setdiff(List7, "C7_SOX9LRG5"),test.use = "MAST", only.pos = TRUE)
DEListH8 <- FindMarkers(D, ident.1 = c("C8_SOX9LRG5"), ident.2 = setdiff(List8, "C8_SOX9LRG5"),test.use = "MAST", only.pos = TRUE)

DEListH1 <- DEListH1[which(DEListH1$p_val_adj<0.001),]
DEListH1 <- DEListH1[order(DEListH1$avg_log2FC,decreasing=TRUE),]
DEListH2 <- DEListH2[which(DEListH2$p_val_adj<0.001),]
DEListH2 <- DEListH2[order(DEListH2$avg_log2FC,decreasing=TRUE),]
DEListH6 <- DEListH6[which(DEListH6$p_val_adj<0.001),]
DEListH6 <- DEListH6[order(DEListH6$avg_log2FC,decreasing=TRUE),]
#DEListH7 <- DEListH7[which(DEListH7$p_val_adj<0.001),]
#DEListH7 <- DEListH7[order(DEListH7$avg_log2FC,decreasing=TRUE),]
DEListH8 <- DEListH8[which(DEListH8$p_val_adj<0.001),]
DEListH8 <- DEListH8[order(DEListH8$avg_log2FC,decreasing=TRUE),]

L8 <- Reduce(intersect, list(rownames(DEListH1),rownames(DEListH2),rownames(DEListH6),rownames(DEListH8)))
L8 <- unique( c(
  intersect(rownames(DEListH1),rownames(DEListH2)),
  intersect(rownames(DEListH1),rownames(DEListH6)),
  intersect(rownames(DEListH1),rownames(DEListH8)),
  intersect(rownames(DEListH2),rownames(DEListH6)),
  intersect(rownames(DEListH2),rownames(DEListH8)),
  intersect(rownames(DEListH6),rownames(DEListH8))
) )
#Nothing

#First do STB_d14
DEListI1 <- FindMarkers(D, ident.1 = c("C1_STB_d14"), ident.2 = setdiff(List1, "C1_STB_d14"),test.use = "MAST", only.pos = TRUE)
DEListI2 <- FindMarkers(D, ident.1 = c("C2_STB_d14"), ident.2 = setdiff(List2, "C2_STB_d14"),test.use = "MAST", only.pos = TRUE)
DEListI3 <- FindMarkers(D, ident.1 = c("C3_STB_d14"), ident.2 = setdiff(List3, "C3_STB_d14"),test.use = "MAST", only.pos = TRUE)
DEListI4 <- FindMarkers(D, ident.1 = c("C4_STB_d14"), ident.2 = setdiff(List4, "C4_STB_d14"),test.use = "MAST", only.pos = TRUE)
DEListI6 <- FindMarkers(D, ident.1 = c("C6_STB_d14"), ident.2 = setdiff(List6, "C6_STB_d14"),test.use = "MAST", only.pos = TRUE)
DEListI7 <- FindMarkers(D, ident.1 = c("C7_STB_d14"), ident.2 = setdiff(List7, "C7_STB_d14"),test.use = "MAST", only.pos = TRUE)
DEListI8 <- FindMarkers(D, ident.1 = c("C8_STB_d14"), ident.2 = setdiff(List8, "C8_STB_d14"),test.use = "MAST", only.pos = TRUE)

DEListI1 <- DEListI1[which(DEListI1$p_val_adj<0.001),]
DEListI1 <- DEListI1[order(DEListI1$avg_log2FC,decreasing=TRUE),]
DEListI2 <- DEListI2[which(DEListI2$p_val_adj<0.001),]
DEListI2 <- DEListI2[order(DEListI2$avg_log2FC,decreasing=TRUE),]
DEListI3 <- DEListI3[which(DEListI3$p_val_adj<0.001),]
DEListI3 <- DEListI3[order(DEListI3$avg_log2FC,decreasing=TRUE),]
DEListI4 <- DEListI4[which(DEListI4$p_val_adj<0.001),]
DEListI4 <- DEListI4[order(DEListI4$avg_log2FC,decreasing=TRUE),]
DEListI6 <- DEListI6[which(DEListI6$p_val_adj<0.001),]
DEListI6 <- DEListI6[order(DEListI6$avg_log2FC,decreasing=TRUE),]
DEListI7 <- DEListI7[which(DEListI7$p_val_adj<0.001),]
DEListI7 <- DEListI7[order(DEListI7$avg_log2FC,decreasing=TRUE),]
DEListI8 <- DEListI8[which(DEListI8$p_val_adj<0.001),]
DEListI8 <- DEListI8[order(DEListI8$avg_log2FC,decreasing=TRUE),]

L9 <- Reduce(intersect, list(rownames(DEListI1),rownames(DEListI2),rownames(DEListI3),rownames(DEListI4),rownames(DEListI6),rownames(DEListI7),rownames(DEListI8)))
L9 <- unique( c(
  intersect(rownames(DEListI1),rownames(DEListI2)),
  intersect(rownames(DEListI1),rownames(DEListI3)),
  intersect(rownames(DEListI1),rownames(DEListI4)),
  intersect(rownames(DEListI1),rownames(DEListI6)),
  intersect(rownames(DEListI1),rownames(DEListI7)),
  intersect(rownames(DEListI1),rownames(DEListI8)),
  intersect(rownames(DEListI2),rownames(DEListI3)),
  intersect(rownames(DEListI2),rownames(DEListI4)),
  intersect(rownames(DEListI2),rownames(DEListI6)),
  intersect(rownames(DEListI2),rownames(DEListI7)),
  intersect(rownames(DEListI2),rownames(DEListI8)),
  intersect(rownames(DEListI3),rownames(DEListI4)),
  intersect(rownames(DEListI3),rownames(DEListI6)),
  intersect(rownames(DEListI3),rownames(DEListI7)),
  intersect(rownames(DEListI3),rownames(DEListI8)),
  intersect(rownames(DEListI4),rownames(DEListI6)),
  intersect(rownames(DEListI4),rownames(DEListI7)),
  intersect(rownames(DEListI4),rownames(DEListI8)),
  intersect(rownames(DEListI6),rownames(DEListI7)),
  intersect(rownames(DEListI6),rownames(DEListI8)),
  intersect(rownames(DEListI7),rownames(DEListI8))
) )
#4

#First do CTB_d14
DEListJ1 <- FindMarkers(D, ident.1 = c("C1_CTB_d14"), ident.2 = setdiff(List1, "C1_CTB_d14"),test.use = "MAST", only.pos = TRUE)
DEListJ2 <- FindMarkers(D, ident.1 = c("C2_CTB_d14"), ident.2 = setdiff(List2, "C2_CTB_d14"),test.use = "MAST", only.pos = TRUE)
#DEListJ3 <- FindMarkers(D, ident.1 = c("C3_CTB_d14"), ident.2 = setdiff(List3, "C3_CTB_d14"),test.use = "MAST", only.pos = TRUE)
DEListJ4 <- FindMarkers(D, ident.1 = c("C4_CTB_d14"), ident.2 = setdiff(List4, "C4_CTB_d14"),test.use = "MAST", only.pos = TRUE)
DEListJ6 <- FindMarkers(D, ident.1 = c("C6_CTB_d14"), ident.2 = setdiff(List6, "C6_CTB_d14"),test.use = "MAST", only.pos = TRUE)
#DEListJ7 <- FindMarkers(D, ident.1 = c("C7_CTB_d14"), ident.2 = setdiff(List7, "C7_CTB_d14"),test.use = "MAST", only.pos = TRUE)
DEListJ8 <- FindMarkers(D, ident.1 = c("C8_CTB_d14"), ident.2 = setdiff(List8, "C8_CTB_d14"),test.use = "MAST", only.pos = TRUE)

DEListJ1 <- DEListJ1[which(DEListJ1$p_val_adj<0.01),]
DEListJ1 <- DEListJ1[order(DEListJ1$avg_log2FC,decreasing=TRUE),]
DEListJ2 <- DEListJ2[which(DEListJ2$p_val_adj<0.01),]
DEListJ2 <- DEListJ2[order(DEListJ2$avg_log2FC,decreasing=TRUE),]
DEListJ4 <- DEListJ4[which(DEListJ4$p_val_adj<0.01),]
DEListJ4 <- DEListJ4[order(DEListJ4$avg_log2FC,decreasing=TRUE),]
DEListJ6 <- DEListJ6[which(DEListJ6$p_val_adj<0.01),]
DEListJ6 <- DEListJ6[order(DEListJ6$avg_log2FC,decreasing=TRUE),]
DEListJ8 <- DEListJ8[which(DEListJ8$p_val_adj<0.01),]
DEListJ8 <- DEListJ8[order(DEListJ8$avg_log2FC,decreasing=TRUE),]

L10 <- Reduce(intersect, list(rownames(DEListJ1),rownames(DEListJ2),rownames(DEListJ4),rownames(DEListJ6),rownames(DEListJ8)))

L10 <- unique( c(
  intersect(rownames(DEListJ1),rownames(DEListJ2)),
  intersect(rownames(DEListJ1),rownames(DEListJ4)),
  intersect(rownames(DEListJ1),rownames(DEListJ6)),
  intersect(rownames(DEListJ1),rownames(DEListJ8)),
  intersect(rownames(DEListJ2),rownames(DEListJ4)),
  intersect(rownames(DEListJ2),rownames(DEListJ6)),
  intersect(rownames(DEListJ2),rownames(DEListJ8)),
  intersect(rownames(DEListJ4),rownames(DEListJ6)),
  intersect(rownames(DEListJ4),rownames(DEListJ8)),
  intersect(rownames(DEListJ6),rownames(DEListJ8))
) )

L10 <- unique( c(rownames(DEListJ1),rownames(DEListJ2),rownames(DEListJ4), rownames(DEListJ6), rownames(DEListJ8) ) )


#None

#First do EVT_d14
#DEListK1 <- FindMarkers(D, ident.1 = c("C1_EVT_d14"), ident.2 = setdiff(List1, "C1_EVT_d14"),test.use = "MAST", only.pos = TRUE)
DEListK2 <- FindMarkers(D, ident.1 = c("C2_EVT_d14"), ident.2 = setdiff(List2, "C2_EVT_d14"),test.use = "MAST", only.pos = TRUE)
#DEListK3 <- FindMarkers(D, ident.1 = c("C3_EVT_d14"), ident.2 = setdiff(List3, "C3_EVT_d14"),test.use = "MAST", only.pos = TRUE)
DEListK4 <- FindMarkers(D, ident.1 = c("C4_EVT_d14"), ident.2 = setdiff(List4, "C4_EVT_d14"),test.use = "MAST", only.pos = TRUE)
DEListK6 <- FindMarkers(D, ident.1 = c("C6_EVT_d14"), ident.2 = setdiff(List6, "C6_EVT_d14"),test.use = "MAST", only.pos = TRUE)
#DEListK8 <- FindMarkers(D, ident.1 = c("C8_EVT_d14"), ident.2 = setdiff(List8, "C8_EVT_d14"),test.use = "MAST", only.pos = TRUE)

DEListK2 <- DEListK2[which(DEListK2$p_val_adj<0.001),]
DEListK2 <- DEListK2[order(DEListK2$avg_log2FC,decreasing=TRUE),]
DEListK4 <- DEListK4[which(DEListK4$p_val_adj<0.001),]
DEListK4 <- DEListK4[order(DEListK4$avg_log2FC,decreasing=TRUE),]
DEListK6 <- DEListK6[which(DEListK6$p_val_adj<0.001),]
DEListK6 <- DEListK6[order(DEListK6$avg_log2FC,decreasing=TRUE),]

L11 <- Reduce(intersect, list(rownames(DEListK2),rownames(DEListK4),rownames(DEListK6)))
L11 <- unique( c(
  intersect(rownames(DEListK2),rownames(DEListK4)),
  intersect(rownames(DEListK2),rownames(DEListK6)),
  intersect(rownames(DEListK4),rownames(DEListK6))
) )
#First do ExMes_d14
DEListL1 <- FindMarkers(D, ident.1 = c("C1_ExMes_d14"), ident.2 = setdiff(List1, "C1_ExMes_d14"),test.use = "MAST", only.pos = TRUE)
DEListL2 <- FindMarkers(D, ident.1 = c("C2_ExMes_d14"), ident.2 = setdiff(List2, "C2_ExMes_d14"),test.use = "MAST", only.pos = TRUE)
#DEListL3 <- FindMarkers(D, ident.1 = c("C3_ExMes_d14"), ident.2 = setdiff(List3, "C3_ExMes_d14"),test.use = "MAST", only.pos = TRUE)
#DEListL4 <- FindMarkers(D, ident.1 = c("C4_ExMes_d14"), ident.2 = setdiff(List4, "C4_ExMes_d14"),test.use = "MAST", only.pos = TRUE)
DEListL6 <- FindMarkers(D, ident.1 = c("C6_ExMes_d14"), ident.2 = setdiff(List6, "C6_ExMes_d14"),test.use = "MAST", only.pos = TRUE)
DEListL7 <- FindMarkers(D, ident.1 = c("C7_ExMes_d14"), ident.2 = setdiff(List7, "C7_ExMes_d14"),test.use = "MAST", only.pos = TRUE)
DEListL8 <- FindMarkers(D, ident.1 = c("C8_ExMes_d14"), ident.2 = setdiff(List8, "C8_ExMes_d14"),test.use = "MAST", only.pos = TRUE)

DEListL1 <- DEListL1[which(DEListL1$p_val_adj<0.05),]
DEListL1 <- DEListL1[order(DEListL1$avg_log2FC,decreasing=TRUE),]
DEListL2 <- DEListL2[which(DEListL2$p_val_adj<0.05),]
DEListL2 <- DEListL2[order(DEListL2$avg_log2FC,decreasing=TRUE),]
DEListL6 <- DEListL6[which(DEListL6$p_val_adj<0.05),]
DEListL6 <- DEListL6[order(DEListL6$avg_log2FC,decreasing=TRUE),]
DEListL7 <- DEListL7[which(DEListL7$p_val_adj<0.05),]
DEListL7 <- DEListL7[order(DEListL7$avg_log2FC,decreasing=TRUE),]
DEListL8 <- DEListL8[which(DEListL8$p_val_adj<0.05),]
DEListL8 <- DEListL8[order(DEListL8$avg_log2FC,decreasing=TRUE),]

L12 <- Reduce(intersect, list(rownames(DEListL1),rownames(DEListL2),rownames(DEListL6),rownames(DEListL7),rownames(DEListL8)))

L12 <- unique( c(
  intersect(rownames(DEListL1),rownames(DEListL2)),
  intersect(rownames(DEListL1),rownames(DEListL6)),
  intersect(rownames(DEListL1),rownames(DEListL7)),
  intersect(rownames(DEListL1),rownames(DEListL8)),
  intersect(rownames(DEListL2),rownames(DEListL6)),
  intersect(rownames(DEListL2),rownames(DEListL7)),
  intersect(rownames(DEListL2),rownames(DEListL8)),
  intersect(rownames(DEListL6),rownames(DEListL7)),
  intersect(rownames(DEListL6),rownames(DEListL8)),
  intersect(rownames(DEListL7),rownames(DEListL8))
) )
#very few
L12 <- unique(c(rownames(DEListL1), rownames(DEListL2), rownames(DEListL6), rownames(DEListL7), rownames(DEListL8) ))

#First do Am_d14
DEListM1 <- FindMarkers(D, ident.1 = c("C1_Am_d14"), ident.2 = setdiff(List1, "C1_Am_d14"),test.use = "MAST", only.pos = TRUE)
DEListM2 <- FindMarkers(D, ident.1 = c("C2_Am_d14"), ident.2 = setdiff(List2, "C2_Am_d14"),test.use = "MAST", only.pos = TRUE)
#DEListM3 <- FindMarkers(D, ident.1 = c("C3_Am_d14"), ident.2 = setdiff(List3, "C3_Am_d14"),test.use = "MAST", only.pos = TRUE)
#DEListM4 <- FindMarkers(D, ident.1 = c("C4_Am_d14"), ident.2 = setdiff(List4, "C4_Am_d14"),test.use = "MAST", only.pos = TRUE)
DEListM6 <- FindMarkers(D, ident.1 = c("C6_Am_d14"), ident.2 = setdiff(List6, "C6_Am_d14"),test.use = "MAST", only.pos = TRUE)
DEListM7 <- FindMarkers(D, ident.1 = c("C7_Am_d14"), ident.2 = setdiff(List7, "C7_Am_d14"),test.use = "MAST", only.pos = TRUE)
#DEListM8 <- FindMarkers(D, ident.1 = c("C8_Am_d14"), ident.2 = setdiff(List8, "C8_Am_d14"),test.use = "MAST", only.pos = TRUE)

DEListM1 <- DEListM1[which(DEListM1$p_val_adj<0.01),]
DEListM1 <- DEListM1[order(DEListM1$avg_log2FC,decreasing=TRUE),]
DEListM2 <- DEListM2[which(DEListM2$p_val_adj<0.01),]
DEListM2 <- DEListM2[order(DEListM2$avg_log2FC,decreasing=TRUE),]
DEListM6 <- DEListM6[which(DEListM6$p_val_adj<0.01),]
DEListM6 <- DEListM6[order(DEListM6$avg_log2FC,decreasing=TRUE),]
DEListM7 <- DEListM7[which(DEListM7$p_val_adj<0.01),]
DEListM7 <- DEListM7[order(DEListM7$avg_log2FC,decreasing=TRUE),]

L13 <- Reduce(intersect, list(rownames(DEListM1),rownames(DEListM2),rownames(DEListM6),rownames(DEListM7)))
L13 <- unique( c(
  intersect(rownames(DEListM1),rownames(DEListM2)),
  intersect(rownames(DEListM1),rownames(DEListM6)),
  intersect(rownames(DEListM1),rownames(DEListM7)),
  intersect(rownames(DEListM2),rownames(DEListM6)),
  intersect(rownames(DEListM2),rownames(DEListM7)),
  intersect(rownames(DEListM6),rownames(DEListM7))
) )
L13 <- unique( c(rownames(DEListM1),rownames(DEListM2),rownames(DEListM6),rownames(DEListM7)) )

#First do Hyp_d14
DEListN1 <- FindMarkers(D, ident.1 = c("C1_Hyp_d14"), ident.2 = setdiff(List1, "C1_Hyp_d14"),test.use = "MAST", only.pos = TRUE)
DEListN2 <- FindMarkers(D, ident.1 = c("C2_Hyp_d14"), ident.2 = setdiff(List2, "C2_Hyp_d14"),test.use = "MAST", only.pos = TRUE)
#DEListN3 <- FindMarkers(D, ident.1 = c("C3_Hyp_d14"), ident.2 = setdiff(List3, "C3_Hyp_d14"),test.use = "MAST", only.pos = TRUE)
#DEListN4 <- FindMarkers(D, ident.1 = c("C4_Hyp_d14"), ident.2 = setdiff(List4, "C4_Hyp_d14"),test.use = "MAST", only.pos = TRUE)
#DEListN6 <- FindMarkers(D, ident.1 = c("C6_Hyp_d14"), ident.2 = setdiff(List6, "C6_Hyp_d14"),test.use = "MAST", only.pos = TRUE)
DEListN7 <- FindMarkers(D, ident.1 = c("C7_Hyp_d14"), ident.2 = setdiff(List7, "C7_Hyp_d14"),test.use = "MAST", only.pos = TRUE)
#DEListN8 <- FindMarkers(D, ident.1 = c("C8_Hyp_d14"), ident.2 = setdiff(List8, "C8_Hyp_d14"),test.use = "MAST", only.pos = TRUE)

DEListN1 <- DEListN1[which(DEListN1$p_val_adj<0.01),]
DEListN1 <- DEListN1[order(DEListN1$avg_log2FC,decreasing=TRUE),]
DEListN2 <- DEListN2[which(DEListN2$p_val_adj<0.01),]
DEListN2 <- DEListN2[order(DEListN2$avg_log2FC,decreasing=TRUE),]
DEListN7 <- DEListN7[which(DEListN7$p_val_adj<0.01),]
DEListN7 <- DEListN7[order(DEListN7$avg_log2FC,decreasing=TRUE),]

L14 <- Reduce(intersect, list(rownames(DEListN1),rownames(DEListN2),rownames(DEListN7)))
L14 <- unique( c(
  intersect(rownames(DEListN1),rownames(DEListN2)),
  intersect(rownames(DEListN1),rownames(DEListN7)),
  intersect(rownames(DEListN2),rownames(DEListN7))
) )

L14 <- unique( c(rownames(DEListN1),rownames(DEListN2),rownames(DEListN7)) )


#First do Am/EmDisc_d14
DEListO1 <- FindMarkers(D, ident.1 = c("C1_Am/EmDisc_d14"), ident.2 = setdiff(List1, "C1_Am/EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
DEListO2 <- FindMarkers(D, ident.1 = c("C2_Am/EmDisc_d14"), ident.2 = setdiff(List2, "C2_Am/EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
#DEListO3 <- FindMarkers(D, ident.1 = c("C3_Am/EmDisc_d14"), ident.2 = setdiff(List3, "C3_Am/EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
DEListO4 <- FindMarkers(D, ident.1 = c("C4_Am/EmDisc_d14"), ident.2 = setdiff(List4, "C4_Am/EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
DEListO6 <- FindMarkers(D, ident.1 = c("C6_Am/EmDisc_d14"), ident.2 = setdiff(List6, "C6_Am/EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
#DEListO7 <- FindMarkers(D, ident.1 = c("C7_Am/EmDisc_d14"), ident.2 = setdiff(List7, "C7_Am/EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
#DEListO8 <- FindMarkers(D, ident.1 = c("C8_Am/EmDisc_d14"), ident.2 = setdiff(List8, "C8_Am/EmDisc_d14"),test.use = "MAST", only.pos = TRUE)

DEListO1 <- DEListO1[which(DEListO1$p_val_adj<0.01),]
DEListO1 <- DEListO1[order(DEListO1$avg_log2FC,decreasing=TRUE),]
DEListO2 <- DEListO2[which(DEListO2$p_val_adj<0.01),]
DEListO2 <- DEListO2[order(DEListO2$avg_log2FC,decreasing=TRUE),]
DEListO4 <- DEListO4[which(DEListO4$p_val_adj<0.01),]
DEListO4 <- DEListO4[order(DEListO4$avg_log2FC,decreasing=TRUE),]
DEListO6 <- DEListO6[which(DEListO6$p_val_adj<0.01),]
DEListO6 <- DEListO6[order(DEListO6$avg_log2FC,decreasing=TRUE),]

L15 <- Reduce(intersect, list(rownames(DEListO1),rownames(DEListO2),rownames(DEListO4),rownames(DEListO6)))
L15 <- unique( c(
  intersect(rownames(DEListO1),rownames(DEListO2)),
  intersect(rownames(DEListO1),rownames(DEListO4)),
  intersect(rownames(DEListO1),rownames(DEListO6)),
  intersect(rownames(DEListO2),rownames(DEListO4)),
  intersect(rownames(DEListO2),rownames(DEListO6)),
  intersect(rownames(DEListO4),rownames(DEListO6))
) )


L15 <- unique( c(rownames(DEListO1),rownames(DEListO2),rownames(DEListO4),rownames(DEListO6) ) )

#First do EmDisc_d14
DEListP1 <- FindMarkers(D, ident.1 = c("C1_EmDisc_d14"), ident.2 = setdiff(List1, "C1_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
DEListP2 <- FindMarkers(D, ident.1 = c("C2_EmDisc_d14"), ident.2 = setdiff(List2, "C2_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
#DEListP3 <- FindMarkers(D, ident.1 = c("C3_EmDisc_d14"), ident.2 = setdiff(List3, "C3_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
#DEListP4 <- FindMarkers(D, ident.1 = c("C4_EmDisc_d14"), ident.2 = setdiff(List4, "C4_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
DEListP6 <- FindMarkers(D, ident.1 = c("C6_EmDisc_d14"), ident.2 = setdiff(List6, "C6_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
DEListP7 <- FindMarkers(D, ident.1 = c("C7_EmDisc_d14"), ident.2 = setdiff(List7, "C7_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
DEListP8 <- FindMarkers(D, ident.1 = c("C8_EmDisc_d14"), ident.2 = setdiff(List8, "C8_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)

DEListP1 <- DEListP1[which(DEListP1$p_val_adj<0.05),]
DEListP1 <- DEListP1[order(DEListP1$avg_log2FC,decreasing=TRUE),]
DEListP2 <- DEListP2[which(DEListP2$p_val_adj<0.05),]
DEListP2 <- DEListP2[order(DEListP2$avg_log2FC,decreasing=TRUE),]
DEListP6 <- DEListP6[which(DEListP6$p_val_adj<0.05),]
DEListP6 <- DEListP6[order(DEListP6$avg_log2FC,decreasing=TRUE),]
DEListP7 <- DEListP7[which(DEListP7$p_val_adj<0.05),]
DEListP7 <- DEListP7[order(DEListP7$avg_log2FC,decreasing=TRUE),]
DEListP8 <- DEListP8[which(DEListP8$p_val_adj<0.05),]
DEListP8 <- DEListP8[order(DEListP8$avg_log2FC,decreasing=TRUE),]

L16 <- Reduce(intersect, list(rownames(DEListP1),rownames(DEListP2),rownames(DEListP6),rownames(DEListP7),rownames(DEListP8)))
L16 <- unique( c(intersect(rownames(DEListP1),rownames(DEListP2)),
intersect(rownames(DEListP1),rownames(DEListP6)),
intersect(rownames(DEListP1),rownames(DEListP7)),
intersect(rownames(DEListP1),rownames(DEListP8)),
intersect(rownames(DEListP2),rownames(DEListP6)),
intersect(rownames(DEListP2),rownames(DEListP7)),
intersect(rownames(DEListP2),rownames(DEListP8)),
intersect(rownames(DEListP6),rownames(DEListP7)),
intersect(rownames(DEListP6),rownames(DEListP8)),
intersect(rownames(DEListP7),rownames(DEListP8)) ) )
L16 <- unique( c(rownames(DEListP1),rownames(DEListP2),rownames(DEListP6),rownames(DEListP7), rownames(DEListP8) ) )



Dat <- Av[L1,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L1_Stromal.pdf",sep=""),scale="row",width=20,height=200)
mat_breaks <- seq(10, 100, length.out = 20)
pheatmap((Dat)*100-1,color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L1_Unscaled_Stromal.pdf",sep=""),width=20,height=200)




Dat <- Av[L2,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L2_Lumenal.pdf",sep=""),scale="row",width=20,height=160)

mat_breaks <- seq(10, 100, length.out = 20)
pheatmap((Dat)*100-1,color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L2_Unscaled_Lumenal.pdf",sep=""),width=20,height=160)



Dat <- Av[L3,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L3_Prolif.pdf",sep=""),scale="row",width=20,height=100)
mat_breaks <- seq(10, 100, length.out = 20)
pheatmap((Dat)*100-1,color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L3_Unscaled_Prolif.pdf",sep=""),width=20,height=100)



Dat <- Av[L4,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L4_Glands.pdf",sep=""),scale="row",width=20,height=80)
mat_breaks <- seq(10, 100, length.out = 20)
pheatmap((Dat)*100-1,color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L4_Unscaled_Glands.pdf",sep=""),width=20,height=80)




Dat <- Av[L5,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L5_Cil.pdf",sep=""),scale="row",width=20,height=30)

mat_breaks <- seq(10, 100, length.out = 20)
pheatmap((Dat)*100-1,color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L5_Unscaled_Cil.pdf",sep=""),width=20,height=30)


#Dat <- Av[L5,exps]
#mat_breaks <- seq(-2, 2, length.out = 20)
#redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
#pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L5_Cil.pdf",sep=""),scale="row",width=20,height=30)


Dat <- Av[L7,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L7_SOX9.pdf",sep=""),scale="row",width=20,height=30)

mat_breaks <- seq(10, 100, length.out = 20)
pheatmap((Dat)*100-1,color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L7_Unscaled_Cil.pdf",sep=""),width=20,height=30)

#Dat <- Av[L8,exps]
#mat_breaks <- seq(-2, 2, length.out = 20)
#redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
#pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L7_SOX9LRG5.pdf",sep=""),scale="row",width=20,height=30)


Dat <- Av[L9,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L9_STB.pdf",sep=""),scale="row",width=20,height=100)

mat_breaks <- seq(10, 100, length.out = 20)
pheatmap((Dat)*100-1,color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L9_Unscaled_STB.pdf",sep=""),width=20,height=100)





Dat <- Av[L10,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L10_CTB.pdf",sep=""),scale="row",width=20,height=100)


mat_breaks <- seq(10, 100, length.out = 20)
pheatmap((Dat)*100-1,color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L10_Unscaled_CTB.pdf",sep=""),width=20,height=100)






Dat <- Av[L11,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L11_EVT.pdf",sep=""),scale="row",width=20,height=30)


mat_breaks <- seq(10, 100, length.out = 20)
pheatmap((Dat)*100-1,color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L11_Unscaled_EVT.pdf",sep=""),width=20,height=30)




Dat <- Av[L12,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L12_ExMes.pdf",sep=""),scale="row",width=20,height=30)

mat_breaks <- seq(10, 100, length.out = 20)
pheatmap((Dat)*100-1,color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L12_Unscaled_ExMes.pdf",sep=""),width=20,height=30)


Dat <- Av[L13,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L13_Am.pdf",sep=""),scale="row",width=20,height=20)
mat_breaks <- seq(10, 100, length.out = 20)
pheatmap((Dat)*100-1,color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L13_Unscaled_Am.pdf",sep=""),width=20,height=20)


Dat <- Av[L14,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L14_Hyp.pdf",sep=""),scale="row",width=20,height=100)
mat_breaks <- seq(10, 100, length.out = 20)
pheatmap((Dat)*100-1,color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L14_Unscaled_Hyp.pdf",sep=""),width=20,height=100)


Dat <- Av[L15,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L15_EmDisc.pdf",sep=""),scale="row",width=20,height=100)

mat_breaks <- seq(10, 100, length.out = 20)
pheatmap((Dat)*100-1,color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L15_Unscaled_EmDisc.pdf",sep=""),width=20,height=100)


write.table(as.data.frame(Av[L15,exps]),file=paste(saveext,"/EmDiscExpression.csv",sep=""))

#First do CTB_d14
DEListQ1 <- FindMarkers(D, ident.1 = c("C1_CTB_d14"), ident.2 = "C1_STB_d14",test.use = "MAST", only.pos = TRUE)
DEListQ2 <- FindMarkers(D, ident.1 = c("C2_CTB_d14"), ident.2 = "C2_STB_d14",test.use = "MAST", only.pos = TRUE)
#DEListJ3 <- FindMarkers(D, ident.1 = c("C3_CTB_d14"), ident.2 = setdiff(List3, "C3_STB_d14"),test.use = "MAST", only.pos = TRUE)
DEListQ4 <- FindMarkers(D, ident.1 = c("C4_CTB_d14"), ident.2 = "C4_STB_d14",test.use = "MAST", only.pos = TRUE)
DEListQ6 <- FindMarkers(D, ident.1 = c("C6_CTB_d14"), ident.2 = "C6_STB_d14",test.use = "MAST", only.pos = TRUE)
#DEListJ7 <- FindMarkers(D, ident.1 = c("C7_CTB_d14"), ident.2 = setdiff(List7, "C7_CTB_d14"),test.use = "MAST", only.pos = TRUE)
DEListQ8 <- FindMarkers(D, ident.1 = c("C8_CTB_d14"), ident.2 = "C8_STB_d14",test.use = "MAST", only.pos = TRUE)

DEListQ1 <- DEListQ1[which(DEListQ1$p_val_adj<0.01),]
DEListQ1 <- DEListQ1[order(DEListQ1$avg_log2FC,decreasing=TRUE),]
DEListQ2 <- DEListQ2[which(DEListQ2$p_val_adj<0.01),]
DEListQ2 <- DEListQ2[order(DEListQ2$avg_log2FC,decreasing=TRUE),]
DEListQ4 <- DEListQ4[which(DEListQ4$p_val_adj<0.01),]
DEListQ4 <- DEListQ4[order(DEListQ4$avg_log2FC,decreasing=TRUE),]
DEListQ6 <- DEListQ6[which(DEListQ6$p_val_adj<0.01),]
DEListQ6 <- DEListQ6[order(DEListQ6$avg_log2FC,decreasing=TRUE),]
DEListQ8 <- DEListQ8[which(DEListQ8$p_val_adj<0.01),]
DEListQ8 <- DEListQ8[order(DEListQ8$avg_log2FC,decreasing=TRUE),]

L17 <- unique( c(
  intersect(rownames(DEListQ1),rownames(DEListQ2)),
  intersect(rownames(DEListQ1),rownames(DEListQ4)),
  intersect(rownames(DEListQ1),rownames(DEListQ6)),
  intersect(rownames(DEListQ1),rownames(DEListQ8)),
  intersect(rownames(DEListQ2),rownames(DEListQ4)),
  intersect(rownames(DEListQ2),rownames(DEListQ6)),
  intersect(rownames(DEListQ2),rownames(DEListQ8)),
  intersect(rownames(DEListQ4),rownames(DEListQ6)),
  intersect(rownames(DEListQ4),rownames(DEListQ8)),
  intersect(rownames(DEListQ6),rownames(DEListQ8))
) )

L17 <- unique( c(rownames(DEListQ1),rownames(DEListQ2),rownames(DEListQ4), rownames(DEListQ6), rownames(DEListQ8) ) )



Dat <- Av[L17,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L17_CTB_STB.pdf",sep=""),scale="row",width=20,height=150)

mat_breaks <- seq(0, 20, length.out = 20)
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L17_Unscaled_CTB_STB.pdf",sep=""),width=20,height=150)




#First do CTB_d14
DEListR1 <- FindMarkers(D, ident.1 = c("C1_EmDisc_d14"), ident.2 = "C1_Am_d14",test.use = "MAST", only.pos = FALSE)
DEListR2 <- FindMarkers(D, ident.1 = c("C2_EmDisc_d14"), ident.2 = "C2_Am_d14",test.use = "MAST", only.pos = FALSE)
DEListR6 <- FindMarkers(D, ident.1 = c("C6_EmDisc_d14"), ident.2 = "C6_Am_d14",test.use = "MAST", only.pos = FALSE)
DEListR7 <- FindMarkers(D, ident.1 = c("C7_EmDisc_d14"), ident.2 = "C7_Am_d14",test.use = "MAST", only.pos = FALSE)

DEListR1 <- DEListR1[which(DEListR1$p_val<0.01),]
DEListR1 <- DEListR1[order(DEListR1$avg_log2FC,decreasing=TRUE),]
DEListR2 <- DEListR2[which(DEListR2$p_val<0.01),]
DEListR2 <- DEListR2[order(DEListR2$avg_log2FC,decreasing=TRUE),]
DEListR6 <- DEListR6[which(DEListR6$p_val<0.01),]
DEListR6 <- DEListR6[order(DEListR6$avg_log2FC,decreasing=TRUE),]
DEListR7 <- DEListR7[which(DEListR7$p_val<0.01),]
DEListR7 <- DEListR7[order(DEListR7$avg_log2FC,decreasing=TRUE),]


L18 <- unique( c(rownames(DEListR1),rownames(DEListR2),rownames(DEListR6), rownames(DEListR7) ) )



Dat <- Av[L18,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L18_Am_EmD.pdf",sep=""),scale="row",width=20,height=60)

mat_breaks <- seq(0, 20, length.out = 20)
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L18_Unscaled_Am_EmD.pdf",sep=""),width=20,height=60)







DEListL1 <- FindMarkers(D, ident.1 = c("C1_ExMes_d14"), ident.2 = c("C1_CTB_d14","C1_STB_d14","C1_Hyp_d14","C1_Am/EmDisc_d14","C1_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
DEListL2 <- FindMarkers(D, ident.1 = c("C2_ExMes_d14"), ident.2 = c("C2_CTB_d14","C2_STB_d14","C2_EVT_d14","C2_Hyp_d14","C2_Am/EmDisc_d14","C2_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
#DEListL3 <- FindMarkers(D, ident.1 = c("C3_ExMes_d14"), ident.2 = setdiff(List3, "C3_ExMes_d14"),test.use = "MAST", only.pos = TRUE)
#DEListL4 <- FindMarkers(D, ident.1 = c("C4_ExMes_d14"), ident.2 = setdiff(List4, "C4_ExMes_d14"),test.use = "MAST", only.pos = TRUE)
DEListL6 <- FindMarkers(D, ident.1 = c("C6_ExMes_d14"), ident.2 = c("C6_CTB_d14","C6_STB_d14","C6_EVT_d14","C6_Am_d14","C6_Am/EmDisc_d14","C6_EmDisc_d14" ),test.use = "MAST", only.pos = TRUE)
DEListL7 <- FindMarkers(D, ident.1 = c("C7_ExMes_d14"), ident.2 = c("C7_STB_d14","C7_EVT_d14","C7_Hyp_d14","C7_Am_d14","C7_Am/EmDisc_d14","C7_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
DEListL8 <- FindMarkers(D, ident.1 = c("C8_ExMes_d14"), ident.2 = c("C8_CTB_d14","C8_STB_d14","C8_EVT_d14","C8_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)

DEListL1 <- DEListL1[which(DEListL1$p_val_adj<0.01),]
DEListL1 <- DEListL1[order(DEListL1$avg_log2FC,decreasing=TRUE),]
DEListL2 <- DEListL2[which(DEListL2$p_val_adj<0.01),]
DEListL2 <- DEListL2[order(DEListL2$avg_log2FC,decreasing=TRUE),]
DEListL6 <- DEListL6[which(DEListL6$p_val_adj<0.01),]
DEListL6 <- DEListL6[order(DEListL6$avg_log2FC,decreasing=TRUE),]
DEListL7 <- DEListL7[which(DEListL7$p_val_adj<0.01),]
DEListL7 <- DEListL7[order(DEListL7$avg_log2FC,decreasing=TRUE),]
DEListL8 <- DEListL8[which(DEListL8$p_val_adj<0.01),]
DEListL8 <- DEListL8[order(DEListL8$avg_log2FC,decreasing=TRUE),]

L19 <- unique( c(
  intersect(rownames(DEListL1),rownames(DEListL2)),
  intersect(rownames(DEListL1),rownames(DEListL6)),
  intersect(rownames(DEListL1),rownames(DEListL7)),
  intersect(rownames(DEListL1),rownames(DEListL8)),
  intersect(rownames(DEListL2),rownames(DEListL6)),
  intersect(rownames(DEListL2),rownames(DEListL7)),
  intersect(rownames(DEListL2),rownames(DEListL8)),
  intersect(rownames(DEListL6),rownames(DEListL7)),
  intersect(rownames(DEListL6),rownames(DEListL8)),
  intersect(rownames(DEListL7),rownames(DEListL8))
) )


Dat <- Av[L19,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L19_ExMes_Emb.pdf",sep=""),scale="row",width=20,height=60)


mat_breaks <- seq(0, 20, length.out = 20)
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L19_Unscaled_ExMes_Emb.pdf",sep=""),width=20,height=60)







#First do EmDisc_d14
DEListP1 <- FindMarkers(D, ident.1 = c("C1_EmDisc_d14"), ident.2 = setdiff(List1, "C1_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
DEListP2 <- FindMarkers(D, ident.1 = c("C2_EmDisc_d14"), ident.2 = setdiff(List2, "C2_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
#DEListP3 <- FindMarkers(D, ident.1 = c("C3_EmDisc_d14"), ident.2 = setdiff(List3, "C3_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
#DEListP4 <- FindMarkers(D, ident.1 = c("C4_EmDisc_d14"), ident.2 = setdiff(List4, "C4_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
DEListP6 <- FindMarkers(D, ident.1 = c("C6_EmDisc_d14"), ident.2 = setdiff(List6, "C6_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
DEListP7 <- FindMarkers(D, ident.1 = c("C7_EmDisc_d14"), ident.2 = setdiff(List7, "C7_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)
DEListP8 <- FindMarkers(D, ident.1 = c("C8_EmDisc_d14"), ident.2 = setdiff(List8, "C8_EmDisc_d14"),test.use = "MAST", only.pos = TRUE)

DEListP1 <- DEListP1[which(DEListP1$p_val<0.05),]
DEListP1 <- DEListP1[order(DEListP1$avg_log2FC,decreasing=TRUE),]
DEListP2 <- DEListP2[which(DEListP2$p_val<0.05),]
DEListP2 <- DEListP2[order(DEListP2$avg_log2FC,decreasing=TRUE),]
DEListP6 <- DEListP6[which(DEListP6$p_val<0.05),]
DEListP6 <- DEListP6[order(DEListP6$avg_log2FC,decreasing=TRUE),]
DEListP7 <- DEListP7[which(DEListP7$p_val<0.05),]
DEListP7 <- DEListP7[order(DEListP7$avg_log2FC,decreasing=TRUE),]
DEListP8 <- DEListP8[which(DEListP8$p_val<0.05),]
DEListP8 <- DEListP8[order(DEListP8$avg_log2FC,decreasing=TRUE),]

L16 <- unique( c(intersect(rownames(DEListP1),rownames(DEListP2)),
                 intersect(rownames(DEListP1),rownames(DEListP6)),
                 intersect(rownames(DEListP1),rownames(DEListP7)),
                 intersect(rownames(DEListP1),rownames(DEListP8)),
                 intersect(rownames(DEListP2),rownames(DEListP6)),
                 intersect(rownames(DEListP2),rownames(DEListP7)),
                 intersect(rownames(DEListP2),rownames(DEListP8)),
                 intersect(rownames(DEListP6),rownames(DEListP7)),
                 intersect(rownames(DEListP6),rownames(DEListP8)),
                 intersect(rownames(DEListP7),rownames(DEListP8)) ) )



Dat <- Av[L16,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L16_EmDisc_Emb.pdf",sep=""),scale="row",width=20,height=60)


NaiveList <- intersect(PrimedNaive$Gene[which(PrimedNaive$logFC>0 & PrimedNaive$FDR<0.05)],L16)
PrimeList <- intersect(PrimedNaive$Gene[which(PrimedNaive$logFC<0 & PrimedNaive$FDR<0.05)],L16)


NaiveList2 <- intersect(PrimedNaiveForm$Gene[which(PrimedNaiveForm$vsNaive>0 & PrimedNaiveForm$P.Value<0.05)],L16)
PrimeList2 <- intersect(PrimedNaiveForm$Gene[which(PrimedNaiveForm$vsPrimed<0 & PrimedNaiveForm$P.Value<0.05)],L16)

Idents(D) <- D$Cells
DotPlot(D, features =  NaiveList )
ggsave(filename=paste(saveext,"/DimRed/Naive1_dot.pdf",sep=""),width = 20, height = 20)
DotPlot(D, features =  PrimeList )
ggsave(filename=paste(saveext,"/DimRed/Primed1_dot.pdf",sep=""),width = 20, height = 20)

DotPlot(D, features =  NaiveList2 )
ggsave(filename=paste(saveext,"/DimRed/Naive2_dot.pdf",sep=""),width = 20, height = 20)
DotPlot(D, features =  PrimeList2 )
ggsave(filename=paste(saveext,"/DimRed/Primed2_dot.pdf",sep=""),width = 20, height = 20)



DotPlot(D, features =  AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E (+ hsTE)")] )
ggsave(filename=paste(saveext,"/DimRed/hsAME-E_hsTE_dot.pdf",sep=""),width = 20, height = 20)


DotPlot(D, features =  AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E (+ cyAME-L)")] )
ggsave(filename=paste(saveext,"/DimRed/hsAME-EAMEL_dot.pdf",sep=""),width = 20, height = 20)


DotPlot(D, features =  AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E + cyAME-L")] )
ggsave(filename=paste(saveext,"/DimRed/hsAME-EAMEL__dot.pdf",sep=""),width = 20, height = 20)

DotPlot(D, features =  AmnionMarkers$Gene[which(AmnionMarkers$marker=="cyAME-L")] )
ggsave(filename=paste(saveext,"/DimRed/cyAME-L_dot.pdf",sep=""),width = 20, height = 20)


mat_breaks <- seq(0, 20, length.out = 20)
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L17_Unscaled_EmbDisc_Emb.pdf",sep=""),width=20,height=60)





#First do Am_d14
DEListM1 <- FindMarkers(D, ident.1 = c("C1_Am_d14"), ident.2 = setdiff(List1, "C1_Am_d14"),test.use = "MAST", only.pos = TRUE)
DEListM2 <- FindMarkers(D, ident.1 = c("C2_Am_d14"), ident.2 = setdiff(List2, "C2_Am_d14"),test.use = "MAST", only.pos = TRUE)
#DEListM3 <- FindMarkers(D, ident.1 = c("C3_Am_d14"), ident.2 = setdiff(List3, "C3_Am_d14"),test.use = "MAST", only.pos = TRUE)
#DEListM4 <- FindMarkers(D, ident.1 = c("C4_Am_d14"), ident.2 = setdiff(List4, "C4_Am_d14"),test.use = "MAST", only.pos = TRUE)
DEListM6 <- FindMarkers(D, ident.1 = c("C6_Am_d14"), ident.2 = setdiff(List6, "C6_Am_d14"),test.use = "MAST", only.pos = TRUE)
DEListM7 <- FindMarkers(D, ident.1 = c("C7_Am_d14"), ident.2 = setdiff(List7, "C7_Am_d14"),test.use = "MAST", only.pos = TRUE)
#DEListM8 <- FindMarkers(D, ident.1 = c("C8_Am_d14"), ident.2 = setdiff(List8, "C8_Am_d14"),test.use = "MAST", only.pos = TRUE)

DEListM1 <- DEListM1[which(DEListM1$p_val<0.01),]
DEListM1 <- DEListM1[order(DEListM1$avg_log2FC,decreasing=TRUE),]
DEListM2 <- DEListM2[which(DEListM2$p_val<0.01),]
DEListM2 <- DEListM2[order(DEListM2$avg_log2FC,decreasing=TRUE),]
DEListM6 <- DEListM6[which(DEListM6$p_val<0.01),]
DEListM6 <- DEListM6[order(DEListM6$avg_log2FC,decreasing=TRUE),]
DEListM7 <- DEListM7[which(DEListM7$p_val<0.01),]
DEListM7 <- DEListM7[order(DEListM7$avg_log2FC,decreasing=TRUE),]

L13 <- Reduce(intersect, list(rownames(DEListM1),rownames(DEListM2),rownames(DEListM6),rownames(DEListM7)))
L13 <- unique( c(
  intersect(rownames(DEListM1),rownames(DEListM2)),
  intersect(rownames(DEListM1),rownames(DEListM6)),
  intersect(rownames(DEListM1),rownames(DEListM7)),
  intersect(rownames(DEListM2),rownames(DEListM6)),
  intersect(rownames(DEListM2),rownames(DEListM7)),
  intersect(rownames(DEListM6),rownames(DEListM7))
) )
L13 <- unique( c(rownames(DEListM1),rownames(DEListM2),rownames(DEListM6),rownames(DEListM7)) )



Dat <- Av[L13,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L13_Am.pdf",sep=""),scale="row",width=20,height=100)


mat_breaks <- seq(0, 20, length.out = 20)
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L13_Unscaled_Am.pdf",sep=""),width=20,height=100)





#First do Am_d14
#First do lumenal
DEListN1 <- FindMarkers(D, ident.1 = c("C1_Lumenal"), ident.2 = c("C1_Prolif","C1_Ciliated"),test.use = "MAST", only.pos = TRUE)
DEListN2 <- FindMarkers(D, ident.1 = c("C2_Lumenal"), ident.2 = c("C2_Prolif","C2_Ciliated"),test.use = "MAST", only.pos = TRUE)
DEListN3 <- FindMarkers(D, ident.1 = c("C3_Lumenal"), ident.2 = c("C3_Prolif","C3_Ciliated"),test.use = "MAST", only.pos = TRUE)
DEListN4 <- FindMarkers(D, ident.1 = c("C4_Lumenal"), ident.2 = c("C4_Prolif","C4_Ciliated"),test.use = "MAST", only.pos = TRUE)
DEListN6 <- FindMarkers(D, ident.1 = c("C6_Lumenal"), ident.2 = c("C6_Prolif","C6_Ciliated"),test.use = "MAST", only.pos = TRUE)
DEListN7 <- FindMarkers(D, ident.1 = c("C7_Lumenal"), ident.2 = c("C7_Prolif","C7_Ciliated"),test.use = "MAST", only.pos = TRUE)
DEListN8 <- FindMarkers(D, ident.1 = c("C8_Lumenal"), ident.2 = c("C8_Prolif","C8_Ciliated"),test.use = "MAST", only.pos = TRUE)

DEListN1 <- DEListN1[which(DEListN1$p_val<0.01),]
DEListN1 <- DEListN1[order(DEListN1$avg_log2FC,decreasing=TRUE),]
DEListN2 <- DEListN2[which(DEListN2$p_val<0.01),]
DEListN2 <- DEListN2[order(DEListN2$avg_log2FC,decreasing=TRUE),]
DEListN6 <- DEListN6[which(DEListN6$p_val<0.01),]
DEListN6 <- DEListN6[order(DEListN6$avg_log2FC,decreasing=TRUE),]
DEListN7 <- DEListN7[which(DEListN7$p_val<0.01),]
DEListN7 <- DEListN7[order(DEListN7$avg_log2FC,decreasing=TRUE),]

L20 <- Reduce(intersect, list(rownames(DEListN1),rownames(DEListN2),rownames(DEListN6),rownames(DEListN7)))
L20 <- unique( c(
  intersect(rownames(DEListN1),rownames(DEListN2)),
  intersect(rownames(DEListN1),rownames(DEListN6)),
  intersect(rownames(DEListN1),rownames(DEListN7)),
  intersect(rownames(DEListN2),rownames(DEListN6)),
  intersect(rownames(DEListN2),rownames(DEListN7)),
  intersect(rownames(DEListN6),rownames(DEListN7))
) )
L20 <- unique( c(rownames(DEListN1),rownames(DEListN2),rownames(DEListN6),rownames(DEListN7)) )



Dat <- Av[L20,exps]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L20_LumenalvsProlfCil.pdf",sep=""),scale="row",width=20,height=100)


mat_breaks <- seq(0, 20, length.out = 20)
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(4,11,17,17,17,23,23,23,26,26,26,41,41,41,88,88,88), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L20_Unscaled_LumenalvsProlfCil.pdf",sep=""),width=20,height=100)




expssub <- c("C1_Lumenal",
          "C2_Lumenal",
          "C3_Lumenal",
          "C4_Lumenal",
          "C6_Lumenal",
          "C7_Lumenal",
          "C8_Lumenal",
          "C1_Prolif",
          "C2_Prolif",
          "C3_Prolif",
          "C4_Prolif",
          "C6_Prolif",
          "C7_Prolif",
          "C8_Prolif",
          "C1_SOX9P",
          "C2_SOX9P",
          "C3_SOX9P",
          "C4_SOX9P",
          "C6_SOX9P",
          "C7_SOX9P",
          "C8_SOX9P",
          "C1_SOX9LRG5",
          "C2_SOX9LRG5",
          "C4_SOX9LRG5",
          "C6_SOX9LRG5",
          "C8_SOX9LRG5",
          "C1_Glandular",
          "C2_Glandular",
          "C3_Glandular",
          "C4_Glandular",
          "C6_Glandular",
          "C7_Glandular",
          "C8_Glandular",
          "C1_Ciliated",
          "C2_Ciliated",
          "C3_Ciliated",
          "C4_Ciliated",
          "C6_Ciliated",
          "C7_Ciliated",
          "C8_Ciliated",
          "C1_Stromal fibroblasts",
          "C2_Stromal fibroblasts",
          "C3_Stromal fibroblasts",
          "C4_Stromal fibroblasts",
          "C6_Stromal fibroblasts",
          "C7_Stromal fibroblasts",
          "C8_Stromal fibroblasts")

L21 <- unique(c(L1,L2,L3,L4,L5,L7,L8,L20))


Dat <- Av[L21,expssub]
mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(7,14,21,26,33,40), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L21_AllMat.pdf",sep=""),scale="row",width=15,height=400)


mat_breaks <- seq(0, 20, length.out = 20)
pheatmap((Dat),color =  redblue1(20), breaks = mat_breaks, border_color = NA, gaps_col=c(7,14,21,26,33,40), cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_L21_Unscaled_AllMat.pdf",sep=""),width=15,height=400)



############

Dsubset1 <- subset(D,idents = c(          "C1_Am_d14",
                                          "C2_Am_d14" ,
                                          "C6_Am_d14",
                                          "C7_Am_d14",
                                          "C1_Am/EmDisc_d14",
                                          "C2_Am/EmDisc_d14",
                                          "C4_Am/EmDisc_d14",
                                          "C6_Am/EmDisc_d14",
                                          "C7_Am/EmDisc_d14",
                                          "C1_EmDisc_d14",
                                          "C2_EmDisc_d14",
                                          "C4_EmDisc_d14",
                                          "C6_EmDisc_d14",
                                          "C7_EmDisc_d14",
                                          "C8_EmDisc_d14"))

DefaultAssay(Dsubset1) <- "RNA"

p <- FeatureScatter(Dsubset1, feature1 = "POU5F1", feature2 = "WNT6")
ggsave(filename=paste(saveext,"/DimRed/Scatter_POU5F1_WNT6.pdf",sep=""),width = 10, height = 10, useDingbats = FALSE, limitsize = FALSE, p)


p <- FeatureScatter(Dsubset1, feature1 = "POU5F1", feature2 = "PDGFA")
ggsave(filename=paste(saveext,"/DimRed/Scatter_POU5F1_PDGFA.pdf",sep=""),width = 10, height = 10, useDingbats = FALSE, limitsize = FALSE, p)



p <- FeatureScatter(Dsubset1, feature1 = "POU5F1", feature2 = "SOX2")
ggsave(filename=paste(saveext,"/DimRed/Scatter_POU5F1_SOX2.pdf",sep=""),width = 10, height = 10, useDingbats = FALSE, limitsize = FALSE, p)


p <- FeatureScatter(Dsubset1, feature1 = "POU5F1", feature2 = "ISL1")
ggsave(filename=paste(saveext,"/DimRed/Scatter_POU5F1_ISL1.pdf",sep=""),width = 10, height = 10, useDingbats = FALSE, limitsize = FALSE, p)



p <- FeatureScatter(Dsubset1, feature1 = "WNT6", feature2 = "PDGFA")
ggsave(filename=paste(saveext,"/DimRed/Scatter_WNT6_PDGFA.pdf",sep=""),width = 10, height = 10, useDingbats = FALSE, limitsize = FALSE, p)




p <- FeatureScatter(Dsubset1, feature1 = "WNT6", feature2 = "VTCN1")
ggsave(filename=paste(saveext,"/DimRed/Scatter_WNT6_VTCN1.pdf",sep=""),width = 10, height = 10, useDingbats = FALSE, limitsize = FALSE, p)



DefaultAssay(Dsubset1) <- "integrated"
#Dsubset1 <- FindVariableFeatures(Dsubset1, nfeatures = 5000) 
#Dsubset1 <- ScaleData(Dsubset1)
Dsubset1 <- RunPCA(Dsubset1, npcs = 20, verbose = FALSE)
Dsubset1 <- RunUMAP(Dsubset1, reduction = "pca", dims = 1:20)
Dsubset1 <- FindNeighbors(Dsubset1, reduction = "pca", dims = 1:2)

#p<- DimPlot(Dsubset1, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/Amsubset_RNAslot_UMAP_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE, p)
#p<-DimPlot(Dsubset1, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/Amsubset_RNA_slot_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)


#Idents(MesotheliaData) <- MesotheliaData$ID5
p<- DimPlot(Dsubset1, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Amsubset_UMAP_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE, p)
p<-DimPlot(Dsubset1, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Amsubset_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)


p<- DimPlot(Dsubset1, pt.size = 4, reduction = "umap", split.by = "ID3",label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Amsubset_UMAP_Allsplit.pdf",sep=""),width = 62, height = 10, useDingbats = FALSE, limitsize = FALSE, p)
p<-DimPlot(Dsubset1, pt.size = 4, reduction = "pca",  split.by = "ID3",label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Amsubset_PCA_Allsplit.pdf",sep=""),width = 62, height = 10, useDingbats = FALSE, limitsize = FALSE,p)


Idents(Dsubset1) <- Dsubset1$Genotype
p<- DimPlot(Dsubset1, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/AmsubsetGT_UMAP_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE, p)
p<-DimPlot(Dsubset1, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/AmsubsetGT_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)


DefaultAssay(Dsubset1) <- "RNA"

FeaturePlot(Dsubset1,  reduction = "umap", features = "WNT6", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset_WNT6.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Dsubset1,  reduction = "umap", features = "VTCN1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset_VTCN1.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Dsubset1,  reduction = "umap", features = "POU5F1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset_POU5F1.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Dsubset1,  reduction = "umap", features = "PDGFA", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset_PDGFA.pdf",sep=""),width = 12, height = 10)


DefaultAssay(Dsubset1) <- "RNA"

FeaturePlot(Dsubset1,  reduction = "umap", split.by = "ID3",features = "WNT6", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset_WNT6split.pdf",sep=""),width = 62, height = 10,limitsize = FALSE)
FeaturePlot(Dsubset1,  reduction = "umap", split.by = "ID3",features = "VTCN1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset_VTCN1split.pdf",sep=""),width = 62, height = 10,limitsize = FALSE)
FeaturePlot(Dsubset1,  reduction = "umap", split.by = "ID3",features = "POU5F1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset_POU5F1split.pdf",sep=""),width = 62, height = 10,limitsize = FALSE)
FeaturePlot(Dsubset1,  reduction = "umap", split.by = "ID3",features = "PDGFA", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset_PDGFAsplit.pdf",sep=""),width = 62, height = 10,limitsize = FALSE)


DefaultAssay(Dsubset1) <- "integrated"

Dsubset1$ExtID <- Idents(Dsubset1)
Dsubset1 <- FindClusters(Dsubset1, resolution = 1)
Dsubset1$Cluster <- Idents(Dsubset1)

Dsubset1 <- FindClusters(Dsubset1, resolution = 2)
Dsubset1$FineCluster <- Idents(Dsubset1)


p<-DimPlot(Dsubset1, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/AmsubsetCl_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-DimPlot(Dsubset1, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/AmsubsetCl_UMAP_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)

DefaultAssay(Dsubset1) <- "RNA"
AllMarkers <- FindAllMarkers(Dsubset1, only.pos = TRUE, test.use = "MAST")
write.table(AllMarkers, file = paste(saveext,"/EmDiscMarkers.csv",sep=""), sep="," )
DefaultAssay(Dsubset1) <- "RNA"
AllMarkers <- FindMarkers(Dsubset1, ident.1 = c(1,2), ident.2 = c(0,3,5), only.pos = TRUE, test.use = "MAST")
write.table(AllMarkers, file = paste(saveext,"/EmDiscMarkers_subset.csv",sep=""), sep="," )



NewG1 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D1_emb.tsv", header=TRUE)
NewG2 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D2_emb.tsv", header=TRUE)
NewG3B <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D3B_emb.tsv", header=TRUE)
NewG6 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D6_emb.tsv", header=TRUE)
NewG8 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D8_emb.tsv", header=TRUE)

Idents(Dsubset1) <- Dsubset1$ID3
Idents(Dsubset1, cells= paste(NewG1$barcode[which(NewG1$assignment=="0")],2,sep="_") ) <- "NC1_1"
Idents(Dsubset1, cells= paste(NewG1$barcode[which(NewG1$assignment=="1")],2,sep="_") ) <- "NC1_2"
Idents(Dsubset1, cells= paste(NewG1$barcode[which(NewG1$assignment=="0/1")],2,sep="_") ) <- "NC1_1"
Idents(Dsubset1, cells= paste(NewG1$barcode[which(NewG1$assignment=="1/0")],2,sep="_") ) <- "NC1_2"

Idents(Dsubset1, cells= paste(NewG2$barcode[which(NewG2$assignment=="0")],3,sep="_") ) <- "NC2_1"
Idents(Dsubset1, cells= paste(NewG2$barcode[which(NewG2$assignment=="1")],3,sep="_") ) <- "NC2_2"
Idents(Dsubset1, cells= paste(NewG2$barcode[which(NewG2$assignment=="0/1")],3,sep="_") ) <- "NC2_1"
Idents(Dsubset1, cells= paste(NewG2$barcode[which(NewG2$assignment=="1/0")],3,sep="_") ) <- "NC2_2"

Idents(Dsubset1, cells= paste(NewG6$barcode[which(NewG6$assignment=="0")],6,sep="_") ) <- "NC6_1"
Idents(Dsubset1, cells= paste(NewG6$barcode[which(NewG6$assignment=="1")],6,sep="_") ) <- "NC6_2"
#Idents(Dsubset1, cells= paste(NewG6$barcode[which(NewG6$assignment=="0/1")],6,sep="_") ) <- "NC6_1"
Idents(Dsubset1, cells= paste(NewG6$barcode[which(NewG6$assignment=="1/0")],6,sep="_") ) <- "NC6_2"

Idents(Dsubset1, cells= paste(NewG8$barcode[which(NewG8$assignment=="0")],8,sep="_") ) <- "NC2_1"
#Idents(Dsubset1, cells= paste(NewG8$barcode[which(NewG8$assignment=="1")],8,sep="_") ) <- "NC2_2"
#Idents(Dsubset1, cells= paste(NewG8$barcode[which(NewG8$assignment=="0/1")],8,sep="_") ) <- "NC2_1"
#Idents(Dsubset1, cells= paste(NewG8$barcode[which(NewG8$assignment=="1/0")],8,sep="_") ) <- "NC2_2"

#Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0")],5,sep="_") ) <- "NC3_1"
#Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/2")],5,sep="_") ) <- "NC3_1"
Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/1")],5,sep="_") ) <- "NC3_1"
Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1")],5,sep="_") ) <- "NC3_2"
#Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1/2")],5,sep="_") ) <- "NC3_2"
#Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1/0")],5,sep="_") ) <- "NC3_2"
Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="2")],5,sep="_") ) <- "NC3_3"
Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="2/0")],5,sep="_") ) <- "NC3_3"
Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/1")],5,sep="_") ) <- "NC2_1"
#Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1/0")],5,sep="_") ) <- "NC2_2"


Dsubset1$Genotype2 <- Idents(Dsubset1)
Idents(Dsubset1) <- Dsubset1$Cells

#p<- DimPlot(Dsubset1, pt.size = 4, reduction = "umap", label = TRUE, split.by = "Genotype2", repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/Amsubset_RNAslot_UMAP_Allsplit.pdf",sep=""),width = 72, height = 10, useDingbats = FALSE, limitsize = FALSE, p)
#p<-DimPlot(Dsubset1, pt.size = 4, reduction = "pca",  split.by = "Genotype2", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/Amsubset_RNA_slot_PCA_Allsplit.pdf",sep=""),width = 72, height = 10, useDingbats = FALSE, limitsize = FALSE,p)


DefaultAssay(Dsubset1) <- "integrated"
Dsubset1 <- FindClusters(Dsubset1, resolution = 0.5)
p<-DimPlot(Dsubset1, pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/AmsubsetCl_PCA_All_more.pdf",sep=""),width = 62, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-DimPlot(Dsubset1, pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/AmsubsetCl_UMAP_All_more.pdf",sep=""),width = 62, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
Dsubset1$Cl <- Idents(Dsubset1)
Embryo1_1 <- subset(Dsubset1,cells= names(which(Dsubset1$ID3=="C1" & Dsubset1$Cl==0)) )
Embryo1_2 <- subset(Dsubset1,cells= names(which(Dsubset1$ID3=="C1" & Dsubset1$Cl%in%c(1,2) )) )
Embryo2_2 <- subset(Dsubset1,cells= names(which(Dsubset1$ID3=="C2" & Dsubset1$Cl==0)) )
Embryo2_2 <- subset(Dsubset1,cells= names(which(Dsubset1$ID3=="C2" & Dsubset1$Cl%in%c(1,2) )) )
Embryo4_1 <- subset(Dsubset1,cells= names(which(Dsubset1$ID3=="C4" & Dsubset1$Cl==0)) )
Embryo4_2 <- subset(Dsubset1,cells= names(which(Dsubset1$ID3=="C4" & Dsubset1$Cl%in%c(1,2) )) )
Embryo6_1 <- subset(Dsubset1,cells= names(which(Dsubset1$ID3=="C6" & Dsubset1$Cl==0)) )
Embryo6_2 <- subset(Dsubset1,cells= names(which(Dsubset1$ID3=="C6" & Dsubset1$Cl%in%c(1,2) )) )
Embryo7_1 <- subset(Dsubset1,cells= names(which(Dsubset1$ID3=="C7" & Dsubset1$Cl==0)) )
Embryo7_2 <- subset(Dsubset1,cells= names(which(Dsubset1$ID3=="C7" & Dsubset1$Cl%in%c(1,2) )) )
Embryo8_1 <- subset(Dsubset1,cells= names(which(Dsubset1$ID3=="C8" & Dsubset1$Cl==0)) )
Embryo8_2 <- subset(Dsubset1,cells= names(which(Dsubset1$ID3=="C8" & Dsubset1$Cl%in%c(1,2) )) )

DefaultAssay(Embryo1_2) <- "RNA"
Embryo1_2 <- FindVariableFeatures(Embryo1_2, nfeatures = 5000) 
Embryo1_2 <- ScaleData(Embryo1_2)
Embryo1_2 <- RunPCA(Embryo1_2, npcs = 15, verbose = FALSE)
#Embryo1_2 <- RunUMAP(Embryo1_2, reduction = "pca", dims = 1:20)
Embryo1_2 <- FindNeighbors(Embryo1_2, reduction = "pca", dims = 1:2)
p<- DimPlot(Embryo1_2, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Embryo1_2_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE, p)
FeaturePlot(Embryo1_2,  reduction = "pca", features = "WNT6", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryo1_2_WNT6.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Embryo1_2,  reduction = "pca", features = "VTCN1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryo1_2_VTCN1.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Embryo1_2,  reduction = "pca", features = "POU5F1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryo1_2_POU5F1.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Embryo1_2,  reduction = "pca", features = "PDGFA", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryo1_2_PDGFA.pdf",sep=""),width = 12, height = 10)


DefaultAssay(Embryo1_1) <- "RNA"
Embryo1_1 <- FindVariableFeatures(Embryo1_1, nfeatures = 5000) 
Embryo1_1 <- ScaleData(Embryo1_1)
Embryo1_1 <- RunPCA(Embryo1_1, npcs = 5, verbose = FALSE)
#Embryo1_2 <- RunUMAP(Embryo1_2, reduction = "pca", dims = 1:20)
Embryo1_2 <- FindNeighbors(Embryo1_1, reduction = "pca", dims = 1:2)
p<- DimPlot(Embryo1_1, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Embryo1_1_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE, p)
FeaturePlot(Embryo1_1,  reduction = "pca", features = "WNT6", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryo1_1_WNT6.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Embryo1_1,  reduction = "pca", features = "VTCN1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryo1_1_VTCN1.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Embryo1_1,  reduction = "pca", features = "POU5F1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryo1_1_POU5F1.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Embryo1_1,  reduction = "pca", features = "PDGFA", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryo1_1_PDGFA.pdf",sep=""),width = 12, height = 10)


EmbryoSubset1 <- subset(Dsubset1,idents=0)
EmbryoSubset2 <- subset(Dsubset1,idents=c(1,2))
DefaultAssay(EmbryoSubset1) <- "RNA"
EmbryoSubset1 <- FindVariableFeatures(EmbryoSubset1, nfeatures = 5000) 
EmbryoSubset1 <- ScaleData(EmbryoSubset1)
EmbryoSubset1 <- RunPCA(EmbryoSubset1, npcs = 5, verbose = FALSE)
#Embryo1_2 <- RunUMAP(Embryo1_2, reduction = "pca", dims = 1:20)
EmbryoSubset1 <- FindNeighbors(EmbryoSubset1, reduction = "pca", dims = 1:2)
p<- DimPlot(EmbryoSubset1, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Embryo1_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE, p)
FeaturePlot(EmbryoSubset1,  reduction = "pca", features = "WNT6", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryo1_WNT6.pdf",sep=""),width = 12, height = 10)
FeaturePlot(EmbryoSubset1,  reduction = "pca", features = "VTCN1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryo1_VTCN1.pdf",sep=""),width = 12, height = 10)
FeaturePlot(EmbryoSubset1,  reduction = "pca", features = "POU5F1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryo1_POU5F1.pdf",sep=""),width = 12, height = 10)
FeaturePlot(EmbryoSubset1,  reduction = "pca", features = "PDGFA", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryo1_PDGFA.pdf",sep=""),width = 12, height = 10)


DefaultAssay(EmbryoSubset2) <- "RNA"
EmbryoSubset2 <- FindVariableFeatures(EmbryoSubset2, nfeatures = 5000) 
EmbryoSubset2 <- ScaleData(EmbryoSubset2)
EmbryoSubset2 <- RunPCA(EmbryoSubset2, npcs = 5, verbose = FALSE)
#Embryo1_2 <- RunUMAP(Embryo1_2, reduction = "pca", dims = 1:20)
EmbryoSubset2 <- FindNeighbors(EmbryoSubset2, reduction = "pca", dims = 1:2)
p<- DimPlot(EmbryoSubset2, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Embryo2_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE, p)
FeaturePlot(EmbryoSubset2,  reduction = "pca", features = "WNT6", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryo2_WNT6.pdf",sep=""),width = 12, height = 10)
FeaturePlot(EmbryoSubset2,  reduction = "pca", features = "VTCN1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryo2_VTCN1.pdf",sep=""),width = 12, height = 10)
FeaturePlot(EmbryoSubset2,  reduction = "pca", features = "POU5F1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryo2_POU5F1.pdf",sep=""),width = 12, height = 10)
FeaturePlot(EmbryoSubset2,  reduction = "pca", features = "PDGFA", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Embryo2_PDGFA.pdf",sep=""),width = 12, height = 10)


DE1 <- FindMarkers(EmbryoSubset2,ident.1 = "1", ident.2 = "2", test.use = "MAST", only.pos = FALSE)
write.table(DE1, file = paste(saveext,"/EmDiscMarkers_cluster1_vs_2.csv",sep=""), sep="," )


DefaultAssay(Dsubset1) <- "RNA"
DE2 <- FindAllMarkers(Dsubset1, test.use = "MAST", only.pos = FALSE)
write.table(DE2, file = paste(saveext,"/EmDiscMarkers_cluster1_vs_2_vs_3.csv",sep=""), sep="," )

#p<-DimPlot(Dsubset2, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/Amsubset2_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)





Idents(Dsubset1) <- Dsubset1$Cluster
p<-DimPlot(Dsubset1, pt.size = 4, reduction = "pca", split.by = "Genotype2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/AmsubsetCluster_PCA_All_more.pdf",sep=""),width = 72, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-DimPlot(Dsubset1, pt.size = 4, reduction = "umap", split.by = "Genotype2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/AmsubsetCluster_UMAP_All_more.pdf",sep=""),width = 72, height = 10, useDingbats = FALSE, limitsize = FALSE,p)


DefaultAssay(Dsubset1) <- "RNA"

FeaturePlot(Dsubset1,  reduction = "umap", split.by = "Genotype2",features = "WNT6", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset_WNT6splitmore.pdf",sep=""),width = 72, height = 10,limitsize = FALSE)
FeaturePlot(Dsubset1,  reduction = "umap", split.by = "Genotype2",features = "VTCN1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset_VTCN1splitmore.pdf",sep=""),width = 72, height = 10,limitsize = FALSE)
FeaturePlot(Dsubset1,  reduction = "umap", split.by = "Genotype2",features = "POU5F1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset_POU5F1splitmore.pdf",sep=""),width = 72, height = 10,limitsize = FALSE)
FeaturePlot(Dsubset1,  reduction = "umap", split.by = "Genotype2",features = "PDGFA", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset_PDGFAsplitmore.pdf",sep=""),width = 72, height = 10,limitsize = FALSE)


uID <- as.character(Dsubset1$Cells)
uID[which(Dsubset1$Genotype2=="NC2_1" & Dsubset1$Cluster%in%c(0))] <- "EmDisc"
uID[which(Dsubset1$Genotype2=="NC2_1" & Dsubset1$Cluster%in%c(1,2,3,4))] <- "Am"

uID[which(Dsubset1$Genotype2=="NC3_3" & Dsubset1$Cluster%in%c(0))] <- "EmDisc"
uID[which(Dsubset1$Genotype2=="NC3_3" & Dsubset1$Cluster%in%c(1,2,3,4))] <- "Am"

uID[which(Dsubset1$Genotype2=="NC3_2" & Dsubset1$Cluster%in%c(0))] <- "EmDisc"
uID[which(Dsubset1$Genotype2=="NC3_2" & Dsubset1$Cluster%in%c(1,2,3,4))] <- "Am"


uID[which(Dsubset1$Genotype2=="NC6_2" & Dsubset1$Cluster%in%c(0))] <- "EmDisc"
uID[which(Dsubset1$Genotype2=="NC6_2" & Dsubset1$Cluster%in%c(1,2,3,4))] <- "Am"


#Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0")],4,sep="_") ) <- "NC4_1"
#Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1")],4,sep="_") ) <- "NC4_2"
#Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/1")],4,sep="_") ) <- "NC4_Mix"

#Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0")],4,sep="_") ) <- "NC4_1"
#Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1")],4,sep="_") ) <- "NC4_2"
#Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/1")],4,sep="_") ) <- "NC4_Mix"



#Dsubset2 <- subset(Dsubset1,idents = c(0,3,5))
#DefaultAssay(Dsubset2) <- "integrated"

#Dsubset2 <- RunPCA(Dsubset2, npcs = 20, verbose = FALSE)
#Dsubset2 <- RunUMAP(Dsubset2, reduction = "pca", dims = 1:20)
#Dsubset2 <- FindNeighbors(Dsubset2, reduction = "pca", dims = 1:2)

#p<- DimPlot(Dsubset2, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/Amsubset2_UMAP_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE, p)
#p<-DimPlot(Dsubset2, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/Amsubset2_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)


#DefaultAssay(Dsubset2) <- "RNA"
#FeaturePlot(Dsubset2,  reduction = "umap", features = "WNT6", cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/Dsubset2_WNT6.pdf",sep=""),width = 12, height = 10)
#FeaturePlot(Dsubset2,  reduction = "umap", features = "VTCN1", cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/Dsubset2_VTCN1.pdf",sep=""),width = 12, height = 10)
#FeaturePlot(Dsubset2,  reduction = "umap", features = "POU5F1", cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/Dsubset2_POU5F1.pdf",sep=""),width = 12, height = 10)
#FeaturePlot(Dsubset2,  reduction = "umap", features = "PDGFA", cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/Dsubset2_PDGFA.pdf",sep=""),width = 12, height = 10)

#DefaultAssay(Dsubset2) <- "integrated"
#Dsubset2 <- FindClusters(Dsubset2, resolution = 1)


#p<- DimPlot(Dsubset2, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/Amsubset2Cl_UMAP_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE, p)
#p<-DimPlot(Dsubset2, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/Amsubset2Cl_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)

#AllMarkers <- FindAllMarkers(Dsubset2, only.pos = TRUE, test.use = "MAST")
#write.table(AllMarkers, file = paste(saveext,"/EmDiscMarkers_subset2.csv",sep=""), sep="," )

#DefaultAssay(Dsubset2) <- "RNA"
#FeaturePlot(Dsubset2,  reduction = "pca", features = "WNT6", cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/Dsubset2_PCA_WNT6.pdf",sep=""),width = 12, height = 10)
#FeaturePlot(Dsubset2,  reduction = "pca", features = "VTCN1", cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/Dsubset2_PCA_VTCN1.pdf",sep=""),width = 12, height = 10)
#FeaturePlot(Dsubset2,  reduction = "pca", features = "POU5F1", cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/Dsubset2_PCA_POU5F1.pdf",sep=""),width = 12, height = 10)
#FeaturePlot(Dsubset2,  reduction = "pca", features = "PDGFA", cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/Dsubset2_PCA_PDGFA.pdf",sep=""),width = 12, height = 10)

#AllMarkers <- FindMarkers(Dsubset2, ident.1 = c(0,3), ident.2 = c(1,2), only.pos = TRUE, test.use = "MAST")
#write.table(AllMarkers, file = paste(saveext,"/EmDiscMarkers_subset_EmD.csv",sep=""), sep="," )


#AllMarkers <- FindMarkers(Dsubset2, ident.2 = c(0,3), ident.1 = c(1,2), only.pos = TRUE, test.use = "MAST")
#write.table(AllMarkers, file = paste(saveext,"/EmDiscMarkers_subset_Am.csv",sep=""), sep="," )

#Idents(Dsubset2,WhichCells(Dsubset2,idents=c(0,3))) <- "EmDisc"
#Idents(Dsubset2,WhichCells(Dsubset2,idents=c(1,2))) <- "Am"


feats <- c("RBFOX2",
           "SMARCA1",
           "FBXO15",
           "SOX15",
           "LIFR",
           "MID1",
           "PML",
           "FRAS1",
           "CD9",
           "AMOTL2",
           "IL6ST",
           "EPHA2",
           "AJUBA",
           "ROR1",
           "L1TD1",
           "ABCG2",
           "RAI14",
           "LEF1",
           "YAP1",
           "SLC46A2",
           "CTNND1",
           "PRDM5",
           "KHDC3L",
           "LRP6",
           "LRP4",
           "FBLN1",
           "FNBP1L",
           "SFRP1",
           "UACA",
           "LMNA",
           "PRTG",
           "ANTXR1",
           "LAPTM4B",
           "ARHGAP42",
           "NACC1",
           "ZSCAN10",
           "ADGRA3",
           "SH3BP4",
           "CDH1",
           "MYC",
           "HES1",
           "PITX2",
           "YES1",
           "SMAD3",
           "ESRRB",
           "FZD1",
           "SMAD2",
           "SMAD5",
           "SMAD4",
           "FZD5",
           "SMAD1",
           "FZD7",
           "RBPMS",
           "GDF3",
           "PARD3",
           "PTK7",
           "PUM2",
           "IL17RD",
           "TEX19",
           "SMAD9",
           "PTK2",
           "GJB1",
           "FERMT2",
           "ALDH7A1",
           "ZFP42",
           "GJB4",
           "FAT1",
           "HOXB5",
           "NANOG",
           "TDGF1",
           "CD276",
           "WWC2",
           "LYPD6",
           "FGF4",
           "SOX2",
           "GJA1",
           "PODXL",
           "SALL4",
           "TSPAN6",
           "SERPINH1",
           "PROM1",
           "PLS3",
           "FOXD3",
           "ITGA4",
           "KCNIP3",
           "DNMT3B",
           "SHROOM3",
           "HMGA2",
           "GPX8",
           "NAV2",
           "LAMB1",
           "PTPN14",
           "GNG12",
           "UTF1",
           "GULP1",
           "PAWR",
           "DCBLD2",
           "ENAH",
           "PTPRF",
           "CDK8",
           "KITLG",
           "MFAP2",
           "DPPA5",
           "ALPP",
           "DPPA4",
           "DPPA3",
           "DPPA2",
           "BCL3",
           "PECAM1",
           "DSG2",
           "ADGRL2",
           "NES",
           "ITGA6",
           "MDFI",
           "CD24",
           "DOCK1",
           "TENM3",
           "FBLIM1",
           "LAMC1",
           "PTPN21",
           "THY1",
           "PTPRG",
           "FUT4",
           "GJC1",
           "EPCAM",
           "SH3PXD2B",
           "TRIM6",
           "ERBB2",
           "TEAD1",
           "TEAD2",
           "TEAD3",
           "TEAD4",
           "MAGI1",
           "DSP",
           "STAT3",
           "TPBG",
           "MYO10",
           "ZFX",
           "KLF4",
           "ANO6",
           "POU5F1",
           "NR6A1",
           "GAL",
           "TJP1",
           "CLDN6",
           "COL4A1",
           "KIT",
           "MYO1B",
           "COL4A2",
           "COL4A5",
           "SPRY4",
           "COL4A6",
           "TAF8",
           "LGR4",
           "ITGB1",
           "PIWIL2",
           "NXN",
           "FKBP10",
           "PIWIL4",
           "NCKAP1",
           "CTNNB1",
           "FSTL1",
           "GLI2",
           "RND3",
           "HHEX",
           "SUMO2",
           "ZIC1",
           "BNIP3",
           "PCGF2",
           "CTNNA1",
           "TRIM28",
           "LRIG3",
           "CD59",
           "SH3D19",
           "PIWIL1",
           "EPHB4")

feats2 <- c("KRT1",
            "KRT5",
            "IL1R1",
            "KNL1",
            "KRT3",
            "KRT2",
            "KRT8",
            "KRT7",
            "DKK3",
            "SFRP5",
            "CD2",
            "SCGB1A1",
            "SPINT1",
            "SLPI",
            "KIF20B",
            "SLC46A2",
            "EHF",
            "CBLC",
            "MST1R",
            "AQP3",
            "CLDN1",
            "BAIAP2L1",
            "MUC1",
            "IFI16",
            "S100A14",
            "PRSS8",
            "MAL2",
            "POF1B",
            "TMEM184A",
            "PAQR5",
            "FZD6",
            "C9ORF152",
            "ST14",
            "ST6GAL1",
            "PLA2G4A",
            "SLC6A14",
            "MYO5B",
            "PRG4",
            "OVOL2",
            "OVOL1",
            "FERMT1",
            "IL22RA1",
            "OCLN",
            "PLEKHS1",
            "GJB3",
            "GJB4",
            "SPRR1B",
            "GRB7",
            "GATA2",
            "ICAM1",
            "LYPD8",
            "LGALS4",
            "C19ORF33",
            "ANPEP",
            "SI",
            "TSPAN8",
            "LIPH",
            "SFN",
            "TSPAN1",
            "PROM1",
            "TNS4",
            "TP63",
            "FGFBP1",
            "RIPK4",
            "CCL20",
            "KRT14",
            "KRT13",
            "KLK10",
            "TMC4",
            "ANLN",
            "CEACAM1",
            "KRT19",
            "INCENP",
            "CEACAM7",
            "FXYD3",
            "KRT16",
            "CEACAM6",
            "CEACAM5",
            "KRT15",
            "CD24",
            "SLC44A4",
            "KLK1",
            "TMPRSS4",
            "KIF15",
            "TMPRSS2",
            "BICDL2",
            "BRCA1",
            "KLK3",
            "EVPL",
            "BMI1",
            "DUOX2",
            "EPCAM",
            "MECOM",
            "STAP2",
            "CTSE",
            "MPZL2",
            "TM4SF20",
            "AP1M2",
            "SERPINB5",
            "DGAT2",
            "CDKN2A",
            "FMO5",
            "F5",
            "GUCA2A",
            "RNF128",
            "B3GNT3",
            "HBEGF",
            "CDKL1",
            "IVL",
            "CRB3",
            "NRP2",
            "FOXA1",
            "LAD1",
            "PIGR",
            "CRYBA1",
            "PRSS22",
            "RBBP8NL",
            "AGR3",
            "GPA33",
            "AGR2",
            "FAM3D",
            "GABRP",
            "GRHL2",
            "CXCL10",
            "SCNN1G",
            "C1ORF116",
            "SCNN1D",
            "SCNN1B",
            "PKP3",
            "CHMP4C",
            "LY6D",
            "PRR15L",
            "RAI14",
            "SCNN1A",
            "ATP2C2",
            "NKD1",
            "DEFB1",
            "C6ORF132",
            "CXCL17",
            "SCEL",
            "KRT6A",
            "TMEM125",
            "CDH3",
            "DUOXA2",
            "CDH1",
            "RAB25",
            "MSX1",
            "SPTSSB",
            "MISP",
            "SLC25A48",
            "KRT6B",
            "GGT6",
            "PHGR1",
            "CKAP2",
            "KIF23",
            "SYT8",
            "ZFP42",
            "SYCN",
            "LCN2",
            "PGR",
            "COL17A1",
            "C1ORF210",
            "TACSTD2",
            "AGER",
            "MUC13",
            "MUC16",
            "GGT1",
            "ITGA1",
            "GPX2",
            "ITGA5",
            "LAMB3",
            "ITGA4",
            "ITGA2",
            "MUC5AC",
            "EPN3",
            "ELF3",
            "ESRP2",
            "APOA1",
            "ESRP1",
            "IRF6",
            "THRSP",
            "CXCL8",
            "DLX5",
            "ALOX15",
            "ITPR2",
            "HMMR",
            "LAMC2",
            "FUT3",
            "EMILIN2",
            "ECT2",
            "BOK",
            "FAM83E",
            "SDR16C5",
            "FAM83B",
            "KRT72",
            "CKAP2L",
            "PLEKHG6",
            "TSTD1",
            "PRC1",
            "KLF5",
            "SELENBP1",
            "CLDN7",
            "P2RX7",
            "CLDN4",
            "LTF",
            "NECTIN4",
            "SEC23B",
            "TRIM31",
            "ITGB6",
            "PSCA",
            "ITGB4",
            "ITGAL",
            "TTR",
            "HIRIP3",
            "EXO1",
            "TRIM29",
            "VSIG2",
            "UGT1A6")
DotPlot(Dsubset1, features =  feats )
ggsave(filename=paste(saveext,"/DimRed/Dsubset2_dot.pdf",sep=""),width = 20, height = 4)


Idents(Dsubset1,WhichCells(Dsubset1,idents=c(1))) <- "EmDisc"
Idents(Dsubset1,WhichCells(Dsubset1,idents=c(2))) <- "Am"
Idents(Dsubset1,WhichCells(Dsubset1,idents=c(0))) <- "Epith"

DefaultAssay(Dsubset1) <- "RNA"
DotPlot(Dsubset1, features =  feats )
ggsave(filename=paste(saveext,"/DimRed/Dsubset1_dot.pdf",sep=""),width = 20, height = 3)


DotPlot(Dsubset1, features =  feats ) 
ggsave(filename=paste(saveext,"/DimRed/Dsubset1_dot_feat.pdf",sep=""),width = 20, height = 3)

DotPlot(Dsubset1, features =  intersect( feats2, rownames(Dsubset1) ) ) 
ggsave(filename=paste(saveext,"/DimRed/Dsubset1_dot_feat2.pdf",sep=""),width = 20, height = 3)





Idents(Dsubset1) <- Dsubset1$Genotype
Dsubset3 <- subset(Dsubset1,idents = c("NAssigned","Embryonic_G1","Embryonic_G2","Embryoic_G3","Embryonic_G4"))
DefaultAssay(Dsubset3) <- "integrated"

Dsubset3 <- RunPCA(Dsubset3, npcs = 20, verbose = FALSE)
Dsubset3 <- RunUMAP(Dsubset3, reduction = "pca", dims = 1:20)
Dsubset3 <- FindNeighbors(Dsubset3, reduction = "pca", dims = 1:2)

p<-DimPlot(Dsubset3, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Amsubset3_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-DimPlot(Dsubset3, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Amsubset3_UMAP_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)



DefaultAssay(Dsubset3) <- "RNA"
FeaturePlot(Dsubset3,  reduction = "umap", features = "WNT6", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset3_WNT6.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Dsubset3,  reduction = "umap", features = "VTCN1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset3_VTCN1.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Dsubset3,  reduction = "umap", features = "POU5F1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset3_POU5F1.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Dsubset3,  reduction = "umap", features = "PDGFA", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset3_PDGFA.pdf",sep=""),width = 12, height = 10)


DefaultAssay(Dsubset3) <- "integrated"
Dsubset3 <- FindClusters(Dsubset3, resolution = 0.3)
p<-DimPlot(Dsubset3, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Amsubset3_Cl_UMAP_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)

Dsubset3$newCl <- Idents(Dsubset3)

uID <- as.character(Dsubset3$Genotype)
uID[which(uID=="NAssigned" & Dsubset3$newCl==1)] <- "Cont"

Dsubset3$Genotype <- uID
Idents(Dsubset3) <- Dsubset3$Genotype
Dsubset3 <- subset(Dsubset3,idents=c("Embryonic_G1", "NAssigned", "Embryonic_G4", "Embryonic_G2", "Embryoic_G3") )
p<-DimPlot(Dsubset3, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Amsubset3_filterd_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-DimPlot(Dsubset3, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Amsubset3_UMAP_filterd_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)




DefaultAssay(Dsubset3) <- "integrated"
Dsubset3 <- RunPCA(Dsubset3, npcs = 20, verbose = FALSE)
Dsubset3 <- RunUMAP(Dsubset3, reduction = "pca", dims = 1:20)
Dsubset3 <- FindNeighbors(Dsubset3, reduction = "pca", dims = 1:2)

Dsubset3 <- FindClusters(Dsubset3, resolution = .3)
p<-DimPlot(Dsubset3, pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Amsubset3filtered_Cl_UMAP_All.pdf",sep=""),width = 70, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
Idents(Dsubset3,WhichCells(Dsubset3,idents=c(2))) <- "EmDisc/Am"
Idents(Dsubset3,WhichCells(Dsubset3,idents=c(1))) <- "Am"
Idents(Dsubset3,WhichCells(Dsubset3,idents=c(0))) <- "EmDisc"

#Idents(Dsubset3) <- Dsubset3$Genotype
#p<-DimPlot(Dsubset3, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/Amsubset3_filterd_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
#p<-DimPlot(Dsubset3, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/Amsubset3_UMAP_filterd_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
DefaultAssay(Dsubset3) <- "RNA"
FeaturePlot(Dsubset3,  reduction = "umap", features = "WNT6", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset3filtered_WNT6.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Dsubset3,  reduction = "umap", features = "VTCN1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset3filtered_VTCN1.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Dsubset3,  reduction = "umap", features = "POU5F1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset3filtered_POU5F1.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Dsubset3,  reduction = "umap", features = "PDGFA", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset3filtered_PDGFA.pdf",sep=""),width = 12, height = 10)

FeaturePlot(Dsubset3,  reduction = "pca", features = "WNT6", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset3filtered_PCA_WNT6.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Dsubset3,  reduction = "pca", features = "VTCN1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset3filtered_PCA_VTCN1.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Dsubset3,  reduction = "pca", features = "POU5F1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset3filtered_PCA_POU5F1.pdf",sep=""),width = 12, height = 10)
FeaturePlot(Dsubset3,  reduction = "pca", features = "PDGFA", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Dsubset3filtered_PCA_PDGFA.pdf",sep=""),width = 12, height = 10)


#Idents(Dsubset3) <- paste(Dsubset3$Genotype,Dsubset3$newCl,sep="_")#
#DotPlot(Dsubset3, features =  feats )
#ggsave(filename=paste(saveext,"/DimRed/Dsubset3filtered_dot.pdf",sep=""),width = 20, height = 4)


#Idents(Dsubset3) <- Dsubset3$newCl
DotPlot(Dsubset3, features =  feats )
ggsave(filename=paste(saveext,"/DimRed/Dsubset3filtered_dot.pdf",sep=""),width = 20, height = 4)

p<-DimPlot(Dsubset3, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Amsubset3filtered_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-DimPlot(Dsubset3, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Amsubset3filtered_UMAP_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)




DsubsetTb <- subset(D,idents = c(         "C2_EVT_d14","C4_CTB_d14","C7_STB_d14",
                                          "C8_CTB_d14","C8_STB_d14","C1_CTB_d14",
                                          "C1_STB_d14","C2_CTB_d14","C6_EVT_d14",
                                          "C2_STB_d14","C3_EVT_d14","C4_EVT_d14",
                                          "C4_STB_d14","C6_CTB_d14","C6_STB_d14",
                                          "C7_EVT_d14","C8_EVT_d14") )

DefaultAssay(DsubsetTb) <- "integrated"
DsubsetTb <- RunPCA(DsubsetTb, npcs = 20, verbose = FALSE)
DsubsetTb <- RunUMAP(DsubsetTb, reduction = "pca", dims = 1:20)
DsubsetTb <- FindNeighbors(DsubsetTb, reduction = "pca", dims = 1:2)

p<-DimPlot(DsubsetTb, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tbsubset3_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-DimPlot(DsubsetTb, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tbsubset3_UMAP_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)

Idents(DsubsetTb) <- DsubsetTb$Cells

p<-DimPlot(DsubsetTb, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tbsubset3_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-DimPlot(DsubsetTb, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tbsubset3_UMAP_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)

Idents(DsubsetTb) <- DsubsetTb$ID3

p<-DimPlot(DsubsetTb, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tbsubset3id_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-DimPlot(DsubsetTb, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tbsubset3id_UMAP_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)


DefaultAssay(DsubsetTb) <- "RNA"
p<-FeaturePlot(DsubsetTb,  reduction = "umap", features = "HLA-G", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Tbsubset3_UMAP_HLAG.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-FeaturePlot(DsubsetTb,  reduction = "pca", features = "HLA-G", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Tbsubset3_PCA_HLAG.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)

p<-FeaturePlot(DsubsetTb,  reduction = "pca", features = "GATA2", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Tbsubset3_PCA_GATA2.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)

p<-FeaturePlot(DsubsetTb,  reduction = "pca", features = "TFAP2A", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Tbsubset3_PCA_TFA2A.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-FeaturePlot(DsubsetTb,  reduction = "pca", features = "TFAP2C", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Tbsubset3_PCA_TFA2C.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)


p<-FeaturePlot(DsubsetTb,  reduction = "pca", features = "GATA3", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Tbsubset3_PCA_GATA3.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)

p<-FeaturePlot(DsubsetTb,  reduction = "pca", features = "CGA", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Tbsubset3_PCA_CGA.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)




NewG1 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D1_emb.tsv", header=TRUE)
NewG2 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D2_emb.tsv", header=TRUE)
NewG3B <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D3B_emb.tsv", header=TRUE)
NewG6 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D6_emb.tsv", header=TRUE)
NewG8 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D8_emb.tsv", header=TRUE)

Idents(DsubsetTb) <- DsubsetTb$ID3
Idents(DsubsetTb, cells= paste(NewG1$barcode[which(NewG1$assignment=="0")],2,sep="_") ) <- "NC1_1"
Idents(DsubsetTb, cells= paste(NewG1$barcode[which(NewG1$assignment=="1")],2,sep="_") ) <- "NC1_2"
Idents(DsubsetTb, cells= paste(NewG1$barcode[which(NewG1$assignment=="0/1")],2,sep="_") ) <- "NC1_1"
Idents(DsubsetTb, cells= paste(NewG1$barcode[which(NewG1$assignment=="1/0")],2,sep="_") ) <- "NC1_2"

Idents(DsubsetTb, cells= paste(NewG2$barcode[which(NewG2$assignment=="0")],3,sep="_") ) <- "NC2_1"
Idents(DsubsetTb, cells= paste(NewG2$barcode[which(NewG2$assignment=="1")],3,sep="_") ) <- "NC2_2"
Idents(DsubsetTb, cells= paste(NewG2$barcode[which(NewG2$assignment=="0/1")],3,sep="_") ) <- "NC2_1"
Idents(DsubsetTb, cells= paste(NewG2$barcode[which(NewG2$assignment=="1/0")],3,sep="_") ) <- "NC2_2"

Idents(DsubsetTb, cells= paste(NewG6$barcode[which(NewG6$assignment=="0")],6,sep="_") ) <- "NC6_1"
Idents(DsubsetTb, cells= paste(NewG6$barcode[which(NewG6$assignment=="1")],6,sep="_") ) <- "NC6_2"
Idents(DsubsetTb, cells= paste(NewG6$barcode[which(NewG6$assignment=="0/1")],6,sep="_") ) <- "NC6_1"
Idents(DsubsetTb, cells= paste(NewG6$barcode[which(NewG6$assignment=="1/0")],6,sep="_") ) <- "NC6_2"

Idents(DsubsetTb, cells= paste(NewG8$barcode[which(NewG8$assignment=="0")],8,sep="_") ) <- "NC2_1"
Idents(DsubsetTb, cells= paste(NewG8$barcode[which(NewG8$assignment=="1")],8,sep="_") ) <- "NC2_2"
Idents(DsubsetTb, cells= paste(NewG8$barcode[which(NewG8$assignment=="0/1")],8,sep="_") ) <- "NC2_1"
Idents(DsubsetTb, cells= paste(NewG8$barcode[which(NewG8$assignment=="1/0")],8,sep="_") ) <- "NC2_2"

Idents(DsubsetTb, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0")],5,sep="_") ) <- "NC3_1"
Idents(DsubsetTb, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/2")],5,sep="_") ) <- "NC3_1"
#Idents(DsubsetTb, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/1")],5,sep="_") ) <- "NC3_1"
Idents(DsubsetTb, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1")],5,sep="_") ) <- "NC3_2"
Idents(DsubsetTb, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1/2")],5,sep="_") ) <- "NC3_2"
Idents(DsubsetTb, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1/0")],5,sep="_") ) <- "NC3_2"
Idents(DsubsetTb, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="2")],5,sep="_") ) <- "NC3_3"
#Idents(DsubsetTb, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="2/0")],5,sep="_") ) <- "NC3_3"
#Idents(DsubsetTb, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/1")],5,sep="_") ) <- "NC2_1"
Idents(DsubsetTb, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1/0")],5,sep="_") ) <- "NC2_2"

DsubsetTb$Genotype2 <- Idents(DsubsetTb)
Idents(DsubsetTb) <- DsubsetTb$Cells

p<-DimPlot(DsubsetTb, pt.size = 4, reduction = "pca", split.by = "Genotype2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/TbsubsetCl_PCA_All_more.pdf",sep=""),width = 72, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-DimPlot(DsubsetTb, pt.size = 4, reduction = "umap", split.by = "Genotype2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/TbsubsetCl_UMAP_All_more.pdf",sep=""),width = 72, height = 10, useDingbats = FALSE, limitsize = FALSE,p)




DefaultAssay(DsubsetTb) <- "integrated"
DsubsetTb <- FindClusters(DsubsetTb, resolution = .3)
p<-DimPlot(DsubsetTb, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tbsubset3filtered_Cl_UMAP_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-DimPlot(DsubsetTb, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tbsubset3filtered_Cl_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
DsubsetTb$newCl <- Idents(DsubsetTb)

DN1 <- as.data.frame(DsubsetTb[["pca"]]@cell.embeddings)[,1:20]
dm <- DiffusionMap(DN1)
DM <- dm@eigenvectors*100
rownames(DM) <- colnames(DsubsetTb)
DsubsetTb[["dm"]] <- CreateDimReducObject(embeddings = DM, key = "DM_", assay = DefaultAssay(DsubsetTb))
p <- DimPlot(DsubsetTb, reduction = "dm", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/DM_Tb",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE,p)
p<-FeaturePlot(DsubsetTb,  reduction = "dm", features = "HLA-G", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Tbsubset3_PCA_DM.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)

Idents(DsubsetTb) <- DsubsetTb$Genotype
p <- DimPlot(DsubsetTb, reduction = "dm", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/DM_GTTb",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE,p)
Idents(DsubsetTb) <- DsubsetTb$newCl
p <- DimPlot(DsubsetTb, reduction = "dm", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/DM_ClTb",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE,p)




DsubsetHyp <- subset(D,idents = c( "C2_Hyp_d14","C2_Hyp_d14","C1_Hyp_d14","C7_Hyp_d14"))
Idents(DsubsetHyp) <- DsubsetHyp$Genotype  
DefaultAssay(DsubsetHyp) <- "integrated"
DsubsetHyp <- RunPCA(DsubsetHyp, npcs = 20, verbose = FALSE)
DsubsetHyp <- RunUMAP(DsubsetHyp, reduction = "pca", dims = 1:20)
DsubsetHyp <- FindNeighbors(DsubsetHyp, reduction = "pca", dims = 1:2)

Idents(DsubsetHyp) <- DsubsetHyp$ID3
#Idents(DsubsetHyp,cells=WhichCells(DsubsetHyp,idents="VE_d14")) <- "Hyp_d14"
p<-DimPlot(DsubsetHyp, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Hypsubset3_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-DimPlot(DsubsetHyp, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Hypsubset3_UMAP_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)

DefaultAssay(DsubsetHyp) <- "RNA"
p<-FeaturePlot(DsubsetHyp,  reduction = "umap", features = "GATA4", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Hypsubset3_UMAP_GATA4.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-FeaturePlot(DsubsetHyp,  reduction = "umap", features = "CER1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Hypsubset3_UMAP_CER1.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)

p<-FeaturePlot(DsubsetHyp,  reduction = "umap", features = "SOX17", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Hypsubset3_UMAP_SOX17.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-FeaturePlot(DsubsetHyp,  reduction = "umap", features = "GATA6", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Hypsubset3_UMAP_GATA6.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-FeaturePlot(DsubsetHyp,  reduction = "umap", features = "PDGFRA", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Hypsubset3_UMAP_PDGFRA.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-FeaturePlot(DsubsetHyp,  reduction = "umap", features = "APOA1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Hypsubset3_UMAP_APOA1.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-FeaturePlot(DsubsetHyp,  reduction = "umap", features = "APOB", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Hypsubset3_UMAP_APOB.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)


p<-FeaturePlot(DsubsetHyp,  reduction = "pca", features = "VIM", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Messubset3_PCA_VIM.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-FeaturePlot(DsubsetHyp,  reduction = "pca", features = "HGF", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Messubset3_PCA_HGF.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-FeaturePlot(DsubsetHyp,  reduction = "pca", features = "PDGFRA", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Messubset3_PCA_PDGFRA.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)


DsubsetMes <- subset(D,idents = c( "C6_ExMes_d14","C3_ExMes_d14","C1_ExMes_d14",
                                   "C8_ExMes_d14","C7_ExMes_d14","C2_ExMes_d14"))
Idents(DsubsetMes) <- DsubsetMes$Genotype  
#DsubsetMes <- subset(DsubsetMes,idents=c("Stromal_G4","Epithelial_G3","Stromal_G3"), invert = TRUE)
DefaultAssay(DsubsetMes) <- "integrated"
DsubsetMes <- RunPCA(DsubsetMes, npcs = 20, verbose = FALSE)
DsubsetMes <- RunUMAP(DsubsetMes, reduction = "pca", dims = 1:20)
DsubsetMes <- FindNeighbors(DsubsetMes, reduction = "pca", dims = 1:2)

Idents(DsubsetMes) <- DsubsetMes$ID3
p<-DimPlot(DsubsetMes, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Messubset3_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-DimPlot(DsubsetMes, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Messubset3_UMAP_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)


p<-FeaturePlot(DsubsetMes,  reduction = "pca", features = "VIM", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Messubset3_PCA_VIM.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-FeaturePlot(DsubsetMes,  reduction = "pca", features = "HGF", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Messubset3_PCA_HGF.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-FeaturePlot(DsubsetMes,  reduction = "pca", features = "PDGFRA", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Messubset3_PCA_PDGFRA.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)


NewG1 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D1_emb.tsv", header=TRUE)
NewG2 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D2_emb.tsv", header=TRUE)
NewG3B <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D3B_emb.tsv", header=TRUE)
NewG6 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D6_emb.tsv", header=TRUE)
NewG8 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D8_emb.tsv", header=TRUE)

Idents(DsubsetMes) <- DsubsetMes$ID3
Idents(DsubsetMes, cells= paste(NewG1$barcode[which(NewG1$assignment=="0")],2,sep="_") ) <- "NC1_1"
#Idents(DsubsetMes, cells= paste(NewG1$barcode[which(NewG1$assignment=="1")],2,sep="_") ) <- "NC1_2"
Idents(DsubsetMes, cells= paste(NewG1$barcode[which(NewG1$assignment=="0/1")],2,sep="_") ) <- "NC1_1"
#Idents(DsubsetMes, cells= paste(NewG1$barcode[which(NewG1$assignment=="1/0")],2,sep="_") ) <- "NC1_2"

#Idents(DsubsetMes, cells= paste(NewG2$barcode[which(NewG2$assignment=="0")],3,sep="_") ) <- "NC2_1"
Idents(DsubsetMes, cells= paste(NewG2$barcode[which(NewG2$assignment=="1")],3,sep="_") ) <- "NC2_2"
#Idents(DsubsetMes, cells= paste(NewG2$barcode[which(NewG2$assignment=="0/1")],3,sep="_") ) <- "NC2_1"
#Idents(DsubsetMes, cells= paste(NewG2$barcode[which(NewG2$assignment=="1/0")],3,sep="_") ) <- "NC2_2"

#Idents(DsubsetMes, cells= paste(NewG6$barcode[which(NewG6$assignment=="0")],6,sep="_") ) <- "NC6_1"
Idents(DsubsetMes, cells= paste(NewG6$barcode[which(NewG6$assignment=="1")],6,sep="_") ) <- "NC6_2"
#Idents(DsubsetMes, cells= paste(NewG6$barcode[which(NewG6$assignment=="0/1")],6,sep="_") ) <- "NC6_1"
#Idents(DsubsetMes, cells= paste(NewG6$barcode[which(NewG6$assignment=="1/0")],6,sep="_") ) <- "NC6_2"

Idents(DsubsetMes, cells= paste(NewG8$barcode[which(NewG8$assignment=="0")],8,sep="_") ) <- "NC2_1"
#Idents(DsubsetMes, cells= paste(NewG8$barcode[which(NewG8$assignment=="1")],8,sep="_") ) <- "NC2_2"
#Idents(DsubsetMes, cells= paste(NewG8$barcode[which(NewG8$assignment=="0/1")],8,sep="_") ) <- "NC2_1"
#Idents(DsubsetMes, cells= paste(NewG8$barcode[which(NewG8$assignment=="1/0")],8,sep="_") ) <- "NC2_2"

#Idents(DsubsetMes, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0")],5,sep="_") ) <- "NC3_1"
#Idents(DsubsetMes, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/2")],5,sep="_") ) <- "NC3_1"
#Idents(DsubsetMes, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/1")],5,sep="_") ) <- "NC3_1"
#Idents(DsubsetMes, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1")],5,sep="_") ) <- "NC3_2"
#Idents(DsubsetMes, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1/2")],5,sep="_") ) <- "NC3_2"
#Idents(DsubsetMes, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1/0")],5,sep="_") ) <- "NC3_2"
#Idents(DsubsetMes, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="2")],5,sep="_") ) <- "NC3_3"
#Idents(DsubsetMes, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="2/0")],5,sep="_") ) <- "NC3_3"
#Idents(DsubsetMes, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/1")],5,sep="_") ) <- "NC2_1"
#Idents(DsubsetMes, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1/0")],5,sep="_") ) <- "NC2_2"

DsubsetMes$Genotype2 <- Idents(DsubsetMes)
Idents(DsubsetMes) <- DsubsetMes$Cells

p<-DimPlot(DsubsetMes, pt.size = 4, reduction = "pca", split.by = "Genotype2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/MessubsetCl_PCA_All_more.pdf",sep=""),width = 72, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-DimPlot(DsubsetMes, pt.size = 4, reduction = "umap", split.by = "Genotype2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/MessubsetCl_UMAP_All_more.pdf",sep=""),width = 72, height = 10, useDingbats = FALSE, limitsize = FALSE,p)




SubsetCombined <- subset(D,idents = c("C2_EVT_d14","C4_CTB_d14","C7_STB_d14",
                                                            "C8_CTB_d14","C8_STB_d14","C1_CTB_d14",
                                                            "C1_STB_d14","C2_CTB_d14","C6_EVT_d14",
                                                            "C2_STB_d14","C3_EVT_d14","C4_EVT_d14",
                                                            "C4_STB_d14","C6_CTB_d14","C6_STB_d14",
                                                            "C7_EVT_d14","C8_EVT_d14",
                                      "C2_Hyp_d14","C2_Hyp_d14","C1_Hyp_d14","C7_Hyp_d14",
                                      "C1_Am_d14",
                                      "C2_Am_d14" ,
                                      "C6_Am_d14",
                                      "C7_Am_d14",
                                      "C1_Am/EmDisc_d14",
                                      "C2_Am/EmDisc_d14",
                                      "C4_Am/EmDisc_d14",
                                      "C6_Am/EmDisc_d14",
                                      "C7_Am/EmDisc_d14",
                                      "C1_EmDisc_d14",
                                      "C2_EmDisc_d14",
                                      "C4_EmDisc_d14",
                                      "C6_EmDisc_d14",
                                      "C7_EmDisc_d14",
                                      "C8_EmDisc_d14","C6_ExMes_d14","C3_ExMes_d14","C1_ExMes_d14",
                                      "C8_ExMes_d14","C7_ExMes_d14","C2_ExMes_d14"
                                      ) )

Idents(SubsetCombined) <- SubsetCombined$Genotype  
#SubsetCombineds <- subset(SubsetCombined,idents=c("Stromal_G4","Epithelial_G3","Stromal_G3","Epithelial_G1"), invert = TRUE)
Idents(SubsetCombined) <- SubsetCombined$Cells 
DefaultAssay(SubsetCombined) <- "integrated"
SubsetCombined <- RunPCA(SubsetCombined, npcs = 20, verbose = FALSE)
SubsetCombined <- RunUMAP(SubsetCombined, reduction = "pca", dims = 1:20)
SubsetCombined <- FindNeighbors(SubsetCombined, reduction = "pca", dims = 1:2)


p<-DimPlot(SubsetCombined, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/AllEmbsubset_PCA_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-DimPlot(SubsetCombined, pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/AllEmbsubset_UMAP_All.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)





NewG1 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D1_emb.tsv", header=TRUE)
NewG2 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D2_emb.tsv", header=TRUE)
NewG3B <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D3B_emb.tsv", header=TRUE)
NewG6 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D6_emb.tsv", header=TRUE)
NewG8 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D8_emb.tsv", header=TRUE)

Idents(SubsetCombined) <- SubsetCombined$ID3
Idents(SubsetCombined, cells= paste(NewG1$barcode[which(NewG1$assignment=="0")],2,sep="_") ) <- "NC1_1"
Idents(SubsetCombined, cells= paste(NewG1$barcode[which(NewG1$assignment=="1")],2,sep="_") ) <- "NC1_2"
Idents(SubsetCombined, cells= paste(NewG1$barcode[which(NewG1$assignment=="0/1")],2,sep="_") ) <- "NC1_1"
Idents(SubsetCombined, cells= paste(NewG1$barcode[which(NewG1$assignment=="1/0")],2,sep="_") ) <- "NC1_2"

Idents(SubsetCombined, cells= paste(NewG2$barcode[which(NewG2$assignment=="0")],3,sep="_") ) <- "NC2_1"
Idents(SubsetCombined, cells= paste(NewG2$barcode[which(NewG2$assignment=="1")],3,sep="_") ) <- "NC2_2"
Idents(SubsetCombined, cells= paste(NewG2$barcode[which(NewG2$assignment=="0/1")],3,sep="_") ) <- "NC2_1"
Idents(SubsetCombined, cells= paste(NewG2$barcode[which(NewG2$assignment=="1/0")],3,sep="_") ) <- "NC2_2"

Idents(SubsetCombined, cells= paste(NewG6$barcode[which(NewG6$assignment=="0")],6,sep="_") ) <- "NC6_1"
Idents(SubsetCombined, cells= paste(NewG6$barcode[which(NewG6$assignment=="1")],6,sep="_") ) <- "NC6_2"
Idents(SubsetCombined, cells= paste(NewG6$barcode[which(NewG6$assignment=="0/1")],6,sep="_") ) <- "NC6_1"
Idents(SubsetCombined, cells= paste(NewG6$barcode[which(NewG6$assignment=="1/0")],6,sep="_") ) <- "NC6_2"

Idents(SubsetCombined, cells= paste(NewG8$barcode[which(NewG8$assignment=="0")],8,sep="_") ) <- "NC2_1"
Idents(SubsetCombined, cells= paste(NewG8$barcode[which(NewG8$assignment=="1")],8,sep="_") ) <- "NC2_2"
Idents(SubsetCombined, cells= paste(NewG8$barcode[which(NewG8$assignment=="0/1")],8,sep="_") ) <- "NC2_1"
Idents(SubsetCombined, cells= paste(NewG8$barcode[which(NewG8$assignment=="1/0")],8,sep="_") ) <- "NC2_2"

Idents(SubsetCombined, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0")],5,sep="_") ) <- "NC3_1"
Idents(SubsetCombined, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/2")],5,sep="_") ) <- "NC3_1"
Idents(SubsetCombined, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/1")],5,sep="_") ) <- "NC3_1"
Idents(SubsetCombined, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1")],5,sep="_") ) <- "NC3_2"
Idents(SubsetCombined, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1/2")],5,sep="_") ) <- "NC3_2"
Idents(SubsetCombined, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1/0")],5,sep="_") ) <- "NC3_2"
Idents(SubsetCombined, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="2")],5,sep="_") ) <- "NC3_3"
Idents(SubsetCombined, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="2/0")],5,sep="_") ) <- "NC3_3"
Idents(SubsetCombined, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/1")],5,sep="_") ) <- "NC2_1"
Idents(SubsetCombined, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1/0")],5,sep="_") ) <- "NC2_2"

SubsetCombined$Genotype2 <- Idents(SubsetCombined)
Idents(SubsetCombined) <- SubsetCombined$Cells

p<-DimPlot(SubsetCombined, pt.size = 4, reduction = "pca", split.by = "Genotype2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/AllCombsubsetCl_PCA_All_more.pdf",sep=""),width = 72, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-DimPlot(SubsetCombined, pt.size = 4, reduction = "umap", split.by = "Genotype2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/AllCombubsetCl_UMAP_All_more.pdf",sep=""),width = 72, height = 10, useDingbats = FALSE, limitsize = FALSE,p)



