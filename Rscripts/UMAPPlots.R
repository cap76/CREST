library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
#library(destiny)
library(ggtern)

set.seed(1)

saveext = "./FinalAlignB/"
dir.create(saveext)
dir.create(paste(saveext,"DimRed/",sep=""))
dir.create(paste(saveext,"Markers/",sep=""))


mammal.combined <- readRDS(file = paste(saveext,"Seurat_combined_FinalAno.rds",sep=""))

mammal.combined1 <- readRDS(file = paste(saveext,"FinalAno.rds",sep=""))


#[1] Unciliated epithelia SOX9P                Hyp_d14             
# [4] Prolif               Cont                 ExMes_d14           
# [7] Stromal fibroblasts  Lumenal              STB_d14             
#[10] SOX9LRG5             Am/EmDisc_d14        EmDisc_d14          
#[13] Am_d14               Ciliated             CTB_d14             
#[16] Glandular            EVT_d14              VE_d14              
Idents(mammal.combined,cells=WhichCells(mammal.combined1,idents="EVT_d14") ) <- "EVT_d14"
Idents(mammal.combined,cells=WhichCells(mammal.combined1,idents="CTB_d14") ) <- "CTB_d14"
Idents(mammal.combined,cells=WhichCells(mammal.combined1,idents="STB_d14") ) <- "STB_d14"
Idents(mammal.combined,cells=WhichCells(mammal.combined1,idents="ExMes_d14") ) <- "ExMes_d14"
Idents(mammal.combined,cells=WhichCells(mammal.combined1,idents="Am_d14") ) <- "Am_d14"
Idents(mammal.combined,cells=WhichCells(mammal.combined1,idents="VE_d14") ) <- "VE_d14"
Idents(mammal.combined,cells=WhichCells(mammal.combined1,idents="EmDisc_d14") ) <- "EmDisc_d14"
Idents(mammal.combined,cells=WhichCells(mammal.combined1,idents="Hyp_d14") ) <- "Hyp_d14"
Idents(mammal.combined,cells=WhichCells(mammal.combined1,idents="Am/EmDisc_d14") ) <- "Am/EmDisc_d14"

Idents(mammal.combined,cells=WhichCells(mammal.combined1,idents="Cont") ) <- "Cont"
Idents(mammal.combined,cells=WhichCells(mammal.combined1,idents="Lumenal") ) <- "Lumenal"
Idents(mammal.combined,cells=WhichCells(mammal.combined1,idents="Ciliated") ) <- "Ciliated"
Idents(mammal.combined,cells=WhichCells(mammal.combined1,idents="SOX9P") ) <- "SOX9P"

Idents(mammal.combined,cells=WhichCells(mammal.combined1,idents="Glandular") ) <- "Glandular"
Idents(mammal.combined,cells=WhichCells(mammal.combined1,idents="SOX9LRG5") ) <- "SOX9LRG5"
Idents(mammal.combined,cells=WhichCells(mammal.combined1,idents="Prolif") ) <- "Prolif"
Idents(mammal.combined,cells=WhichCells(mammal.combined1,idggplot
                                        ents="Stromal fibroblasts") ) <- "Stromal fibroblasts"


p<-DimPlot(mammal.combined, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_R1",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE,p)

mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20,n.neighbors = 10)
p<-DimPlot(mammal.combined, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_R2",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE,p)

mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20,n.neighbors = 20)
p<-DimPlot(mammal.combined, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_R3",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE,p)

mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20,n.neighbors = 25)
p<-DimPlot(mammal.combined, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_R4",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE,p)

mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20,n.neighbors = 50)
p<-DimPlot(mammal.combined, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_R5",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE,p)

mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20,n.neighbors = 50)
p<-DimPlot(mammal.combined, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_R6",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE,p)

mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20,n.neighbors = 100)
p<-DimPlot(mammal.combined, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_R7",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE,p)

mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20,n.neighbors = 200)
p<-DimPlot(mammal.combined, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_R8",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE,p)



wwdasdasdasdasd

DefaultAssay(mammal.combined) <- "RNA"

Idents(mammal.combined) <- factor(mammal.combined$Tissue,levels=c("Embryonic","Epithelial","Stromal"))

VlnPlot(mammal.combined,feature="GATA2") 
ggsave(filename=paste(saveext,"Vln_GATA2",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
VlnPlot(mammal.combined,feature="GATA3")
ggsave(filename=paste(saveext,"Vln_GATA3",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
VlnPlot(mammal.combined,feature="TFAP2C")
ggsave(filename=paste(saveext,"Vln_TFAP2C",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
VlnPlot(mammal.combined,feature="TFAP2A")
ggsave(filename=paste(saveext,"Vln_TFAP2A",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
VlnPlot(mammal.combined,feature="CGA")
ggsave(filename=paste(saveext,"Vln_CGA",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
VlnPlot(mammal.combined,feature="KRT7")
ggsave(filename=paste(saveext,"Vln_KRT7",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
VlnPlot(mammal.combined,feature="MMP7")
ggsave(filename=paste(saveext,"Vln_MMP7",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
VlnPlot(mammal.combined,feature="DCN")
ggsave(filename=paste(saveext,"Vln_DCN",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE)

VlnPlot(mammal.combined,feature="VIM")
ggsave(filename=paste(saveext,"Vln_VIM",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE)



VlnPlot(mammal.combined,feature="MSLN")
ggsave(filename=paste(saveext,"Vln_MSLN",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
VlnPlot(mammal.combined,feature="ESR1")
ggsave(filename=paste(saveext,"Vln_ESR1",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
VlnPlot(mammal.combined,feature="PAEP")
ggsave(filename=paste(saveext,"Vln_PAEP",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
VlnPlot(mammal.combined,feature="SORL1")
ggsave(filename=paste(saveext,"Vln_SORL1",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
VlnPlot(mammal.combined,feature="ERBB4")
ggsave(filename=paste(saveext,"Vln_ERBB4",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
VlnPlot(mammal.combined,feature="FGF14")
ggsave(filename=paste(saveext,"Vln_FGF14",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
VlnPlot(mammal.combined,feature="SPON1")
ggsave(filename=paste(saveext,"Vln_SPON1",".pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
v

s`dddsadsad
#Blastoid2_Tr
DefaultAssay(mammal.combined) <- "RNA"
M1 <- FindMarkers(mammal.combined, ident.1=c("Human (Ara EBs)_Am"),ident.2=c("Blastoid2_Tr"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1=c("Cyno (Lanner)_Am"),ident.2 = c("Cyno (Lanner)_Tr"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M3 <- FindMarkers(mammal.combined, ident.1=c("HumanCS5-7_Am"),ident.2 = c("HumanCS5-7_Tr"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M4 <- FindMarkers(mammal.combined, ident.1=c("Cyno (in vitro)_Am"),ident.2 = c("Cyno (in vitro)_Tr"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M5 <- FindMarkers(mammal.combined, ident.1=c("Marmoset_Am"),ident.2 = c("Marmoset_Tr"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
cwrite.csv(as.data.frame(M1), file=paste(saveext,"/Am_Tr_Hum2.csv",sep=""))
write.csv(as.data.frame(M2), file=paste(saveext,"/Am_Tr_Cyno2.csv",sep=""))
write.csv(as.data.frame(M3), file=paste(saveext,"/Am_Tr_Hum1.csv",sep=""))
#write.csv(as.data.frame(M4), file=paste(saveext,"/Am_Tr_Cyno1.csv",sep=""))
write.csv(as.data.frame(M5), file=paste(saveext,"/Am_Tr_Marm.csv",sep=""))

M6 <- FindMarkers(mammal.combined, ident.2=c("Human (Ara EBs)_Am"),ident.1=c("Human (Ara EBs)_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M7 <- FindMarkers(mammal.combined, ident.2=c("Cyno (Lanner)_Am"),ident.1 = c("Cyno (Lanner)_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M8 <- FindMarkers(mammal.combined, ident.2=c("HumanCS5-7_Am"),ident.1 = c("HumanCS5-7_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M9 <- FindMarkers(mammal.combined, ident.2=c("Cyno (in vitro)_Am"),ident.1 = c("Cyno (in vitro)_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M10 <- FindMarkers(mammal.combined, ident.2=c("Marmoset_Am"),ident.1 = c("Marmoset_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(M6), file=paste(saveext,"/STB_Am_Hum2.csv",sep=""))
write.csv(as.data.frame(M7), file=paste(saveext,"/STB_Am_Cyno2.csv",sep=""))
write.csv(as.data.frame(M8), file=paste(saveext,"/STB_Am_Hum1.csv",sep=""))
write.csv(as.data.frame(M9), file=paste(saveext,"/STB_Am_Cyno1.csv",sep=""))
write.csv(as.data.frame(M10), file=paste(saveext,"/STB_Am_Marm.csv",sep=""))

M11 <- FindMarkers(mammal.combined, ident.2=c("Human (Ara EBs)_Am"),ident.1=c("Human (Ara EBs)_CTB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M12 <- FindMarkers(mammal.combined, ident.2=c("Cyno (Lanner)_Am"),ident.1 = c("Cyno (Lanner)_CTB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M13 <- FindMarkers(mammal.combined, ident.2=c("HumanCS5-7_Am"),ident.1 = c("HumanCS5-7_CTB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M14 <- FindMarkers(mammal.combined, ident.2=c("Cyno (in vitro)_Am"),ident.1 = c("Cyno (in vitro)_CTB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M15 <- FindMarkers(mammal.combined, ident.2=c("Marmoset_Am"),ident.1 = c("Marmoset_CTB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(M11), file=paste(saveext,"/CTB_Am_Hum2.csv",sep=""))
write.csv(as.data.frame(M12), file=paste(saveext,"/CTB_Am_Cyno2.csv",sep=""))
write.csv(as.data.frame(M13), file=paste(saveext,"/CTB_Am_Hum1.csv",sep=""))
write.csv(as.data.frame(M14), file=paste(saveext,"/CTB_Am_Cyno1.csv",sep=""))
write.csv(as.data.frame(M15), file=paste(saveext,"/CTB_Am_Marm.csv",sep=""))

M16 <- FindMarkers(mammal.combined, ident.2=c("Human (Ara EBs)_Am"),ident.1=c("Human (Ara EBs)_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M17 <- FindMarkers(mammal.combined, ident.2=c("Cyno (Lanner)_Am"),ident.1 = c("Cyno (Lanner)_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M18 <- FindMarkers(mammal.combined, ident.2=c("HumanCS5-7_Am"),ident.1 = c("HumanCS5-7_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M19 <- FindMarkers(mammal.combined, ident.2=c("Cyno (in vitro)_Am"),ident.1 = c("Cyno (in vitro)_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(M16), file=paste(saveext,"/EVT_Am_Hum2.csv",sep=""))
write.csv(as.data.frame(M17), file=paste(saveext,"/EVT_Am_Cyno2.csv",sep=""))
write.csv(as.data.frame(M18), file=paste(saveext,"/EVT_Am_Hum1.csv",sep=""))
write.csv(as.data.frame(M19), file=paste(saveext,"/EVT_Am_Cyno1.csv",sep=""))

M20 <- FindMarkers(mammal.combined, ident.2=c("Blastoid2_Tr"),ident.1=c("Human (Ara EBs)_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M21 <- FindMarkers(mammal.combined, ident.2=c("Cyno (Lanner)_Tr"),ident.1 = c("Cyno (Lanner)_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M22 <- FindMarkers(mammal.combined, ident.2=c("HumanCS5-7_Tr"),ident.1 = c("HumanCS5-7_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M23 <- FindMarkers(mammal.combined, ident.2=c("Cyno (in vitro)_Tr"),ident.1 = c("Cyno (in vitro)_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M24 <- FindMarkers(mammal.combined, ident.2=c("Marmoset_Tr"),ident.1 = c("Marmoset_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(M20), file=paste(saveext,"/STB_Tr_Hum2.csv",sep=""))
write.csv(as.data.frame(M21), file=paste(saveext,"/STB_Tr_Cyno2.csv",sep=""))
write.csv(as.data.frame(M22), file=paste(saveext,"/STB_Tr_Hum1.csv",sep=""))
#write.csv(as.data.frame(M23), file=paste(saveext,"/STB_Tr_Cyno1.csv",sep=""))
write.csv(as.data.frame(M24), file=paste(saveext,"/STB_Tr_Marm.csv",sep=""))

#STB_Tr_Cyno2.csv
#Cyno (in vitro)_Tr

M25 <- FindMarkers(mammal.combined, ident.2=c("Blastoid2_Tr"),ident.1=c("Human (Ara EBs)_CTB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M26 <- FindMarkers(mammal.combined, ident.2=c("Cyno (Lanner)_Tr"),ident.1 = c("Cyno (Lanner)_CTB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M27 <- FindMarkers(mammal.combined, ident.2=c("HumanCS5-7_Tr"),ident.1 = c("HumanCS5-7_CTB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M28 <- FindMarkers(mammal.combined, ident.2=c("Cyno (in vitro)_Tr"),ident.1 = c("Cyno (in vitro)_CTB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M29 <- FindMarkers(mammal.combined, ident.2=c("Marmoset_Tr"),ident.1 = c("Marmoset_CTB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(M25), file=paste(saveext,"/CTB_Tr_Hum2.csv",sep=""))
write.csv(as.data.frame(M26), file=paste(saveext,"/CTB_Tr_Cyno2.csv",sep=""))
write.csv(as.data.frame(M27), file=paste(saveext,"/CTB_Tr_Hum1.csv",sep=""))
##write.csv(as.data.frame(M28), file=paste(saveext,"/CTB_Tr_Cyno1.csv",sep=""))
write.csv(as.data.frame(M29), file=paste(saveext,"/CTB_Tr_Marm.csv",sep=""))


M30 <- FindMarkers(mammal.combined, ident.2=c("Blastoid2_Tr"),ident.1=c("Human (Ara EBs)_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M31 <- FindMarkers(mammal.combined, ident.2=c("Cyno (Lanner)_Tr"),ident.1 = c("Cyno (Lanner)_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M32 <- FindMarkers(mammal.combined, ident.2=c("HumanCS5-7_Tr"),ident.1 = c("HumanCS5-7_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M33 <- FindMarkers(mammal.combined, ident.2=c("Cyno (in vitro)_Tr"),ident.1 = c("Cyno (in vitro)_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(M30), file=paste(saveext,"/Tr_EVT_Hum2.csv",sep=""))
write.csv(as.data.frame(M31), file=paste(saveext,"/Tr_EVT_Cyno2.csv",sep=""))
write.csv(as.data.frame(M32), file=paste(saveext,"/Tr_EVT_Hum1.csv",sep=""))
#write.csv(as.data.frame(M33), file=paste(saveext,"/Tr_EVT_Cyno1.csv",sep=""))

M34 <- FindMarkers(mammal.combined, ident.2=c("Human (Ara EBs)_CTB"),ident.1=c("Human (Ara EBs)_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M35 <- FindMarkers(mammal.combined, ident.2=c("Cyno (Lanner)_CTB"),ident.1 = c("Cyno (Lanner)_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M36 <- FindMarkers(mammal.combined, ident.2=c("HumanCS5-7_CTB"),ident.1 = c("HumanCS5-7_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M37 <- FindMarkers(mammal.combined, ident.2=c("Cyno (in vitro)_CTB"),ident.1 = c("Cyno (in vitro)_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M38 <- FindMarkers(mammal.combined, ident.2=c("Marmoset_CTB"),ident.1 = c("Marmoset_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(M34), file=paste(saveext,"/STB_CTB_Hum2.csv",sep=""))
write.csv(as.data.frame(M35), file=paste(saveext,"/STB_CTB_Cyno2.csv",sep=""))
write.csv(as.data.frame(M36), file=paste(saveext,"/STB_CTB_Hum1.csv",sep=""))
write.csv(as.data.frame(M37), file=paste(saveext,"/STB_CTB_Cyno1.csv",sep=""))
write.csv(as.data.frame(M38), file=paste(saveext,"/STB_CTB_Marm.csv",sep=""))

M39 <- FindMarkers(mammal.combined, ident.2=c("Human (Ara EBs)_CTB"),ident.1=c("Human (Ara EBs)_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M40 <- FindMarkers(mammal.combined, ident.2=c("Cyno (Lanner)_CTB"),ident.1 = c("Cyno (Lanner)_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M41 <- FindMarkers(mammal.combined, ident.2=c("HumanCS5-7_CTB"),ident.1 = c("HumanCS5-7_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M42 <- FindMarkers(mammal.combined, ident.2=c("Cyno (in vitro)_CTB"),ident.1 = c("Cyno (in vitro)_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(M39), file=paste(saveext,"/EVT_CTB_Hum2.csv",sep=""))
write.csv(as.data.frame(M40), file=paste(saveext,"/EVT_CTB_Cyno2.csv",sep=""))
write.csv(as.data.frame(M41), file=paste(saveext,"/EVT_CTB_Hum1.csv",sep=""))
write.csv(as.data.frame(M42), file=paste(saveext,"/EVT_CTB_Cyno1.csv",sep=""))

M43 <- FindMarkers(mammal.combined, ident.2=c("Human (Ara EBs)_STB"),ident.1=c("Human (Ara EBs)_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M44 <- FindMarkers(mammal.combined, ident.2=c("Cyno (Lanner)_STB"),ident.1 = c("Cyno (Lanner)_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M45 <- FindMarkers(mammal.combined, ident.2=c("HumanCS5-7_STB"),ident.1 = c("HumanCS5-7_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M46 <- FindMarkers(mammal.combined, ident.2=c("Cyno (in vitro)_STB"),ident.1 = c("Cyno (in vitro)_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(M43), file=paste(saveext,"/EVT_STB_Hum2.csv",sep=""))
write.csv(as.data.frame(M44), file=paste(saveext,"/EVT_STB_Cyno2.csv",sep=""))
write.csv(as.data.frame(M45), file=paste(saveext,"/EVT_STB_Hum1.csv",sep=""))
write.csv(as.data.frame(M46), file=paste(saveext,"/EVT_STB_Cyno1.csv",sep=""))

#STB_Tr_Cyno2.csv
#Cyno (in vitro)_Tr

saveRDS(M1,file = paste(saveext,"/M1.rds",sep=""))
saveRDS(M2,file = paste(saveext,"/M2.rds",sep=""))
saveRDS(M3,file = paste(saveext,"/M3.rds",sep=""))
#saveRDS(M4,file = paste(saveext,"/M4.rds",sep=""))
saveRDS(M5,file = paste(saveext,"/M5.rds",sep=""))
saveRDS(M6,file = paste(saveext,"/M6.rds",sep=""))
saveRDS(M7,file = paste(saveext,"/M7.rds",sep=""))
saveRDS(M8,file = paste(saveext,"/M8.rds",sep=""))
saveRDS(M9,file = paste(saveext,"/M9.rds",sep=""))
saveRDS(M10,file = paste(saveext,"/M10.rds",sep=""))
saveRDS(M11,file = paste(saveext,"/M11.rds",sep=""))
saveRDS(M12,file = paste(saveext,"/M12.rds",sep=""))
saveRDS(M13,file = paste(saveext,"/M13.rds",sep=""))
saveRDS(M14,file = paste(saveext,"/M14.rds",sep=""))
saveRDS(M15,file = paste(saveext,"/M15.rds",sep=""))
saveRDS(M16,file = paste(saveext,"/M16.rds",sep=""))
saveRDS(M17,file = paste(saveext,"/M17.rds",sep=""))
saveRDS(M18,file = paste(saveext,"/M18.rds",sep=""))
saveRDS(M19,file = paste(saveext,"/M19.rds",sep=""))
saveRDS(M20,file = paste(saveext,"/M20.rds",sep=""))
saveRDS(M21,file = paste(saveext,"/M21.rds",sep=""))
saveRDS(M22,file = paste(saveext,"/M22.rds",sep=""))
#saveRDS(M23,file = paste(saveext,"/M23.rds",sep=""))
saveRDS(M24,file = paste(saveext,"/M24.rds",sep=""))
saveRDS(M25,file = paste(saveext,"/M25.rds",sep=""))
saveRDS(M26,file = paste(saveext,"/M26.rds",sep=""))
saveRDS(M27,file = paste(saveext,"/M27.rds",sep=""))
#saveRDS(M28,file = paste(saveext,"/M28.rds",sep=""))
saveRDS(M29,file = paste(saveext,"/M29.rds",sep=""))
saveRDS(M30,file = paste(saveext,"/M30.rds",sep=""))
saveRDS(M31,file = paste(saveext,"/M31.rds",sep=""))
saveRDS(M32,file = paste(saveext,"/M32.rds",sep=""))
#saveRDS(M33,file = paste(saveext,"/M33.rds",sep=""))
saveRDS(M34,file = paste(saveext,"/M34.rds",sep=""))
saveRDS(M35,file = paste(saveext,"/M35.rds",sep=""))
saveRDS(M36,file = paste(saveext,"/M36.rds",sep=""))
saveRDS(M37,file = paste(saveext,"/M37.rds",sep=""))
saveRDS(M38,file = paste(saveext,"/M38.rds",sep=""))
saveRDS(M39,file = paste(saveext,"/M39.rds",sep=""))
saveRDS(M40,file = paste(saveext,"/M40.rds",sep=""))
saveRDS(M41,file = paste(saveext,"/M41.rds",sep=""))
saveRDS(M42,file = paste(saveext,"/M42.rds",sep=""))
saveRDS(M43,file = paste(saveext,"/M43.rds",sep=""))
saveRDS(M44,file = paste(saveext,"/M44.rds",sep=""))
saveRDS(M45,file = paste(saveext,"/M45.rds",sep=""))
saveRDS(M46,file = paste(saveext,"/M46.rds",sep=""))

M1 <- readRDS(file = paste(saveext,"/M1.rds",sep=""))
M2 <- readRDS(file = paste(saveext,"/M2.rds",sep=""))
M3 <- readRDS(file = paste(saveext,"/M3.rds",sep=""))
#M4 <- readRDS(file = paste(saveext,"/M4.rds",sep=""))
M5 <- readRDS(file = paste(saveext,"/M5.rds",sep=""))
M6 <- readRDS(file = paste(saveext,"/M6.rds",sep=""))
M7 <- readRDS(file = paste(saveext,"/M7.rds",sep=""))
M8 <- readRDS(file = paste(saveext,"/M8.rds",sep=""))
M9 <- readRDS(file = paste(saveext,"/M9.rds",sep=""))
M10 <- readRDS(file = paste(saveext,"/M10.rds",sep=""))
M11 <- readRDS(file = paste(saveext,"/M11.rds",sep=""))
M12 <- readRDS(file = paste(saveext,"/M12.rds",sep=""))
M13 <- readRDS(file = paste(saveext,"/M13.rds",sep=""))
M14 <- readRDS(file = paste(saveext,"/M14.rds",sep=""))
M15 <- readRDS(file = paste(saveext,"/M15.rds",sep=""))
M16 <- readRDS(file = paste(saveext,"/M16.rds",sep=""))
M17 <- readRDS(file = paste(saveext,"/M17.rds",sep=""))
M18 <- readRDS(file = paste(saveext,"/M18.rds",sep=""))
M19 <- readRDS(file = paste(saveext,"/M19.rds",sep=""))
M20 <- readRDS(file = paste(saveext,"/M20.rds",sep=""))
M21 <- readRDS(file = paste(saveext,"/M21.rds",sep=""))
M22 <- readRDS(file = paste(saveext,"/M22.rds",sep=""))
#M23 <- readRDS(file = paste(saveext,"/M23.rds",sep=""))
M24 <- readRDS(file = paste(saveext,"/M24.rds",sep=""))
M25 <- readRDS(file = paste(saveext,"/M25.rds",sep=""))
M26 <- readRDS(file = paste(saveext,"/M26.rds",sep=""))
M27 <- readRDS(file = paste(saveext,"/M27.rds",sep=""))
#M28 <- readRDS(file = paste(saveext,"/M28.rds",sep=""))
M29 <- readRDS(file = paste(saveext,"/M29.rds",sep=""))
M30 <- readRDS(file = paste(saveext,"/M30.rds",sep=""))
M31 <- readRDS(file = paste(saveext,"/M31.rds",sep=""))
M32 <- readRDS(file = paste(saveext,"/M32.rds",sep=""))
#M33 <- readRDS(file = paste(saveext,"/M33.rds",sep=""))
M34 <- readRDS(file = paste(saveext,"/M34.rds",sep=""))
M35 <- readRDS(file = paste(saveext,"/M35.rds",sep=""))
M36 <- readRDS(file = paste(saveext,"/M36.rds",sep=""))
M37 <- readRDS(file = paste(saveext,"/M37.rds",sep=""))
M38 <- readRDS(file = paste(saveext,"/M38.rds",sep=""))
M39 <- readRDS(file = paste(saveext,"/M39.rds",sep=""))
M40 <- readRDS(file = paste(saveext,"/M40.rds",sep=""))
M41 <- readRDS(file = paste(saveext,"/M41.rds",sep=""))
M42 <- readRDS(file = paste(saveext,"/M42.rds",sep=""))
M43 <- readRDS(file = paste(saveext,"/M43.rds",sep=""))
M44 <- readRDS(file = paste(saveext,"/M44.rds",sep=""))
M45 <- readRDS(file = paste(saveext,"/M45.rds",sep=""))
M46 <- readRDS(file = paste(saveext,"/M46.rds",sep=""))

#Now do some plotting
Ae <- AverageExpression(mammal.combined)
Ae <- Ae$RNA
Ae$gene <- rownames(Ae)

SIGNAL<-read.table("./LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]
TF<- read.table("./epifactors.txt",header = F)
TF <- TF$V1

saveRDS(Ae,file = paste(saveext,"/Ae.rds",sep=""))
Ae <- readRDS(file = paste(saveext,"/Ae.rds",sep=""))


cells <- c("HumanCS5-7_Am","Cyno (in vitro)_Am","Human (Ara EBs)_Am","Cyno (Lanner)_Am",
"HumanCS5-7_Tr","Blastoid2_Tr","Cyno (Lanner)_Tr",
"HumanCS5-7_CTB","Cyno (in vitro)_CTB","Human (Ara EBs)_CTB","Cyno (Lanner)_CTB",
"HumanCS5-7_STB","Cyno (in vitro)_STB","Human (Ara EBs)_STB","Cyno (Lanner)_STB",
"HumanCS5-7_EVT","Cyno (in vitro)_EVT","Human (Ara EBs)_EVT","Cyno (Lanner)_EVT")

cells1 <- c("HumanCS5-7_Am","HumanCS5-7_Tr","HumanCS5-7_CTB","HumanCS5-7_STB","HumanCS5-7_EVT")

x <- as.data.frame(Ae[TF,cells1])
#x[,cells]
row_sub = apply(x, 1, function(row) all(row !=0 ))
x <- x[row_sub,]
library(pheatmap)
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(-2, 2, length.out = 60)
pheatmap(x, breaks = mat_breaks, color =  redblue1(60), border_color = NA, cluster_rows=TRUE,cluster_cols=FALSE, scale ="row", filename = paste(saveext,"/HM1_norm",".pdf",sep=""),width=6,height=15)

cells1 <- c("Cyno (in vitro)_Am","Cyno (in vitro)_CTB","Cyno (in vitro)_STB","Cyno (in vitro)_EVT")
x <- as.data.frame(Ae[TF,cells1])
row_sub = apply(x, 1, function(row) all(row !=0 ))
x <- x[row_sub,]
library(pheatmap)
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(-2, 2, length.out = 60)
pheatmap(x, breaks = mat_breaks, color =  redblue1(60), border_color = NA, cluster_rows=TRUE,cluster_cols=FALSE, scale ="row", filename = paste(saveext,"/HM2_norm",".pdf",sep=""),width=6,height=15)



#sadasd

###1) Tr Amm
Cl1 <- M5
Cl2 <- M3
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Marmoset_Am'>0.1 | Ae5$'Marmoset_Tr'>0.1 | Ae5$'HumanCS5-7_Tr'>0.1 | Ae5$'HumanCS5-7_Am'>0.1)  ) ,]
#Create a set of genes we 
§want to label e.g., TFs
genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1Marm_AmTr_TF_fin.pdf",sep=""),width = 13, height = 13, plot = p1)



dev.off()


####Cyno1 Huma1 Am Tr
Cl1 <- M5
Cl2 <- M3
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Marmoset_Am'>0.1 | Ae5$'Marmoset_Tr'>0.1 | Ae5$'HumanCS5-7_Tr'>0.1 | Ae5$'HumanCS5-7_Am'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1Marm_AmTr_TF_seg1__fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
#Create a set of genes we want to label e.g., TFs
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1Marm_AmTr_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#Cyno1 Huma1 Am Tr
Cl1 <- M5
Cl2 <- M1
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Marmoset_Am'>0.1 | Ae5$'Marmoset_Tr'>0.1 | Ae5$'Blastoid2_Tr'>0.1 | Ae5$'Human (Ara EBs)_Am'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2Marm_AmTr_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2Marm_AmTr_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#2)   ###################
#Next  Am STB
#Cyno2 Human2 Am Tr
Cl1 <- M6
Cl2 <- M7
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Human (Ara EBs)_Am'>0.1 | Ae5$'Human (Ara EBs)_STB'>0.1 | Ae5$'Cyno (Lanner)_STB'>0.1 | Ae5$'Cyno (Lanner)_Am'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2Human2_AmSTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2Human2_AmSTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




#Cyno1 Huma1
Cl1 <- M8
Cl2 <- M9
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'HumanCS5-7_Am'>0.1 | Ae5$'HumanCS5-7_STB'>0.1 | Ae5$'Cyno (in vitro)_STB'>0.1 | Ae5$'Cyno (in vitro)_Am'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Human1_AmSTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Human1_AmSTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#Cyno1 Huma1 Am Tr
Cl1 <- M10
Cl2 <- M9
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Marmoset_Am'>0.1 | Ae5$'Marmoset_STB'>0.1 | Ae5$'Cyno (in vitro)_Am'>0.1 | Ae5$'Cyno (in vitro)_STB'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Hum1_AmSTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Hum1_AmSTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#Cyno1 Huma1 Am Tr
Cl1 <- M10
Cl2 <- M8
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Marmoset_Am'>0.1 | Ae5$'Marmoset_STB'>0.1 | Ae5$'HumanCS5-7_STB'>0.1 | Ae5$'HumanCS5-7_Am'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Hum1Marm_AmSTB_TF_seg1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Hum1Marm_AmSTB_TF_seg2.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



#3)   ###################
#Next  Am CTB
#Cyno2 Human2 Am Tr
Cl1 <- M11
Cl2 <- M12
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Human (Ara EBs)_Am'>0.1 | Ae5$'Human (Ara EBs)_CTB'>0.1 | Ae5$'Cyno (Lanner)_CTB'>0.1 | Ae5$'Cyno (Lanner)_Am'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))#
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))


p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2Human2_AmCTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2Human2_AmCTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Cyno1 Huma1 Am Tr
Cl1 <- M13
Cl2 <- M14
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'HumanCS5-7_Am'>0.1 | Ae5$'HumanCS5-7_CTB'>0.1 | Ae5$'Cyno (in vitro)_CTB'>0.1 | Ae5$'Cyno (in vitro)_Am'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))

#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Human1_AmCTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Human1_AmCTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Cyno1 Huma1 Am Tr
Cl1 <- M15
Cl2 <- M14
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Marmoset_Am'>0.1 | Ae5$'Marmoset_CTB'>0.1 | Ae5$'Cyno (in vitro)_CTB'>0.1 | Ae5$'Cyno (in vitro)_Am'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))

#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Marm_AmCTB_TF_seg1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Marm_AmCTB_TF_seg2.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Cyno1 Huma1 Am Tr
Cl1 <- M15
Cl2 <- M13
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Marmoset_Am'>0.1 | Ae5$'Marmoset_CTB'>0.1 | Ae5$'HumanCS5-7_CTB'>0.1 | Ae5$'HumanCS5-7_Am'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Hum1Marm_AmCTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Hum1Marm_AmCTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#4)   ###################
#Next  Am EVT
#Cyno2 Human2 Am Tr
Cl1 <- M16
Cl2 <- M17
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Human (Ara EBs)_Am'>0.1 | Ae5$'Human (Ara EBs)_EVT'>0.1 | Ae5$'Cyno (Lanner)_EVT'>0.1 | Ae5$'Cyno (Lanner)_Am'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2Human2_AmEVT_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2Human2_AmEVT_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




#Cyno1Human1_STBEVT_TF_seg2_fin.pdf


#Cyno1 Huma1 Am Tr
Cl1 <- M18
Cl2 <- M19
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'HumanCS5-7_Am'>0.1 | Ae5$'HumanCS5-7_EVT'>0.1 | Ae5$'Cyno (in vitro)_EVT'>0.1 | Ae5$'Cyno (in vitro)_Am'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Human1_AmEVT_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Human1_AmEVT_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Cyno1Human1_STBEVT_TF_seg2_fin.pdf
#5)   ###################
#Next  Tr STB
Cl1 <- M24
Cl2 <- M22
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Marmoset_Tr'>0.1 | Ae5$'Marmoset_STB'>0.1 | Ae5$'HumanCS5-7_STB'>0.1 | Ae5$'HumanCS5-7_Tr'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Hum1Marm_TrSTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Hum1Marm_TrSTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- M24
Cl2 <- M20
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Marmoset_Tr'>0.1 | Ae5$'Marmoset_STB'>0.1 | Ae5$'Human (Ara EBs)_STB'>0.1 | Ae5$'Blastoid2_Tr'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Hum2Marm_TrSTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Hum2Marm_TrSTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#6)   ###################
#Next  Tr CTB
#Cyno2 Human2 Am Tr
#Cyno1 Huma1 Am Tr
Cl1 <- M29
Cl2 <- M27
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Marmoset_Tr'>0.1 | Ae5$'Marmoset_CTB'>0.1 | Ae5$'HumanCS5-7_CTB'>0.1 | Ae5$'HumanCS5-7_Tr'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Hum1Marm_TrCTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Hum1Marm_TrCTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- M29
Cl2 <- M25
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Marmoset_Tr'>0.1 | Ae5$'Marmoset_CTB'>0.1 | Ae5$'Human (Ara EBs)_CTB'>0.1 | Ae5$'Blastoid2_Tr'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Hum2Marm_TrCTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Hum2Marm_TrCTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#8)   ###################
#Next  stB CTB
#Cyno2 Human2 Am Tr
Cl1 <- M34
Cl2 <- M35
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Human (Ara EBs)_STB'>0.1 | Ae5$'Human (Ara EBs)_CTB'>0.1 | Ae5$'Cyno (Lanner)_CTB'>0.1 | Ae5$'Cyno (Lanner)_STB'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2Human2_STBCTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2Human2_STBCTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#Cyno1 Huma1 Am Tr
Cl1 <- M36
Cl2 <- M37
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'HumanCS5-7_STB'>0.1 | Ae5$'HumanCS5-7_CTB'>0.1 | Ae5$'Cyno (in vitro)_CTB'>0.1 | Ae5$'Cyno (in vitro)_STB'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Human1_STBCTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Human1_STBCTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Cyno1 Huma1 Am Tr
Cl1 <- M38
Cl2 <- M37
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Marmoset_STB'>0.1 | Ae5$'Marmoset_CTB'>0.1 | Ae5$'Cyno (in vitro)_CTB'>0.1 | Ae5$'Cyno (in vitro)_STB'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Marm_STBCTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Marm_STBCTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



#Cyno1 Huma1 Am Tr
Cl1 <- M38
Cl2 <- M36
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'HumanCS5-7_STB'>0.1 | Ae5$'HumanCS5-7_CTB'>0.1 | Ae5$'Marmoset_STB'>0.1 | Ae5$'Marmoset_CTB'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1Marm_STBCTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1Marm_STBCTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#9)   ###################
#Next  CTB EVT
#Cyno2 Human2 Am Tr
Cl1 <- M39
Cl2 <- M40
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Human (Ara EBs)_EVT'>0.1 | Ae5$'Human (Ara EBs)_CTB'>0.1 | Ae5$'Cyno (Lanner)_CTB'>0.1 | Ae5$'Cyno (Lanner)_EVT'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2Human2_CTBEVT_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2Human2_CTBEVT_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Cyno1 Huma1 Am Tr
Cl1 <- M41
Cl2 <- M42
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'HumanCS5-7_EVT'>0.1 | Ae5$'HumanCS5-7_CTB'>0.1 | Ae5$'Cyno (in vitro)_CTB'>0.1 | Ae5$'Cyno (in vitro)_EVT'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Human1_CTBEVT_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Human1_CTBEVT_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#10)   ###################
#Next  STB EVT
#Cyno2 Human2 Am Tr
Cl1 <- M43
Cl2 <- M44
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Human (Ara EBs)_EVT'>0.1 | Ae5$'Human (Ara EBs)_CTB'>0.1 | Ae5$'Cyno (Lanner)_CTB'>0.1 | Ae5$'Cyno (Lanner)_EVT'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2Human2_STBEVT_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2Human2_STBEVT_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
#Cyno1Human1_STBEVT_TF_seg2_fin.pdf
#Cyno1 Huma1 Am Tr
Cl1 <- M45
Cl2 <- M46
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'HumanCS5-7_EVT'>0.1 | Ae5$'HumanCS5-7_CTB'>0.1 | Ae5$'Cyno (in vitro)_CTB'>0.1 | Ae5$'Cyno (in vitro)_EVT'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Human1_STBEVT_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),TF))
genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
genes.to.label1 = setdiff(genes.to.label5,genes.to.label3)
genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1Human1_STBEVT_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Cyno1Human1_STBEVT_TF_seg2_fin.pdf


#M34 <- FindMarkers(mammal.combined, ident.2=c("Human (Ara EBs)_CTB"),ident.1=c("Human (Ara EBs)_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M36 <- FindMarkers(mammal.combined, ident.2=c("HumanCS5-7_CTB"),ident.1 = c("HumanCS5-7_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#write.csv(as.data.frame(M34), file=paste(saveext,"/STB_CTB_Hum2.csv",sep=""))
#write.csv(as.data.frame(M35), file=paste(saveext,"/STB_CTB_Cyno2.csv",sep=""))
#write.csv(as.data.frame(M36), file=paste(saveext,"/STB_CTB_Hum1.csv",sep=""))
#write.csv(as.data.frame(M37), file=paste(saveext,"/STB_CTB_Cyno1.csv",sep=""))
#write.csv(as.data.frame(M38), file=paste(saveext,"/STB_CTB_Marm.csv",sep=""))

#M39 <- FindMarkers(mammal.combined, ident.2=c("Human (Ara EBs)_CTB"),ident.1=c("Human (Ara EBs)_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M41 <- FindMarkers(mammal.combined, ident.2=c("HumanCS5-7_CTB"),ident.1 = c("HumanCS5-7_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#write.csv(as.data.frame(M39), file=paste(saveext,"/EVT_CTB_Hum2.csv",sep=""))
#write.csv(as.data.frame(M40), file=paste(saveext,"/EVT_CTB_Cyno2.csv",sep=""))
#write.csv(as.data.frame(M41), file=paste(saveext,"/EVT_CTB_Hum1.csv",sep=""))
#write.csv(as.data.frame(M42), file=paste(saveext,"/EVT_CTB_Cyno1.csv",sep=""))



#Cyno1 Huma1 Am Tr
Cl1 <- M34
Cl2 <- M39
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'Human (Ara EBs)_EVT'>0.1 | Ae5$'Human (Ara EBs)_STB'>0.1 | Ae5$'Human (Ara EBs)_CTB'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))]),TF),Ae5$gene)
#intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF)
genes.to.label2 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF),Ae5$gene)
genes.to.label3 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)> log(1.2))]),TF),Ae5$gene)
genes.to.label4 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)> log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF),Ae5$gene)
genes.to.label5 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label6 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label7 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)
genes.to.label8 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)


saveRDS(genes.to.label1,file=paste(saveext,"Human2_CTBSTBEVT_seg1.rds",sep=""))
saveRDS(genes.to.label2,file=paste(saveext,"Human2_CTBSTBEVT_seg2.rds",sep=""))
saveRDS(genes.to.label3,file=paste(saveext,"Human2_CTBSTBEVT_seg3.rds",sep=""))
saveRDS(genes.to.label4,file=paste(saveext,"Human2_CTBSTBEVT_seg4.rds",sep=""))
saveRDS(genes.to.label5,file=paste(saveext,"Human2_CTBSTBEVT_seg5.rds",sep=""))
saveRDS(genes.to.label6,file=paste(saveext,"Human2_CTBSTBEVT_seg6.rds",sep=""))
saveRDS(genes.to.label7,file=paste(saveext,"Human2_CTBSTBEVT_seg7.rds",sep=""))
saveRDS(genes.to.label8,file=paste(saveext,"Human2_CTBSTBEVT_seg8.rds",sep=""))

#genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
#intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_CTBSTBEVT_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_CTBSTBEVT_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
#intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
#genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
#intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
#genes.to.label5 = intersect(intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF),Ae5$gene)
#genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
#genes.to.label3 = setdiff(genes.to.label5,c(genes.to.label1,genes.to.label2) )
#genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label3, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_CTBSTBEVT_TF_seg3_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label4, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_CTBSTBEVT_TF_seg4_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_CTBSTBEVT_TF_seg5_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label6, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_CTBSTBEVT_TF_seg6_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label7, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_CTBSTBEVT_TF_seg7_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label8, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_CTBSTBEVT_TF_seg8_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Cyno1Human1_STBEVT_TF_seg2_fin.pdf


#Cyno1 Huma1 Am Tr
Cl1 <- M36
Cl2 <- M41
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'HumanCS5-7_EVT'>0.1 | Ae5$'HumanCS5-7_CTB'>0.1 | Ae5$'Cyno (in vitro)_STB'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))]),TF)
genes.to.label2 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF)
genes.to.label3 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)> log(1.2))]),TF),Ae5$gene)
genes.to.label4 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)> log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF),Ae5$gene)
genes.to.label5 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label6 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label7 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)
genes.to.label8 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)

saveRDS(genes.to.label1,file=paste(saveext,"Human1_CTBSTBEVT_seg1.rds",sep=""))
saveRDS(genes.to.label2,file=paste(saveext,"Human1_CTBSTBEVT_seg2.rds",sep=""))
saveRDS(genes.to.label3,file=paste(saveext,"Human1_CTBSTBEVT_seg3.rds",sep=""))
saveRDS(genes.to.label4,file=paste(saveext,"Human1_CTBSTBEVT_seg4.rds",sep=""))
saveRDS(genes.to.label5,file=paste(saveext,"Human1_CTBSTBEVT_seg5.rds",sep=""))
saveRDS(genes.to.label6,file=paste(saveext,"Human1_CTBSTBEVT_seg6.rds",sep=""))
saveRDS(genes.to.label7,file=paste(saveext,"Human1_CTBSTBEVT_seg7.rds",sep=""))
saveRDS(genes.to.label8,file=paste(saveext,"Human1_CTBSTBEVT_seg8.rds",sep=""))

#genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
#intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_CTBSTBEVT_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_CTBSTBEVT_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
#genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
#intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
#genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
#intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
#genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
#genes.to.label3 = setdiff(genes.to.label5,c(genes.to.label1,genes.to.label2) )
#genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label3, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_CTBSTBEVT_TF_seg3_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label4, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_CTBSTBEVT_TF_seg4_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_CTBSTBEVT_TF_seg5_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label6, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_CTBSTBEVT_TF_seg6_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label7, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_CTBSTBEVT_TF_seg7_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label8, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_CTBSTBEVT_TF_seg8_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Cyno1 Huma1 Am Tr

#M35 <- FindMarkers(mammal.combined, ident.2=c("Cyno (Lanner)_CTB"),ident.1 = c("Cyno (Lanner)_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M36 <- FindMarkers(mammal.combined, ident.2=c("HumanCS5-7_CTB"),ident.1 = c("HumanCS5-7_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M37 <- FindMarkers(mammal.combined, ident.2=c("Cyno (in vitro)_CTB"),ident.1 = c("Cyno (in vitro)_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M38 <- FindMarkers(mammal.combined, ident.2=c("Marmoset_CTB"),ident.1 = c("Marmoset_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#write.csv(as.data.frame(M34), file=paste(saveext,"/STB_CTB_Hum2.csv",sep=""))
#write.csv(as.data.frame(M35), file=paste(saveext,"/STB_CTB_Cyno2.csv",sep=""))
#write.csv(as.data.frame(M36), file=paste(saveext,"/STB_CTB_Hum1.csv",sep=""))
#write.csv(as.data.frame(M37), file=paste(saveext,"/STB_CTB_Cyno1.csv",sep=""))
#write.csv(as.data.frame(M38), file=paste(saveext,"/STB_CTB_Marm.csv",sep=""))

#M39 <- FindMarkers(mammal.combined, ident.2=c("Human (Ara EBs)_CTB"),ident.1=c("Human (Ara EBs)_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M40 <- FindMarkers(mammal.combined, ident.2=c("Cyno (Lanner)_CTB"),ident.1 = c("Cyno (Lanner)_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M41 <- FindMarkers(mammal.combined, ident.2=c("HumanCS5-7_CTB"),ident.1 = c("HumanCS5-7_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M42 <- FindMarkers(mammal.combined, ident.2=c("Cyno (in vitro)_CTB"),ident.1 = c("Cyno (in vitro)_EVT"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))#M35 <- FindMarkers(mammal.combined, ident.2=c("Cyno (Lanner)_CTB"),ident.1 = c("Cyno (Lanner)_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Cl1 <- M35
Cl2 <- M40
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$'Cyno (Lanner)_CTB'>0.1 | Ae5$'Cyno (Lanner)_EVT'>0.1 | Ae5$'Cyno (Lanner)_STB'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))]),TF)
genes.to.label2 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF)
genes.to.label3 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)> log(1.2))]),TF),Ae5$gene)
genes.to.label4 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)> log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF),Ae5$gene)
genes.to.label5 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label6 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label7 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)
genes.to.label8 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)

saveRDS(genes.to.label1,file=paste(saveext,"Cyno1_CTBSTBEVT_seg1.rds",sep=""))
saveRDS(genes.to.label2,file=paste(saveext,"Cyno1_CTBSTBEVT_seg2.rds",sep=""))
saveRDS(genes.to.label3,file=paste(saveext,"Cyno1_CTBSTBEVT_seg3.rds",sep=""))
saveRDS(genes.to.label4,file=paste(saveext,"Cyno1_CTBSTBEVT_seg4.rds",sep=""))
saveRDS(genes.to.label5,file=paste(saveext,"Cyno1_CTBSTBEVT_seg5.rds",sep=""))
saveRDS(genes.to.label6,file=paste(saveext,"Cyno1_CTBSTBEVT_seg6.rds",sep=""))
saveRDS(genes.to.label7,file=paste(saveext,"Cyno1_CTBSTBEVT_seg7.rds",sep=""))
saveRDS(genes.to.label8,file=paste(saveext,"Cyno1_CTBSTBEVT_seg8.rds",sep=""))

#genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
#intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_CTBSTBEVT_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_CTBSTBEVT_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
#genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
#intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
#genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
#intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
#genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
#genes.to.label3 = setdiff(genes.to.label5,c(genes.to.label1,genes.to.label2) )
#genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label3, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_CTBSTBEVT_TF_seg3_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label4, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_CTBSTBEVT_TF_seg4_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_CTBSTBEVT_TF_seg5_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label6, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_CTBSTBEVT_TF_seg6_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label7, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_CTBSTBEVT_TF_seg7_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label8, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_CTBSTBEVT_TF_seg8_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- M37
Cl2 <- M42
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$'Cyno (in vitro)_CTB'>0.1 | Ae5$'Cyno (in vitro)_CTB'>0.1 | Ae5$'Cyno (in vitro)_CTB'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
#genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),SIGNAL$V1))
genes.to.label1 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))]),TF)
genes.to.label2 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF)
genes.to.label3 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)> log(1.2))]),TF),Ae5$gene)
genes.to.label4 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)> log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF),Ae5$gene)


genes.to.label5 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label6 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label7 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)
genes.to.label8 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)


saveRDS(genes.to.label1,file=paste(saveext,"Cyno2_CTBSTBEVT_seg1.rds",sep=""))
saveRDS(genes.to.label2,file=paste(saveext,"Cyno2_CTBSTBEVT_seg2.rds",sep=""))
saveRDS(genes.to.label3,file=paste(saveext,"Cyno2_CTBSTBEVT_seg3.rds",sep=""))
saveRDS(genes.to.label4,file=paste(saveext,"Cyno2_CTBSTBEVT_seg4.rds",sep=""))
saveRDS(genes.to.label5,file=paste(saveext,"Cyno2_CTBSTBEVT_seg5.rds",sep=""))
saveRDS(genes.to.label6,file=paste(saveext,"Cyno2_CTBSTBEVT_seg6.rds",sep=""))
saveRDS(genes.to.label7,file=paste(saveext,"Cyno2_CTBSTBEVT_seg7.rds",sep=""))
saveRDS(genes.to.label8,file=paste(saveext,"Cyno2_CTBSTBEVT_seg8.rds",sep=""))

#genes.to.label2 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
#intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_CTBSTBEVT_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_CTBSTBEVT_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
#genes.to.label3 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]),TF),
#intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]),TF))
#genes.to.label4 = c(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.5))]), SIGNAL$V1),
#intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.5))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.5))]), SIGNAL$V1))
#genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),TF))
#genes.to.label6 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.5))]),SIGNAL$V1))
#genes.to.label3 = setdiff(genes.to.label5,c(genes.to.label1,genes.to.label2) )
#genes.to.label2 = setdiff(genes.to.label6,genes.to.label4)
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label3, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_CTBSTBEVT_TF_seg3_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
#p1 <- LabelPoints(plot = p1, points = genes.to.label4, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_CTBSTBEVT_TF_seg4_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_CTBSTBEVT_TF_seg5_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label6, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_CTBSTBEVT_TF_seg6_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label7, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_CTBSTBEVT_TF_seg7_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label8, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_CTBSTBEVT_TF_seg8_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
#Cyno2_CTBSTBEVT_TF_seg3_fin.pdf


#M11 <- FindMarkers(mammal.combined, ident.2=c("Human (Ara EBs)_Am"),ident.1=c("Human (Ara EBs)_CTB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M12 <- FindMarkers(mammal.combined, ident.2=c("Cyno (Lanner)_Am"),ident.1 = c("Cyno (Lanner)_CTB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M13 <- FindMarkers(mammal.combined, ident.2=c("HumanCS5-7_Am"),ident.1 = c("HumanCS5-7_CTB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M14 <- FindMarkers(mammal.combined, ident.2=c("Cyno (in vitro)_Am"),ident.1 = c("Cyno (in vitro)_CTB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M15 <- FindMarkers(mammal.combined, ident.2=c("Marmoset_Am"),ident.1 = c("Marmoset_CTB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#Huma2 Am CTB STB
Cl1 <- M11
Cl2 <- M34
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$'Human (Ara EBs)_Am'>0.1 | Ae5$'Human (Ara EBs)_STB'>0.1 | Ae5$'Human (Ara EBs)_CTB'>0.1)  ) ,]
genes.to.label1 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))]),TF),Ae5$gene)
genes.to.label2 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF),Ae5$gene)
genes.to.label3 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)> log(1.2))]),TF),Ae5$gene)
genes.to.label4 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)> log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF),Ae5$gene)
genes.to.label5 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label6 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label7 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)
genes.to.label8 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)

saveRDS(genes.to.label1,file=paste(saveext,"Human2_AmCTBSTB_seg1.rds",sep=""))
saveRDS(genes.to.label2,file=paste(saveext,"Human2_AmCTBSTB_seg2.rds",sep=""))
saveRDS(genes.to.label3,file=paste(saveext,"Human2_AmCTBSTB_seg3.rds",sep=""))
saveRDS(genes.to.label4,file=paste(saveext,"Human2_AmCTBSTB_seg4.rds",sep=""))
saveRDS(genes.to.label5,file=paste(saveext,"Human2_AmCTBSTB_seg5.rds",sep=""))
saveRDS(genes.to.label6,file=paste(saveext,"Human2_AmCTBSTB_seg6.rds",sep=""))
saveRDS(genes.to.label7,file=paste(saveext,"Human2_AmCTBSTB_seg7.rds",sep=""))
saveRDS(genes.to.label8,file=paste(saveext,"Human2_AmCTBSTB_seg8.rds",sep=""))

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_AmCTBSTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_AmCTBSTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label3, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_AmCTBSTB_TF_seg3_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label4, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_AmCTBSTB_TF_seg4_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_AmCTBSTB_TF_seg5_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label6, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_AmCTBSTB_TF_seg6_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label7, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_AmCTBSTB_TF_seg7_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label8, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_AmCTBSTB_TF_seg8_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

#Cyno2_CTBSTBEVT_TF_seg3_fin.pdf

Cl1 <- M11
Cl2 <- M39
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$'Human (Ara EBs)_Am'>0.1 | Ae5$'Human (Ara EBs)_EVT'>0.1 | Ae5$'Human (Ara EBs)_CTB'>0.1)  ) ,]
genes.to.label1 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))]),TF),Ae5$gene)
genes.to.label2 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF),Ae5$gene)
genes.to.label3 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)> log(1.2))]),TF),Ae5$gene)
genes.to.label4 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)> log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF),Ae5$gene)
genes.to.label5 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label6 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label7 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)
genes.to.label8 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)

saveRDS(genes.to.label1,file=paste(saveext,"Human2_AmCTBEVT_seg1.rds",sep=""))
saveRDS(genes.to.label2,file=paste(saveext,"Human2_AmCTBEVT_seg2.rds",sep=""))
saveRDS(genes.to.label3,file=paste(saveext,"Human2_AmCTBEVT_seg3.rds",sep=""))
saveRDS(genes.to.label4,file=paste(saveext,"Human2_AmCTBEVT_seg4.rds",sep=""))
saveRDS(genes.to.label5,file=paste(saveext,"Human2_AmCTBEVT_seg5.rds",sep=""))
saveRDS(genes.to.label6,file=paste(saveext,"Human2_AmCTBEVT_seg6.rds",sep=""))
saveRDS(genes.to.label7,file=paste(saveext,"Human2_AmCTBEVT_seg7.rds",sep=""))
saveRDS(genes.to.label8,file=paste(saveext,"Human2_AmCTBEVT_seg8.rds",sep=""))


p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_AmCTBEVT_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_AmCTBEVT_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label3, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_AmCTBEVT_TF_seg3_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label4, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_AmCTBEVT_TF_seg4_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_AmCTBEVT_TF_seg5_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label6, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_AmCTBEVVT_TF_seg6_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label7, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_AmCTBEVT_TF_seg7_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label8, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human2_AmCTBEVT_TF_seg8_fin.pdf",sep=""),width = 13, height = 13, plot = p1)


#Cyno1 Huma1 Am Tr
Cl1 <- M13
Cl2 <- M36
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$'HumanCS5-7_Am'>0.1 | Ae5$'HumanCS5-7_CTB'>0.1 | Ae5$'Cyno (in vitro)_STB'>0.1)  ) ,]
genes.to.label1 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))]),TF)
genes.to.label2 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF)
genes.to.label3 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)> log(1.2))]),TF),Ae5$gene)
genes.to.label4 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)> log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF),Ae5$gene)
genes.to.label5 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label6 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label7 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)
genes.to.label8 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)


saveRDS(genes.to.label1,file=paste(saveext,"Human1_AmCTBSTB_seg1.rds",sep=""))
saveRDS(genes.to.label2,file=paste(saveext,"Human1_AmCTBSTB_seg2.rds",sep=""))
saveRDS(genes.to.label3,file=paste(saveext,"Human1_AmCTBSTB_seg3.rds",sep=""))
saveRDS(genes.to.label4,file=paste(saveext,"Human1_AmCTBSTB_seg4.rds",sep=""))
saveRDS(genes.to.label5,file=paste(saveext,"Human1_AmCTBSTB_seg5.rds",sep=""))
saveRDS(genes.to.label6,file=paste(saveext,"Human1_AmCTBSTB_seg6.rds",sep=""))
saveRDS(genes.to.label7,file=paste(saveext,"Human1_AmCTBSTB_seg7.rds",sep=""))
saveRDS(genes.to.label8,file=paste(saveext,"Human1_AmCTBSTB_seg8.rds",sep=""))

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_AmCTBSTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_AmCTBSTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label3, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_AmCTBSTB_TF_seg3_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label4, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_AmCTBSTB_TF_seg4_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_AmCTBSTB_TF_seg5_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label6, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_AmCTBSTB_TF_seg6_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label7, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_AmCTBSTB_TF_seg7_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label8, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_AmCTBSTB_TF_seg8_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


#Cyno1 Huma1 Am Tr
Cl1 <- M13
Cl2 <- M41
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$'HumanCS5-7_Am'>0.1 | Ae5$'HumanCS5-7_CTB'>0.1 | Ae5$'Cyno (in vitro)_EVT'>0.1)  ) ,]
genes.to.label1 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))]),TF)
genes.to.label2 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF)
genes.to.label3 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)> log(1.2))]),TF),Ae5$gene)
genes.to.label4 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)> log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF),Ae5$gene)
genes.to.label5 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label6 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label7 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)
genes.to.label8 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)


saveRDS(genes.to.label1,file=paste(saveext,"Human1_AmCTBEVT_seg1.rds",sep=""))
saveRDS(genes.to.label2,file=paste(saveext,"Human1_AmCTBEVT_seg2.rds",sep=""))
saveRDS(genes.to.label3,file=paste(saveext,"Human1_AmCTBEVT_seg3.rds",sep=""))
saveRDS(genes.to.label4,file=paste(saveext,"Human1_AmCTBEVT_seg4.rds",sep=""))
saveRDS(genes.to.label5,file=paste(saveext,"Human1_AmCTBEVT_seg5.rds",sep=""))
saveRDS(genes.to.label6,file=paste(saveext,"Human1_AmCTBEVT_seg6.rds",sep=""))
saveRDS(genes.to.label7,file=paste(saveext,"Human1_AmCTBEVT_seg7.rds",sep=""))
saveRDS(genes.to.label8,file=paste(saveext,"Human1_AmCTBEVT_seg8.rds",sep=""))

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_AmCTBEVT_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_AmCTBEVT_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label3, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_AmCTBEVT_TF_seg3_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label4, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_AmCTBEVT_TF_seg4_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_AmCTBEVT_TF_seg5_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label6, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_AmCTBEVT_TF_seg6_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label7, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_AmCTBEVT_TF_seg7_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label8, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Human1_AmCTBEVT_TF_seg8_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



Cl1 <- M12
Cl2 <- M35
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$'Cyno (Lanner)_CTB'>0.1 | Ae5$'Cyno (Lanner)_Am'>0.1 | Ae5$'Cyno (Lanner)_STB'>0.1)  ) ,]
genes.to.label1 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))]),TF)
genes.to.label2 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF)
genes.to.label3 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)> log(1.2))]),TF),Ae5$gene)
genes.to.label4 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)> log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF),Ae5$gene)
genes.to.label5 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label6 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label7 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)
genes.to.label8 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)

saveRDS(genes.to.label1,file=paste(saveext,"Cyno1_AmCTBSTB_seg1.rds",sep=""))
saveRDS(genes.to.label2,file=paste(saveext,"Cyno1_AmCTBSTB_seg2.rds",sep=""))
saveRDS(genes.to.label3,file=paste(saveext,"Cyno1_AmCTBSTB_seg3.rds",sep=""))
saveRDS(genes.to.label4,file=paste(saveext,"Cyno1_AmCTBSTB_seg4.rds",sep=""))
saveRDS(genes.to.label5,file=paste(saveext,"Cyno1_AmCTBSTB_seg5.rds",sep=""))
saveRDS(genes.to.label6,file=paste(saveext,"Cyno1_AmCTBSTB_seg6.rds",sep=""))
saveRDS(genes.to.label7,file=paste(saveext,"Cyno1_AmCTBSTB_seg7.rds",sep=""))
saveRDS(genes.to.label8,file=paste(saveext,"Cyno1_AmCTBSTB_seg8.rds",sep=""))

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_AmCTBSTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_AmCTBSTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label3, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_AmCTBSTB_TF_seg3_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label4, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_AmCTBSTB_TF_seg4_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_AmCTBSTB_TF_seg5_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label6, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_AmCTBSTB_TF_seg6_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label7, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_AmCTBSTB_TF_seg7_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label8, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_AmCTBSTB_TF_seg8_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- M12
Cl2 <- M40
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$'Cyno (Lanner)_CTB'>0.1 | Ae5$'Cyno (Lanner)_EVT'>0.1 | Ae5$'Cyno (Lanner)_Am'>0.1)  ) ,]
genes.to.label1 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))]),TF)
genes.to.label2 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF)
genes.to.label3 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)> log(1.2))]),TF),Ae5$gene)
genes.to.label4 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)> log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF),Ae5$gene)
genes.to.label5 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label6 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label7 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)
genes.to.label8 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)


saveRDS(genes.to.label1,file=paste(saveext,"Cyno1_AmCTBEVT_seg1.rds",sep=""))
saveRDS(genes.to.label2,file=paste(saveext,"Cyno1_AmCTBEVT_seg2.rds",sep=""))
saveRDS(genes.to.label3,file=paste(saveext,"Cyno1_AmCTBEVT_seg3.rds",sep=""))
saveRDS(genes.to.label4,file=paste(saveext,"Cyno1_AmCTBEVT_seg4.rds",sep=""))
saveRDS(genes.to.label5,file=paste(saveext,"Cyno1_AmCTBEVT_seg5.rds",sep=""))
saveRDS(genes.to.label6,file=paste(saveext,"Cyno1_AmCTBEVT_seg6.rds",sep=""))
saveRDS(genes.to.label7,file=paste(saveext,"Cyno1_AmCTBEVT_seg7.rds",sep=""))
saveRDS(genes.to.label8,file=paste(saveext,"Cyno1_AmCTBEVT_seg8.rds",sep=""))


p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_AmCTBEVT_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_AmCTBEVT_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label3, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_AmCTBEVT_TF_seg3_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label4, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_AmCTBEVT_TF_seg4_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_AmCTBEVT_TF_seg5_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label6, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_AmCTBEVT_TF_seg6_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label7, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_AmCTBEVT_TF_seg7_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label8, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno1_AmCTBEVT_TF_seg8_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- M14
Cl2 <- M37
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$'Cyno (in vitro)_STB'>0.1 | Ae5$'Cyno (in vitro)_CTB'>0.1 | Ae5$'Cyno (in vitro)_Am'>0.1)  ) ,]
genes.to.label1 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))]),TF)
genes.to.label2 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF)
genes.to.label3 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)> log(1.2))]),TF),Ae5$gene)
genes.to.label4 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)> log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF),Ae5$gene)
genes.to.label5 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label6 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label7 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)
genes.to.label8 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)

saveRDS(genes.to.label1,file=paste(saveext,"Cyno2_AmCTBSTB_seg1.rds",sep=""))
saveRDS(genes.to.label2,file=paste(saveext,"Cyno2_AmCTBSTB_seg2.rds",sep=""))
saveRDS(genes.to.label3,file=paste(saveext,"Cyno2_AmCTBSTB_seg3.rds",sep=""))
saveRDS(genes.to.label4,file=paste(saveext,"Cyno2_AmCTBSTB_seg4.rds",sep=""))
saveRDS(genes.to.label5,file=paste(saveext,"Cyno2_AmCTBSTB_seg5.rds",sep=""))
saveRDS(genes.to.label6,file=paste(saveext,"Cyno2_AmCTBSTB_seg6.rds",sep=""))
saveRDS(genes.to.label7,file=paste(saveext,"Cyno2_AmCTBSTB_seg7.rds",sep=""))
saveRDS(genes.to.label8,file=paste(saveext,"Cyno2_AmCTBSTB_seg8.rds",sep=""))


p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_AmCTBSTB_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_AmCTBSTB_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label3, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_AmCTBSTB_TF_seg3_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label4, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_AmCTBSTB_TF_seg4_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_AmCTBSTB_TF_seg5_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label6, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_AmCTBSTB_TF_seg6_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label7, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_AmCTBSTB_TF_seg7_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label8, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_AmCTBSTB_TF_seg8_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- M14
Cl2 <- M42
Ae5 <- Ae
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_logFC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
Ae5 <- Ae5[which((Ae5$'Cyno (in vitro)_CTB'>0.1 | Ae5$'Cyno (in vitro)_EVT'>0.1 | Ae5$'Cyno (in vitro)_Am'>0.1)  ) ,]
genes.to.label1 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))]),TF)
genes.to.label2 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF)
genes.to.label3 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)> log(1.2))]),TF),Ae5$gene)
genes.to.label4 = intersect(intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)> log(1.2))],rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))]),TF),Ae5$gene)
genes.to.label5 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)>log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label6 = intersect(intersect(rownames(Cl1)[which((Cl1$avg_logFC)< -log(1.2))],rownames(Cl2)[which(abs(Cl2$avg_logFC)<log(1.2))]),TF)
genes.to.label7 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)>log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)
genes.to.label8 = intersect(intersect(rownames(Cl2)[which((Cl2$avg_logFC)< -log(1.2))],rownames(Cl1)[which(abs(Cl1$avg_logFC)<log(1.2))]),TF)

saveRDS(genes.to.label1,file=paste(saveext,"Cyno2_AmCTBEVT_seg1.rds",sep=""))
saveRDS(genes.to.label2,file=paste(saveext,"Cyno2_AmCTBEVT_seg2.rds",sep=""))
saveRDS(genes.to.label3,file=paste(saveext,"Cyno2_AmCTBEVT_seg3.rds",sep=""))
saveRDS(genes.to.label4,file=paste(saveext,"Cyno2_AmCTBEVT_seg4.rds",sep=""))
saveRDS(genes.to.label5,file=paste(saveext,"Cyno2_AmCTBEVT_seg5.rds",sep=""))
saveRDS(genes.to.label6,file=paste(saveext,"Cyno2_AmCTBEVT_seg6.rds",sep=""))
saveRDS(genes.to.label7,file=paste(saveext,"Cyno2_AmCTBEVT_seg7.rds",sep=""))
saveRDS(genes.to.label8,file=paste(saveext,"Cyno2_AmCTBEVT_seg8.rds",sep=""))


p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_AmCTBEVT_TF_seg1_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_AmCTBEVT_TF_seg2_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label3, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_AmCTBEVT_TF_seg3_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label4, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_AmCTBEVT_TF_seg4_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_AmCTBEVT_TF_seg5_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label6, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_AmCTBEVT_TF_seg6_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label7, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_AmCTBEVT_TF_seg7_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label8, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression datasete1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cyno2_AmCTBEVT_TF_seg8_fin.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

