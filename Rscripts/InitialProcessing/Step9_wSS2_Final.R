library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library("Matrix")

set.seed(1)
#This reruns the previous analysis extracting out the embryonic lineags

#Save folder
saveext = "./FinalAlign/"
dir.create(saveext)
dir.create(paste(saveext,"DimRed/",sep=""))
dir.create(paste(saveext,"Markers/",sep=""))

mammal.combined <- readRDS(paste(saveext,"Seurat_combined3.rds",sep=""))

D1 <- readRDS("../Trasferv6/Mole_experimentRv6_1.rds")
D2 <- readRDS("../Trasferv6/Mole_experimentRv6_2.rds")
D3 <- readRDS("../Trasferv6/Mole_experimentRv6_3.rds")
D3B <- readRDS("../Trasferv6/Mole_experimentRv6_3B.rds")
D6 <- readRDS("../Trasferv6/Mole_experimentRv6_6.rds")
D7 <- readRDS("../Trasferv6/Mole_experimentRv6_7.rds")
D8 <- readRDS("../Trasferv6/Mole_experimentRv6_8.rds")

#Load SS2 data
SSD <- readRDS('../Trasferv6/MergedDataSS2.rds') #Smart Seq
#Simplify anotations for SS2 data
uID <- as.character(Idents(SSD))
uID[which(uID%in%c("CTb_CS4","EPI_CS4","TB_CS3"))] <- "Tr"
uID[which(uID%in%c("CTb_CS5A/B","CTb_CS5C","CTb_CS6"))] <- "CTB"
uID[which(uID%in%c("EmDisc_CS5A/B","EmDisc_CS5C","EmDisc_CS6","EmDiscPS_CS6"))] <- "EmDisc"
uID[which(uID%in%c("EVTb_CS5C","EVTb_CS6"))] <- "EVT"
uID[which(uID%in%c("STb_CS5A/B","STb_CS5C","STb_CS6"))] <- "STB"
uID[which(uID%in%c("VE_CS4","VE_CS5","VE_CS5A/B","VE_CS5C","SYS_CS6"))] <- "VE"
Idents(SSD) <- uID
SSD <- subset(SSD,idents="Unknown",invert=TRUE)
SSD <- subset(SSD,idents=c("Tr","CTB","EmDisc","EVT","STB","VE","ICM_CS3"))
SSD$ID3 <- "D1"
SSD$species1 <- "SS2"

#Load genotype anotations
BC1 <- read.table("Genotype/D1_4Cl.tsv",sep="\t",header=TRUE)
BC2 <- read.table("Genotype/D2_4Cl.tsv",sep="\t",header=TRUE)
BC3 <- read.table("Genotype/D3_6Cl.tsv",sep="\t",header=TRUE)
BC3B <- read.table("Genotype/D3B_6Cl.tsv",sep="\t",header=TRUE)
#BC4 <- read.table("Genotype/D5_cl.tsv",sep="\t",header=TRUE)
BC6 <- read.table("Genotype/D6_4Cl.tsv",sep="\t",header=TRUE)
BC7 <- read.table("Genotype/D7_4Cl.tsv",sep="\t",header=TRUE)
BC8 <- read.table("Genotype/D8_4Cl.tsv",sep="\t",header=TRUE)

d1 <- as.data.frame(colnames(D1))
d2 <- as.data.frame(colnames(D2))
d3 <- as.data.frame(colnames(D3))
d3B <- as.data.frame(colnames(D3B))
d6 <- as.data.frame(colnames(D6))
d7 <- as.data.frame(colnames(D7))
d8 <- as.data.frame(colnames(D8))

colnames(d1) <- "barcode"
colnames(d2) <- "barcode"
colnames(d3) <- "barcode"
colnames(d3B) <- "barcode"
colnames(d6) <- "barcode"
colnames(d7) <- "barcode"
colnames(d8) <- "barcode"

m1 <- merge(x = d1, y = BC1, by = "barcode", all = TRUE)
m2 <- merge(x = d2, y = BC2, by = "barcode", all = TRUE)
m3 <- merge(x = d3, y = BC3, by = "barcode", all = TRUE)
m3B <- merge(x = d3B, y = BC3B, by = "barcode", all = TRUE)
m6 <- merge(x = d6, y = BC6, by = "barcode", all = TRUE)
m7 <- merge(x = d7, y = BC7, by = "barcode", all = TRUE)
m8 <- merge(x = d8, y = BC8, by = "barcode", all = TRUE)

#Load in the embryonic data
Emb <- readRDS("../Trasferv6/MagdaHumanDataHarmonyhg38_labelled.rds")
DefaultAssay(Emb) <- "RNA"
#Load inn the endometrial data
#Dep1 <- readRDS("../../Endometrial/Endom1.rds")
#Dep2 <- readRDS("../../Endometrial/Endom2.rds")
#Dep3 <- readRDS("../../Endometrial/Endom3.rds")
#Dep4 <- readRDS("../../Endometrial/Endom4.rds")
#Dep5 <- readRDS("../../Endometrial/Endom5.rds")
#Dep6 <- readRDS("../../Endometrial/Endom6.rds")

#Assign endometrial data a matching "batch" to the embryo (there are 6 in total). This is somewhat arbitrary.
#Dep1$ID3 <- "B1"
#Dep2$ID3 <- "B2"
#Dep3$ID3 <- "B3"
#Dep4$ID3 <- "B4"
#Dep5$ID3 <- "B5"#
#Dep6$ID3 <- "B6"
#Merge data togetehr
#D <- merge(Emb, y = c(Dep1,Dep2,Dep3,Dep4,Dep5,Dep6), project = "merged1")
D <- Emb #merge(Emb, y = c(Dep1,Dep2,Dep3,Dep4,Dep5,Dep6), project = "merged1")
#Save cell type anotation to a new slot
D$Cells <- Idents(D)

#Set identity to the "batch" and split
Idents(D) <- D$ID3
MergedData1 <- subset(D,idents="B1")
MergedData2 <- subset(D,idents="B2")
MergedData3 <- subset(D,idents="B3")
MergedData4 <- subset(D,idents="B4")
MergedData5 <- subset(D,idents="B5")
MergedData6 <- subset(D,idents="B6")

#Reset the cell anotations
Idents(MergedData1) <- MergedData1$Cells
Idents(MergedData2) <- MergedData2$Cells
Idents(MergedData3) <- MergedData3$Cells
Idents(MergedData4) <- MergedData4$Cells
Idents(MergedData5) <- MergedData5$Cells
Idents(MergedData6) <- MergedData6$Cells

#Save individual files
#saveRDS(MergedData1,file="Merged1.rds")
#saveRDS(MergedData2,file="Merged2.rds")
#saveRDS(MergedData3,file="Merged3.rds")
#saveRDS(MergedData4,file="Merged4.rds")
#saveRDS(MergedData5,file="Merged5.rds")
#saveRDS(MergedData6,file="Merged6.rds")

#Give old dataset a unique ID
MergedData1$species1 <- "Old"
MergedData2$species1 <- "Old"
MergedData3$species1 <- "Old"
MergedData4$species1 <- "Old"
MergedData5$species1 <- "Old"
MergedData6$species1 <- "Old"

#Load in new data
#D1 <- readRDS("../Mole_experimentRv5_1.rds")
Idents(D1) <- "NAssigned"
Idents(D1,cells=BC1$barcode[which(BC1$assignment==0)]) <- "Stromal_G1"
Idents(D1,cells=BC1$barcode[which(BC1$assignment==1)]) <- "Epithelial_G2"
Idents(D1,cells=BC1$barcode[which(BC1$assignment==2)]) <- "Embryoni_G3"
Idents(D1,cells=BC1$barcode[which(BC1$assignment==3)]) <- "Epithelial_G4"
D1$PCl1 <- m1[,"cluster3"]
D1$PCl2 <- m1[,"cluster0"]
D1$PCl3 <- m1[,"cluster1"]
D1$PCl4 <- m1[,"cluster2"]
#D2 <- readRDS("../Mole_experimentRv5_2.rds")

Idents(D2) <- "NAssigned"
Idents(D2,cells=BC2$barcode[which(BC2$assignment==0)]) <- "Epithelial_G1"
Idents(D2,cells=BC2$barcode[which(BC2$assignment==1)]) <- "Epithelial_G2"
Idents(D2,cells=BC2$barcode[which(BC2$assignment==2)]) <- "Stromal_G3"
Idents(D2,cells=BC2$barcode[which(BC2$assignment==3)]) <- "Embryonic_G4"
D2$PCl1 <- m2[,"cluster3"]
D2$PCl2 <- m2[,"cluster2"]
D2$PCl3 <- m2[,"cluster0"]
D2$PCl4 <- m2[,"cluster1"]

#D3 <- readRDS("../Mole_experimentRv5_3.rds")
Idents(D3) <- "NAssigned"
Idents(D3,cells=BC3$barcode[which(BC3$assignment==0)]) <- "Epithelial_G1"
Idents(D3,cells=BC3$barcode[which(BC3$assignment==1)]) <- "Embryonic_G2"
#Idents(D3,cells=BC3$barcode[which(BC3$assignment==2)]) <- "G3"
Idents(D3,cells=BC3$barcode[which(BC3$assignment==3)]) <- "Stromal_G4"
Idents(D3,cells=BC3$barcode[which(BC3$assignment==4)]) <- "Embryonic_G5"
Idents(D3,cells=BC3$barcode[which(BC3$assignment==5)]) <- "Epithelial_G6"
D3$PCl1 <- m3[,"cluster1"]
D3$PCl2 <- m3[,"cluster4"]
D3$PCl3 <- m3[,"cluster3"]
D3$PCl4 <- m3[,"cluster0"]
D3$PCl5 <- m3[,"cluster5"]
D3$PCl6 <- m3[,"cluster2"]

#D3B <- readRDS("../Mole_experimentRv5_3B.rds")
Idents(D3B) <- "NAssigned"
Idents(D3B,cells=BC3B$barcode[which(BC3B$assignment==0)]) <- "Epithelial_G1"
Idents(D3B,cells=BC3B$barcode[which(BC3B$assignment==1)]) <- "Epithelial_G2"
Idents(D3B,cells=BC3B$barcode[which(BC3B$assignment==2)]) <- "Stromal_G3"
Idents(D3B,cells=BC3B$barcode[which(BC3B$assignment==3)]) <- "Embryonic_G4"
Idents(D3B,cells=BC3B$barcode[which(BC3B$assignment==4)]) <- "Embryonic_G5"
Idents(D3B,cells=BC3B$barcode[which(BC3B$assignment==5)]) <- "Epithelial_G6"

D3B$PCl1 <- m3B[,"cluster3"]
D3B$PCl2 <- m3B[,"cluster4"]
D3B$PCl3 <- m3B[,"cluster2"]
D3B$PCl4 <- m3B[,"cluster0"]
D3B$PCl5 <- m3B[,"cluster1"]
D3B$PCl6 <- m3B[,"cluster5"]

#D6 <- readRDS("../Mole_experimentRv5_6.rds")
Idents(D6) <- "NAssigned"
Idents(D6,cells=BC6$barcode[which(BC6$assignment==0)]) <- "Epithelial_G1"
Idents(D6,cells=BC6$barcode[which(BC6$assignment==1)]) <- "Embryonic_G2"
Idents(D6,cells=BC6$barcode[which(BC6$assignment==2)]) <- "Epithelial_G3"
Idents(D6,cells=BC6$barcode[which(BC6$assignment==3)]) <- "Stromal_G4"

D6$PCl1 <- m6[,"cluster2"]
D6$PCl2 <- m6[,"cluster3"]
D6$PCl3 <- m6[,"cluster0"]
D6$PCl4 <- m6[,"cluster1"]

D7 <- readRDS("../Mole_experimentRv5_7.rds")
Idents(D7) <- "NAssigned"
Idents(D7,cells=BC7$barcode[which(BC7$assignment==0)]) <- "Stromal_G1"
Idents(D7,cells=BC7$barcode[which(BC7$assignment==1)]) <- "Epithelial_G2"
Idents(D7,cells=BC7$barcode[which(BC7$assignment==2)]) <- "Epithelial_G3"
Idents(D7,cells=BC7$barcode[which(BC7$assignment==3)]) <- "Embryonic_G4"


D7$PCl1 <- m7[,"cluster3"]
D7$PCl2 <- m7[,"cluster0"]
D7$PCl3 <- m7[,"cluster1"]
D7$PCl4 <- m7[,"cluster2"]

D8 <- readRDS("../Mole_experimentRv5_8.rds")
Idents(D8) <- "NAssigned"
Idents(D8,cells=BC8$barcode[which(BC8$assignment==0)]) <- "G1"
Idents(D8,cells=BC8$barcode[which(BC8$assignment==1)]) <- "Stromal_G2"
Idents(D8,cells=BC8$barcode[which(BC8$assignment==2)]) <- "G3"
Idents(D8,cells=BC8$barcode[which(BC8$assignment==3)]) <- "G4"


#D8$PCl1 <- m8[,"cluster0"]
D8$PCl2 <- m8[,"cluster1"]
D8$PCl3 <- m8[,"cluster0"]
D8$PCl4 <- m8[,"cluster2"]
D8$PCl5 <- m8[,"cluster3"]


#Give new dataset an ID
D1$species1 <- "New"
D2$species1 <- "New"
D3$species1 <- "New"
D3B$species1 <- "New"
D6$species1 <- "New"
D7$species1 <- "New"
D8$species1 <- "New"

D1$ID3 <- "C1"
D2$ID3 <- "C2"
D3$ID3 <- "C3"
D3B$ID3 <- "C4"
D6$ID3 <- "C6"
D7$ID3 <- "C7"
D8$ID3 <- "C8"

D1$Genotype <- Idents(D1)
D2$Genotype <- Idents(D2)
D3$Genotype <- Idents(D3)
D3B$Genotype <- Idents(D3B)
D6$Genotype <- Idents(D6)
D7$Genotype <- Idents(D7)
D8$Genotype <- Idents(D8)

D1$Dataset <- "10X Ours"
D2$Dataset <- "10X Ours"
D3$Dataset <- "10X Ours"
D3B$Dataset <- "10X Ours"
D6$Dataset <- "10X Ours"
D7$Dataset <- "10X Ours"
D8$Dataset <- "10X Ours"

# [1] "D1" "B1" "B2" "B3" "B4" "B5" "B6" "C1" "C2" "C3" "C4" "C6" "C7" "C8"
#> unique(mammal.combined$Dataset)
#[1] "SS2 Reference 1" "10X Reference 1" "10X Ours"       
#> unique(mammal.combined$Tissue)
#[1] "Embryonic"  "VE"         "Tr"         "Maternal"   "Hyp_d12"   
#[6] "EmDisc_d12" "Cil"       
#> unique(mammal.combined$Tissue2)
#[1] "Tb"         "Embryonic"  "VE"         "Tr"         "Epithelial"
#[6] "Stromal"    "Hyp_d12"    "EmDisc_d12" "Cil"     

Idents(D1,cells=gsub("-1_8","-1",names(which(mammal.combined$Tissue2%in%c("EmDisc_d12","Embryonic","VE","Hyp_d12") & mammal.combined$ID3=="C1"))) ) <- "EmDisc"
Idents(D2,cells=gsub("-1_9","-1",names(which(mammal.combined$Tissue2%in%c("EmDisc_d12","Embryonic","VE","Hyp_d12") & mammal.combined$ID3=="C2"))) ) <- "EmDisc"
#Idents(D3,cells=gsub("-1_10","-1",names(which(mammal.combined$Tissue2%in%c("EmDisc_d12","Embryonic","VE","Hyp_d12") & mammal.combined$ID3=="C3"))) ) <- "EmDisc"
Idents(D3B,cells=gsub("-1_11","-1",names(which(mammal.combined$Tissue2%in%c("EmDisc_d12","Embryonic","VE","Hyp_d12") & mammal.combined$ID3=="C4"))) ) <- "EmDisc"
Idents(D6,cells=gsub("-1_12","-1",names(which(mammal.combined$Tissue2%in%c("EmDisc_d12","Embryonic","VE","Hyp_d12") & mammal.combined$ID3=="C6"))) ) <- "EmDisc"
Idents(D7,cells=gsub("-1_13","-1",names(which(mammal.combined$Tissue2%in%c("EmDisc_d12","Embryonic","VE","Hyp_d12") & mammal.combined$ID3=="C7"))) ) <- "EmDisc"
Idents(D8,cells=gsub("-1_14","-1",names(which(mammal.combined$Tissue2%in%c("EmDisc_d12","Embryonic","VE","Hyp_d12") & mammal.combined$ID3=="C8"))) ) <- "EmDisc"

Idents(D1,cells=gsub("-1_8","-1",names(which(mammal.combined$Tissue2%in%c("Tb","Tr") & mammal.combined$ID3=="C1"))) ) <- "Tb"
Idents(D2,cells=gsub("-1_9","-1",names(which(mammal.combined$Tissue2%in%c("Tb","Tr") & mammal.combined$ID3=="C2"))) ) <- "Tb"
#Idents(D3,cells=gsub("-1_10","-1",names(which(mammal.combined$Tissue2%in%c("Tb","Tr") & mammal.combined$ID3=="C3"))) ) <- "Tb"
Idents(D3B,cells=gsub("-1_11","-1",names(which(mammal.combined$Tissue2%in%c("Tb","Tr") & mammal.combined$ID3=="C4"))) ) <- "Tb"
Idents(D6,cells=gsub("-1_12","-1",names(which(mammal.combined$Tissue2%in%c("Tb","Tr") & mammal.combined$ID3=="C6"))) ) <- "Tb"
Idents(D7,cells=gsub("-1_13","-1",names(which(mammal.combined$Tissue2%in%c("Tb","Tr") & mammal.combined$ID3=="C7"))) ) <- "Tb"
Idents(D8,cells=gsub("-1_14","-1",names(which(mammal.combined$Tissue2%in%c("Tb","Tr") & mammal.combined$ID3=="C8"))) ) <- "Tb"

List1 <-  which(mammal.combined$Tissue2%in%c("EmDisc_d12","Embryonic","VE","Hyp_d12","Tb","Tr") & mammal.combined$ID3=="C1")
List2 <-  which(mammal.combined$Tissue2%in%c("EmDisc_d12","Embryonic","VE","Hyp_d12","Tb","Tr") & mammal.combined$ID3=="C2")
#List3 <-  which(mammal.combined$Tissue2%in%c("EmDisc_d12","Embryonic","VE","Hyp_d12","Tb","Tr") & mammal.combined$ID3=="C3")
List4 <-  which(mammal.combined$Tissue2%in%c("EmDisc_d12","Embryonic","VE","Hyp_d12","Tb","Tr") & mammal.combined$ID3=="C4")
List6 <-  which(mammal.combined$Tissue2%in%c("EmDisc_d12","Embryonic","VE","Hyp_d12","Tb","Tr") & mammal.combined$ID3=="C6")
List7 <-  which(mammal.combined$Tissue2%in%c("EmDisc_d12","Embryonic","VE","Hyp_d12","Tb","Tr") & mammal.combined$ID3=="C7")
List8 <-  which(mammal.combined$Tissue2%in%c("EmDisc_d12","Embryonic","VE","Hyp_d12","Tb","Tr") & mammal.combined$ID3=="C8")

D1 <- subset(D1,cells=gsub("-1_8","-1",names(List1)))
D2 <- subset(D2,cells=gsub("-1_9","-1",names(List2)))
#D3 <- subset(D3,cells=gsub("-1_10","-1",names(List3)))
D3B <- subset(D3B,cells=gsub("-1_11","-1",names(List4)))
D6 <- subset(D6,cells=gsub("-1_12","-1",names(List6)))
D7 <- subset(D7,cells=gsub("-1_13","-1",names(List7)))
D8 <- subset(D8,cells=gsub("-1_14","-1",names(List8)))

MergedData1$Dataset <- "10X Reference 1"
MergedData2$Dataset <- "10X Reference 1"
MergedData3$Dataset <- "10X Reference 1"
MergedData4$Dataset <- "10X Reference 1"
MergedData5$Dataset <- "10X Reference 1"
MergedData6$Dataset <- "10X Reference 1"
SSD$Dataset <- "SS2 Reference 1"


saveRDS(SSD,file="prD1.rds")
saveRDS(MergedData1,file="prD2.rds")
#saveRDS(MergedData2,file="prD3.rds")
saveRDS(MergedData3,file="prD4.rds")
saveRDS(MergedData4,file="prD5.rds")
saveRDS(MergedData5,file="prD6.rds")
saveRDS(MergedData6,file="prD7.rds")

saveRDS(D1,file="prD8.rds")
saveRDS(D2,file="prD9.rds")
saveRDS(D3,file="prD10.rds")
saveRDS(D3B,file="prD11.rds")
saveRDS(D6,file="prD12.rds")
saveRDS(D7,file="prD13.rds")
saveRDS(D8,file="prD14.rds")


#Do integration each dataset seperately
mammal.anchors <- FindIntegrationAnchors(object.list = list(SSD,MergedData1,MergedData2,MergedData3,MergedData4,MergedData5,MergedData6,D1,D2,D3B,D6,D8), dims = 1:20, anchor.features = 4000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

uID <- as.character(Idents(mammal.combined))
uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & Idents(mammal.combined)%in%c("CTB_d9"))] <- "STB_d9"
uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & Idents(mammal.combined)%in%c("CTB_d11"))] <- "STB_d11"
uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & Idents(mammal.combined)%in%c("CTB_d12"))] <- "STB_d12"

uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & Idents(mammal.combined)%in%c("STB_d9"))] <- "CTB_d9"
uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & Idents(mammal.combined)%in%c("STB_d11"))] <- "CTB_d11"
uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & Idents(mammal.combined)%in%c("STB_d12"))] <- "CTB_d12"

Idents(mammal.combined) <- uID

mammal.combined$Tissue2 <- Idents(mammal.combined) 

#Save this
#saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined_justEmbryonic.rds",sep=""))

#Do initial plots to take a look
DimPlot(mammal.combined,  reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat_splitII_justEmbryonic",".pdf",sep=""),width = 140, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined,  reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat_splitII_justEmbryonic",".pdf",sep=""),width = 140, height = 10, limitsize = FALSE, useDingbats=FALSE)

#Do initial plots to take a look
DimPlot(mammal.combined,  reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat_splitII_justEmbryonic_bydataset",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE, useDingbats=FALSE)

DimPlot(mammal.combined,  reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat_splitII_justEmbryonic_byddataset",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE, useDingbats=FALSE)

#Do some clustering at different resolutions 
mammal.combined <- FindClusters(mammal.combined, resolution = 1.0)
p<-DimPlot(mammal.combined, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-Cl10_justEmbryonic",".pdf",sep=""),width = 60, height = 10, limitsize = FALSE, useDingbats=FALSE,p)
p<-DimPlot(mammal.combined, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-Cl10_justEmbryonic",".pdf",sep=""),width = 60, height = 10, limitsize = FALSE, useDingbats=FALSE,p)
p<-DimPlot(mammal.combined, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-Cl10_justEmbryonic_bydataset",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE, useDingbats=FALSE,p)
p<-DimPlot(mammal.combined, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-Cl10_justEmbryonic_bydataset",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE, useDingbats=FALSE,p)
mammal.combined$Cl10 <- Idents(mammal.combined)
mammal.combined <- FindClusters(mammal.combined, resolution = 1.5)
p<-DimPlot(mammal.combined, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-Cl15_justEmbryonic",".pdf",sep=""),width = 60, height = 10, limitsize = FALSE, useDingbats=FALSE,p)
p<-DimPlot(mammal.combined, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-Cl15_justEmbryonic",".pdf",sep=""),width = 60, height = 10, limitsize = FALSE, useDingbats=FALSE,p)
mammal.combined$Cl15 <- Idents(mammal.combined)
p<-DimPlot(mammal.combined, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-Cl15_justEmbryonic_bydataset",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE, useDingbats=FALSE,p)
p<-DimPlot(mammal.combined, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-Cl15_justEmbryonic_bydataset",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE, useDingbats=FALSE,p)
mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
p<-DimPlot(mammal.combined, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-Cl05_justEmbryonic",".pdf",sep=""),width = 60, height = 10, limitsize = FALSE, useDingbats=FALSE,p)
p<-DimPlot(mammal.combined, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-Cl05_justEmbryonic",".pdf",sep=""),width = 60, height = 10, limitsize = FALSE, useDingbats=FALSE,p)
mammal.combined$Cl05 <- Idents(mammal.combined)
p<-DimPlot(mammal.combined, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-Cl05_justEmbryonic_bydataset",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE, useDingbats=FALSE,p)
p<-DimPlot(mammal.combined, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-Cl05_justEmbryonic_bydataset",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE, useDingbats=FALSE,p)
saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined_justEmbryonic.rds",sep=""))


sadsadadas

newID <- as.character(mammal.combined$Cells2)
newID[which(mammal.combined$Cl05%in%c(0,12,14) & mammal.combined$species2=="D")] <- "Stromal fibroblasts"
newID[which(mammal.combined$Cl05%in%c(7) & mammal.combined$species2=="D")] <- "Ciliated"
newID[which(mammal.combined$Cl05%in%c(1,3) & mammal.combined$species2=="D")] <- "Unciliated epithelia 2"
newID[which(mammal.combined$Cl05%in%c(2,5) & mammal.combined$species2=="D")] <- "Unciliated epithelia 1"
newID[which(mammal.combined$Cl05%in%c(10) )] <- "Unciliated epithelia 3"
newID[which(mammal.combined$Cl05%in%c(17) & mammal.combined$species2=="D")] <- "Hyp_d14"
newID[which(mammal.combined$Cl05%in%c(4,9,13,15) & mammal.combined$species2=="D")] <- "CTB_d14"
newID[which(mammal.combined$Cl05%in%c(6,8,11) & mammal.combined$species2=="D")] <- "STB_d14"
#newID[which(mammal.combined$Cl05%in%c(9) & mammal.combined$species2=="D")] <- "EVT_d14"
newID[which(mammal.combined$Cl05%in%c(18) & mammal.combined$species2=="D")] <- "EmDisc_d14"

newID[which(mammal.combined$Cl05%in%c(16) & mammal.combined$species2=="D")] <- "Unciliated epithelia 4"
newID[which(mammal.combined$Cl15%in%c(17) & mammal.combined$species2=="D")] <- "EVT"


Idents(mammal.combined) <- newID

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- cols[colind]

o


xzcxzcxzcx

#> unique(Dep1$ID5)#
#[1] 17
#> unique(Dep2$ID5)
#[1] 22
#> unique(Dep3$ID5)
#[1] 20
#> unique(Dep4$ID5)
#[1] 23
#> unique(Dep5$ID5)
#[1] 26
#> unique(Dep6$ID5)
#[1] 16



Dep1$species1 <- "Ref1"
Dep2$species1 <- "Ref2"
Dep3$species1 <- "Ref3"
Dep4$species1 <- "Ref4"
Dep5$species1 <- "Ref5"
Dep6$species1 <- "Ref6"

Dp1 <- readRDS("Endo1.rds") #22
Dp2 <- readRDS("Endo2.rds") #17
Dp3 <- readRDS("Endo3.rds") #22
Dp4 <- readRDS("Endo4.rds") #26
Dp5 <- readRDS("Endo5.rds") #20
Dp6 <- readRDS("Endo6.rds") #26
Dp7 <- readRDS("Endo7.rds") #19
Dp8 <- readRDS("Endo8.rds") #23
Dp9 <- readRDS("Endo9.rds") #16
#20

Dp1$species1 <- "Proc1"
Dp2$species1 <- "Proc2"
Dp3$species1 <- "Proc3"
Dp4$species1 <- "Proc4"
Dp5$species1 <- "Proc5"
Dp6$species1 <- "Proc6"
Dp7$species1 <- "Proc7"
Dp8$species1 <- "Proc8"
Dp9$species1 <- "Proc9"


saveRDS(Dep1,'prD1.rds')
saveRDS(Dep2,'prD2.rds')
saveRDS(Dep3,'prD3.rds')
saveRDS(Dep4,'prD4.rds')
saveRDS(Dep5,'prD5.rds')
saveRDS(Dep6,'prD6.rds')

saveRDS(De1,'prD6.rds')
saveRDS(Dp2,'prD6.rds')
saveRDS(Dp3,'prD6.rds')
saveRDS(Dp4,'prD6.rds')
saveRDS(Dp5,'prD6.rds')


mammal.anchors <- FindIntegrationAnchors(object.list = list(Dep1,Dep2,Dep3,Dep3,Dep4,Dep5,Dep6,Dp1,Dp3,Dp4,Dp5,Dp6,Dp7,Dp8,Dp9), dims = 1:20, anchor.features = 4000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined.rds",sep=""))


DimPlot(mammal.combined,  reduction = "umap", split.by = "species1", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat_splitII",".pdf",sep=""),width = 140, height = 10, limitsize = FALSE, useDingbats=FALSE)

DimPlot(mammal.combined,  reduction = "pca", split.by = "species1", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat_splitII",".pdf",sep=""),width = 140, height = 10, limitsize = FALSE, useDingbats=FALSE)

#saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined.rds",sep=""))


sadsadasdas

mammal.combined <- FindClusters(mammal.combined, resolution = 1.0)
DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-Cl10",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-Cl10",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
mammal.combined$Cl10 <- Idents(mammal.combined)

mammal.combined <- FindClusters(mammal.combined, resolution = 1.5)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-Cl15",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
mammal.combined$Cl15 <- Idents(mammal.combined)

mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-Cl05",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-Cl05",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
mammal.combined$Cl05 <- Idents(mammal.combined)


saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined.rds",sep=""))

