library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library("Matrix")
library(e1071)
library(pracma)
set.seed(1)
library(pheatmap)
#Same folder as before
saveext = "./FinalAlignB/"
dir.create(saveext)
dir.create(paste(saveext,"DimRed/",sep=""))
dir.create(paste(saveext,"Markers/",sep=""))

#Cell anotation and colours
cols <- c("#0767DA",
"#BF9600",
"#E68600",
"#E6B500",
"#F04C04",
"#D74404",
"#D74404",
"#BF3C04",
"#00BFBF",
"#0c9cf5",
"#0767DA",
"#0233BF",
"#0767DA",
"#877BD6",
"#7368CF",
"#5F54C7",
"#5F54C7",
"#5F54C7",
"#BF9600",
"#E68600",
"#AA8600",
"#921FE6",
"#8717D1",
"#7C10BB",
"#7108A6",
"#7C10BB",
"#877BD6",
"#6355B5",
"#3E2E94",
"#1A0873",
"#3E2E94",
"#BF9600",
"#BF9600",
"#BF0489",
"#00E6E6",
"#BF9600",
"#BF9600",
"#BF9600",
"#BF9600",
"#BF9600",
"#BF9600",
"#E68600",
"#E68600",
"#E68600",
"#E68600",
"#E68600",
"#E68600",
"#0767DA",
"#0767DA",
"#0767DA",
"#0767DA",
"#0767DA",
"#0767DA",
"lightgrey",
"#00BFBF",
"red",
"#5F54C7",
"black")


cType <-c("Embryonic",
"Stromal",
"Epithelial",
"Hyp_d9",
"Hyp_d11",
"Hyp_d12",
"VE",
"VE_d14",
"EmDisc_d9",
"EmDisc_d11",
"EmDisc_d12",
"EmDisc_d14",
"EmDisc",
"Am_d9",
"Am_d11",
"Am_d12",
"Am_d14",
"Am",
"Stromal",
"Uncilliated epithelial",
"Cilliated epithelial",
"CTB_d9",
"CTB_d11",
"CTB_d12",
"CTB_d14",
"CTB",
"STB_d9",
"STB_d11",
"STB_d12",
"STB_d14",
"STB",
"EVT_d14",
"EVT",
"Tr",
"ICM",
"Stromal_G1",
"Stromal_G2",
"Stromal_G3",
"Stromal_G4",
"Stromal_G5",
"Stromal_G6",
"Epithelial_G1",
"Epithelial_G2",
"Epithelial_G3",
"Epithelial_G4",
"Epithelial_G5",
"Epithelial_G6",
"Embryonic_G1",
"Embryonic_G2",
"Embryonic_G3",
"Embryonic_G4",
"Embryonic_G5",
"Embryonic_G6",
"NGT",
"Epi",
"Cont",
"Am/EmDisc_d14",
"Unkown_Emb")

cType <- c("9","15","14","11","21","20","23",
"7","17","6","10","8","16","13",
"2","1","5","3","12","4","0","19","18","22")

cols <- c("#00BFBF",
"#00ABBF",
"#0197BF",
"#0183BF",
"#016FBF",
"#015BBF",
"#0247BF",
"#E68600",
"#D9AD0D",
"#C49C0A",
"#AE8A06",
"#997903", 
"#826701",
"#6e5701", 
"#634e00",
"#f7ad45",
"#FAA734",
"#F8A32E",
"#F6A028",
"#F39C23",
"#F1981D",
"#EF9517",
"#ED9111",
"#EA8D0C",
"#E88A06")



#Load in the previous run
mammal.combined <- readRDS(file = paste(saveext,"Seurat_combined_FinalAno.rds",sep=""))


Idents(mammal.combined) <- mammal.combined$Dataset
mammal.combined <- subset(mammal.combined,idents=c("10X Ours"))
mito.genes <- grep(pattern = "^MT-", x = rownames(mammal.combined@assays[["RNA"]]), value = TRUE)
percent.mito2 <- Matrix::colSums(mammal.combined@assays[["RNA"]][mito.genes, ])/Matrix::colSums(mammal.combined@assays[["RNA"]])
mammal.combined <- AddMetaData(object = mammal.combined, metadata = percent.mito2, col.name = "percent.mito2")

FeaturePlot(mammal.combined, reduction = "pca", features = "nFeature_RNA", pt.size = 2, cols = c("lightgrey", "black")) # + xlim(c(-13, 13)) + ylim(c(-13, 13))#
ggsave(filename=paste(saveext,"/Markers/","AraMarkers", "_nFeature_RNA_nosplit.pdf",sep=""),width = 10, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, reduction = "pca", features = "nCount_RNA", pt.size = 2, cols = c("lightgrey", "black")) # + xlim(c(-13, 13)) + ylim(c(-13, 13))#
ggsave(filename=paste(saveext,"/Markers/","AraMarkers", "_nCount_RNA_nosplit.pdf",sep=""),width = 10, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, reduction = "pca", features = "percent.mito2", pt.size = 2, cols = c("lightgrey", "black")) # + xlim(c(-13, 13)) + ylim(c(-13, 13))#
ggsave(filename=paste(saveext,"/Markers/","AraMarkers", "_mito_nosplit.pdf",sep=""),width = 10, height = 10,limitsize = FALSE)

VlnPlot(mammal.combined,features=c("nCount_RNA"),split.by="ID3",group.by="ID3")
ggsave(filename=paste(saveext,"/Markers/","Violin_", "_nCount_RNA.pdf",sep=""),width = 10, height = 15,limitsize = FALSE)
VlnPlot(mammal.combined,features=c("nFeature_RNA"),split.by="ID3",group.by="ID3")
ggsave(filename=paste(saveext,"/Markers/","Violin_", "nFeature_RNA.pdf",sep=""),width = 10, height = 15,limitsize = FALSE)

VlnPlot(mammal.combined,features=c("percent.mito2"),split.by="ID3",group.by="ID3")
ggsave(filename=paste(saveext,"/Markers/","Violin_", "mito_RNA.pdf",sep=""),width = 10, height = 15,limitsize = FALSE)



#colind <- integer( length( levels(Idents(mammal.combined)) )  )
#for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
#  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
#}
#coluse <- cols[colind]
#
#p1 <- DimPlot(mammal.combined, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
#ggsave(filename=paste(saveext,"/DimRed/ClusterBasedAnotation",".pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)



require(tidyverse)


mammal.combined3 <- readRDS(file = paste(saveext,"Seurat_combined_FinalAno.rds",sep=""))
X <- mammal.combined3$Dataset
X2 <- mammal.combined3$ID3
Y <- factor(mammal.combined3$Cl05,levels=c("9","15","14","11","21","20","23","7","17","6","10","8","16","13","2","1","5","3","12","4","0","19","18","22"))

Z <- as.character(Idents(mammal.combined3))
Z2 <- mammal.combined3$Genotype
Z2[which(Z2=="NAssigned")] <- "NGT"
Z2[which(Z2=="Epithelial_G1")] <- "Epithelail"
Z2[which(Z2=="Epithelial_G2")] <- "Epithelail"
Z2[which(Z2=="Epithelial_G3")] <- "Epithelail"
Z2[which(Z2=="Epithelial_G4")] <- "Epithelail"
Z2[which(Z2=="Epithelial_G5")] <- "Epithelail"
Z2[which(Z2=="Epithelial_G6")] <- "Epithelail"

Z2[which(Z2=="Stromal_G1")] <- "Stromal"
Z2[which(Z2=="Stromal_G2")] <- "Stromal"
Z2[which(Z2=="Stromal_G3")] <- "Stromal"
Z2[which(Z2=="Stromal_G4")] <- "Stromal"
Z2[which(Z2=="Stromal_G5")] <- "Stromal"
Z2[which(Z2=="Stromal_G6")] <- "Stromal"

Z2[which(Z2=="Embryonic_G1")] <- "Embryonic"
Z2[which(Z2=="Embryonic_G2")] <- "Embryonic"
Z2[which(Z2=="Embryonic_G3")] <- "Embryonic"
Z2[which(Z2=="Embryonic_G4")] <- "Embryonic"
Z2[which(Z2=="Embryonic_G5")] <- "Embryonic"
Z2[which(Z2=="Embryonic_G6")] <- "Embryonic"
Z2[which(Z2=="Embryoic_G3")] <- "Embryonic"


#Embryoic_G3

Dat <- data.frame(x=Y[which(X=="SS2 Reference 1")],y=Z[which(X=="SS2 Reference 1")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"Table1","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)


Dat <- data.frame(x=Y[which(X=="10X Reference 1")],y=Z[which(X=="10X Reference 1")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"Table2","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)


Dat <- data.frame(x=Y[which(X2=="C1")],y=Z2[which(X2=="C1")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","lightgrey","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"TableC1","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)

Dat <- data.frame(x=Y[which(X2=="C2")],y=Z2[which(X2=="C2")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","lightgrey","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"TableC2","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
Dat <- data.frame(x=Y[which(X2=="C3")],y=Z2[which(X2=="C3")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","lightgrey","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"TableC3","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
Dat <- data.frame(x=Y[which(X2=="C4")],y=Z2[which(X2=="C4")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","lightgrey","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"TableC4","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)

Dat <- data.frame(x=Y[which(X2=="C6")],y=Z2[which(X2=="C6")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","lightgrey","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"TableC6","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)
Dat <- data.frame(x=Y[which(X2=="C7")],y=Z2[which(X2=="C7")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","lightgrey","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"TableC7","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)

Dat <- data.frame(x=Y[which(X2=="C8")],y=Z2[which(X2=="C8")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#E68600","lightgrey","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"TableC8","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)


Dat <- data.frame(x=Y[which(X2%in%c("C1","C2","C3","C4","C6","C7","C8"))],y=Z2[which(X2%in%c("C1","C2","C3","C4","C6","C7","C8"))] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","lightgrey","#BF9600"))+ theme_classic()
ggsave(filename=paste(saveext,"TableC9","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)


sadsdadsadadsadsdas
#Idents(SSD) <- uID

#uID[which(uID=="EPI_CS4")] <- "TB_CS4"

uID[which(mammal.combined1$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined1$Cl05%in%c(13))] <- "EmDisc_d14"
uID[which(mammal.combined1$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined1$Cl05%in%c(16))] <- "Hyp_d14"
uID[which(mammal.combined1$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined1$Cl05%in%c(3,5,12,17))] <- "CTB_d14"
uID[which(mammal.combined1$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined1$Cl05%in%c(0,1,2,4,6,7,8,9,10,11,14,15,18,19,20))] <- "STB_d14"
uID[which(mammal.combined1$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined1$Cl15%in%c(5))] <- "EVT_d14"


Idents(mammal.combined1) <- uID #mammal.combined1$Tissue

saveRDS(mammal.combined1,file = paste(saveext,"Seurat_combined_justEmbryonic_ano.rds",sep=""))


colind <- integer( length( levels(Idents(mammal.combined1)) )  )
for (i in 1:length( levels(Idents(mammal.combined1)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined1))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined1, cols=coluse,  pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/PCA_1_2_Embryoniconly_ano","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)
p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,3))
ggsave(filename=paste(saveext,"/DimRed/PCA_1_3_Embryoniconly_ano","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)
p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(2,3))
ggsave(filename=paste(saveext,"/DimRed/PCA_2_3_Embryoniconly_ano","boundary.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE)

p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "umap", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/UMAP_1_2Embryoniconly_ano","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)

p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "umap", label = TRUE, split.by="ID3", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/UMAP_1_2Embryoniconly_ano","boundary_byID3.pdf",sep=""),width = 120, height = 10, useDingbats = FALSE, limitsize = FALSE)


mammal.combined1$ano <- Idents(mammal.combined1)
Idents(mammal.combined1) <- mammal.combined1$Dataset

mammal.combined3 <- subset(mammal.combined1,idents=c("10X Ours"))
mammal.combined2 <- subset(mammal.combined1,idents=c("SS2 Reference 1"))
mammal.combined1 <- subset(mammal.combined1,idents=c("10X Reference 1"))


#sadsadsad

uID <- as.character(mammal.combined1$ano)
uID[which(uID%in%c("EmDisc_d9","EmDisc_d11","EmDisc_d12"))] <- "1) EmDisc"
uID[which(uID%in%c("CTB_d9","CTB_d11","CTB_d12","STB_d9","STB_d11","STB_d12","Hyp_d9","Hyp_d11","Hyp_d12"))] <- "2) Not EmDisc"
Idents(mammal.combined1) <- uID

#Grid of points for later
make.grid = function(x, n = 1000) {
  grange = apply(x, 2, range)
  x1 = seq(from = grange[1,1], to = grange[2,1], length = n)
  x2 = seq(from = grange[1,2], to = grange[2,2], length = n)
  expand.grid(X1 = x1, X2 = x2)
}

dat <- as.data.frame(mammal.combined1 [["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined1)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "radial", cost = .5, scale = TRUE)
xgrid = make.grid(as.matrix(dat[,1:2]))
ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
grange = apply(dat[,1:2], 2, range)

m1 <- seq(grange[1,1],grange[2,1],length.out=1000)
m2 <- seq(grange[1,2],grange[2,2],length.out=1000)
m3 <- matrix(func, 1000, 1000)

XX <- t(repmat(seq(grange[1,1],grange[2,1],length.out=1000),1000,1))
YY <- (repmat(seq(grange[1,2],grange[2,2],length.out=1000),1000,1))
res <- matrix(m3, dim(m3)[1]*dim(m3)[2],1)

DATA <- data.frame(x=matrix(XX, dim(m3)[1]*dim(m3)[2],1),y=matrix(YY, dim(m3)[1]*dim(m3)[2],1),z=matrix(m3, dim(m3)[1]*dim(m3)[2],1))

Idents(mammal.combined1) <- mammal.combined1$ano
colind <- integer( length( levels(Idents(mammal.combined1)) ) )  
for (i in 1:length( levels(Idents(mammal.combined1)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined1))[i])
}
coluse <- cols[colind]
p1 <- DimPlot(mammal.combined1, cols=coluse,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_B1-B6_Embryonic","boundary_v_maternal_embryoniconly.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)


Idents(mammal.combined3) <- mammal.combined3$ano
colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]
p1 <- DimPlot(mammal.combined3, cols=coluse,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z))
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_C1-C8_Embryonic_v_maternal","boundary_embryoniconly.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

uID <- as.character(mammal.combined3$ano) 
uID[which(mammal.combined3$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined3$Cl05%in%c(10))] <- "Novel1"
uID[which(mammal.combined3$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined3$Cl05%in%c(19))] <- "Novel2"
uID[which(mammal.combined3$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined3$Cl05%in%c(20))] <- "Novel3"
Idents(mammal.combined3) <- uID

DefaultAssay(mammal.combined3) <- "RNA"
#Cl <- FindAllMarkers(mammal.combined3, test.use = "MAST")

saveRDS(uID,file="IDs.rds")
saveRDS(GetAssayData(mammal.combined3,assay="RNA"),file="RNAData.rds")



Idents(mammal.combined3) <- mammal.combined3$Genotype
p1 <- DimPlot(mammal.combined3,  pt.size = 4, reduction = "umap", label = TRUE, split.by="ID3", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/UMAP_1_2Embryoniconly_ano","genotype_byID3.pdf",sep=""),width = 120, height = 10, useDingbats = FALSE, limitsize = FALSE)

mammal.combined3B <- subset(mammal.combined3,idents=c("Epithelial_G2","Epithelial_G4","Stromal_G1","Stromal_G3","Epithelial_G1","Epithelial_G6","Stromal_G4","Epithelial_G3","Stromal_G2"),invert=TRUE)

Idents(mammal.combined3B) <- mammal.combined3B$ano

colind <- integer( length( levels(Idents(mammal.combined3B)) )  )
for (i in 1:length( levels(Idents(mammal.combined3B)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3B))[i])
}
coluse <- cols[colind]
p1 <- DimPlot(mammal.combined3B, cols=coluse,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_C1-C8_Embryonic_testfilter","boundary_embryoniconly.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

saveRDS(mammal.combined3B,file=paste(saveext,"Seurat_combined_justEmbryonic_filtered.rds",sep=""))



saveRDS(mammal.combined3B$ano,file="IDs.rds")
saveRDS(GetAssayData(mammal.combined3B,assay="RNA"),file="RNAData.rds")

#write.table(as.data.frame(Cl),file="ClustersEmbryonicOnly.csv",quote = FALSE, sep = ",")

DefaultAssay(mammal.combined3) <- "integrated"
genelist <- rownames(GetAssayData(mammal.combined3))

Idents(mammal.combined3) <- mammal.combined3$Cl05
DefaultAssay(mammal.combined3) <- "RNA"
Ae <- AverageExpression(mammal.combined3)
Ae <- Ae$RNA
#Ae$gene <- rownames(Ae)

library(pheatmap)

redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
mat_breaks <- seq(0, 6, length.out = 20)

C1 <- cor(  log2(Ae +1),log2(Ae +1) )
C2 <- cor(  log2(Ae[genelist,] +1),log2(Ae[genelist,] +1) )

pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 3), border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/LineageHeatmap_C1_RNAAll.pdf",sep=""))
pheatmap(C2,color =  redblue1(20), display_numbers = round(C1, digits = 3), border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/LineageHeatmap_C1_RNAIntgenes.pdf",sep=""))



#sadsadsadsad


mammal.combined <- readRDS(file = paste(saveext,"Seurat_combined_justEmbryonic.rds",sep=""))

uID <- as.character(mammal.combined$Cells)
#uID[which(uID=="EPI_CS4")] <- "TB_CS4"

uID[which(mammal.combined$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined$Cl05%in%c(13))] <- "EmDisc_d14"
uID[which(mammal.combined$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined$Cl05%in%c(16))] <- "Hyp_d14"
uID[which(mammal.combined$ID3%in%c("C1","C2","C3","C4","C6","C7","C8") & mammal.combined$Cl05%in%c(0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,17,18,19,20))] <- "Tb_d14"

#uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & mammal.combined$Cells%in%c("CTB_d9"))] <- "STB_d9"
#uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & mammal.combined$Cells%in%c("CTB_d11"))] <- "STB_d11"
#uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & mammal.combined$Cells%in%c("CTB_d12"))] <- "STB_d12"

#uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & mammal.combined$Cells%in%c("STB_d9"))] <- "CTB_d9"
#uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & mammal.combined$Cells%in%c("STB_d11"))] <- "CTB_d11"
#uID[which(mammal.combined$ID3%in%c("B1","B2","B3","B4","B5","B6") & mammal.combined$Cells%in%c("STB_d12"))] <- "CTB_d12"




mammal.combined$anno <- uID
Idents(mammal.combined) <- mammal.combined$Dataset

mammal.combined1 <- subset(mammal.combined,idents=c("10X Reference 1"))
mammal.combined2 <- subset(mammal.combined,idents=c("SS2 Reference 1"))
mammal.combined3 <- subset(mammal.combined,idents=c("10X Ours"))


Idents(mammal.combined1) <- mammal.combined1$anno
colind <- integer( length( levels(Idents(mammal.combined1)) )  )
for (i in 1:length( levels(Idents(mammal.combined1)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined1))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/PCA1_Embryoniconly","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)


#Idents(mammal.combined2) <- mammal.combined2$anno 
#colind <- integer( length( levels(Idents(mammal.combined2)) )  )
#for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
#  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
#}
#coluse <- cols[colind]

#p1 <- DimPlot(mammal.combined2, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, dims=c(1,3))
#ggsave(filename=paste(saveext,"/DimRed/PCA2_Embryoniconly","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)


Idents(mammal.combined3) <- mammal.combined3$anno 
colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined3, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/PCA3_Embryoniconly","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

#p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE, dims=c(1,2))
#ggsave(filename=paste(saveext,"/DimRed/UMAP_B1-B6_EmDHyp","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)



sdsadsaadadad


uID <- as.character(mammal.combined1$Tissue2)
uID[which(uID%in%c("Tb"))] <- "1) Tb"
uID[which(uID%in%c("EmDisc_d12","Embryonic","Hyp_d12","Stromal","Epithelial"))] <- "2) Not Tb"
Idents(mammal.combined1) <- uID

dat <- as.data.frame(mammal.combined1 [["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined1)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "radial", cost = 1, scale = TRUE)

xgrid = make.grid(as.matrix(dat[,1:2]))

ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
grange = apply(dat[,1:2], 2, range)

m1 <- seq(grange[1,1],grange[2,1],length.out=1000)
m2 <- seq(grange[1,2],grange[2,2],length.out=1000)
m3 <- matrix(func, 1000, 1000)

XX <- t(repmat(seq(grange[1,1],grange[2,1],length.out=1000),1000,1))
YY <- (repmat(seq(grange[1,2],grange[2,2],length.out=1000),1000,1))
res <- matrix(m3, dim(m3)[1]*dim(m3)[2],1)

DATA <- data.frame(x=matrix(XX, dim(m3)[1]*dim(m3)[2],1),y=matrix(YY, dim(m3)[1]*dim(m3)[2],1),z=matrix(m3, dim(m3)[1]*dim(m3)[2],1))

Idents(mammal.combined1) <- mammal.combined1$Tissue2
p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z))
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_B1-B6_Tb","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

uID <- as.character(mammal.combined1$Tissue2)
uID[which(uID%in%c("Epithelial"))] <- "1) Epithelial"
uID[which(uID%in%c("Tb","EmDisc_d12","Embryonic","Hyp_d12","Stromal"))] <- "2) Noot Epithelial"
Idents(mammal.combined1) <- uID

dat <- as.data.frame(mammal.combined1 [["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined1)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "radial", cost = 1, scale = TRUE)

ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
grange = apply(dat[,1:2], 2, range)

m1 <- seq(grange[1,1],grange[2,1],length.out=1000)
m2 <- seq(grange[1,2],grange[2,2],length.out=1000)
m3 <- matrix(func, 1000, 1000)

XX <- t(repmat(seq(grange[1,1],grange[2,1],length.out=1000),1000,1))
YY <- (repmat(seq(grange[1,2],grange[2,2],length.out=1000),1000,1))
res <- matrix(m3, dim(m3)[1]*dim(m3)[2],1)
DATA <- data.frame(x=matrix(XX, dim(m3)[1]*dim(m3)[2],1),y=matrix(YY, dim(m3)[1]*dim(m3)[2],1),z=matrix(m3, dim(m3)[1]*dim(m3)[2],1))
Idents(mammal.combined1) <- mammal.combined1$Tissue2
p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z))
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_B1-B6_Epithelial","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

uID <- as.character(mammal.combined1$Tissue2)
uID[which(uID%in%c("Stromal"))] <- "1) Stromal"
uID[which(uID%in%c("Tb","EmDisc_d12","Embryonic","Hyp_d12","Epithelial"))] <- "2) Not Stromal"
Idents(mammal.combined1) <- uID

dat <- as.data.frame(mammal.combined1 [["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined1)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "radial", cost = 1, scale = TRUE)

ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
grange = apply(dat[,1:2], 2, range)

m1 <- seq(grange[1,1],grange[2,1],length.out=1000)
m2 <- seq(grange[1,2],grange[2,2],length.out=1000)
m3 <- matrix(func, 1000, 1000)

XX <- t(repmat(seq(grange[1,1],grange[2,1],length.out=1000),1000,1))
YY <- (repmat(seq(grange[1,2],grange[2,2],length.out=1000),1000,1))
res <- matrix(m3, dim(m3)[1]*dim(m3)[2],1)
DATA <- data.frame(x=matrix(XX, dim(m3)[1]*dim(m3)[2],1),y=matrix(YY, dim(m3)[1]*dim(m3)[2],1),z=matrix(m3, dim(m3)[1]*dim(m3)[2],1))
Idents(mammal.combined1) <- mammal.combined1$Tissue2
p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z))
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_B1-B6_Stromal","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)



uID <- as.character(mammal.combined1$Tissue2)
uID[which(uID%in%c("Stromal","Epithelial"))] <- "1) Maternal"
uID[which(uID%in%c("Tb","EmDisc_d12","Embryonic","Hyp_d12"))] <- "2) Not Maternal"
Idents(mammal.combined1) <- uID

dat <- as.data.frame(mammal.combined1 [["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined1)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "radial", cost = 1, scale = TRUE)

ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
grange = apply(dat[,1:2], 2, range)

m1 <- seq(grange[1,1],grange[2,1],length.out=1000)
m2 <- seq(grange[1,2],grange[2,2],length.out=1000)
m3 <- matrix(func, 1000, 1000)

XX <- t(repmat(seq(grange[1,1],grange[2,1],length.out=1000),1000,1))
YY <- (repmat(seq(grange[1,2],grange[2,2],length.out=1000),1000,1))
res <- matrix(m3, dim(m3)[1]*dim(m3)[2],1)
DATA <- data.frame(x=matrix(XX, dim(m3)[1]*dim(m3)[2],1),y=matrix(YY, dim(m3)[1]*dim(m3)[2],1),z=matrix(m3, dim(m3)[1]*dim(m3)[2],1))
Idents(mammal.combined1) <- mammal.combined1$Tissue2
p1 <- DimPlot(mammal.combined1, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z))
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_B1-B6_Mat","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)



Idents(mammal.combined3) <- mammal.combined3$Tissue2
colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined3, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/PCA_1_2_C1-C6_EmDHyp","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
p1 <- DimPlot(mammal.combined3, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, dims=c(1,3))
ggsave(filename=paste(saveext,"/DimRed/PCA_1_3_C1-C6_EmDHyp","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
p1 <- DimPlot(mammal.combined3, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, dims=c(2,3))
ggsave(filename=paste(saveext,"/DimRed/PCA_2_3_C1-C6_EmDHyp","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)
p1 <- DimPlot(mammal.combined3, cols=coluse, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/UMAP_C1-C6_EmDHyp","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)

uID <- as.character(mammal.combined3$Tissue2)
uID[which(uID%in%c("Tb","Stromal","Epithelial"))] <- "2) Not embryonic"
uID[which(uID%in%c("EmDisc_d12","Embryonic","Hyp_d12"))] <- "1) EmDisc/Hyp"
Idents(mammal.combined3) <- uID

#Grid of points for later
make.grid = function(x, n = 1000) {
  grange = apply(x, 2, range)
  x1 = seq(from = grange[1,1], to = grange[2,1], length = n)
  x2 = seq(from = grange[1,2], to = grange[2,2], length = n)
  expand.grid(X1 = x1, X2 = x2)
}

dat <- as.data.frame(mammal.combined3[["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined3)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "radial", cost = .5, scale = TRUE)

xgrid = make.grid(as.matrix(dat[,1:2]))
ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
grange = apply(dat[,1:2], 2, range)

m1 <- seq(grange[1,1],grange[2,1],length.out=1000)
m2 <- seq(grange[1,2],grange[2,2],length.out=1000)
m3 <- matrix(func, 1000, 1000)

XX <- t(repmat(seq(grange[1,1],grange[2,1],length.out=1000),1000,1))
YY <- (repmat(seq(grange[1,2],grange[2,2],length.out=1000),1000,1))
res <- matrix(m3, dim(m3)[1]*dim(m3)[2],1)

DATA <- data.frame(x=matrix(XX, dim(m3)[1]*dim(m3)[2],1),y=matrix(YY, dim(m3)[1]*dim(m3)[2],1),z=matrix(m3, dim(m3)[1]*dim(m3)[2],1))

Idents(mammal.combined3) <- mammal.combined3$Tissue2
colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]
p1 <- DimPlot(mammal.combined3, cols=coluse,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z))
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_C1-C8_EmDHyp","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)


uID <- as.character(mammal.combined3$Tissue2)
uID[which(uID%in%c("Tb"))] <- "1) Tb"
uID[which(uID%in%c("EmDisc_d12","Embryonic","Hyp_d12","Stromal","Epithelial"))] <- "2) Not Tb"
Idents(mammal.combined3) <- uID

dat <- as.data.frame(mammal.combined3[["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined3)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "radial", cost = 1, scale = TRUE)

xgrid = make.grid(as.matrix(dat[,1:2]))

ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
grange = apply(dat[,1:2], 2, range)

m1 <- seq(grange[1,1],grange[2,1],length.out=1000)
m2 <- seq(grange[1,2],grange[2,2],length.out=1000)
m3 <- matrix(func, 1000, 1000)

XX <- t(repmat(seq(grange[1,1],grange[2,1],length.out=1000),1000,1))
YY <- (repmat(seq(grange[1,2],grange[2,2],length.out=1000),1000,1))
res <- matrix(m3, dim(m3)[1]*dim(m3)[2],1)

DATA <- data.frame(x=matrix(XX, dim(m3)[1]*dim(m3)[2],1),y=matrix(YY, dim(m3)[1]*dim(m3)[2],1),z=matrix(m3, dim(m3)[1]*dim(m3)[2],1))

Idents(mammal.combined3) <- mammal.combined3$Tissue2
p1 <- DimPlot(mammal.combined3, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z))
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_C1-C6_Tb","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)




uID <- as.character(mammal.combined3$Tissue2)
uID[which(uID%in%c("Stromal","Epithelial"))] <- "1) Maternal"
uID[which(uID%in%c("Tb","EmDisc_d12","Embryonic","Hyp_d12"))] <- "2) Not Maternal"
Idents(mammal.combined3) <- uID

dat <- as.data.frame(mammal.combined3[["pca"]]@cell.embeddings)
dat$ID <- Idents(mammal.combined3)
svmfit = svm(x = as.matrix(dat[,1:2]),  y = dat$ID, kernel = "radial", cost = 1, scale = TRUE)

ygrid = predict(svmfit, xgrid)
func = predict(svmfit, xgrid, decision.values = TRUE)
func = attributes(func)$decision
grange = apply(dat[,1:2], 2, range)

m1 <- seq(grange[1,1],grange[2,1],length.out=1000)
m2 <- seq(grange[1,2],grange[2,2],length.out=1000)
m3 <- matrix(func, 1000, 1000)

XX <- t(repmat(seq(grange[1,1],grange[2,1],length.out=1000),1000,1))
YY <- (repmat(seq(grange[1,2],grange[2,2],length.out=1000),1000,1))
res <- matrix(m3, dim(m3)[1]*dim(m3)[2],1)
DATA <- data.frame(x=matrix(XX, dim(m3)[1]*dim(m3)[2],1),y=matrix(YY, dim(m3)[1]*dim(m3)[2],1),z=matrix(m3, dim(m3)[1]*dim(m3)[2],1))
Idents(mammal.combined3) <- mammal.combined3$Tissue2
p1 <- DimPlot(mammal.combined3, cols=coluse, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE)
p1 <- p1 + geom_contour(data = DATA, aes(x=x,y=y,z = z))
ggsave(filename=paste(saveext,"/DimRed/PCA_decision_C1-C6_Mat","boundary.pdf",sep=""),width = 12, height = 8, useDingbats = FALSE)


x


sdsadas

newID <- as.character(mammal.combined$Cells)


newID[which(mammal.combined$Cl05%in%c(0,10,12) & mammal.combined$species2=="New")] <- "Stromal fibroblasts"
newID[which(mammal.combined$Cl05%in%c(8) & mammal.combined$species2=="New")] <- "Ciliated"

newID[which(mammal.combined$Cl05%in%c(3,4,5,7) & mammal.combined$species2=="New")] <- "Unciliated epithelia 2"
newID[which(mammal.combined$Cl05%in%c(1) & mammal.combined$species2=="New")] <- "Unciliated epithelia 1"
newID[which(mammal.combined$Cl05%in%c(13) )] <- "Unciliated epithelia 3"

newID[which(mammal.combined$Cl05%in%c(17) & mammal.combined$species2=="New")] <- "Hyp_d14"

newID[which(mammal.combined$Cl05%in%c(2,15) & mammal.combined$species2=="New")] <- "CTB_d14"
newID[which(mammal.combined$Cl05%in%c(6,9,11,14) & mammal.combined$species2=="New")] <- "STB_d14"
newID[which(mammal.combined$Cl05%in%c(19) & mammal.combined$species2=="New")] <- "EmDisc_d14"

#newID[which(mammal.combined$Cl05%in%c(16) & mammal.combined$species2=="New")] <- "Unciliated epithelia 4"
newID[which(mammal.combined$Cl05%in%c(21) & mammal.combined$species2=="New")] <- "EVT"

Idents(mammal.combined) <- newID

#colind <- integer( length( levels(Idents(mammal.combined)) )  )
#for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
#  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
#}
#coluse <- cols[colind]

DimPlot(mammal.combined, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat_splitIICol",".pdf",sep=""),width = 60, height = 10, limitsize = FALSE, useDingbats=FALSE)
#DimPlot(mammal.combined, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"UMAP_SeuratIICol",".pdf",sep=""),width = 60, height = 10, limitsize = FALSE, useDingbats=FALSE)

sadsadsa

DimPlot(mammal.combined, cols = coluse,  reduction = "umap", split.by = "species1", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat_splitIICol",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, cols = coluse, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_SeuratIICol",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE)



DimPlot(mammal.combined,cols = coluse,  reduction = "pca", split.by = "species1", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat_splitIICol",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, cols = coluse, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_SeuratIICol",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE)

saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined_ano1.rds",sep=""))

