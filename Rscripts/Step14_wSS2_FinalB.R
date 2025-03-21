library(Seurat) 
library(ggplot2) 
library(cowplot) 
library(Matrix) 
library(dplyr) 
library("Matrix") 
library(pheatmap)
#library(destiny)
set.seed(1)

#Save folder
saveext = "./FinalAlignB/"
dir.create(saveext)
dir.create(paste(saveext,"DimRed/",sep=""))
dir.create(paste(saveext,"Markers/",sep=""))


other_endometrium_data <- readRDS(file="../../Endometrial/OtherEndo/other_endometrium_data.rds")
other_endometrium_data <- subset(other_endometrium_data,idents=c("eS","dS","SOX9","Ciliated","Lumenal","Fibroblast C7","Glandular"))

#DE1 <- FindMarkers(other_endometrium_data,ident.1 = "eS", ident.2 = "dS",only.pos=TRUE,test.use="MAST")
#DE1A <- FindMarkers(other_endometrium_data,ident.2 = "eS", ident.1 = "dS", only.pos=TRUE,test.use="MAST")

#DE2 <- FindMarkers(other_endometrium_data,ident.1 = "Ciliated", ident.2 = "Lumenal",only.pos=TRUE,test.use="MAST")
#DE2A <- FindMarkers(other_endometrium_data,ident.2 = "Ciliated", ident.1 = "Lumenal",only.pos=TRUE,test.use="MAST")

#DE3 <- FindMarkers(other_endometrium_data,ident.1 = "Lumenal", ident.2 = "Glandular",only.pos=TRUE,test.use="MAST")
#DE3A <- FindMarkers(other_endometrium_data,ident.2 = "Lumenal", ident.1 = "Glandular",only.pos=TRUE,test.use="MAST")

DE4 <- FindMarkers(other_endometrium_data,ident.1 = "Ciliated", ident.2 = "Glandular",only.pos=TRUE,test.use="MAST")
DE4A <- FindMarkers(other_endometrium_data,ident.2 = "Ciliated", ident.1 = "Glandular",only.pos=TRUE,test.use="MAST")

#saveRDS(DE1,file="DE_es_dS1.rds")
#saveRDS(DE1A,file="DE_es_dS2.rds")

#saveRDS(DE2,file="DE_Cil_Lum1.rds")
#saveRDS(DE2A,file="DE_Cil_Lum2.rds")

#saveRDS(DE3,file="DE_Lum_Gla1.rds")
#saveRDS(DE3A,file="DE_Lum_Gla2.rds")

saveRDS(DE4,file="DE_Cil_Gla1.rds")
saveRDS(DE4A,file="DE_Cil_Gla2.rds")

#DE1<-readRDS(DE1,file="DE_es_dS1.rds")
#DE1<-DE1<-readRDS(DE1A,file="DE_es_dS2.rds")

DE2<-readRDS(file="DE_Cil_Lum1.rds")
DE2A<-readRDS(file="DE_Cil_Lum2.rds")

DE3<-readRDS(file="DE_Lum_Gla1.rds")
DE3A<-readRDS(file="DE_Lum_Gla2.rds")

DE4<-readRDS(file="DE_Cil_Gla1.rds")
DE4A<-readRDS(file="DE_Cil_Gla2.rds")


dsfsfsfdfsfsd

mammal.combined2 <- readRDS(file = paste(saveext,"Seurat_combined_FinalAno.rds",sep=""))
mammal.combined2$Cells <- Idents(mammal.combined2)

#X1 <- GetAssayData(D1)
#saveRDS(X1,file="Tb_expession.rds")
#saveRDS(Idents(D1),file="Tb_expession_label.rds")
#saveRDS(D1$Dataset,file="Tb_expession_Dataset.rds")

#mammal.combined2 <- subset(mammal.combined2,idents=c("Unciliated epithelia","Ciliated"))
#mammal.combined2$Cells <- Idents(mammal.combined2)
#Idents(mammal.combined2) <- mammal.combined2$Dataset
#mammal.combined2 <- subset(mammal.combined2,idents="10X Ours")

Idents(mammal.combined2) <- mammal.combined2$Cells
DefaultAssay(mammal.combined2) <- "RNA"
Type0<-WhichCells(mammal.combined2,idents=c("Unciliated epithelia","Ciliated"))

Type1<-WhichCells(mammal.combined2,expression=SPP1>0.250961574,slot="data")
Type2<-WhichCells(mammal.combined2,expression=SCGB2A2>0.250961574,slot="data")
Gland <- intersect(intersect(Type1,Type2),Type0)
Idents(mammal.combined2,cells=Gland) <- "Glandular"

Idents(mammal.combined2) <- mammal.combined2$Cl05
Type3<-WhichCells(mammal.combined2,idents=c(1,5))
Type3B<-WhichCells(mammal.combined2,idents=c(0,19))
Type3C<-WhichCells(mammal.combined2,idents=c(22,13,12,3,8,4,2))

Idents(mammal.combined2) <- mammal.combined2$Dataset
Type4<-WhichCells(mammal.combined2,idents=c("10X Ours"))

Idents(mammal.combined2) <- mammal.combined2$Cells
Idents(mammal.combined2,cells=intersect(Type3,Type4)) <- "Glandular"

Idents(mammal.combined2) <- mammal.combined2$Cells
Idents(mammal.combined2,cells=intersect(Type3C,Type4)) <- "Other"


EVT <- readRDS("newEVT.rds")
Idents(mammal.combined2,cells=Gland) <- "Gland"
Idents(mammal.combined2,cells=EVT) <- "EVT_d14"

Idents(mammal.combined2,WhichCells(mammal.combined2,idents="EmDisc1_d14")) <- "EmDisc_d14"
Idents(mammal.combined2,WhichCells(mammal.combined2,idents="EmDisc2_d14")) <- "EmDisc_d14"

DefaultAssay(mammal.combined2) <- "RNA"

Markers <- c("DPPA3",
"JAM3",
"GATA2",
"GATA3",
"TFAP2A",
"TFAP2C",
"CGB8",
"CGB5",
"PPARG",
"CGA",
"ERVW-1",
"PRDM6",
"DIO2",
"TBX3",
"NOTUM",
"HLA-G",
"ASCL2",
"STAT1",
"HGF",
"SNAI1",
"SNAI2",
"HAND1",
"HAND2",
"BST2",
"GATA6",
"GATA4",
"CER1",
"NODAL",
"PDGFRA",
"LEFTY1",
"LEFTY1",
"APOA1",
"WNT6",
"BAMBI",
"GABRP",
"VTCN1",
"AKAP12",
"PODXL",
"SFRP1",
"SOX15",
"PDGFA",
"POU5F1")


#[1] EVT                  Am                   STB                 
# [4] EmDisc               VE                   CTB                 
# [7] Epi                  Tr                   ICM                 
#[10] Unciliated epithelia Stromal fibroblasts  Ciliated            
#[13] Hyp_d14              ExMes_d14            STB_d14             
#[16] Am/EmDisc_d14        Am_d14               EmDisc1_d14         
#[19] CTB_d14              EmDisc2_d14          EVT_d14             
#21 Levels: Epi VE STB Tr EmDisc Am ICM EVT CTB ExMes_d14 Hyp_d14 ... STB_d14
#> unique((mammal.combined2$Dataset))
#[1] "SS2 Reference 1" "10X Ours"  

Av <- AverageExpression(mammal.combined2)
Av <- Av$RNA

D <- Av[Markers,c("CTB_d14","STB_d14","EVT_d14","ExMes_d14","Hyp_d14","Am_d14","Am/EmDisc_d14","EmDisc_d14")]

mat_breaks <- seq(-2, 2, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((D),color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_selectmarkers.pdf",sep=""),scale="row",width=10,height=12)

Idents(mammal.combined2) <- paste(mammal.combined2$Dataset,Idents(mammal.combined2),sep="_")
Av <- AverageExpression(mammal.combined2)
Av <- Av$RNA

#M1 <- FindMarkers(mammal.combined2, ident.1=c("10X Ours_EmDisc_d14"),ident.2=c("10X Ours_STB_d14"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
#M2 <- FindMarkers(mammal.combined2, ident.1=c("SS2 Reference 1_EmDisc"),ident.2 = c("SS2 Reference 1_STB"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Cl1 <- M1
#Cl2 <- M2



#write.table(as.data.frame(Cl),file="ClustersEmbryonicOnly.csv",quote = FALSE, sep = ",")
AvExp3 <- as.data.frame(Av)
AvExp3$Sum <- AvExp3[,"10X Ours_Unciliated epithelia"]+AvExp3[,"10X Ours_Glandular"]+AvExp3[,"10X Ours_Ciliated"]
AvExp3$Delta <- log(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, max)+1)-log(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, min) +1)
AvExp3$Delta <- 3*(AvExp3$Delta/max(AvExp3$Delta))

AvExp3$MAX <- 100*(apply(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")], 1, max))
AvExp3$MAX[which(AvExp3$MAX<10)] <- 0.1
AvExp3$MAX[which(AvExp3$MAX>10)] <- 5
#AVE <- readRDS(file=paste(saveext,"/","ExpDataForTern.rds",sep=""))
AVEs <- Av #$RNA
list2 <- read.table("./epifactors.txt",header = F)
list2 <- read.table("./TF.txt",header = F)
SIGNAL<-read.table("./LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]
list2 <- data.frame(V1= c(SIGNAL1, SIGNAL2, SIGNAL3))
list2 <- c("SOX9","PGR","ESR1","MMP7","CPM","MKI67","HMGB2","PLAU","IL32","TNF","WNT7A","KRT7","KRT5","LGR5","IHH","EMID1","PPARG","MUC12",
"CDC20B","CCNO","HES6","FOXJ1","PIFO","TP73","HEY1","ABCG1","SCGB2A2","C2CD4A","SLC18A2","PAEP","CXCL14","SPP1","DPP4","PtSG1",
"CLD22","PAX2","VTNC1","SLC26A7","MSLN")
AVE2 <-  as.data.frame(AvExp3[,c("10X Ours_Unciliated epithelia","10X Ours_Glandular","10X Ours_Ciliated")])
colnames(AVE2) <- c("Epit","Gland","Cil")
AVE2$Sum <- log2(AvExp3[rownames(AVE2),"Sum"]+1)
AVE2$Delta <- AvExp3[rownames(AVE2),"Delta"]
AVE2$MAX <- AvExp3[rownames(AVE2),"MAX"]
AVE2$Sum[which(is.na(AVE2$Sum)==1)] <- 0
AVE2$Delta[which(is.na(AVE2$Delta)==1)] <- 0
AVE2$Delta <- 5*(AVE2$Delta/max(AVE2$Delta))
AVE2$Sum <- 5*(AVE2$Sum/max(AVE2$Sum))
AVE2<- AVE2[which(AVE2$MAX>1),]

AVE2 <- AVE2[intersect(list2,rownames(AVE2)),]
#AVE2 <- AVE2[intersect(list2$V1,rownames(AVE2)),]

library(ggtern)
p1 <- ggtern(data=AVE2, aes(x =Epit,z=Gland,y=Cil)) + geom_point(size=AVE2$MAX, color = "black", position= position_jitter_tern(x=0.01, y=0.01, z=0.01))
p1 <- p1+ theme_showgrid()
#theme_custom(col.grid.minor = black")+theme_showgrid() ; p1
p1 <- p1 + annotate(geom  = 'text', x = (AVE2$Epit)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),
z= (AVE2$Gland)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),
y= (AVE2$Cil)/(AVE2$Epit+AVE2$Gland+AVE2$Cil),label = rownames(AVE2),color = c("black"))
ggsave(filename=paste(saveext,"GGTern_Marklers.pdf",sep=""),width = 20, height = 20, plot = p1)

