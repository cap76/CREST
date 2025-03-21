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
library("ggtern")

saveext = "./FinalAlignB/"
dir.create(saveext)


mammal.combined <- readRDS(file = paste(saveext,"Seurat_combined_FinalAno.rds",sep=""))
#mammal.combined4C <- subset(mammal.combined,idents=c("EVT_d14","CTB_d14","STB_d14","ExMes_d14"),invert=TRUE)

mammal.combined4C <- subset(mammal.combined,idents=c("EmD/Hyp","Hyp_d14","Am_d14","EmDisc1_d14","EmDisc2_d14","Am/EmDisc_d14"))


D <- GetAssayData(mammal.combined4C,assay="RNA")
Tern <- as.data.frame(t(exp(D)-1))
saveRDS(Tern,'Data4Tern.rds')


Tern$POU5F1 <- (Tern$POU5F1-min(Tern$POU5F1))/max(Tern$POU5F1-min(Tern$POU5F1))
Tern$WNT6 <- (Tern$WNT6-min(Tern$WNT6))/max(Tern$WNT6-min(Tern$WNT6))
Tern$GATA4 <- (Tern$GATA4-min(Tern$GATA4))/max(Tern$GATA4-min(Tern$GATA4))

Tern$SOX17 <- (Tern$SOX17-min(Tern$SOX17))/max(Tern$SOX17-min(Tern$SOX17))
Tern$PDGFA <- (Tern$PDGFA-min(Tern$PDGFA))/max(Tern$PDGFA-min(Tern$PDGFA))
Tern$PDGFRA <- (Tern$PDGFRA-min(Tern$PDGFRA))/max(Tern$PDGFRA-min(Tern$PDGFRA))

#Tern[] <- lapply(Tern, jitter, 3)

p1 <- ggtern(data=Tern, aes(x = POU5F1,z= WNT6,y= GATA4)) + geom_point(color = "black",size=4, position= position_jitter_tern(x=0.1, y=0.1, z=0.1)) 
p1 <- p1+ theme_showgrid()
#theme_custom(col.grid.minor = black")+theme_showgrid() ; p1
#p1 <- p1 + annotate(geom  = 'text', x = (Tern$POU5F1+0.1)/(Tern$POU5F1+Tern$WNT6+Tern$GATA4),
#z= (Tern$WNT6+0.1)/(Tern$WNT6+Tern$GATA4+Tern$POU5F1),
#y= (Tern$GATA4+0.1)/(Tern$POU5F1+Tern$GATA4+Tern$WNT6),label = rownames(Tern),color = c("green")) 
ggsave(filename=paste(saveext,"GGTern_WNT6_GATA4_POU5F1.pdf",sep=""),width = 20, height = 20, plot = p1)

#p1 <- ggtern(data=D1, aes(x=Newworld,y=Oldworld, z=Ape)) + geom_point(aes(color=colsi)) + scale_color_manual(values=c("grey", "black")) 

p1 <- ggtern(data=Tern,aes(x = POU5F1,z= WNT6,y= GATA4)) + geom_point(color = "black",size=4,position= position_jitter_tern(x=0.1, y=0.1, z=0.1)) #+ geom_segment(size=.5,arrow=arrow(),color="black")
p1 <- p1+theme_showgrid()
#p1 <- p1 + annotate(geom  = 'text', x = Tern$POU5F1/(Tern$POU5F1+Tern$WNT6+Tern$GATA4),
#z= Tern$WNT6/(Tern$WNT6+Tern$GATA4+Tern$POU5F1),
#y= Tern$GATA4/(Tern$POU5F1+Tern$GATA4+Tern$WNT6),label = rownames(Tern),color = c("red"))
ggsave(filename=paste(saveext,"GGTern_WNT6_GATA4_POU5F1.pdf",sep=""),width = 20, height = 20, plot = p1)

p1 <- ggtern(data=Tern,aes(x = POU5F1,z= WNT6,y= GATA4)) + geom_point(color = "black",size=4,position= position_jitter_tern(x=0.1, y=0.1, z=0.1)) #+ geom_segment(size=.5,arrow=arrow(),color="black")
p1 <- p1+theme_showgrid()
#p1 <- p1 + annotate(geom  = 'text', x = Tern$POU5F1/(Tern$POU5F1+Tern$WNT6+Tern$GATA4),
#z= Tern$WNT6/(Tern$WNT6+Tern$GATA4+Tern$POU5F1),
#y= Tern$GATA4/(Tern$POU5F1+Tern$GATA4+Tern$WNT6),label = rownames(Tern),color = c("red"))
ggsave(filename=paste(saveext,"GGTern_WNT6_GATA4_POU5F1.pdf",sep=""),width = 20, height = 20, plot = p1)



p1 <- ggtern(data=Tern,aes(x = PDGFA,z= WNT6,y= PDGFRA)) + geom_point(color = "black",size=4,position= position_jitter_tern(x=0.1, y=0.1, z=0.1)) #+ geom_segment(size=.5,arrow=arrow(),color="black")
p1 <- p1+theme_showgrid()

#p1 <- p1 + annotate(geom  = 'text', x = Tern$PDGFA/(Tern$PDGFA+Tern$WNT6+Tern$PDGFRA),
#z= Tern$WNT6/(Tern$WNT6+Tern$PDGFRA+Tern$PDGFA),
#y= Tern$PDGFRA/(Tern$PDGFA+Tern$PDGFRA+Tern$WNT6),label = rownames(Tern),color = c("red"))
ggsave(filename=paste(saveext,"GGTern_WNT6_PDGFRA_PDGFA.pdf",sep=""),width = 20, height = 20, plot = p1)



p1 <- ggtern(data=Tern,aes(x = POU5F1,z= WNT6,y= SOX17)) + geom_point(color = "black",size=4,position= position_jitter_tern(x=0.1, y=0.1, z=0.1)) #+ geom_segment(size=.5,arrow=arrow(),color="black")
p1 <- p1+theme_showgrid()

#p1 <- p1 + annotate(geom  = 'text', x = Tern$POU5F1/(Tern$POU5F1+Tern$WNT6+Tern$SOX17),
#z= Tern$WNT6/(Tern$WNT6+Tern$SOX17+Tern$POU5F1),
#y= Tern$SOX17/(Tern$POU5F1+Tern$SOX17+Tern$WNT6),label = rownames(Tern),color = c("red"))
ggsave(filename=paste(saveext,"GGTern_WNT6_SOX17_POU5F1.pdf",sep=""),width = 20, height = 20, plot = p1)



p1 <- ggtern(data=Tern,aes(x = PDGFA,z= WNT6,y= SOX17)) + geom_point(color = "black",size=4,position= position_jitter_tern(x=0.1, y=0.1, z=0.1)) #+ geom_segment(size=.5,arrow=arrow(),color="black")
p1 <- p1+theme_showgrid()

#p1 <- p1 + annotate(geom  = 'text', x = Tern$PDGFA/(Tern$PDGFA+Tern$WNT6+Tern$SOX17),
#z= Tern$WNT6/(Tern$WNT6+Tern$SOX17+Tern$PDGFA),
#y= Tern$SOX17/(Tern$PDGFA+Tern$SOX17+Tern$WNT6),label = rownames(Tern),color = c("red"))
ggsave(filename=paste(saveext,"GGTern_WNT6_SOX17_PDGFA.pdf",sep=""),width = 20, height = 20, plot = p1)


FeaturePlot(mammal.combined4C,feature="WNT6",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_WNT6_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)



FeaturePlot(mammal.combined4C,feature="FOXF1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_FOXF1_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="MESP1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_MESP1_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)








DefaultAssay(mammal.combined2) <- "RNA"
FeaturePlot(mammal.combined2,feature="WNT6",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_WNT6_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="WNT6",split.by = "Dataset",reduction = "umap",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalUMAP_WNT6_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="LGALS9",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_LGALS9_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="HAND2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_HAND2_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="CD44",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_CD44_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="SFRP1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SFRP1_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="WNT5A",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_WNT5A_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="RARRES2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_RARRES2_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="GPR1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_GPR1_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="C5AR1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_C5AR1_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="NRG1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_NRG1_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="ERBB4",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_ERBB4_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="EGFR",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_EGFR_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="COPA",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_COPA_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="SFRP1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SFRP1_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="HLA-C",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_HLA-C_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="ERVW-1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_ERVW1_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="EPCAM",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_EPCAM_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="MCAM",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_MCAM_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="MAOA",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_MAOA_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="CCL7",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_CCL7_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="NETO1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_NETO1_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="NRG1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_NRG1_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="CCL3L1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_CCL3L1_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="AGR3",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_AGR3_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="ERBB4",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_ERBB4_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="SCGB2A1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SCGB2A1_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="PADI1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_PADI1_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="SCGB2A1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SCGB2A1_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="MUC5B",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_MUC5B_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="ESR1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_ESR1_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="FGF14",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_FGF14_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="IL24",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_IL24_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="GPR32",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_GPR32_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="CGB8",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_CGB8_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="CGB5",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_CGB5_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="PSG3",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_PSG3_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="FGF2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_FGF2_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="FGF7",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_FGF7_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)




asdasdsadsaas
s
dasdsadadsadsadsa
FeaturePlot(mammal.combined, reduction = "pca", features = "percent.mito2", pt.size = 2, cols = c("lightgrey", "black")) # + xlim(c(-13, 13)) + ylim(c(-13, 13))#
ggsave(filename=paste(saveext,"/Markers/","AraMarkers", "_mito_nosplit.pdf",sep=""),width = 10, height = 10,limitsize = FALSE)

VlnPlot(mammal.combined,features=c("nCount_RNA"),split.by="ID3",group.by="ID3")
ggsave(filename=paste(saveext,"/Markers/","Violin_", "_nCount_RNA.pdf",sep=""),width = 10, height = 15,limitsize = FALSE)
VlnPlot(mammal.combined,features=c("nFeature_RNA"),split.by="ID3",group.by="ID3")
ggsave(filename=paste(saveext,"/Markers/","Violin_", "nFeature_RNA.pdf",sep=""),width = 10, height = 15,limitsize = FALSE)

VlnPlot(mammal.combined,features=c("percent.mito2"),split.by="ID3",group.by="ID3")
ggsave(filename=paste(saveext,"/Markers/","Violin_", "mito_RNA.pdf",sep=""),width = 10, height = 15,limitsize = FALSE)


#sadsadsad

mammal.combined1 <- readRDS(file = paste(saveext,"Seurat_combined_justEmbryonic_filtered.rds",sep=""))
#mammal.combined2 <- readRDS(file = paste(saveext,"Seurat_combined3.rds",sep=""))

uID <- Idents(mammal.combined1) #as.character(mammal.combined1$Cells)
uID[which(uID%in%c("CTb_CS4","EPI_CS4","TB_CS3"))] <- "Tr"
uID[which(uID%in%c("CTb_CS5A/B","CTb_CS5C","CTb_CS6"))] <- "CTB"
uID[which(uID%in%c("EmDisc_CS5A/B","EmDisc_CS5C","EmDisc_CS6","EmDiscPS_CS6"))] <- "EmDisc"
uID[which(uID%in%c("EVTb_CS5C","EVTb_CS6"))] <- "EVT"
uID[which(uID%in%c("STb_CS5A/B","STb_CS5C","STb_CS6"))] <- "STB"
uID[which(uID%in%c("VE_CS4","VE_CS5","VE_CS5A/B","VE_CS5C","SYS_CS6"))] <- "VE"
uID[which(uID%in%c("ICM_CS3"))] <- "ICM"

uID[which(Idents(mammal.combined1)%in%c("STB_d9","STB_d11","STB_d12"))] <- "CTB"
uID[which(Idents(mammal.combined1)%in%c("CTB_d9","CTB_d11","CTB_d12"))] <- "STB"
uID[which(uID%in%c("Hyp_d9","Hyp_d11","Hyp_d12"))] <- "VE"
uID[which(uID%in%c("EmDisc_d9","EmDisc_d11","EmDisc_d12"))] <- "EmDisc"

#uID[which(uID%in%c("Unciliated epithelia 1","Unciliated epithelia 2","Cil"))] <- "Unciliated epithelia"
#uID[which(uID%in%c("Cil"))] <- "Ciliated epithelia"
#uID[which(uID%in%c("ICM_CS3"))] <- "ICM"
#uID[which(uID%in%c("ICM_CS3"))] <- "ICM"
Idents(mammal.combined1) <- uID

AmList <- c("SRR10039958Aligned.sortedByCoord.out.bam_1",
"SRR10039959Aligned.sortedByCoord.out.bam_1",
"SRR10039969Aligned.sortedByCoord.out.bam_1",
"SRR10039972Aligned.sortedByCoord.out.bam_1",
"SRR10039973Aligned.sortedByCoord.out.bam_1",
"SRR10039974Aligned.sortedByCoord.out.bam_1",
"SRR10039975Aligned.sortedByCoord.out.bam_1",
"SRR10039980Aligned.sortedByCoord.out.bam_1",
"SRR10039983Aligned.sortedByCoord.out.bam_1",
"SRR10040003Aligned.sortedByCoord.out.bam_1",
"SRR10040011Aligned.sortedByCoord.out.bam_1")

CTBList <- c("SRR10039789Aligned.sortedByCoord.out.bam_1",
"SRR10039560Aligned.sortedByCoord.out.bam_1",
"SRR10039561Aligned.sortedByCoord.out.bam_1",
"SRR10039564Aligned.sortedByCoord.out.bam_1",
"SRR10039572Aligned.sortedByCoord.out.bam_1",
"SRR10039574Aligned.sortedByCoord.out.bam_1",
"SRR10039578Aligned.sortedByCoord.out.bam_1",
"SRR10039580Aligned.sortedByCoord.out.bam_1",
"SRR10039582Aligned.sortedByCoord.out.bam_1",
"SRR10039586Aligned.sortedByCoord.out.bam_1",
"SRR10039587Aligned.sortedByCoord.out.bam_1",
"SRR10039588Aligned.sortedByCoord.out.bam_1",
"SRR10039589Aligned.sortedByCoord.out.bam_1",
"SRR10039591Aligned.sortedByCoord.out.bam_1",
"SRR10039592Aligned.sortedByCoord.out.bam_1",
"SRR10039593Aligned.sortedByCoord.out.bam_1",
"SRR10039594Aligned.sortedByCoord.out.bam_1",
"SRR10039595Aligned.sortedByCoord.out.bam_1",
"SRR10039596Aligned.sortedByCoord.out.bam_1",
"SRR10039600Aligned.sortedByCoord.out.bam_1",
"SRR10039601Aligned.sortedByCoord.out.bam_1",
"SRR10039602Aligned.sortedByCoord.out.bam_1",
"SRR10039603Aligned.sortedByCoord.out.bam_1",
"SRR10039605Aligned.sortedByCoord.out.bam_1",
"SRR10039607Aligned.sortedByCoord.out.bam_1",
"SRR10039608Aligned.sortedByCoord.out.bam_1",
"SRR10039628Aligned.sortedByCoord.out.bam_1",
"SRR10039630Aligned.sortedByCoord.out.bam_1",
"SRR10039632Aligned.sortedByCoord.out.bam_1",
"SRR10039633Aligned.sortedByCoord.out.bam_1",
"SRR10039635Aligned.sortedByCoord.out.bam_1",
"SRR10039637Aligned.sortedByCoord.out.bam_1",
"SRR10039638Aligned.sortedByCoord.out.bam_1",
"SRR10039640Aligned.sortedByCoord.out.bam_1",
"SRR10039642Aligned.sortedByCoord.out.bam_1",
"SRR10039645Aligned.sortedByCoord.out.bam_1",
"SRR10039646Aligned.sortedByCoord.out.bam_1",
"SRR10039655Aligned.sortedByCoord.out.bam_1",
"SRR10039658Aligned.sortedByCoord.out.bam_1",
"SRR10039659Aligned.sortedByCoord.out.bam_1",
"SRR10039664Aligned.sortedByCoord.out.bam_1",
"SRR10039670Aligned.sortedByCoord.out.bam_1",
"SRR10039672Aligned.sortedByCoord.out.bam_1",
"SRR10039673Aligned.sortedByCoord.out.bam_1",
"SRR10039674Aligned.sortedByCoord.out.bam_1",
"SRR10039675Aligned.sortedByCoord.out.bam_1",
"SRR10039676Aligned.sortedByCoord.out.bam_1",
"SRR10039677Aligned.sortedByCoord.out.bam_1",
"SRR10039678Aligned.sortedByCoord.out.bam_1",
"SRR10039679Aligned.sortedByCoord.out.bam_1",
"SRR10039680Aligned.sortedByCoord.out.bam_1",
"SRR10039681Aligned.sortedByCoord.out.bam_1",
"SRR10039683Aligned.sortedByCoord.out.bam_1",
"SRR10039685Aligned.sortedByCoord.out.bam_1",
"SRR10039691Aligned.sortedByCoord.out.bam_1",
"SRR10039692Aligned.sortedByCoord.out.bam_1",
"SRR10039693Aligned.sortedByCoord.out.bam_1",
"SRR10039695Aligned.sortedByCoord.out.bam_1",
"SRR10039696Aligned.sortedByCoord.out.bam_1",
"SRR10039697Aligned.sortedByCoord.out.bam_1",
"SRR10039700Aligned.sortedByCoord.out.bam_1",
"SRR10039701Aligned.sortedByCoord.out.bam_1",
"SRR10039702Aligned.sortedByCoord.out.bam_1",
"SRR10039703Aligned.sortedByCoord.out.bam_1",
"SRR10039704Aligned.sortedByCoord.out.bam_1",
"SRR10039705Aligned.sortedByCoord.out.bam_1",
"SRR10039706Aligned.sortedByCoord.out.bam_1",
"SRR10039712Aligned.sortedByCoord.out.bam_1",
"SRR10039713Aligned.sortedByCoord.out.bam_1",
"SRR10039718Aligned.sortedByCoord.out.bam_1",
"SRR10039719Aligned.sortedByCoord.out.bam_1",
"SRR10039721Aligned.sortedByCoord.out.bam_1",
"SRR10039722Aligned.sortedByCoord.out.bam_1",
"SRR10039724Aligned.sortedByCoord.out.bam_1",
"SRR10039725Aligned.sortedByCoord.out.bam_1",
"SRR10039726Aligned.sortedByCoord.out.bam_1",
"SRR10039727Aligned.sortedByCoord.out.bam_1",
"SRR10039728Aligned.sortedByCoord.out.bam_1",
"SRR10039731Aligned.sortedByCoord.out.bam_1",
"SRR10039733Aligned.sortedByCoord.out.bam_1",
"SRR10039738Aligned.sortedByCoord.out.bam_1",
"SRR10039740Aligned.sortedByCoord.out.bam_1",
"SRR10039762Aligned.sortedByCoord.out.bam_1",
"SRR10039763Aligned.sortedByCoord.out.bam_1",
"SRR10039765Aligned.sortedByCoord.out.bam_1",
"SRR10039799Aligned.sortedByCoord.out.bam_1",
"SRR10039801Aligned.sortedByCoord.out.bam_1",
"SRR10039806Aligned.sortedByCoord.out.bam_1",
"SRR10039810Aligned.sortedByCoord.out.bam_1",
"SRR10039812Aligned.sortedByCoord.out.bam_1")


eEVTList <- c("SRR10039652Aligned.sortedByCoord.out.bam_1",
"SRR10039708Aligned.sortedByCoord.out.bam_1",
"SRR10039720Aligned.sortedByCoord.out.bam_1",
"SRR10039729Aligned.sortedByCoord.out.bam_1",
"SRR10039736Aligned.sortedByCoord.out.bam_1",
"SRR10039739Aligned.sortedByCoord.out.bam_1",
"SRR10039741Aligned.sortedByCoord.out.bam_1",
"SRR10039749Aligned.sortedByCoord.out.bam_1",
"SRR10039753Aligned.sortedByCoord.out.bam_1",
"SRR10039761Aligned.sortedByCoord.out.bam_1",
"SRR10039772Aligned.sortedByCoord.out.bam_1",
"SRR10039777Aligned.sortedByCoord.out.bam_1",
"SRR10039779Aligned.sortedByCoord.out.bam_1",
"SRR10039784Aligned.sortedByCoord.out.bam_1",
"SRR10039785Aligned.sortedByCoord.out.bam_1",
"SRR10039790Aligned.sortedByCoord.out.bam_1",
"SRR10039791Aligned.sortedByCoord.out.bam_1",
"SRR10039793Aligned.sortedByCoord.out.bam_1",
"SRR10039794Aligned.sortedByCoord.out.bam_1",
"SRR10039796Aligned.sortedByCoord.out.bam_1",
"SRR10039798Aligned.sortedByCoord.out.bam_1",
"SRR10039802Aligned.sortedByCoord.out.bam_1",
"SRR10039803Aligned.sortedByCoord.out.bam_1",
"SRR10039804Aligned.sortedByCoord.out.bam_1",
"SRR10039805Aligned.sortedByCoord.out.bam_1",
"SRR10039811Aligned.sortedByCoord.out.bam_1",
"SRR10039820Aligned.sortedByCoord.out.bam_1",
"SRR10039821Aligned.sortedByCoord.out.bam_1",
"SRR10039838Aligned.sortedByCoord.out.bam_1",
"SRR10039845Aligned.sortedByCoord.out.bam_1",
"SRR10039846Aligned.sortedByCoord.out.bam_1",
"SRR10039851Aligned.sortedByCoord.out.bam_1",
"SRR10039852Aligned.sortedByCoord.out.bam_1",
"SRR10039853Aligned.sortedByCoord.out.bam_1",
"SRR10039855Aligned.sortedByCoord.out.bam_1",
"SRR10039856Aligned.sortedByCoord.out.bam_1",
"SRR10039860Aligned.sortedByCoord.out.bam_1",
"SRR10039863Aligned.sortedByCoord.out.bam_1",
"SRR10039867Aligned.sortedByCoord.out.bam_1",
"SRR10039870Aligned.sortedByCoord.out.bam_1",
"SRR10039871Aligned.sortedByCoord.out.bam_1",
"SRR10039872Aligned.sortedByCoord.out.bam_1",
"SRR10039874Aligned.sortedByCoord.out.bam_1",
"SRR10039875Aligned.sortedByCoord.out.bam_1",
"SRR10039881Aligned.sortedByCoord.out.bam_1",
"SRR10039882Aligned.sortedByCoord.out.bam_1",
"SRR10039883Aligned.sortedByCoord.out.bam_1",
"SRR10039888Aligned.sortedByCoord.out.bam_1",
"SRR10039891Aligned.sortedByCoord.out.bam_1",
"SRR10039894Aligned.sortedByCoord.out.bam_1",
"SRR10039895Aligned.sortedByCoord.out.bam_1",
"SRR10039897Aligned.sortedByCoord.out.bam_1",
"SRR10039900Aligned.sortedByCoord.out.bam_1",
"SRR10039927Aligned.sortedByCoord.out.bam_1",
"SRR10040022Aligned.sortedByCoord.out.bam_1",
"SRR10040034Aligned.sortedByCoord.out.bam_1",
"SRR10040043Aligned.sortedByCoord.out.bam_1",
"SRR10040045Aligned.sortedByCoord.out.bam_1",
"SRR10040049Aligned.sortedByCoord.out.bam_1")


eSTBList<- c("SRR10039954Aligned.sortedByCoord.out.bam_1",
"SRR10039636Aligned.sortedByCoord.out.bam_1",
"SRR10039639Aligned.sortedByCoord.out.bam_1",
"SRR10039649Aligned.sortedByCoord.out.bam_1",
"SRR10039650Aligned.sortedByCoord.out.bam_1",
"SRR10039654Aligned.sortedByCoord.out.bam_1",
"SRR10039660Aligned.sortedByCoord.out.bam_1",
"SRR10039662Aligned.sortedByCoord.out.bam_1",
"SRR10039663Aligned.sortedByCoord.out.bam_1",
"SRR10039665Aligned.sortedByCoord.out.bam_1",
"SRR10039666Aligned.sortedByCoord.out.bam_1",
"SRR10039667Aligned.sortedByCoord.out.bam_1",
"SRR10039668Aligned.sortedByCoord.out.bam_1",
"SRR10039671Aligned.sortedByCoord.out.bam_1",
"SRR10039684Aligned.sortedByCoord.out.bam_1",
"SRR10039687Aligned.sortedByCoord.out.bam_1",
"SRR10039688Aligned.sortedByCoord.out.bam_1",
"SRR10039698Aligned.sortedByCoord.out.bam_1",
"SRR10039699Aligned.sortedByCoord.out.bam_1",
"SRR10039709Aligned.sortedByCoord.out.bam_1",
"SRR10039710Aligned.sortedByCoord.out.bam_1",
"SRR10039711Aligned.sortedByCoord.out.bam_1",
"SRR10039714Aligned.sortedByCoord.out.bam_1",
"SRR10039716Aligned.sortedByCoord.out.bam_1",
"SRR10039732Aligned.sortedByCoord.out.bam_1",
"SRR10039735Aligned.sortedByCoord.out.bam_1",
"SRR10039737Aligned.sortedByCoord.out.bam_1",
"SRR10039754Aligned.sortedByCoord.out.bam_1",
"SRR10039756Aligned.sortedByCoord.out.bam_1",
"SRR10039757Aligned.sortedByCoord.out.bam_1",
"SRR10039758Aligned.sortedByCoord.out.bam_1",
"SRR10039759Aligned.sortedByCoord.out.bam_1",
"SRR10039760Aligned.sortedByCoord.out.bam_1",
"SRR10039764Aligned.sortedByCoord.out.bam_1",
"SRR10039768Aligned.sortedByCoord.out.bam_1",
"SRR10039774Aligned.sortedByCoord.out.bam_1",
"SRR10039776Aligned.sortedByCoord.out.bam_1",
"SRR10039782Aligned.sortedByCoord.out.bam_1",
"SRR10039787Aligned.sortedByCoord.out.bam_1",
"SRR10039792Aligned.sortedByCoord.out.bam_1",
"SRR10039795Aligned.sortedByCoord.out.bam_1",
"SRR10039797Aligned.sortedByCoord.out.bam_1",
"SRR10039807Aligned.sortedByCoord.out.bam_1",
"SRR10039808Aligned.sortedByCoord.out.bam_1",
"SRR10039809Aligned.sortedByCoord.out.bam_1",
"SRR10039813Aligned.sortedByCoord.out.bam_1",
"SRR10039814Aligned.sortedByCoord.out.bam_1",
"SRR10039815Aligned.sortedByCoord.out.bam_1",
"SRR10039817Aligned.sortedByCoord.out.bam_1",
"SRR10039819Aligned.sortedByCoord.out.bam_1",
"SRR10039830Aligned.sortedByCoord.out.bam_1",
"SRR10039854Aligned.sortedByCoord.out.bam_1",
"SRR10039866Aligned.sortedByCoord.out.bam_1",
"SRR10039873Aligned.sortedByCoord.out.bam_1",
"SRR10039911Aligned.sortedByCoord.out.bam_1",
"SRR10039948Aligned.sortedByCoord.out.bam_1",
"SRR10039985Aligned.sortedByCoord.out.bam_1",
"SRR10039986Aligned.sortedByCoord.out.bam_1")

EVTList <- c("SRR10039952Aligned.sortedByCoord.out.bam_1",
"SRR10039818Aligned.sortedByCoord.out.bam_1",
"SRR10039910Aligned.sortedByCoord.out.bam_1",
"SRR10039912Aligned.sortedByCoord.out.bam_1",
"SRR10039913Aligned.sortedByCoord.out.bam_1",
"SRR10039916Aligned.sortedByCoord.out.bam_1",
"SRR10039918Aligned.sortedByCoord.out.bam_1",
"SRR10039920Aligned.sortedByCoord.out.bam_1",
"SRR10039922Aligned.sortedByCoord.out.bam_1",
"SRR10039923Aligned.sortedByCoord.out.bam_1",
"SRR10039929Aligned.sortedByCoord.out.bam_1",
"SRR10039930Aligned.sortedByCoord.out.bam_1",
"SRR10039938Aligned.sortedByCoord.out.bam_1",
"SRR10039940Aligned.sortedByCoord.out.bam_1",
"SRR10040013Aligned.sortedByCoord.out.bam_1",
"SRR10040014Aligned.sortedByCoord.out.bam_1",
"SRR10040016Aligned.sortedByCoord.out.bam_1",
"SRR10040017Aligned.sortedByCoord.out.bam_1",
"SRR10040020Aligned.sortedByCoord.out.bam_1",
"SRR10040021Aligned.sortedByCoord.out.bam_1",
"SRR10040023Aligned.sortedByCoord.out.bam_1",
"SRR10040025Aligned.sortedByCoord.out.bam_1",
"SRR10040026Aligned.sortedByCoord.out.bam_1",
"SRR10040027Aligned.sortedByCoord.out.bam_1",
"SRR10040028Aligned.sortedByCoord.out.bam_1",
"SRR10040029Aligned.sortedByCoord.out.bam_1",
"SRR10040031Aligned.sortedByCoord.out.bam_1",
"SRR10040032Aligned.sortedByCoord.out.bam_1",
"SRR10040033Aligned.sortedByCoord.out.bam_1",
"SRR10040036Aligned.sortedByCoord.out.bam_1",
"SRR10040037Aligned.sortedByCoord.out.bam_1",
"SRR10040038Aligned.sortedByCoord.out.bam_1",
"SRR10040040Aligned.sortedByCoord.out.bam_1",
"SRR10040042Aligned.sortedByCoord.out.bam_1",
"SRR10040044Aligned.sortedByCoord.out.bam_1",
"SRR10040046Aligned.sortedByCoord.out.bam_1",
"SRR10040048Aligned.sortedByCoord.out.bam_1",
"SRR10040050Aligned.sortedByCoord.out.bam_1",
"SRR10040051Aligned.sortedByCoord.out.bam_1")

ICMList <- c("SRR10039497Aligned.sortedByCoord.out.bam_1",
"SRR10039501Aligned.sortedByCoord.out.bam_1",
"SRR10039503Aligned.sortedByCoord.out.bam_1",
"SRR10039504Aligned.sortedByCoord.out.bam_1",
"SRR10039505Aligned.sortedByCoord.out.bam_1",
"SRR10039506Aligned.sortedByCoord.out.bam_1",
"SRR10039508Aligned.sortedByCoord.out.bam_1",
"SRR10039509Aligned.sortedByCoord.out.bam_1",
"SRR10039512Aligned.sortedByCoord.out.bam_1",
"SRR10039513Aligned.sortedByCoord.out.bam_1",
"SRR10039514Aligned.sortedByCoord.out.bam_1",
"SRR10039519Aligned.sortedByCoord.out.bam_1",
"SRR10039526Aligned.sortedByCoord.out.bam_1",
"SRR10039527Aligned.sortedByCoord.out.bam_1",
"SRR10039530Aligned.sortedByCoord.out.bam_1",
"SRR10039539Aligned.sortedByCoord.out.bam_1",
"SRR10039540Aligned.sortedByCoord.out.bam_1",
"SRR10039544Aligned.sortedByCoord.out.bam_1",
"SRR10039545Aligned.sortedByCoord.out.bam_1")

PostEpiAME <- c("SRR10039869Aligned.sortedByCoord.out.bam_1",
"SRR10039887Aligned.sortedByCoord.out.bam_1",
"SRR10039890Aligned.sortedByCoord.out.bam_1",
"SRR10039987Aligned.sortedByCoord.out.bam_1",
"SRR10039991Aligned.sortedByCoord.out.bam_1",
"SRR10040000Aligned.sortedByCoord.out.bam_1",
"SRR10040001Aligned.sortedByCoord.out.bam_1",
"SRR10040004Aligned.sortedByCoord.out.bam_1",
"SRR10040041Aligned.sortedByCoord.out.bam_1",
"SRR10040047Aligned.sortedByCoord.out.bam_1")

PostEPIE1 <- c("SRR10039641Aligned.sortedByCoord.out.bam_1",
"SRR10039644Aligned.sortedByCoord.out.bam_1",
"SRR10039647Aligned.sortedByCoord.out.bam_1",
"SRR10039651Aligned.sortedByCoord.out.bam_1",
"SRR10039653Aligned.sortedByCoord.out.bam_1",
"SRR10039656Aligned.sortedByCoord.out.bam_1",
"SRR10039657Aligned.sortedByCoord.out.bam_1",
"SRR10039682Aligned.sortedByCoord.out.bam_1",
"SRR10039690Aligned.sortedByCoord.out.bam_1",
"SRR10039694Aligned.sortedByCoord.out.bam_1",
"SRR10039715Aligned.sortedByCoord.out.bam_1",
"SRR10039723Aligned.sortedByCoord.out.bam_1",
"SRR10039730Aligned.sortedByCoord.out.bam_1",
"SRR10039745Aligned.sortedByCoord.out.bam_1",
"SRR10039746Aligned.sortedByCoord.out.bam_1",
"SRR10039748Aligned.sortedByCoord.out.bam_1",
"SRR10039766Aligned.sortedByCoord.out.bam_1",
"SRR10039767Aligned.sortedByCoord.out.bam_1",
"SRR10039769Aligned.sortedByCoord.out.bam_1",
"SRR10039770Aligned.sortedByCoord.out.bam_1",
"SRR10039771Aligned.sortedByCoord.out.bam_1",
"SRR10039773Aligned.sortedByCoord.out.bam_1",
"SRR10039788Aligned.sortedByCoord.out.bam_1",
"SRR10039824Aligned.sortedByCoord.out.bam_1",
"SRR10039825Aligned.sortedByCoord.out.bam_1",
"SRR10039826Aligned.sortedByCoord.out.bam_1",
"SRR10039827Aligned.sortedByCoord.out.bam_1",
"SRR10039844Aligned.sortedByCoord.out.bam_1",
"SRR10039868Aligned.sortedByCoord.out.bam_1",
"SRR10039909Aligned.sortedByCoord.out.bam_1")

PostEPIE2 <- c("SRR10039832Aligned.sortedByCoord.out.bam_1",
"SRR10039847Aligned.sortedByCoord.out.bam_1",
"SRR10039849Aligned.sortedByCoord.out.bam_1",
"SRR10039850Aligned.sortedByCoord.out.bam_1",
"SRR10039858Aligned.sortedByCoord.out.bam_1",
"SRR10039859Aligned.sortedByCoord.out.bam_1",
"SRR10039861Aligned.sortedByCoord.out.bam_1",
"SRR10039862Aligned.sortedByCoord.out.bam_1",
"SRR10039865Aligned.sortedByCoord.out.bam_1",
"SRR10039886Aligned.sortedByCoord.out.bam_1",
"SRR10039898Aligned.sortedByCoord.out.bam_1",
"SRR10039901Aligned.sortedByCoord.out.bam_1",
"SRR10039908Aligned.sortedByCoord.out.bam_1",
"SRR10039915Aligned.sortedByCoord.out.bam_1",
"SRR10039924Aligned.sortedByCoord.out.bam_1",
"SRR10039925Aligned.sortedByCoord.out.bam_1",
"SRR10039949Aligned.sortedByCoord.out.bam_1",
"SRR10039984Aligned.sortedByCoord.out.bam_1",
"SRR10039989Aligned.sortedByCoord.out.bam_1",
"SRR10039993Aligned.sortedByCoord.out.bam_1",
"SRR10039994Aligned.sortedByCoord.out.bam_1",
"SRR10039995Aligned.sortedByCoord.out.bam_1",
"SRR10039996Aligned.sortedByCoord.out.bam_1",
"SRR10039997Aligned.sortedByCoord.out.bam_1",
"SRR10039999Aligned.sortedByCoord.out.bam_1",
"SRR10040005Aligned.sortedByCoord.out.bam_1",
"SRR10040006Aligned.sortedByCoord.out.bam_1",
"SRR10040008Aligned.sortedByCoord.out.bam_1")

PostEpiGast <- c("SRR10039934Aligned.sortedByCoord.out.bam_1",
"SRR10039936Aligned.sortedByCoord.out.bam_1",
"SRR10039937Aligned.sortedByCoord.out.bam_1",
"SRR10039951Aligned.sortedByCoord.out.bam_1",
"SRR10039960Aligned.sortedByCoord.out.bam_1",
"SRR10039962Aligned.sortedByCoord.out.bam_1",
"SRR10039963Aligned.sortedByCoord.out.bam_1",
"SRR10039966Aligned.sortedByCoord.out.bam_1",
"SRR10039968Aligned.sortedByCoord.out.bam_1",
"SRR10039971Aligned.sortedByCoord.out.bam_1",
"SRR10039976Aligned.sortedByCoord.out.bam_1",
"SRR10039978Aligned.sortedByCoord.out.bam_1",
"SRR10039981Aligned.sortedByCoord.out.bam_1",
"SRR10039982Aligned.sortedByCoord.out.bam_1",
"SRR10040009Aligned.sortedByCoord.out.bam_1",
"SRR10040012Aligned.sortedByCoord.out.bam_1")

PrEn <- c("SRR10039932Aligned.sortedByCoord.out.bam_1",
"SRR10039953Aligned.sortedByCoord.out.bam_1",
"SRR10039957Aligned.sortedByCoord.out.bam_1",
"SRR10039961Aligned.sortedByCoord.out.bam_1",
"SRR10040010Aligned.sortedByCoord.out.bam_1",
"SRR10039529Aligned.sortedByCoord.out.bam_1",
"SRR10039510Aligned.sortedByCoord.out.bam_1",
"SRR10039535Aligned.sortedByCoord.out.bam_1",
"SRR10039943Aligned.sortedByCoord.out.bam_1",
"SRR10039573Aligned.sortedByCoord.out.bam_1",
"SRR10039599Aligned.sortedByCoord.out.bam_1",
"SRR10039604Aligned.sortedByCoord.out.bam_1",
"SRR10039609Aligned.sortedByCoord.out.bam_1",
"SRR10039626Aligned.sortedByCoord.out.bam_1",
"SRR10039661Aligned.sortedByCoord.out.bam_1",
"SRR10039742Aligned.sortedByCoord.out.bam_1",
"SRR10039743Aligned.sortedByCoord.out.bam_1",
"SRR10039786Aligned.sortedByCoord.out.bam_1",
"SRR10039822Aligned.sortedByCoord.out.bam_1",
"SRR10039823Aligned.sortedByCoord.out.bam_1",
"SRR10039835Aligned.sortedByCoord.out.bam_1",
"SRR10039892Aligned.sortedByCoord.out.bam_1",
"SRR10039899Aligned.sortedByCoord.out.bam_1",
"SRR10039902Aligned.sortedByCoord.out.bam_1",
"SRR10039903Aligned.sortedByCoord.out.bam_1",
"SRR10039904Aligned.sortedByCoord.out.bam_1",
"SRR10039906Aligned.sortedByCoord.out.bam_1",
"SRR10039907Aligned.sortedByCoord.out.bam_1")


Epi <- c("SRR10039624Aligned.sortedByCoord.out.bam_1",
"SRR10039629Aligned.sortedByCoord.out.bam_1",
"SRR10039750Aligned.sortedByCoord.out.bam_1",
"SRR10039531Aligned.sortedByCoord.out.bam_1",
"SRR10039532Aligned.sortedByCoord.out.bam_1",
"SRR10039534Aligned.sortedByCoord.out.bam_1",
"SRR10039536Aligned.sortedByCoord.out.bam_1",
"SRR10039537Aligned.sortedByCoord.out.bam_1",
"SRR10039538Aligned.sortedByCoord.out.bam_1",
"SRR10039550Aligned.sortedByCoord.out.bam_1",
"SRR10039565Aligned.sortedByCoord.out.bam_1",
"SRR10039568Aligned.sortedByCoord.out.bam_1",
"SRR10039598Aligned.sortedByCoord.out.bam_1",
"SRR10039611Aligned.sortedByCoord.out.bam_1",
"SRR10039612Aligned.sortedByCoord.out.bam_1",
"SRR10039613Aligned.sortedByCoord.out.bam_1",
"SRR10039615Aligned.sortedByCoord.out.bam_1",
"SRR10039617Aligned.sortedByCoord.out.bam_1",
"SRR10039620Aligned.sortedByCoord.out.bam_1",
"SRR10039622Aligned.sortedByCoord.out.bam_1",
"SRR10039623Aligned.sortedByCoord.out.bam_1")

EpiPrEn <- c("SRR10039498Aligned.sortedByCoord.out.bam_1",
"SRR10039523Aligned.sortedByCoord.out.bam_1",
"SRR10039528Aligned.sortedByCoord.out.bam_1",
"SRR10039542Aligned.sortedByCoord.out.bam_1",
"SRR10039543Aligned.sortedByCoord.out.bam_1",
"SRR10039546Aligned.sortedByCoord.out.bam_1",
"SRR10039547Aligned.sortedByCoord.out.bam_1",
"SRR10039548Aligned.sortedByCoord.out.bam_1",
"SRR10039549Aligned.sortedByCoord.out.bam_1",
"SRR10039552Aligned.sortedByCoord.out.bam_1",
"SRR10039553Aligned.sortedByCoord.out.bam_1",
"SRR10039554Aligned.sortedByCoord.out.bam_1",
"SRR10039557Aligned.sortedByCoord.out.bam_1")


STB <- c("SRR10039778Aligned.sortedByCoord.out.bam_1",
"SRR10039780Aligned.sortedByCoord.out.bam_1",
"SRR10039800Aligned.sortedByCoord.out.bam_1",
"SRR10039816Aligned.sortedByCoord.out.bam_1",
"SRR10039828Aligned.sortedByCoord.out.bam_1",
"SRR10039829Aligned.sortedByCoord.out.bam_1",
"SRR10039831Aligned.sortedByCoord.out.bam_1",
"SRR10039833Aligned.sortedByCoord.out.bam_1",
"SRR10039834Aligned.sortedByCoord.out.bam_1",
"SRR10039836Aligned.sortedByCoord.out.bam_1",
"SRR10039837Aligned.sortedByCoord.out.bam_1",
"SRR10039839Aligned.sortedByCoord.out.bam_1",
"SRR10039841Aligned.sortedByCoord.out.bam_1",
"SRR10039842Aligned.sortedByCoord.out.bam_1",
"SRR10039848Aligned.sortedByCoord.out.bam_1",
"SRR10039876Aligned.sortedByCoord.out.bam_1",
"SRR10039877Aligned.sortedByCoord.out.bam_1",
"SRR10039878Aligned.sortedByCoord.out.bam_1",
"SRR10039880Aligned.sortedByCoord.out.bam_1",
"SRR10039884Aligned.sortedByCoord.out.bam_1",
"SRR10039893Aligned.sortedByCoord.out.bam_1",
"SRR10039896Aligned.sortedByCoord.out.bam_1",
"SRR10039914Aligned.sortedByCoord.out.bam_1",
"SRR10039917Aligned.sortedByCoord.out.bam_1",
"SRR10039928Aligned.sortedByCoord.out.bam_1",
"SRR10039939Aligned.sortedByCoord.out.bam_1",
"SRR10039944Aligned.sortedByCoord.out.bam_1",
"SRR10039945Aligned.sortedByCoord.out.bam_1",
"SRR10039946Aligned.sortedByCoord.out.bam_1",
"SRR10039947Aligned.sortedByCoord.out.bam_1",
"SRR10039998Aligned.sortedByCoord.out.bam_1",
"SRR10040007Aligned.sortedByCoord.out.bam_1",
"SRR10040015Aligned.sortedByCoord.out.bam_1",
"SRR10040018Aligned.sortedByCoord.out.bam_1",
"SRR10040024Aligned.sortedByCoord.out.bam_1",
"SRR10040030Aligned.sortedByCoord.out.bam_1",
"SRR10040035Aligned.sortedByCoord.out.bam_1",
"SRR10040039Aligned.sortedByCoord.out.bam_1")

Tr <- c("SRR10039744Aligned.sortedByCoord.out.bam_1",
"SRR10039747Aligned.sortedByCoord.out.bam_1",
"SRR10039559Aligned.sortedByCoord.out.bam_1",
"SRR10039562Aligned.sortedByCoord.out.bam_1",
"SRR10039563Aligned.sortedByCoord.out.bam_1",
"SRR10039566Aligned.sortedByCoord.out.bam_1",
"SRR10039567Aligned.sortedByCoord.out.bam_1",
"SRR10039569Aligned.sortedByCoord.out.bam_1",
"SRR10039570Aligned.sortedByCoord.out.bam_1",
"SRR10039571Aligned.sortedByCoord.out.bam_1",
"SRR10039575Aligned.sortedByCoord.out.bam_1",
"SRR10039576Aligned.sortedByCoord.out.bam_1",
"SRR10039577Aligned.sortedByCoord.out.bam_1",
"SRR10039581Aligned.sortedByCoord.out.bam_1",
"SRR10039583Aligned.sortedByCoord.out.bam_1",
"SRR10039584Aligned.sortedByCoord.out.bam_1",
"SRR10039585Aligned.sortedByCoord.out.bam_1",
"SRR10039618Aligned.sortedByCoord.out.bam_1",
"SRR10039515Aligned.sortedByCoord.out.bam_1",
"SRR10039516Aligned.sortedByCoord.out.bam_1",
"SRR10039517Aligned.sortedByCoord.out.bam_1",
"SRR10039521Aligned.sortedByCoord.out.bam_1",
"SRR10039522Aligned.sortedByCoord.out.bam_1",
"SRR10039524Aligned.sortedByCoord.out.bam_1",
"SRR10039525Aligned.sortedByCoord.out.bam_1",
"SRR10039579Aligned.sortedByCoord.out.bam_1",
"SRR10039597Aligned.sortedByCoord.out.bam_1",
"SRR10039606Aligned.sortedByCoord.out.bam_1",
"SRR10039625Aligned.sortedByCoord.out.bam_1",
"SRR10039627Aligned.sortedByCoord.out.bam_1",
"SRR10039631Aligned.sortedByCoord.out.bam_1",
"SRR10039634Aligned.sortedByCoord.out.bam_1")


pICM <- c("SRR10039520Aligned.sortedByCoord.out.bam_1",
"SRR10039533Aligned.sortedByCoord.out.bam_1")
pEVT <- c("SRR10039921Aligned.sortedByCoord.out.bam_1")
pSTB <- c("SRR10039648Aligned.sortedByCoord.out.bam_1",
"SRR10039669Aligned.sortedByCoord.out.bam_1",
"SRR10039686Aligned.sortedByCoord.out.bam_1",
"SRR10039689Aligned.sortedByCoord.out.bam_1",
"SRR10039751Aligned.sortedByCoord.out.bam_1",
"SRR10039752Aligned.sortedByCoord.out.bam_1",
"SRR10039755Aligned.sortedByCoord.out.bam_1",
"SRR10039781Aligned.sortedByCoord.out.bam_1",
"SRR10039840Aligned.sortedByCoord.out.bam_1",
"SRR10039843Aligned.sortedByCoord.out.bam_1",
"SRR10039857Aligned.sortedByCoord.out.bam_1",
"SRR10039864Aligned.sortedByCoord.out.bam_1",
"SRR10039879Aligned.sortedByCoord.out.bam_1",
"SRR10039919Aligned.sortedByCoord.out.bam_1",
"SRR10039931Aligned.sortedByCoord.out.bam_1",
"SRR10039933Aligned.sortedByCoord.out.bam_1",
"SRR10039942Aligned.sortedByCoord.out.bam_1",
"SRR10039992Aligned.sortedByCoord.out.bam_1",
"SRR10040019Aligned.sortedByCoord.out.bam_1")

pVE <- c("SRR10039610Aligned.sortedByCoord.out.bam_1",
"SRR10039707Aligned.sortedByCoord.out.bam_1",
"SRR10039905Aligned.sortedByCoord.out.bam_1")
pEpi <- c("SRR10039511Aligned.sortedByCoord.out.bam_1",
"SRR10039541Aligned.sortedByCoord.out.bam_1",
"SRR10039551Aligned.sortedByCoord.out.bam_1",
"SRR10039555Aligned.sortedByCoord.out.bam_1",
"SRR10039556Aligned.sortedByCoord.out.bam_1",
"SRR10039558Aligned.sortedByCoord.out.bam_1",
"SRR10039590Aligned.sortedByCoord.out.bam_1",
"SRR10039614Aligned.sortedByCoord.out.bam_1",
"SRR10039616Aligned.sortedByCoord.out.bam_1",
"SRR10039619Aligned.sortedByCoord.out.bam_1",
"SRR10039621Aligned.sortedByCoord.out.bam_1")
pPS <- c("SRR10039935Aligned.sortedByCoord.out.bam_1",
"SRR10039941Aligned.sortedByCoord.out.bam_1",
"SRR10039955Aligned.sortedByCoord.out.bam_1",
"SRR10039956Aligned.sortedByCoord.out.bam_1",
"SRR10039964Aligned.sortedByCoord.out.bam_1",
"SRR10039965Aligned.sortedByCoord.out.bam_1",
"SRR10039967Aligned.sortedByCoord.out.bam_1",
"SRR10039977Aligned.sortedByCoord.out.bam_1",
"SRR10039979Aligned.sortedByCoord.out.bam_1",
"SRR10039988Aligned.sortedByCoord.out.bam_1")
pEmD <- c("SRR10039643Aligned.sortedByCoord.out.bam_1",
"SRR10039734Aligned.sortedByCoord.out.bam_1",
"SRR10039775Aligned.sortedByCoord.out.bam_1",
"SRR10039783Aligned.sortedByCoord.out.bam_1",
"SRR10039885Aligned.sortedByCoord.out.bam_1",
"SRR10039889Aligned.sortedByCoord.out.bam_1",
"SRR10039926Aligned.sortedByCoord.out.bam_1",
"SRR10039950Aligned.sortedByCoord.out.bam_1",
"SRR10039970Aligned.sortedByCoord.out.bam_1",
"SRR10040002Aligned.sortedByCoord.out.bam_1")


#Idents(mammal.combined1,cells=AmList) <- "Am"
#Idents(mammal.combined1,cells=CTBList) <- "CTB"
#Idents(mammal.combined1,cells=eEVTList) <- "eEVT"
#Idents(mammal.combined1,cells=eSTBList) <- "eSTB"
#Idents(mammal.combined1,cells=EVTList) <- "EVT"
#Idents(mammal.combined1,cells=ICMList) <- "ICM"
#Idents(mammal.combined1,cells=PostEpiAME) <- "LAm"
#Idents(mammal.combined1,cells=PostEPIE1) <- "PostEPIE1"
#Idents(mammal.combined1,cells=PostEPIE2) <- "PostEPIE2"
#Idents(mammal.combined1,cells=PostEpiGast) <- "PostEpiGast"
#Idents(mammal.combined1,cells=PrEn) <- "PrEn"
#Idents(mammal.combined1,cells=Epi) <- "Epi"
#Idents(mammal.combined1,cells=EpiPrEn) <- "EpiPrEn"
#Idents(mammal.combined1,cells=STB) <- "STB "
#Idents(mammal.combined1,cells=Tr) <- "Tr"
#Idents(mammal.combined1,cells=pICM) <- "Unknown"
#Idents(mammal.combined1,cells=pEVT) <- "Unknown"
#Idents(mammal.combined1,cells=pSTB) <- "pSTB"
#Idents(mammal.combined1,cells=pVE) <- "pVE"
#Idents(mammal.combined1,cells=pEpi) <- "pEpi"
#Idents(mammal.combined1,cells=pPS) <- "Unknown"
#Idents(mammal.combined1,cells=pEmD) <- "Unknown"

Idents(mammal.combined1) <- mammal.combined1$ano

Idents(mammal.combined1,cells=AmList) <- "Am"
Idents(mammal.combined1,cells=CTBList) <- "CTB"
Idents(mammal.combined1,cells=eEVTList) <- "CTB"
Idents(mammal.combined1,cells=eSTBList) <- "STB"
Idents(mammal.combined1,cells=EVTList) <- "EVT"
Idents(mammal.combined1,cells=ICMList) <- "ICM"
Idents(mammal.combined1,cells=PostEpiAME) <- "Am"
Idents(mammal.combined1,cells=PostEPIE1) <- "EmDisc"
Idents(mammal.combined1,cells=PostEPIE2) <- "EmDisc"
Idents(mammal.combined1,cells=PostEpiGast) <- "EmDisc"
Idents(mammal.combined1,cells=PrEn) <- "VE"
Idents(mammal.combined1,cells=Epi) <- "Epi"
Idents(mammal.combined1,cells=EpiPrEn) <- "Epi"
Idents(mammal.combined1,cells=STB) <- "STB "
Idents(mammal.combined1,cells=Tr) <- "Tr"
Idents(mammal.combined1,cells=pICM) <- "Unknown"
Idents(mammal.combined1,cells=pEVT) <- "Unknown"
Idents(mammal.combined1,cells=pSTB) <- "STB"
Idents(mammal.combined1,cells=pVE) <- "VE"
Idents(mammal.combined1,cells=pEpi) <- "Epi"
Idents(mammal.combined1,cells=pPS) <- "Unknown"
Idents(mammal.combined1,cells=pEmD) <- "Unknown"

mammal.combined1 <- subset(mammal.combined1,idents=c("Unknown"),invert=TRUE)

uID <- as.character(Idents(mammal.combined1))
uID[which(Idents(mammal.combined1)=="STB_d9")] <- "CTB_d9"
uID[which(Idents(mammal.combined1)=="STB_d11")] <- "CTB_d11"
uID[which(Idents(mammal.combined1)=="STB_d12")] <- "CTB_d12"
uID[which(Idents(mammal.combined1)=="CTB_d9")] <- "STB_d9"
uID[which(Idents(mammal.combined1)=="CTB_d11")] <- "STB_d11"
uID[which(Idents(mammal.combined1)=="CTB_d12")] <- "STB_d12"

Idents(mammal.combined1) <- uID
Idents(mammal.combined1,cells=WhichCells(mammal.combined1,idents="STB ")) <- "STB"
Idents(mammal.combined1,cells=WhichCells(mammal.combined1,idents="Hyp_d14")) <- "VE_d14"

p1 <- DimPlot(mammal.combined1,  pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_1_2_Embryonicoly_ano","boundaryrp.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)

#Idents(mammal.combined1) <- mammal.combined1$Anotations
colind <- integer( length( levels(Idents(mammal.combined1)) )  )
for (i in 1:length( levels(Idents(mammal.combined1)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined1))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined1,  cols = coluse, pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_1_2_Embryonicoly_ano","boundaryrp.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)




DefaultAssay(mammal.combined1) <- "RNA"
FeaturePlot(mammal.combined1,feature="SOX15",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_SOX15_justEmbryonic",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined1,feature="SOX15",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_SOX15_justEmbryonic",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined1,feature="PDGFRA",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_PDGFRA_justEmbryonic",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined1,feature="BAMBI",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_BAMBI_justEmbryonic",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined1,feature="DIO2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_DIO2_justEmbryonic",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined1,feature="HLA-G",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_HLA-G_justEmbryonic",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined1,feature="JAM3",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_JAM3_justEmbryonic",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined1,feature="APOA1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_APOA1_justEmbryonic",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined1,feature="GATA4",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_GATA4_justEmbryonic",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined1,feature="CER1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_CER1_justEmbryonic",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined1,feature="CGA",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_CGA_justEmbryonic",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined1,feature="GATA3",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_GATA3_justEmbryonic",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)


#asdsakjbdsakjdsajkbd

#sadad
#p1 <- DimPlot(mammal.combined2,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
#ggsave(filename=paste(saveext,"/DimRed/CompletePCA_1_2_EmbryonicMat_ano","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)

#Idents(mammal.combined2) <- mammal.combined2$Genotype

#p1 <- DimPlot(mammal.combined2,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
#ggsave(filename=paste(saveext,"/DimRed/CompletePCA_1_2_EmbryonicMat_GT","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)

#Idents(mammal.combined2) <- mammal.combined2$Cl15

#p1 <- DimPlot(mammal.combined2,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
#ggsave(filename=paste(saveext,"/DimRed/CompletePCA_1_2_EmbryonicMat_Cl15","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)



ACTB <- c("ACATTTCAGCTTTCTT-1_8",
"ATGACCACACTACTTT-1_8",
"ACCAAACAGAATCGTA-1_9",
"ACTGCAATCGCTCTCA-1_9",
"CACTGGGCAGATACTC-1_9",
"CAGATCAGTGGAATGC-1_9",
"CATAGACAGGATATAC-1_9",
"CTAGACAAGAAAGTCT-1_9",
"CTTCGGTAGTCTTGGT-1_9",
"GGAATCTTCGTTCTAT-1_9",
"GTAGATCTCGGCTATA-1_9",
"GTTAGACTCGAGTTGT-1_9",
"TATCTTGCACATAACC-1_9",
"TCAATCTGTGAAAGTT-1_9",
"TCCTTCTAGGCCACCT-1_9",
"TCTATACGTTAACCTG-1_9",
"TGAGACTCACTCCGGA-1_9",
"TGCGGGTAGACTACGG-1_9",
"TTCGCTGTCTCTCAAT-1_9",
"TTTCCTCTCTCCGATC-1_9",
"TTTGATCGTAGCTTGT-1_9",
"AAGACTCCACTTCAAG-1_11",
"ACATTTCCAATATCCG-1_11",
"AGCTACACACTCTGCT-1_11",
"ATTCCTAAGGAAAGGT-1_11",
"CCTTCAGGTTAAACCC-1_11",
"CTGTCGTGTGAGACGT-1_11",
"GTGAGGAGTGGTAATA-1_11",
"TGTGATGCATTCCTCG-1_11",
"TTCATGTAGCACTTTG-1_11",
"GTTGCTCAGTAGGGTC-1_12",
"AATGACCGTGGCTAGA-1_14",
"ATCCGTCTCCCTCGTA-1_14",
"CATTCCGTCCATCGTC-1_14")

AEPI <- c("AACAACCGTGCCCAGT-1_8",
"AACAAGAGTCATAGTC-1_8",
"AATTCCTTCTGAATGC-1_8",
"ACACTGAAGGCCTTCG-1_8",
"ACTTTCAAGTACCCTA-1_8",
"AGACCCGAGAAAGCGA-1_8",
"AGCTCAACATCGGATT-1_8",
"AGTTAGCTCCAGTTCC-1_8",
"ATCCTATGTAGGCAAC-1_8",
"ATGCGATCATGACAGG-1_8",
"ATGGGAGCACTCCTGT-1_8",
"CAACCTCGTCACTTCC-1_8",
"CACAGATTCTAGTACG-1_8",
"CGGGACTCAGCTTCCT-1_8",
"CTAACCCGTGTCTTCC-1_8",
"CTGCTCACACTCACTC-1_8",
"CTGTACCGTAGGGAGG-1_8",
"GACCCAGAGTTTGGCT-1_8",
"GACCCTTAGCTCACTA-1_8",
"GCAGCCAAGGTGCAGT-1_8",
"GCAGTTATCTGGCCTT-1_8",
"GGGACTCCACCGTGAC-1_8",
"GGGTTATGTGACCGTC-1_8",
"GTATTGGAGGTTGGAC-1_8",
"TACCCGTCAAACCATC-1_8",
"TCGCACTAGCCAAGGT-1_8",
"TCTTGCGCAGCGTACC-1_8",
"TGTTACTGTTTGTTCT-1_8",
"TTACCATTCTCCAATT-1_8",
"TTCATGTGTTGCAACT-1_8",
"TTCTTCCCAGACCAGA-1_8",
"TTTGATCAGTAATCCC-1_8",
"AAACCCAAGTAGCAAT-1_9",
"AACAGGGGTCCAGCCA-1_9",
"AACCAACCACGTCATA-1_9",
"AACTTCTAGGACATCG-1_9",
"AAGGAATTCACTTGTT-1_9",
"AAGGTAACAGTTAGGG-1_9",
"AATCGACAGACGAAGA-1_9",
"AATCGTGAGCTCGACC-1_9",
"AATGACCGTTCGGACC-1_9",
"ACACTGAGTCTGGTTA-1_9",
"ACAGAAAAGACTCGAG-1_9",
"ACCTGTCAGTGAGTGC-1_9",
"ACGGTTAAGGACTATA-1_9",
"ACGTAGTGTGACGTCC-1_9",
"ACTATTCAGCACCGTC-1_9",
"ACTGATGAGTCATTGC-1_9",
"AGACCCGCAGAGGAAA-1_9",
"AGACCCGCATCGCTAA-1_9",
"ATACCGAAGGCCCAAA-1_9",
"ATCCCTGGTACTCGTA-1_9",
"ATCGGCGGTCCTGAAT-1_9",
"ATGCATGTCTCTTAAC-1_9",
"ATTATCCGTCCATACA-1_9",
"ATTCACTGTAATTGGA-1_9",
"ATTCGTTCAGGAACCA-1_9",
"ATTCGTTTCCTTATGT-1_9",
"ATTTCACTCTGCACCT-1_9",
"CAAGACTAGGACGCTA-1_9",
"CAAGAGGAGTTGGAGC-1_9",
"CACAACACAGGACATG-1_9",
"CAGATACCATGACAGG-1_9",
"CAGCAATAGTTCTACG-1_9",
"CAGGTATGTTCACGAT-1_9",
"CATAGACGTCTACGTA-1_9",
"CATCCGTTCTTGGCTC-1_9",
"CATGCCTGTCACAGAG-1_9",
"CATTGCCTCCGATAAC-1_9",
"CCACCATGTTCAGGTT-1_9",
"CCCATTGCATCTGCGG-1_9",
"CCCATTGGTTGCATAC-1_9",
"CCCGAAGGTCACAGTT-1_9",
"CCCTAACCACCGTGAC-1_9",
"CCGGGTACAATCGCAT-1_9",
"CCTAAGACAGGTCCGT-1_9",
"CCTCAACTCGTCACCT-1_9",
"CCTCAGTGTTCAGCTA-1_9",
"CCTCTAGCACCGCTAG-1_9",
"CCTTGTGAGCCTGGAA-1_9",
"CGGGCATCAACCAGAG-1_9",
"CGGGTCATCAGTCACA-1_9",
"CGGTCAGGTTGTAAAG-1_9",
"CGTCCATGTACAGTCT-1_9",
"CTAACCCAGTAGAGTT-1_9",
"CTAACTTAGCTGTTCA-1_9",
"CTATAGGGTGCCTACG-1_9",
"CTATCCGTCCTGTTGC-1_9",
"CTATCTAAGCGGACAT-1_9",
"CTCAGAAAGTGATAGT-1_9",
"CTCCAACTCACCCTCA-1_9",
"CTGAATGTCGCTTGAA-1_9",
"CTGATCCGTTGCGTAT-1_9",
"CTGCTCATCAGTGATC-1_9",
"CTGTGGGTCAGCATTG-1_9",
"CTTCAATAGGATACGC-1_9",
"CTTCGGTTCCGAAATC-1_9",
"CTTCTAACAAGATCCT-1_9",
"CTTTCAACAACAGTGG-1_9",
"GAACTGTGTATGGAGC-1_9",
"GAAGAATAGGAGAGGC-1_9",
"GACCAATGTGTCATGT-1_9",
"GACCTTCCACTTGGCG-1_9",
"GAGTCATCACGGGCTT-1_9",
"GAGTTACTCTCCGCAT-1_9",
"GATCAGTTCTGGGATT-1_9",
"GATCCCTAGACATACA-1_9",
"GATCCCTCACTACGGC-1_9",
"GATGAGGAGTTCTACG-1_9",
"GCCAGCATCACTACGA-1_9",
"GCGGAAAAGTTGCGCC-1_9",
"GCGTTTCAGTATAACG-1_9",
"GCTGAATCAATGTCTG-1_9",
"GCTGGGTTCATAGCAC-1_9",
"GGATCTACAAGTATAG-1_9",
"GGATGTTAGATCGCTT-1_9",
"GGGACAACAATCGCGC-1_9",
"GGGATCCTCTGGACTA-1_9",
"GGGTAGAAGTAGGGTC-1_9",
"GGGTGAAAGTATGAGT-1_9",
"GGTCACGTCTGAGTCA-1_9",
"GGTGAAGAGGACTTCT-1_9",
"GGTGTCGGTACCCGAC-1_9",
"GGTTAACTCCGGCAGT-1_9",
"GTCGAATCAGCTACAT-1_9",
"GTCTACCAGAGCCTGA-1_9",
"GTCTCACGTATACAGA-1_9",
"GTCTGTCAGCGACATG-1_9",
"GTGCACGCACAGCATT-1_9",
"GTTGCGGTCATGCGGC-1_9",
"TAACACGTCACGTAGT-1_9",
"TAACCAGTCCCAGTGG-1_9",
"TACGGTATCCCACAGG-1_9",
"TATACCTTCCGTGGGT-1_9",
"TATATCCGTGGGTCAA-1_9",
"TCAAGACTCTCCAATT-1_9",
"TCAAGTGGTCACGTGC-1_9",
"TCACATTAGGAAGTAG-1_9",
"TCAGCCTCATCCGGTG-1_9",
"TCATCATGTACCCGCA-1_9",
"TCATCCGGTGCTATTG-1_9",
"TCATGGAAGTCTGGAG-1_9",
"TCATGTTTCCTTCGAC-1_9",
"TCATGTTTCGAGTCCG-1_9",
"TCATTACCAAGTATAG-1_9",
"TCCTCCCCAAACGAGC-1_9",
"TCCTTTCCAACGGCCT-1_9",
"TCGGTCTTCCACCTCA-1_9",
"TCTACCGGTATTCCGA-1_9",
"TCTCCGACATTGCAAC-1_9",
"TCTTAGTGTGGACTAG-1_9",
"TCTTCCTAGTCTAACC-1_9",
"TGAGACTAGACCACGA-1_9",
"TGATGGTAGTATAACG-1_9",
"TGCATCCGTAGTGGCA-1_9",
"TGCGACGCACATTGTG-1_9",
"TGGATGTAGGGTCAAC-1_9",
"TGTTACTAGCCTCTCT-1_9",
"TTACCATGTGGTCTCG-1_9",
"TTAGGGTGTTTCCATT-1_9",
"TTCACCGCACACCTAA-1_9",
"TTCATGTGTTGAATCC-1_9",
"TTCCAATGTGCGGTAA-1_9",
"TTCTCTCAGACTACGG-1_9",
"TTGCATTAGGGCATGT-1_9",
"AATCGACCAGTAACCT-1_11",
"ACGTAGTAGTCTAGCT-1_11",
"AGCGCTGTCGCGGACT-1_11",
"ATCGTAGCAGCTTTCC-1_11",
"CAATTTCTCGAGCCTG-1_11",
"CCTCTAGTCATTGCTT-1_11",
"CCTGCATAGTTACGAA-1_11",
"CTAAGTGTCCGTAATG-1_11",
"CTCATTACAAACACGG-1_11",
"CTGCCTACAACAAGTA-1_11",
"GATGGAGCATCTTAGG-1_11",
"GTAATGCTCGGTGTAT-1_11",
"TATTTCGCACCTGAAT-1_11",
"TCGGGCAGTAGCTCGC-1_11",
"TGCATCCCAAACTAAG-1_11",
"AACAAGATCTCCAAGA-1_12",
"ACCTACCAGTTAGAAC-1_12",
"AGACAGGCAGTGGTGA-1_12",
"AGATAGAAGTCGAAGC-1_12",
"AGGGAGTTCAACCCGG-1_12",
"ATCCGTCGTACAAGTA-1_12",
"CAACGGCTCAAGAAAC-1_12",
"CACACAACAGTCTGGC-1_12",
"CTCCCAAAGGACATCG-1_12",
"CTCCCAAGTGGATCGA-1_12",
"CTGAGCGCAGCTTTCC-1_12",
"CTGTATTCAAGGCGTA-1_12",
"GACGCTGGTTCATCGA-1_12",
"GATTGGTTCGCTAGCG-1_12",
"GCGGAAATCATCCTGC-1_12",
"GGAGGATAGTCGTCTA-1_12",
"GGTGTCGCAGCCGTTG-1_12",
"GTCAAACGTTAGGAGC-1_12",
"GTCTTTAGTCTTGCTC-1_12",
"GTGCTGGCACCGGTCA-1_12",
"GTGGTTATCGCCAATA-1_12",
"TACCCGTTCTACCCAC-1_12",
"TAGTGCACAACTTGCA-1_12",
"TCGCTCATCCACGGAC-1_12",
"TCGTAGACACTGATTG-1_12",
"TCTACATTCGAGAAAT-1_12",
"TCTGCCACATGCCGAC-1_12",
"TCTTTGAAGCAATAGT-1_12",
"TGGATGTGTTAACAGA-1_12",
"TTCATTGCACGAAAGC-1_12",
"TTGGATGAGTTTGTCG-1_12",
"TTTCATGTCCGCAACG-1_12",
"AAAGGGCCACTGGAAG-1_14",
"AACGAAAGTGAACCGA-1_14",
"ACATCCCCAGCTATAC-1_14",
"CATCAAGTCATACGAC-1_14",
"GATGCTACAGCAAGAC-1_14",
"GCCCAGACAACCGCCA-1_14",
"GCTACCTAGGCCACCT-1_14",
"GCTGAATCAGCCATTA-1_14",
"GTTCATTCAGTGACCC-1_14",
"TCAGCCTGTCATACCA-1_14",
"TTTGTTGAGCGTGCCT-1_14")

AEVT <- c("AGCGTATGTGGTACAG-1_8",
"GTGTCCTTCGAACACT-1_8",
"TGTCCACGTGGACCTC-1_8",
"TTCCGGTAGGGAGTTC-1_8",
"TTCCTCTAGTATGCAA-1_8",
"ACAGAAAGTTGGACTT-1_9",
"ACCGTTCCATGATAGA-1_9",
"AGACCATGTACTAACC-1_9",
"ATCGTGACAGCAGGAT-1_9",
"ATCTTCACACTAACCA-1_9",
"ATTCTACAGACCATTC-1_9",
"ATTTACCTCACCACAA-1_9",
"CACGTTCTCAAAGACA-1_9",
"GACCCTTTCTAACACG-1_9",
"GACTTCCTCGCTCATC-1_9",
"GAGTCTATCAGCTAGT-1_9",
"GATGTTGTCTCGACCT-1_9",
"GTAAGTCAGAGTAACT-1_9",
"GTCCACTTCTTCACGC-1_9",
"TACTTACTCAAACCCA-1_9",
"TCAGGTAGTTTCTATC-1_9",
"TCATTGTCAGTTTCAG-1_9",
"TGCGGCACACCTGCAG-1_9",
"TGCTCCAGTCATCCCT-1_9",
"ACACAGTTCCGCTTAC-1_11",
"ACTTTCATCACACGAT-1_11",
"AGAAGCGGTCCTACAA-1_11",
"AGACAAACATGTCGTA-1_11",
"ATTCCTACAACAAGAT-1_11",
"GACGTTATCAGGTAAA-1_11",
"GCAGCTGCAGTCTCTC-1_11",
"GCGGAAACATGACGAG-1_11",
"GGTAATCTCCGTGTCT-1_11",
"GTTCGCTCAGGCATTT-1_11",
"TATCGCCAGGGCGAGA-1_11",
"TTGGGCGCAGTCGAGA-1_11",
"ATCGTGAGTTATCTTC-1_12",
"CAGCCAGTCAGTGTTG-1_12",
"CAGCGTGAGTGGCGAT-1_12",
"CGTTCTGGTTCTTGCC-1_12",
"GACTGATCAGCAGGAT-1_12",
"GCTCAAATCTCACGAA-1_12",
"GGCGTCAAGCTCGAAG-1_12",
"TTGGGATAGACCTGGA-1_12",
"AACACACCAAGCGCAA-1_14",
"ACGATCATCCACCTCA-1_14",
"ACGTAACTCCATTTAC-1_14",
"AGACAGGCAAGAGGTC-1_14",
"ATTTCACAGAATGTTG-1_14",
"CAACGGCTCTTGGCTC-1_14",
"TTTACTGGTCGCCTAG-1_14")

AHYP <- c("ACTTTGTGTGGCTGCT-1_8",
"AGAGCAGGTAACGGTG-1_8",
"ATTCACTCATACTTTC-1_8",
"CACTAAGCAATCTAGC-1_8",
"CCAATGAAGTTCATCG-1_8",
"CCCTCAACACCCTTAC-1_8",
"CCTCACATCCCGAACG-1_8",
"CCTGCATCATTAAGCC-1_8",
"CGAGGCTAGATGCGAC-1_8",
"CTCACTGGTAATGATG-1_8",
"GACCCAGGTCTCCCTA-1_8",
"GACGTTAAGCAGTAAT-1_8",
"GATCACATCCTCTCTT-1_8",
"GTTGCTCTCACCTCGT-1_8",
"TATCAGGCACAAGTTC-1_8",
"TCACTATCATACTGAC-1_8",
"TCATTACTCCGATTAG-1_8",
"TCTTCCTTCAACCCGG-1_8",
"TGTGAGTCAAGAGGCT-1_8",
"TTCCTAACACTATCGA-1_8",
"AGGTAGGTCCCTCTTT-1_11",
"CATACTTGTCCGGTGT-1_12")

ASTB <- c("AAGCCATGTCCCTAAA-1_8",
"ACTCTCGCAGTCGAGA-1_8",
"ACTCTCGGTTCTCTCG-1_8",
"ACTTTGTTCCGAGCTG-1_8",
"AGGCCACCATGTGGTT-1_8",
"ATTCATCGTCTTACAG-1_8",
"CCACACTAGTTGCCCG-1_8",
"CCCTAACGTAGACAGC-1_8",
"CGTAAGTAGCTATCCA-1_8",
"CTTTCAATCAGCCCAG-1_8",
"GTTACAGTCTTTCCGG-1_8",
"TAAGCCAGTTGCGGAA-1_8",
"TACGTCCGTCTGTGCG-1_8",
"TAGCACACATCGATAC-1_8",
"TCATTACGTCTTTATC-1_8",
"AAAGGGCCACCAATTG-1_9",
"AAAGTCCCACTGAGGA-1_9",
"AACAAAGCACAGAGCA-1_9",
"AACCATGGTTATGACC-1_9",
"AACCTTTTCAACTCTT-1_9",
"ACAGAAACAGCAGAAC-1_9",
"ACCAACATCATCACTT-1_9",
"ACCCTTGCAAGTGGCA-1_9",
"ACCTACCTCCCATTTA-1_9",
"ACCTGTCTCGCTTGCT-1_9",
"ACTATTCGTCTGTGAT-1_9",
"ACTGATGCACTTCAGA-1_9",
"ACTTATCAGCGTGCCT-1_9",
"ACTTTCACAGCTAACT-1_9",
"AGAACCTTCGGACTTA-1_9",
"AGACACTGTCCGGTGT-1_9",
"AGACCATTCTAAACGC-1_9",
"AGCCACGCAAGCCTGC-1_9",
"AGCTTCCTCGCAGTTA-1_9",
"AGGTAGGGTTAGCTAC-1_9",
"AGGTCATAGGCCTTGC-1_9",
"AGTACCACAATCCAGT-1_9",
"AGTCAACTCCATCTGC-1_9",
"ATACTTCTCTTGCAGA-1_9",
"ATAGGCTCACACGGAA-1_9",
"ATAGGCTGTCCAATCA-1_9",
"ATCCACCAGGCATCAG-1_9",
"ATCCCTGGTGGCAGAT-1_9",
"ATCGCCTCACTAAACC-1_9",
"ATCTCTAAGTAGGATT-1_9",
"ATGCGATTCGCATTAG-1_9",
"ATTATCCAGGCTCACC-1_9",
"ATTATCCAGTCTCTGA-1_9",
"ATTCATCCAAGTGACG-1_9",
"ATTCCATCACACAGCC-1_9",
"ATTTACCCACCATTCC-1_9",
"CAACAACTCCTACGAA-1_9",
"CAACCAATCAAGAGGC-1_9",
"CAACCTCAGTCACAGG-1_9",
"CAACGATTCAACTGGT-1_9",
"CAACGGCTCAGTCCGG-1_9",
"CAACGGCTCGACACCG-1_9",
"CAATGACCACTACCCT-1_9",
"CACAACAGTCAGTCGC-1_9",
"CACGGGTTCATTTACC-1_9",
"CATGCGGTCACAACCA-1_9",
"CATGGTACACACACTA-1_9",
"CATTCTACAAAGCTCT-1_9",
"CATTGTTTCCGATTAG-1_9",
"CATTTCACAGCAATTC-1_9",
"CCCTCAACAAGTCCAT-1_9",
"CCGGACATCGAGAGCA-1_9",
"CCGTGAGCAAGTCCCG-1_9",
"CCTCAACTCTTGCAAG-1_9",
"CCTCACATCTAACACG-1_9",
"CCTTTGGTCTATGCCC-1_9",
"CGAGTGCCAATACGAA-1_9",
"CGCCAGAAGGCGTTAG-1_9",
"CGGAACCCATGGGATG-1_9",
"CTAGACATCTGCTCTG-1_9",
"CTCAATTGTGTCCAAT-1_9",
"CTCCACAAGTCCCGGT-1_9",
"CTCCCTCCAGCACAAG-1_9",
"CTCGAGGTCTATACTC-1_9",
"CTCTCGACAATGCAGG-1_9",
"CTGCCATTCTGTCTCG-1_9",
"CTTCTCTAGCCATTCA-1_9",
"GAATAGAAGAACCGCA-1_9",
"GACCTTCGTCACCGAC-1_9",
"GACTCAACACTGGACC-1_9",
"GAGGGTAAGAAACCAT-1_9",
"GATAGCTCACTCCGGA-1_9",
"GATGGAGCATTCGGGC-1_9",
"GATTGGTTCGCCCAGA-1_9",
"GCACATAGTCATCCCT-1_9",
"GCAGGCTTCTGCATGA-1_9",
"GCATCGGGTGGTCTAT-1_9",
"GCCGATGTCTGAGGCC-1_9",
"GCGGATCTCCTAGCCT-1_9",
"GCTCAAACATTGCCGG-1_9",
"GCTTCACTCTTCGATT-1_9",
"GCTTTCGAGTTTCGGT-1_9",
"GGAATCTAGCCTCTCT-1_9",
"GGATGTTCATGCACTA-1_9",
"GGCTTTCCAACATACC-1_9",
"GGGCCATTCTTCCACG-1_9",
"GGGTCACCACCGTACG-1_9",
"GGTGATTGTTATGTGC-1_9",
"GGTGTTAGTCCCAAAT-1_9",
"GTAACCATCCGACATA-1_9",
"GTAGAAAAGGTACTGG-1_9",
"GTAGTACTCAAGAGGC-1_9",
"GTATTTCTCAGCAGAG-1_9",
"GTCTACCCAACACGTT-1_9",
"GTCTAGAAGCAACAAT-1_9",
"GTGAGGATCTTGGCTC-1_9",
"TAACGACCAAAGCACG-1_9",
"TAATCTCGTATCTCTT-1_9",
"TAATCTCTCGCAATGT-1_9",
"TAGGGTTTCCTCATAT-1_9",
"TAGGTTGTCATTACTC-1_9",
"TATATCCTCCATTTGT-1_9",
"TATCTGTTCGGTCGAC-1_9",
"TATTGCTGTCTGCAAT-1_9",
"TCACGCTAGCTAAGTA-1_9",
"TCACTATTCACAATGC-1_9",
"TCACTCGTCCACACCT-1_9",
"TCAGGGCCATACCATG-1_9",
"TCAGTGACACTCTCGT-1_9",
"TCATATCAGTAGTCTC-1_9",
"TCATCCGGTATTGACC-1_9",
"TCATTCACACTTACAG-1_9",
"TCCACCATCATAGACC-1_9",
"TCCACCATCGATAACC-1_9",
"TCCCACAGTAATTAGG-1_9",
"TCCTAATTCATTCATC-1_9",
"TCGGTCTCAGCAGGAT-1_9",
"TCTATACAGATACATG-1_9",
"TCTATACCAAGCTACT-1_9",
"TCTATCATCTCCAAGA-1_9",
"TCTCAGCTCGACTCCT-1_9",
"TCTGGCTTCATAGGCT-1_9",
"TCTTAGTGTCCTCCTA-1_9",
"TGAGCATTCGCGCTGA-1_9",
"TGAGGTTCAGTGTGGA-1_9",
"TGATGCATCACACGAT-1_9",
"TGCAGATTCAATCTCT-1_9",
"TGCAGATTCGATGGAG-1_9",
"TGCATCCGTGGGTTGA-1_9",
"TGCGGGTTCCTCTTTC-1_9",
"TGTGATGCAGCCTACG-1_9",
"TGTGGCGCAACCTAAC-1_9",
"TGTTCCGGTGATATAG-1_9",
"TTACAGGAGACCGCCT-1_9",
"TTACAGGGTTGCTCCT-1_9",
"TTACGCCAGTATAACG-1_9",
"TTACGTTCACGACAGA-1_9",
"TTCCTTCTCCAGTGCG-1_9",
"TTCGGTCAGGGAACAA-1_9",
"TTGCCTGCAAATGATG-1_9",
"TTGGGATCATGGATCT-1_9",
"TTGGGCGAGTGGTGAC-1_9",
"TTGTTTGAGATTACCC-1_9",
"TTTCATGGTCACGACC-1_9",
"TTTCCTCAGCATTGAA-1_9",
"TTTCCTCGTTTGAACC-1_9",
"TTTGACTTCCCTCATG-1_9",
"TTTGATCTCCTACCGT-1_9",
"TTTGGAGGTGCTCTTC-1_9",
"AAGTGAAAGTACAGAT-1_11",
"AGTAACCAGGTAGACC-1_11",
"AGTTCCCTCATTTCCA-1_11",
"ATCACAGAGGGTTAAT-1_11",
"ATCCTATTCAGCATTG-1_11",
"ATCTCTAGTCTACAGT-1_11",
"ATGGGTTCAGACCGCT-1_11",
"CATGCAAAGAAGCTGC-1_11",
"CGAGTTAAGCAGGGAG-1_11",
"CGTGATACAGATTCGT-1_11",
"CGTTGGGAGGGAGGGT-1_11",
"CTCAGTCAGTAGCAAT-1_11",
"CTTCTCTCACTCAGAT-1_11",
"GACCTTCTCATATGGC-1_11",
"GACGTTATCCGTGGGT-1_11",
"GCAGCCAAGCGGTAAC-1_11",
"GCTGAATAGAAGTCAT-1_11",
"GGTAGAGGTATGGAGC-1_11",
"GGTTAACAGGTATAGT-1_11",
"GTAACACCAAGCTACT-1_11",
"GTAAGTCAGCCTCCAG-1_11",
"TACTGCCGTGCTATTG-1_11",
"TCATGAGTCGAACTCA-1_11",
"TCATTGTCATTGCCGG-1_11",
"TCGCACTTCATCACCC-1_11",
"TCTTAGTTCGGAACTT-1_11",
"TGGTTAGTCAAAGACA-1_11",
"TTCATTGAGGTTCATC-1_11",
"TTCCGTGAGTGTACAA-1_11",
"TTCTTGACAGTATACC-1_11",
"ACCTACCCAATCGCCG-1_12",
"ACGCACGGTAGCGCCT-1_12",
"ACTTTCAGTTAACAGA-1_12",
"ACTTTGTAGCCGTTGC-1_12",
"AGCATCATCCACGTAA-1_12",
"AGGACGAGTCAGACTT-1_12",
"ATACCGACATCTAACG-1_12",
"CAACAACGTAATGCGG-1_12",
"CAACAGTCATAAGCGG-1_12",
"CAACGATAGAGGATGA-1_12",
"CAGGCCATCAGGAACG-1_12",
"CATACTTAGACCAAGC-1_12",
"CATTCTATCATTTCGT-1_12",
"CCAAGCGTCATTCACT-1_12",
"CCGCAAGGTGCACAAG-1_12",
"CCTCATGTCTCCTGTG-1_12",
"CGAGTTAGTCGGCCTA-1_12",
"CTACCCACACAAAGTA-1_12",
"CTACTATGTAAGCGGT-1_12",
"CTCAGGGCAGGCACTC-1_12",
"CTCCCAATCGTCAAAC-1_12",
"GAACACTGTTCCGCAG-1_12",
"GGCTTGGTCAGACTGT-1_12",
"GGGTGAATCTCGGGAC-1_12",
"GTGGAGAAGGCTCTCG-1_12",
"GTTAGACTCTTTCGAT-1_12",
"GTTGAACCACTGTGAT-1_12",
"GTTGTCCCACTTGAAC-1_12",
"TATCTGTGTTGCGAAG-1_12",
"TCCGGGATCGCCTCTA-1_12",
"TCCTCCCCAAATGGTA-1_12",
"TTCACGCTCGACCAAT-1_12",
"TTGAGTGAGTATGCAA-1_12",
"AAATGGATCAACCTCC-1_14",
"ACAAGCTTCCACACCT-1_14",
"ACAAGCTTCCCTTGGT-1_14",
"CACGGGTGTGTAGCAG-1_14",
"CACTAAGGTGAAGCTG-1_14",
"CATTCATGTCGATTCA-1_14",
"GAGGCCTAGATGACAT-1_14",
"GGCAGTCAGGTCATCT-1_14",
"GTACAGTTCTGCGGAC-1_14",
"TTGATGGCAGGCGATA-1_14")

mammal.combined2 <- readRDS(file = paste(saveext,"Seurat_combined3.rds",sep=""))

uID <- as.character(Idents(mammal.combined2))
uID[which(uID%in%c("CTb_CS4","EPI_CS4","TB_CS3"))] <- "Tr"
uID[which(uID%in%c("CTb_CS5A/B","CTb_CS5C","CTb_CS6"))] <- "CTB"
uID[which(uID%in%c("EmDisc_CS5A/B","EmDisc_CS5C","EmDisc_CS6","EmDiscPS_CS6"))] <- "EmDisc"
uID[which(uID%in%c("EVTb_CS5C","EVTb_CS6"))] <- "EVT"
uID[which(uID%in%c("STb_CS5A/B","STb_CS5C","STb_CS6"))] <- "STB"
uID[which(uID%in%c("VE_CS4","VE_CS5","VE_CS5A/B","VE_CS5C","SYS_CS6"))] <- "VE"
uID[which(uID%in%c("ICM_CS3"))] <- "ICM"
uID[which(uID%in%c("Unciliated epithelia 1","Unciliated epithelia 2","Epith"))] <- "Unciliated epithelia"
uID[which(uID%in%c("Cil"))] <- "Ciliated epithelia"
uID[which(uID%in%c("Ciliated"))] <- "Ciliated epithelia"
uID[which(uID%in%c("Stroma"))] <- "Stroma fibroblasts"
uID[which(Idents(mammal.combined2)%in%c("STB_d9"))] <- "CTB_d9"
uID[which(Idents(mammal.combined2)%in%c("STB_d11"))] <- "CTB_d11"
uID[which(Idents(mammal.combined2)%in%c("STB_d12"))] <- "CTB_d12"
uID[which(Idents(mammal.combined2)%in%c("CTB_d9"))] <- "STB_d9"
uID[which(Idents(mammal.combined2)%in%c("CTB_d11"))] <- "STB_d11"
uID[which(Idents(mammal.combined2)%in%c("CTB_d12"))] <- "STB_d12"
#uID[which(uID%in%c("Hyp_d9","Hyp_d11","Hyp_d12"))] <- "VE"
#uID[which(uID%in%c("EmDisc_d9","EmDisc_d11","EmDisc_d12"))] <- "EmDisc"
uID[which(Idents(mammal.combined2)%in%c("CTB") & mammal.combined2$Dataset==c("10X Ours"))] <- "STB"
uID[which(Idents(mammal.combined2)%in%c("STB") & mammal.combined2$Dataset==c("10X Ours"))] <- "CTB"
uID[which(mammal.combined2$Cl15%in%c(8,18,19,23,4,12,15,20,22,28,29,31,32,37,38) & mammal.combined2$Genotype%in%c("Epithelial_G1","Epithelial_G2","Epithelial_G3","Epithelial_G4"))] <- "Tb_Epith_Mix"
uID[which(mammal.combined2$Cl15%in%c(8,18,19,23,4,12,15,20,22,28,29,31,32,37,38) & mammal.combined2$Genotype%in%c("Stromal_G1","Stromal_G2","Stromal_G3","Stromal_G4","Stromal_G5","Stromal_G6"))] <- "Tb_Stroma_Mix"
Idents(mammal.combined2) <- uID

Idents(mammal.combined2,cells=AmList) <- "Am"
Idents(mammal.combined2,cells=CTBList) <- "CTB"
Idents(mammal.combined2,cells=eEVTList) <- "CTB"
Idents(mammal.combined2,cells=eSTBList) <- "STB"
Idents(mammal.combined2,cells=EVTList) <- "EVT"
Idents(mammal.combined2,cells=ICMList) <- "ICM"
Idents(mammal.combined2,cells=PostEpiAME) <- "Am"
Idents(mammal.combined2,cells=PostEPIE1) <- "EmDisc"
Idents(mammal.combined2,cells=PostEPIE2) <- "EmDisc"

Idents(mammal.combined2,cells=PostEpiGast) <- "EmDisc"
Idents(mammal.combined2,cells=PrEn) <- "VE"
Idents(mammal.combined2,cells=Epi) <- "Epi"
Idents(mammal.combined2,cells=EpiPrEn) <- "Epi"
Idents(mammal.combined2,cells=STB) <- "STB"
Idents(mammal.combined2,cells=Tr) <- "Tr"

Idents(mammal.combined2,cells=pICM) <- "Unknown"
Idents(mammal.combined2,cells=pEVT) <- "Unknown"
Idents(mammal.combined2,cells=pSTB) <- "STB"
Idents(mammal.combined2,cells=pVE) <- "VE"
Idents(mammal.combined2,cells=pEpi) <- "Epi"
Idents(mammal.combined2,cells=pPS) <- "Unknown"
Idents(mammal.combined2,cells=pEmD) <- "Unknown"

mammal.combined2$Cells <- Idents(mammal.combined2)

mammal.combined2 <- subset(mammal.combined2,idents=c("Tb_Epith_Mix","Tb_Stroma_Mix","Unknown"),invert=TRUE)

saveRDS(mammal.combined2,file = paste(saveext,"Seurat_combined3_filtermix.rds",sep=""))
#osddsad

mammal.combined2test <- mammal.combined2
Idents(mammal.combined2test) <- mammal.combined2test$Cl05

p1 <- DimPlot(mammal.combined2test,  cols = coluse, pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_cltest","boundaryrp.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)


DefaultAssay(mammal.combined2) <- "RNA"
FeaturePlot(mammal.combined2,feature="WNT6",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_WNT6_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="WNT6",split.by = "Dataset",reduction = "umap",cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalUMAP_WNT6_EmbryonicMat_alldata",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="LGALS9",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_LGALS9_EmbryonicMat_alldata",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="HAND2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_HAND2_EmbryonicMat_alldata",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="CD44",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_CD44_EmbryonicMat_alldata",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="SFRP1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_SFRP1_EmbryonicMat_alldata",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="WNT5A",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_WNT5A_EmbryonicMat_alldata",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="RARRES2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_RARRES2_EmbryonicMat_alldata",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="GPR1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_GPR1_EmbryonicMat_alldata",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="C5AR1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_C5AR1_EmbryonicMat_alldata",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="NRG1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_NRG1_EmbryonicMat_alldata",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="ERBB4",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_ERBB4_EmbryonicMat_alldata",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="EGFR",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_EGFR_EmbryonicMat_alldata",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="COPA",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
ggsave(filename=paste(saveext,"FinalPCA_COPA_EmbryonicMat_alldata",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

DefaultAssay(mammal.combined2) <- "integrated"
#p1 <- DimPlot(mammal.combined2,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
#ggsave(filename=paste(saveext,"/DimRed/CompletePCA_1_2_EmbryonicMat_ano","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)
#
Idents(mammal.combined2,cells=WhichCells(mammal.combined2,
idents=c("CTB_d9","CTB_d11","CTB_d12","Hyp_d9","Hyp_d11","Hyp_d12","EmD/Hyp","EmDisc","Unkown_Emb","EmDisc_d9","EmDisc_d11","EmDisc_d12","STB_d9","STB_d11","STB_d12","Tr","STB","VE","ICM","EVT","Epi","STB","EVT","CTB","Am") )) <- "Embryonic"

Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents=c("Unciliated epithelia","Ciliated epithelia") )) <- "Epithelial"
Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents=c("Stroma fibroblasts","Stromal fibroblasts") )) <- "Stromal"

mammal.combined2$Cells1 <- Idents(mammal.combined2)

#Idents(mammal.combined2,cells=ACTB) <- "CTB_d14"
#Idents(mammal.combined2,cells=AHYP) <- "VE_d14"
#Idents(mammal.combined2,cells=ASTB) <- "STB_d14"
#Idents(mammal.combined2,cells=AEVT) <- "EVT_d14"
#Idents(mammal.combined2,cells=AEPI) <- "EmDisc_d14"

#Filter edometrial
#mammal.combined2$Cells2 <- Idents(mammal.combined2)
#Idents(mammal.combined2) <- mammal.combined2$Cells1 
mammal.combined2 <- subset(mammal.combined2,idents="Embryonic")

p1 <- DimPlot(mammal.combined2, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_1_2_EmbryonicMat_ano","boundary2.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)



mammal.combinedcl <- mammal.combined2

Idents(mammal.combinedcl) <- mammal.combinedcl$Cl05
p1 <- DimPlot(mammal.combinedcl,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_1_2_EmbryonicCl_ano","boundary2.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)


#Now what's left?
#Idents(mammal.combined2) <- mammal.combined2$Cells

uID <- as.character(mammal.combined2$Cells)
uID[which(mammal.combined2$Dataset=="10X Ours")] <- "STB_d14"
Idents(mammal.combined2) <- uID


#p1 <- DimPlot(mammal.combined2, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
#ggsave(filename=paste(saveext,"/DimRed/CompletePCA_1_2_EmbryonicMat_ano","boundary3.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)

#sadada
#uID[1:length(uID)] <- "STB_d14" #as.character(mammal.combined2$Cells)
uID[which(mammal.combined2$Dataset=="10X Ours" & mammal.combined2$Cl15%in%c("22","35","36"))] <- "EmDisc_d14"
uID[which(mammal.combined2$Dataset=="10X Ours" & mammal.combined2$Cl15%in%c("20"))] <- "CTB_d14"
uID[which(mammal.combined2$Dataset=="10X Ours" & mammal.combined2$Cl15%in%c("23"))] <- "Tr"

Idents(mammal.combined2) <- uID
Idents(mammal.combined2,cells=AHYP) <- "VE_d14"
Idents(mammal.combined2,cells=ACTB) <- "CTB_d14"
Idents(mammal.combined2,cells=AEVT) <- "EVT_d14"
#Idents(mammal.combined2,cells=VEcells) <- "VE_d14"

p1 <- DimPlot(mammal.combined2, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_1_2_EmbryonicMat_ano","boundary3.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)

#EmD <- WhichCells(mammal.combined3,idents=c("EmDisc_d14","STB_d14","VE_d14"))


list1<-c("GATA4",
"PDGFRA",
"CER1",
"SOX17",
"APOA1",
"APOB",
"POU5F1",
"SOX15",
"SFRP1",
"DNMT3B",
"PDGFA",
"TFAP2A",
"TFAP2C",
"WNT6",
"VTCN1",
"GATA2",
"GATA3",
"CGA",
"HLA-G",
"JAM3",
"DPPA3",
"DCN",
"MSLN")





uID <- Idents(mammal.combined2)
Idents(mammal.combined2) <- mammal.combined2$Dataset
ourcells <- WhichCells(mammal.combined2,idents="10X Ours")
mammal.combined2$iid <- uID
dat1 <- subset(mammal.combined2,idents="10X Ours")
Idents(mammal.combined2) <- uID


Idents(dat1) <- dat1$Cl05
Trcells <- colnames(dat1$Cl05[which(dat1$iid=="Tr")])


#Include Tr cells as a reference for some markers
#Idents(dat1,cells=Trcells) <- 20
dat1 <- subset(dat1,idents=c("20","19","14","13","12"))
Idents(dat1) <- dat1$iid
#
#D <- GetAssayData(dat1,assay="RNA")
#D <- D[list1,]

#redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
#pheatmap(D,color =  redblue1(20), border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes.pdf",sep=""),scale="row",width=10,height=10)


uID <- Idents(mammal.combined2)
Idents(mammal.combined2) <- mammal.combined2$Cl05
ourcells1 <- WhichCells(mammal.combined2,idents=c("20","19","14","13","12") )
Idents(mammal.combined2) <- uID

ourcells2 <- intersect(ourcells1, ourcells )


#dat1 <- mammal.combined2
D <- GetAssayData(dat1,assay="RNA")

saveRDS(D,file="subsetdata.rds")
D <- as.data.frame(D)[list1,]

mat_breaks <- seq(-2, 2, length.out = 20)

redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap((D),color =  redblue1(20), breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes.pdf",sep=""),scale="row",width=20,height=10)

#pheatmap(log2(D+1),color =  redblue1(20),kmeans_k=5, breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenesKmea.pdf",sep=""),scale="row",width=10,height=10)

res5 <- pheatmap((D),color =  redblue1(20), cutree_cols = 25,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes.pdf",sep=""),scale="row",width=20,height=10)#
#clust <- cbind(log2(D+1), cluster = cutree(res$tree_col, k = 25))
#saveRDS(clust, file = "Cluster.rds")

res1 <- pheatmap((D),color =  redblue1(20), cutree_cols = 5,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes5.pdf",sep=""),scale="row",width=20,height=10)
res2 <- pheatmap((D),color =  redblue1(20), cutree_cols = 10,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes10.pdf",sep=""),scale="row",width=20,height=10)
res3 <- pheatmap((D),color =  redblue1(20), cutree_cols = 15,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes15.pdf",sep=""),scale="row",width=20,height=10)
res4 <- pheatmap((D),color =  redblue1(20), cutree_cols = 20,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes20.pdf",sep=""),scale="row",width=20,height=10)
res6 <- pheatmap((D),color =  redblue1(20), cutree_cols = 30,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes30.pdf",sep=""),scale="row",width=20,height=10)
res7 <- pheatmap((D),color =  redblue1(20), cutree_cols = 35,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes35.pdf",sep=""),scale="row",width=20,height=10)


#clust <- cbind(log2(D+1), cluster = cutree(res$tree_col, k = 25))
saveRDS(D, file = "ClusterData.rds")

saveRDS(res1, file = "Cluster1.rds")
saveRDS(res2, file = "Cluster2.rds")
saveRDS(res3, file = "Cluster3.rds")
saveRDS(res4, file = "Cluster4.rds")
saveRDS(res5, file = "Cluster5.rds")
saveRDS(res6, file = "Cluster6.rds")
saveRDS(res7, file = "Cluster7.rds")




ct1 <- sort(cutree(res1$tree_col, k=5))
ct2 <- sort(cutree(res2$tree_col, k=10))
ct3 <- sort(cutree(res3$tree_col, k=15))
ct4 <- sort(cutree(res4$tree_col, k=20))
ct5 <- sort(cutree(res5$tree_col, k=25))
ct6 <- sort(cutree(res6$tree_col, k=30))
ct7 <- sort(cutree(res7$tree_col, k=35))

write.table(ct1,file="Cluster1.csv",quote = FALSE)
write.table(ct2,file="Cluster2.csv",quote = FALSE)
write.table(ct3,file="Cluster3.csv",quote = FALSE)
write.table(ct4,file="Cluster4.csv",quote = FALSE)
write.table(ct5,file="Cluster5.csv",quote = FALSE)
write.table(ct6,file="Cluster6.csv",quote = FALSE)
write.table(ct7,file="Cluster7.csv",quote = FALSE)


cl1 <- names(which(ct4==1))
cl2 <- names(which(ct4==2))
cl3 <- names(which(ct4==3))
cl4 <- names(which(ct4==4))
cl5 <- names(which(ct4==5))
cl6 <- names(which(ct4==6))
cl7 <- names(which(ct4==7))
cl8 <- names(which(ct4==8))
cl9 <- names(which(ct4==9))
cl10 <- names(which(ct4==10))
cl11 <- names(which(ct4==11))
cl12 <- names(which(ct4==12))
cl13 <- names(which(ct4==13))
cl14 <- names(which(ct4==14))
cl15 <- names(which(ct4==15))
cl16 <- names(which(ct4==16))
cl17 <- names(which(ct4==17))
cl18 <- names(which(ct4==18))
cl19 <- names(which(ct4==19))
cl20 <- names(which(ct4==20))

saveRDS(dat1,file="dat1.rds")
mammal.combined4B <- dat1
Idents(mammal.combined4B,cells=cl1) <- "Clu1"
Idents(mammal.combined4B,cells=cl2) <- "Clu2"
Idents(mammal.combined4B,cells=cl3) <- "Clu3"
Idents(mammal.combined4B,cells=cl4) <- "Clu4"
Idents(mammal.combined4B,cells=cl5) <- "Clu5"
Idents(mammal.combined4B,cells=cl6) <- "Clu6"
Idents(mammal.combined4B,cells=cl7) <- "Clu7"
Idents(mammal.combined4B,cells=cl8) <- "Clu8"
Idents(mammal.combined4B,cells=cl9) <- "Clu9"
Idents(mammal.combined4B,cells=cl10) <- "Clu10"
Idents(mammal.combined4B,cells=cl11) <- "Clu11"
Idents(mammal.combined4B,cells=cl12) <- "Clu12"
Idents(mammal.combined4B,cells=cl13) <- "Clu13"
Idents(mammal.combined4B,cells=cl14) <- "Clu14"
Idents(mammal.combined4B,cells=cl15) <- "Clu15"
Idents(mammal.combined4B,cells=cl16) <- "Clu16"
Idents(mammal.combined4B,cells=cl17) <- "Clu17"
Idents(mammal.combined4B,cells=cl18) <- "Clu18"
Idents(mammal.combined4B,cells=cl19) <- "Clu19"
Idents(mammal.combined4B,cells=cl20) <- "Clu20"

mammal.combined2test <- mammal.combined2



Idents(mammal.combined2test,cells=cl1) <- "VE_d14"
Idents(mammal.combined2test,cells=cl2) <- "Am/EmDisc_d14"
Idents(mammal.combined2test,cells=cl3) <- "Am_d14"
Idents(mammal.combined2test,cells=cl4) <- "Cont"
Idents(mammal.combined2test,cells=cl5) <- "Am/EmDisc_d14"
Idents(mammal.combined2test,cells=cl6) <- "VE_d14"
Idents(mammal.combined2test,cells=cl7) <- "VE_d14"
Idents(mammal.combined2test,cells=cl8) <- "Cont"
Idents(mammal.combined2test,cells=cl9) <- "Am_d14"
Idents(mammal.combined2test,cells=cl10) <- "Cont"
Idents(mammal.combined2test,cells=cl11) <- "Am/EmDisc_d14"
Idents(mammal.combined2test,cells=cl12) <- "EmDisc_d14"
Idents(mammal.combined2test,cells=cl13) <- "EmDisc_d14"
Idents(mammal.combined2test,cells=cl14) <- "Cont"
Idents(mammal.combined2test,cells=cl15) <- "CTB_d14"
Idents(mammal.combined2test,cells=cl16) <- "CTB_d14"
Idents(mammal.combined2test,cells=cl17) <- "CTB_d14"
Idents(mammal.combined2test,cells=cl18) <- "VE_d14"
Idents(mammal.combined2test,cells=cl19) <- "VE_d14"
Idents(mammal.combined2test,cells=cl20) <- "VE_d14"

colind <- integer( length( levels(Idents(mammal.combined2test)) )  )
for (i in 1:length( levels(Idents(mammal.combined2test)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2test))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined2test,  cols = coluse, pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_fullano",".pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)



#mammal.combined4B <- subset(mammal.combined4B,idents=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))
AvE <- AverageExpression(mammal.combined4B)
saveRDS(AvE,file="AvExp.rds")
AvE <- AvE$RNA
AvE <- AvE[list1,]

PB <- pheatmap(log2(AvE+1),color =  redblue1(20), breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenesPB.pdf",sep=""),scale="row",width=10,height=10)


mammal.combined4B <- subset(mammal.combined4B,idents=c("Clu1","Clu2","Clu3","Clu4","Clu5","Clu6","Clu7","Clu8","Clu9","Clu10","Clu11","Clu12","Clu13","Clu14","Clu15","Clu16","Clu17","Clu18","Clu19","Clu20"))



PB <- pheatmap(log2(AvE+1),color =  redblue1(20),  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenesPB2.pdf",sep=""),scale="row",width=10,height=10)





list1<-c(
"PDGFRA",
"SOX17",
"POU5F1",
"SOX15",
"SFRP1",
"DNMT3B",
"PDGFA",
"TFAP2A",
"TFAP2C",
"WNT6",
"VTCN1",
"GATA2",
"GATA3",
"CGA",
"HLA-G",
"JAM3",
"DPPA3",
"DCN",
"MSLN")

cluster1cells <- c(cl5,cl19)


res1 <- pheatmap((D[list1,cl5]),color =  redblue1(20), cutree_cols = 5,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes5_sse.pdf",sep=""),scale="row",width=20,height=10)
res2 <- pheatmap((D[list1,cl5]),color =  redblue1(20), cutree_cols = 10,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes10_sse.pdf",sep=""),scale="row",width=20,height=10)
res3 <- pheatmap((D[list1,cl5]),color =  redblue1(20), cutree_cols = 15,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes15_sse.pdf",sep=""),scale="row",width=20,height=10)
res4 <- pheatmap((D[list1,cl5]),color =  redblue1(20), cutree_cols = 20,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes20_sse.pdf",sep=""),scale="row",width=20,height=10)
res5 <- pheatmap((D[list1,cl5]),color =  redblue1(20), cutree_cols = 25,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes_sse.pdf",sep=""),scale="row",width=20,height=10)#
res6 <- pheatmap((D[list1,cl5]),color =  redblue1(20), cutree_cols = 30,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes30_sse.pdf",sep=""),scale="row",width=20,height=10)
res7 <- pheatmap((D[list1,cl5]),color =  redblue1(20), cutree_cols = 35,  breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenes35_sse.pdf",sep=""),scale="row",width=20,height=10)

ct4 <- sort(cutree(res6$tree_col, k=30))


ncl1 <- names(which(ct4==1))
ncl2 <- names(which(ct4==2))
ncl3 <- names(which(ct4==3))
ncl4 <- names(which(ct4==4))
ncl5 <- names(which(ct4==5))
ncl6 <- names(which(ct4==6))
ncl7 <- names(which(ct4==7))
ncl8 <- names(which(ct4==8))
ncl9 <- names(which(ct4==9))
ncl10 <- names(which(ct4==10))
ncl11 <- names(which(ct4==11))
ncl12 <- names(which(ct4==12))
ncl13 <- names(which(ct4==13))
ncl14 <- names(which(ct4==14))
ncl15 <- names(which(ct4==15))
ncl16 <- names(which(ct4==16))
ncl17 <- names(which(ct4==17))
ncl18 <- names(which(ct4==18))
ncl19 <- names(which(ct4==19))
ncl20 <- names(which(ct4==20))
ncl21 <- names(which(ct4==21))
ncl22 <- names(which(ct4==22))
ncl23 <- names(which(ct4==23))
ncl24 <- names(which(ct4==24))
ncl25 <- names(which(ct4==25))
ncl26 <- names(which(ct4==26))
ncl27 <- names(which(ct4==27))
ncl28 <- names(which(ct4==28))
ncl29 <- names(which(ct4==29))
ncl30 <- names(which(ct4==30))




Idents(mammal.combined4B,cells=ncl1) <- "Clu1B"
Idents(mammal.combined4B,cells=ncl2) <- "Clu2B"
Idents(mammal.combined4B,cells=ncl3) <- "Clu3B"
Idents(mammal.combined4B,cells=ncl4) <- "Clu4B"
Idents(mammal.combined4B,cells=ncl5) <- "Clu5B"
Idents(mammal.combined4B,cells=ncl6) <- "Clu6B"
Idents(mammal.combined4B,cells=ncl7) <- "Clu7B"
Idents(mammal.combined4B,cells=ncl8) <- "Clu8B"
Idents(mammal.combined4B,cells=ncl9) <- "Clu9B"
Idents(mammal.combined4B,cells=ncl10) <- "Clu10B"
Idents(mammal.combined4B,cells=ncl11) <- "Clu11B"
Idents(mammal.combined4B,cells=ncl12) <- "Clu12B"
Idents(mammal.combined4B,cells=ncl13) <- "Clu13B"
Idents(mammal.combined4B,cells=ncl14) <- "Clu14B"
Idents(mammal.combined4B,cells=ncl15) <- "Clu15B"
Idents(mammal.combined4B,cells=ncl16) <- "Clu16B"
Idents(mammal.combined4B,cells=ncl17) <- "Clu17B"
Idents(mammal.combined4B,cells=ncl18) <- "Clu18B"
Idents(mammal.combined4B,cells=ncl19) <- "Clu19B"
Idents(mammal.combined4B,cells=ncl20) <- "Clu20B"
Idents(mammal.combined4B,cells=ncl21) <- "Clu21B"
Idents(mammal.combined4B,cells=ncl22) <- "Clu22B"
Idents(mammal.combined4B,cells=ncl23) <- "Clu23B"
Idents(mammal.combined4B,cells=ncl24) <- "Clu24B"
Idents(mammal.combined4B,cells=ncl25) <- "Clu25B"
Idents(mammal.combined4B,cells=ncl26) <- "Clu26B"
Idents(mammal.combined4B,cells=ncl27) <- "Clu27B"
Idents(mammal.combined4B,cells=ncl28) <- "Clu28B"
Idents(mammal.combined4B,cells=ncl29) <- "Clu29B"
Idents(mammal.combined4B,cells=ncl30) <- "Clu30B"


AvE <- AverageExpression(mammal.combined4B)
saveRDS(AvE,file="AvExp.rds")
AvE <- AvE$RNA
AvE <- AvE[list1,]

PB <- pheatmap(log2(AvE+1),color =  redblue1(20), breaks=mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, filename = paste(saveext,"/Markers/LineageHeatmap_allgenesPBsse.pdf",sep=""),scale="row",width=10,height=10)

saveRDS(res1, file = "Cluster1sse.rds")
saveRDS(res2, file = "Cluster2sse.rds")
saveRDS(res3, file = "Cluster3sse.rds")
saveRDS(res4, file = "Cluster4sse.rds")
saveRDS(res5, file = "Cluster5sse.rds")
saveRDS(res6, file = "Cluster6sse.rds")
saveRDS(res7, file = "Cluster7sse.rds")

mammal.combined4B <- mammal.combined2
#Idents(mammal.combined4B) <- mammal.combined4B$Cells

Idents(mammal.combined4B,cells=cl1) <- "VE_d14"
Idents(mammal.combined4B,cells=cl2) <- "Am_d14"
Idents(mammal.combined4B,cells=cl3) <- "Am_d14"
Idents(mammal.combined4B,cells=cl4) <- "Cont"
Idents(mammal.combined4B,cells=cl5) <- "Am/EmDisc_d14"
Idents(mammal.combined4B,cells=cl6) <- "VE_d14"
Idents(mammal.combined4B,cells=cl7) <- "VE_d14"
Idents(mammal.combined4B,cells=cl8) <- "Cont"
Idents(mammal.combined4B,cells=cl9) <- "Am_d14"
Idents(mammal.combined4B,cells=cl10) <- "Cont"
Idents(mammal.combined4B,cells=cl11) <- "EmDisc_d14"
Idents(mammal.combined4B,cells=cl12) <- "EmDisc_d14"
Idents(mammal.combined4B,cells=cl13) <- "EmDisc_d14"
Idents(mammal.combined4B,cells=cl14) <- "Cont"
Idents(mammal.combined4B,cells=cl15) <- "CTB_d14"
Idents(mammal.combined4B,cells=cl16) <- "CTB_d14"
Idents(mammal.combined4B,cells=cl17) <- "CTB_d14"
Idents(mammal.combined4B,cells=cl18) <- "VE_d14"
Idents(mammal.combined4B,cells=cl19) <- "VE_d14"
Idents(mammal.combined4B,cells=cl20) <- "VE_d14"

Idents(mammal.combined4B,cells=ncl1) <- "Am_d14"
Idents(mammal.combined4B,cells=ncl2) <- "Am_d14"
Idents(mammal.combined4B,cells=ncl3) <- "EmDisc_d14"
Idents(mammal.combined4B,cells=ncl4) <- "Am_d14"
Idents(mammal.combined4B,cells=ncl5) <- "Cont"
Idents(mammal.combined4B,cells=ncl6) <- "Cont"
Idents(mammal.combined4B,cells=ncl7) <- "Cont"
Idents(mammal.combined4B,cells=ncl8) <- "Am_d14"
Idents(mammal.combined4B,cells=ncl9) <- "Cont"
Idents(mammal.combined4B,cells=ncl10) <- "EmDisc_d14"
Idents(mammal.combined4B,cells=ncl11) <- "Am_d14"
Idents(mammal.combined4B,cells=ncl12) <- "EmDisc_d14"
Idents(mammal.combined4B,cells=ncl13) <- "EmDisc_d14"
Idents(mammal.combined4B,cells=ncl14) <- "EmDisc_d14"
Idents(mammal.combined4B,cells=ncl15) <- "CTB_d14"
Idents(mammal.combined4B,cells=ncl16) <- "EmDisc_d14"
Idents(mammal.combined4B,cells=ncl17) <- "Cont"
Idents(mammal.combined4B,cells=ncl18) <- "EmDisc_d14"
Idents(mammal.combined4B,cells=ncl19) <- "Cont"
Idents(mammal.combined4B,cells=ncl20) <- "Cont"
Idents(mammal.combined4B,cells=ncl21) <- "Am_d14"
Idents(mammal.combined4B,cells=ncl22) <- "Am_d14"
Idents(mammal.combined4B,cells=ncl23) <- "CTB_d14"
Idents(mammal.combined4B,cells=ncl24) <- "Cont"
Idents(mammal.combined4B,cells=ncl25) <- "Am_d14"
Idents(mammal.combined4B,cells=ncl26) <- "Am_d14"
Idents(mammal.combined4B,cells=ncl27) <- "VE_d14"
Idents(mammal.combined4B,cells=ncl28) <- "CTB_d14"
Idents(mammal.combined4B,cells=ncl29) <- "Cont"
Idents(mammal.combined4B,cells=ncl30) <- "Cont"


saveRDS(mammal.combined4B,file="OurDataAnotated.rds")



mammal.combined4B <- subset(mammal.combined4B,idents=c("Cont"),invert=TRUE)

amcells1 <- colnames(mammal.combined4B)[which(mammal.combined4B$Dataset=="10X Reference 1" & mammal.combined4B$Cl05=="14" & Idents(mammal.combined4B)%in%c("CTB_d9","STB_d9") )]
amcells2 <- colnames(mammal.combined4B)[which(mammal.combined4B$Dataset=="10X Reference 1" & mammal.combined4B$Cl05=="14" & Idents(mammal.combined4B)%in%c("CTB_d11","STB_d11"))]
amcells3 <- colnames(mammal.combined4B)[which(mammal.combined4B$Dataset=="10X Reference 1" & mammal.combined4B$Cl05=="14" & Idents(mammal.combined4B)%in%c("CTB_d12","STB_d12"))]

Idents(mammal.combined4B,cells=amcells1) <- "Am_d9"
Idents(mammal.combined4B,cells=amcells2) <- "Am_d11"
Idents(mammal.combined4B,cells=amcells3) <- "Am_d12"


colind <- integer( length( levels(Idents(mammal.combined4B)) )  )
for (i in 1:length( levels(Idents(mammal.combined4B)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined4B))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined4B, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/ClusterBasedAnotation",".pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)


p1 <- DimPlot(mammal.combined4B, cols = coluse,   pt.size = 4, reduction = "umap", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/ClusterBasedAnotationUMAP",".pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)


saveRDS(GetAssayData(mammal.combined4B,assay="RNA"),file="LabelledExp.rds")
saveRDS(mammal.combined4B$Dataset,file="LabelledDataset.rds")
saveRDS(Idents(mammal.combined4B),file="LabelledIDs.rds")






#sadadasdasdd
#
#mammal.combined3 <- mammal.combined2
#Concells <- WhichCells(mammal.combined3,expression=MSLN>0.2)
#Idents(mammal.combined3,cells=intersect(Concells,ourcells)) <- "Cont"
#colind <- integer( length( levels(Idents(mammal.combined3)) )  )
#for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
#  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
#}
#coluse <- cols[colind]

#p1 <- DimPlot(mammal.combined3, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
#ggsave(filename=paste(saveext,"/DimRed/MarkerFilter_Cont1","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)
#
#mammal.combined2 <- subset(mammal.combined3,idents="Cont",invert=TRUE)


#mammal.combined3 <- mammal.combined2
#Concells <- WhichCells(mammal.combined3,expression=DCN>0.1)
#Idents(mammal.combined3,cells=intersect(Concells,ourcells)) <- "Cont"
#colind <- integer( length( levels(Idents(mammal.combined3)) )  )
#for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
#  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
#}
#coluse <- cols[colind]

#p1 <- DimPlot(mammal.combined3, cols = coluse,   pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
#ggsave(filename=paste(saveext,"/DimRed/MarkerFilter_Cont2","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE)

#mammal.combined2 <- subset(mammal.combined3,idents="Cont",invert=TRUE)



DefaultAssay(mammal.combined2) <- "RNA"


#Idents(mammal.combined2) <- mammal.combined2$Dataset
#mammal.combined2 <- subset(mammal.combined2,idents="10X Ours")

FeaturePlot(mammal.combined2,feature="LGALS9",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_LGALS9_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="HAND2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_HAND2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="CD44",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_CD44_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="NOTUM",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_NOTUM_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="MMP2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_MMP2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="FSTL3",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_FSTL3_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="TOP2A",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_TOP2A_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="PARP1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_PARP1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="FSTL3",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_FSTL3_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="FN1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_FN1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="WNT6",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_WNT6_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="IFI6",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_IFI6_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="FSTL3",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_FSTL3_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="FSTL3",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_FSTL3_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="SFRP1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SFRP1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="WNT5A",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_WNT5A_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="RARRES2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_RARRES2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="GPR1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_GPR1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="C5AR1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_C5AR1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="NRG1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_NRG1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="ERBB4",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_ERBB4_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="EGFR",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_EGFR_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="COPA",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_COPA_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)



FeaturePlot(mammal.combined2,feature="ZIC2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_ZIC2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="VIM",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_VIM_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="ZEB2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_ZEB2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="FOXF1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_FOXF1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="BST2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_BST2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)




#FeaturePlot(mammal.combined2,feature="COPA",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4)
#ggsave(filename=paste(saveext,"FinalPCA_COPA_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
#

FeaturePlot(mammal.combined2,feature="ROR2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_ROR2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

#
#FeaturePlot(mammal.combined2,feature="TJP1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
#gsave(filename=paste(saveext,"FinalPCA_ZO-1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

#FeaturePlot(mammal.combined2,feature="TJP2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
#ggsave(filename=paste(saveext,"FinalPCA_ZO-2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)




FeaturePlot(mammal.combined2,feature="FOXJ1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_FOXJ1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="TFAP2C",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_TFAP2C_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="PDGFA",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_PDGFA_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="WNT5A",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_WNT5A_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="MUC4",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_MUC4_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="MSLN",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_MSLN_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="FGF2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_FGF2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="VTCN1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_VTCN1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="TFAP2A",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_TFAP2A_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="ISL1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_ISL1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="POU5F1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_POU5F1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="OTX2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_OTX2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="FSTL3",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_FSTL3_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="EOMES",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_EOMES_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="S100P",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_S100P_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="SDC1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SDC1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="PGF",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_PGF_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="DKK1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_DKK1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="SFRP2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SFRP2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="TDGF1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_TDGF1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="DPPA5",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_DPPA5_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)



#FinalPCA_TIMP2_EmbryonicMat.pdf

FeaturePlot(mammal.combined2,feature="TIMP2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_TIMP2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
#FeaturePlot(mammal.combined2,feature="ID2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
#ggsave(filename=paste(saveext,"FinalPCA_ID2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="NRP1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_NRP1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="DMKN",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_DMKN_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)



FeaturePlot(mammal.combined2,feature="SOX15",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SOX15_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="PDGFRA",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_PDGFRA_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="BAMBI",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_BAMBI_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="DIO2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_DIO2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="HLA-G",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_HLA-G_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="JAM3",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_JAM3_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="APOA1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_APOA1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="GATA4",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_GATA4_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="CER1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_CER1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="CGA",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_CGA_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="GATA3",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_GATA3_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="GATA3",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_GATA3_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="DIO2",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_DIO2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="BAMBI",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_BAMBI_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="FOXJ1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_FOXJ1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="DCN",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_DCN_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="WNT5A",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_WNT5A_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="WNT7A",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_WNT7A_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="LGR5",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_LGR5_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="SOX9",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SOX9_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)



FeaturePlot(mammal.combined2,feature="PAEP",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_PAEP_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="SPP1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SPP1_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="SCGB2A2",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SCGB2A2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="WNT6",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_WNT6_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined2,feature="WNT6",split.by = "Dataset",reduction = "umap",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalUMAP_WNT6_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="GABRP",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_GABRP_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="DNMT3B",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_DNMT3B_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="SFRP2",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SFRP2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="STAT3",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_STAT3_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined2,feature="HAND2",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_HAND2_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined2,feature="TBX4",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_TBX4_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)

#FinalPCA_ROR2_EmbryonicMat.pdf

mammal.combined4C <- subset(mammal.combined4B,idents=c("EmDisc_d14","VE_d14","Am_d14"))

#asdsadada

DefaultAssay(mammal.combined4C) <- "RNA"

#Idents(mammal.combined2) <- mammal.combined2$Cells 
#mammal.combined2a <- mammal.combined2 #subset(mammal.combined2,idents="Embryonic")

FeatureScatter(mammal.combined4C,feature1="CGA",feature2="GATA4")
ggsave(filename=paste(saveext,"Scatter_CGA_GATA4_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined4C,feature1="GATA2",feature2="GATA4")
ggsave(filename=paste(saveext,"Scatter_GATA2_GATA4_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined4C,feature1="CGA",feature2="PDGFA")
ggsave(filename=paste(saveext,"Scatter_CGA_PDGFA_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined4C,feature1="GATA4",feature2="PDGFA")
ggsave(filename=paste(saveext,"Scatter_GATA4_PDGFA_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeatureScatter(mammal.combined4C,feature1="POU5F1",feature2="PDGFA")
ggsave(filename=paste(saveext,"Scatter_POU5F1_PDGFA_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined4C,feature1="POU5F1",feature2="GATA4")
ggsave(filename=paste(saveext,"Scatter_POU5F1_GATA4_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined4C,feature1="CGA",feature2="POU5F1")
ggsave(filename=paste(saveext,"Scatter_CGA_POU5F1_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined4C,feature1="GATA2",feature2="POU5F1")
ggsave(filename=paste(saveext,"Scatter_POU5F1_GATA2_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeatureScatter(mammal.combined4C,feature1="SOX15",feature2="POU5F1")
ggsave(filename=paste(saveext,"Scatter_POU5F1_SOX15_EmbryonicMat",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


#FeaturePlot(mammal.combined2,feature="WNT5A",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "black"))
#ggsave(filename=paste(saveext,"FinalPCA_WNT5A_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
#FeaturePlot(mammal.combined2,feature="WNT5A",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "black"))
#ggsave(filename=paste(saveext,"FinalPCA_WNT5A_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)#
#FeaturePlot(mammal.combined2,feature="WNT5A",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "black"))
#ggsave(filename=paste(saveext,"FinalPCA_WNT5A_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
#FeaturePlot(mammal.combined2,feature="WNT5A",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "black"))
#ggsave(filename=paste(saveext,"FinalPCA_WNT5A_EmbryonicMat",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)



#mammal.combined4C <- subset(mammal.combined4B,idents=c("EmDisc_d14","VE_d14","Am_d14"))

FeaturePlot(mammal.combined4C,feature="WNT6",split.by = "Dataset",reduction = "umap",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalUMAP_WNT6_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="GABRP",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_GABRP_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined4C,feature="DNMT3B",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_DNMT3B_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="SFRP2",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SFRP2_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="STAT3",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_STAT3_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="HAND2",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_HAND2_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="TBX4",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_TBX4_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined4C,feature="EOMES",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_EOMES_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="T",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_T_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="SOX2",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SOX2_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="NANOG",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_NANOG_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="TDGF1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_TDGF1_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="PDGFRA",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_PDGFRA_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="PDGFA",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_PDGFA_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="WNT8A",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_WNT8A_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="VTCN1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_VTCN1_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="GABRP",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_GABRP_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="MIXL1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_MIXL1_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined4C,feature="TDGF1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_GATA3_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="TDGF1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_GATA3_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined4C,feature="SOX17",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SOX17_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="PRDM1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_PRDM1_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="TFAP2A",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_TFAP2A_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="TFAP2C",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_TFAP2C_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="SNAI1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SNAI1_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="SNAI2",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_SNAI2_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)



FeaturePlot(mammal.combined4C,feature="NODAL",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_NODAL_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)



#FeaturePlot(mammal.combined4C,feature="NODAL",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF000,sep=""),width = 32, height = 10, limitsize = FALSE)0"),pt.size = 4, keep.scale="all")
#ggsave(filename=paste(saveext,"FinalPCA_NODAL_Embryonic",".pdf"


FeaturePlot(mammal.combined4C,feature="POU5F1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_POU5F1_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="ZIC2",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_ZIC2_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="VIM",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_VIM_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="ZEB2",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_ZEB3_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined4C,feature="FOXF1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_FOXF1_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="APOB",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_APOB_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="APOA1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_APOA1_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="BST2",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_BST2_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined4C,feature="BAMBI",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_BAMBI_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="APOA1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_APOA1_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="GATA4",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_GATA4_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="CER1",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_CER1_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="OTX2",split.by = "Dataset",reduction = "pca", cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_OTX2_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="LEFTY1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_LEFTY1_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="LEFTY2",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_LEFTY2_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)


FeaturePlot(mammal.combined4C,feature="HGF",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_HGF_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="HAND1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_HAND1_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="GC",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_GC_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="VCAN",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_VCAN_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

FeaturePlot(mammal.combined4C,feature="KLF4",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_KLF4_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="VCAN",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_VCAN_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)
FeaturePlot(mammal.combined4C,feature="KDR",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_KDR_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)



FeaturePlot(mammal.combined4C,feature="COL6A1",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
ggsave(filename=paste(saveext,"FinalPCA_COL6A1_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)

#FeaturePlot(mammal.combined4C,feature="SP5",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
#ggsave(filename=paste(saveext,"FinalPCA_SP5_Embryonic",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)
#FeaturePlot(mammal.combined4C,feature="VCAN",split.by = "Dataset",reduction = "pca",cols =  c("lightgrey", "#FF0000"),pt.size = 4, keep.scale="all")
#ggsave(filename=paste(saveext,"FinalPCA_VCAN_Embryonic",".pdf",sep=""),width = 32, height = 10, limitsize = FALSE)


saveRDS(mammal.combined4C,file="EmDiscAmVE.rds")

pgc_signature_gene_list <-list(c("SOX17","PRDM1", "PRDM14", "TFAP2A", "TFAP2C", "NANOS3","POU5F1"))


mammal.combined4C <- AddModuleScore(object = mammal.combined4C, features = pgc_signature_gene_list, name = "pgc_signature_gene_list")
FeaturePlot(object = mammal.combined4C, features = "pgc_signature_gene_list") #+scale_color_viridis(discrete = FALSE, option="turbo")
ggsave(filename=paste(saveext,"FinalPCA_PGC_Embryonic",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE)



#require(tidyverse)
#X <- mammal.combined2$Dataset
#Y <- mammal.combined2$Cl15   
#Z <- as.character(Idents(mammal.combined2))

#Z2 <- mammal.combined2$Genotype

#Dat <- data.frame(x=Y[which(X=="SS2 Reference 1")],y=Z[which(X=="SS2 Reference 1")] )
#Dat1 <- Dat %>%  count(y, x)
#ggplot(Dat1, aes(fill=y, x=x, y=n)) + geom_bar(position="fill", stat="identity")
#ggsave(filename=paste(saveext,"Table1","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)

#Dat <- data.frame(x=Y[which(X=="10X Reference 1")],y=Z[which(X=="10X Reference 1")] )
#Dat1 <- Dat %>%  count(y, x)
#ggplot(Dat1, aes(fill=y, x=x, y=n)) + geom_bar(position="fill", stat="identity")
#ggsave(filename=paste(saveext,"Table2","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)


#Dat <- data.frame(x=Y[which(X=="10X Ours")],y=Z2[which(X=="10X Ours")] )
#Dat1 <- Dat %>%  count(y, x)
#ggplot(Dat1, aes(fill=y, x=x, y=n)) + geom_bar(position="fill", stat="identity")
#ggsave(filename=paste(saveext,"Table3","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)


saddsadsdasds


mammal.combined3 <- readRDS(file = paste(saveext,"Seurat_combined3.rds",sep=""))
Idents(mammal.combined3) <- mammal.combined3$Dataset
mammal.combined3 <- subset(mammal.combined3,idents=c("10X Ours"))


uID <- as.character(mammal.combined3$Genotype)
uID[which(is.na(uID)==1)] <- "NGT"
Idents(mammal.combined3) <- uID

Idents(mammal.combined3,cells=WhichCells(mammal.combined3,idents=c("NAssigned"))) <- "NGT"
Idents(mammal.combined3,cells=WhichCells(mammal.combined3,idents=c("G1"))) <- "Stromal_G1"
Idents(mammal.combined3,cells=WhichCells(mammal.combined3,idents=c("G3"))) <- "Epithelial_G3"
Idents(mammal.combined3,cells=WhichCells(mammal.combined3,idents=c("G4"))) <- "Epithelial_G4"
Idents(mammal.combined3,cells=WhichCells(mammal.combined3,idents=c("Embryoni_G3"))) <- "Embryonic_G3"

colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined3, cols=coluse,  pt.size = 4, reduction = "pca", label = TRUE, split.by="ID3", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_1_2_Genotype_ano","boundary.pdf",sep=""),width = 72, height = 10, useDingbats = FALSE, limitsize = FALSE)

mammal.combined3 <- readRDS(file = paste(saveext,"Seurat_combined3.rds",sep=""))

Idents(mammal.combined3) <- mammal.combined3$Cl05

cType <- c("4",
"8",
"12",
"14",
"16",
"19",
"20",
"21",
"0",
"2",
"3",
"5",
"6",
"7",
"10",
"13",
"15",
"18",
"1",
"9",
"11",
"17")

cols <- c("#00BFBF",
"#00ABBF",
"#0197BF",
"#0183BF",
"#016FBF",
"#015BBF",
"#0247BF",
"#0233BF",
"#FAA734",
"#F8A32E",
"#F6A028",
"#F39C23",
"#F1981D",
"#EF9517",
"#ED9111",
"#EA8D0C",
"#E88A06",
"#E68600",
"#D9AD0D",
"#C49C0A",
"#AE8A06",
"#997903")

colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- cols[colind]

p1 <- DimPlot(mammal.combined3, cols=coluse,  pt.size = 4, reduction = "pca", label = TRUE, split.by="Dataset", repel = TRUE, dims=c(1,2))
ggsave(filename=paste(saveext,"/DimRed/CompletePCA_1_2_Cl05_ano","boundary.pdf",sep=""),width = 32, height = 10, useDingbats = FALSE, limitsize = FALSE)






mammal.combined3 <- readRDS(file = paste(saveext,"Seurat_combined3.rds",sep=""))

uID <- as.character(Idents(mammal.combined3))
uID[which(uID%in%c("CTb_CS4","EPI_CS4","TB_CS3"))] <- "Tr"
uID[which(uID%in%c("CTb_CS5A/B","CTb_CS5C","CTb_CS6"))] <- "CTB"
uID[which(uID%in%c("EmDisc_CS5A/B","EmDisc_CS5C","EmDisc_CS6","EmDiscPS_CS6"))] <- "EmDisc"
uID[which(uID%in%c("EVTb_CS5C","EVTb_CS6"))] <- "EVT"
uID[which(uID%in%c("STb_CS5A/B","STb_CS5C","STb_CS6"))] <- "STB"
uID[which(uID%in%c("VE_CS4","VE_CS5","VE_CS5A/B","VE_CS5C","SYS_CS6"))] <- "VE"
uID[which(uID%in%c("ICM_CS3"))] <- "ICM"
uID[which(uID%in%c("Unciliated epithelia 1","Unciliated epithelia 2","Epith"))] <- "Unciliated epithelia"
uID[which(uID%in%c("Cil"))] <- "Ciliated epithelia"
uID[which(uID%in%c("Ciliated"))] <- "Ciliated epithelia"
uID[which(uID%in%c("Stroma"))] <- "Stroma fibroblasts"
uID[which(Idents(mammal.combined3)%in%c("STB_d9"))] <- "CTB_d9"
uID[which(Idents(mammal.combined3)%in%c("STB_d11"))] <- "CTB_d11"
uID[which(Idents(mammal.combined3)%in%c("STB_d12"))] <- "CTB_d12"
uID[which(Idents(mammal.combined3)%in%c("CTB_d9"))] <- "STB_d9"
uID[which(Idents(mammal.combined3)%in%c("CTB_d11"))] <- "STB_d11"
uID[which(Idents(mammal.combined3)%in%c("CTB_d12"))] <- "STB_d12"
#uID[which(uID%in%c("Hyp_d9","Hyp_d11","Hyp_d12"))] <- "VE"
#uID[which(uID%in%c("EmDisc_d9","EmDisc_d11","EmDisc_d12"))] <- "EmDisc"
uID[which(Idents(mammal.combined3)%in%c("CTB") & mammal.combined3$Dataset==c("10X Ours"))] <- "STB"
uID[which(Idents(mammal.combined3)%in%c("STB") & mammal.combined3$Dataset==c("10X Ours"))] <- "CTB"
uID[which(mammal.combined3$Cl15%in%c(8,18,19,23,4,12,15,20,22,28,29,31,32,37,38) & mammal.combined3$Genotype%in%c("Epithelial_G1","Epithelial_G2","Epithelial_G3","Epithelial_G4"))] <- "Tb_Epith_Mix"
uID[which(mammal.combined3$Cl15%in%c(8,18,19,23,4,12,15,20,22,28,29,31,32,37,38) & mammal.combined3$Genotype%in%c("Stromal_G1","Stromal_G2","Stromal_G3","Stromal_G4","Stromal_G5","Stromal_G6"))] <- "Tb_Stroma_Mix"
Idents(mammal.combined3) <- uID

Idents(mammal.combined3,cells=AmList) <- "Am"
Idents(mammal.combined3,cells=CTBList) <- "CTB"
Idents(mammal.combined3,cells=eEVTList) <- "CTB"
Idents(mammal.combined3,cells=eSTBList) <- "STB"
Idents(mammal.combined3,cells=EVTList) <- "EVT"
Idents(mammal.combined3,cells=ICMList) <- "ICM"
Idents(mammal.combined3,cells=PostEpiAME) <- "Am"
Idents(mammal.combined3,cells=PostEPIE1) <- "EmDisc"
Idents(mammal.combined3,cells=PostEPIE2) <- "EmDisc"

Idents(mammal.combined3,cells=PostEpiGast) <- "EmDisc"
Idents(mammal.combined3,cells=PrEn) <- "VE"
Idents(mammal.combined3,cells=Epi) <- "Epi"
Idents(mammal.combined3,cells=EpiPrEn) <- "Epi"
Idents(mammal.combined3,cells=STB) <- "STB"
Idents(mammal.combined3,cells=Tr) <- "Tr"

Idents(mammal.combined3,cells=pICM) <- "Unknown"
Idents(mammal.combined3,cells=pEVT) <- "Unknown"
Idents(mammal.combined3,cells=pSTB) <- "STB"
Idents(mammal.combined3,cells=pVE) <- "VE"
Idents(mammal.combined3,cells=pEpi) <- "Epi"
Idents(mammal.combined3,cells=pPS) <- "Unknown"
Idents(mammal.combined3,cells=pEmD) <- "Unknown"

Idents(mammal.combined3,cells=WhichCells(mammal.combined3,
idents=c("CTB_d9","CTB_d11","CTB_d12","Hyp_d9","Hyp_d11","Hyp_d12","EmD/Hyp","EmDisc","Unkown_Emb","EmDisc_d9","EmDisc_d11","EmDisc_d12","STB_d9","STB_d11","STB_d12","Tr","STB","VE","ICM","EVT","Epi","STB","EVT","CTB","Am") )) <- "Embryonic"
Idents(mammal.combined3,cells=WhichCells(mammal.combined3,idents=c("Unciliated epithelia","Ciliated epithelia") )) <- "Epithelial"
Idents(mammal.combined3,cells=WhichCells(mammal.combined3,idents=c("Stroma fibroblasts","Stromal fibroblasts") )) <- "Stromal"

mammal.combined3 <- subset(mammal.combined3,idents="Unknown",invert=TRUE)


require(tidyverse)
X <- mammal.combined2$Dataset
X2 <- mammal.combined2$ID3
Y <- factor(mammal.combined2$Cl05,levels=c(4,8,12,14,16,19,20,15,5,7,21,0,2,3,6,10,13,18,1,9,11,17))
Z <- as.character(Idents(mammal.combined2))

Z2 <- mammal.combined2$Genotype
Z2[which(Z2=="NAssigned")] <- "NGT"

Z2[which(Z2=="G1")] <- "Stromal_G1"
Z2[which(Z2=="G3")] <- "Epithelial_G3"
Z2[which(Z2=="G4")] <- "Epithelial_G4"
Z2[which(Z2=="Embryoni_G3")] <- "Embryonic_G3"


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



Z2[which(Z2=="Epithelial_G1" & mammal.combined3$ID3=="C3")] <- "Embryonic"
Z2[which(Z2=="Epithelial_G6" & mammal.combined3$ID3=="C3")] <- "Embryonic"
Z2[which(Z2=="Epithelial_G3" & mammal.combined3$ID3=="C4")] <- "Embryonic"





Dat <- data.frame(x=Y[which(X=="SS2 Reference 1")],y=Z[which(X=="SS2 Reference 1")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"Table1","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)


Dat <- data.frame(x=Y[which(X=="10X Reference 1")],y=Z[which(X=="10X Reference 1")] )
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n, col=y)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#0767DA","#E68600","#BF9600"))+ theme_classic() 
ggsave(filename=paste(saveext,"Table2","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)



#c("#0767DA","#E68600","#BF9600")

#Dat <- data.frame(x=Y[which(X=="10X Reference 1")],y=Z[which(X=="10X Reference 1")] )
#Dat1 <- Dat %>%  count(y, x)
#ggplot(Dat1, aes(fill=y, x=x, y=n)) + geom_bar(position="fill", stat="identity")
#ggsave(filename=paste(saveext,"Table2","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)


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



sadsadsadsad


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

