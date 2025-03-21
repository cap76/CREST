Idents(Endo3) <- Endo3$ID5
Endo4 <- subset(Endo3,idents=c("SOX9","Ciliated","Lumenal","Glandular","eS","dS","Fibroblast C7"))
D_s <- subset(D,idents=c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular","Stromal fibroblasts"))

strom <- WhichCells(D_s,idents="Stromal fibroblasts")
cil <- WhichCells(D_s,idents="Ciliated")

ours <- colnames(D_s)

newIDs <- as.character(Idents(D_s))
newIDs[which(D_s$ID3 %in% c("C1","C2","C3","C4","C6","C7","C8") )] <- "Lumenal"
SOX <- WhichCells(D_s,expression=SOX9> log(1) )
LRG <- WhichCells(D_s,expression=LGR5>log(1)  )
SPP <- WhichCells(D_s,expression=SPP1>log(1)  )
SCGB <- WhichCells(D_s,expression=SCGB2A2> log(1)  )

Idents(D_s) <- newIDs
Idents(D_s, cells= intersect(SPP,ours)) <- "Glandular"
Idents(D_s, cells=intersect(SCGB,ours)) <- "Glandular"
Idents(D_s, cells=intersect(SOX,ours)) <- "SOX9"
Idents(D_s, cells=intersect(LRG,ours)) <- "LRG5"

Idents(D_s, cells=strom) <- "Stroma"
Idents(D_s, cells=cil) <- "Ciliated"

mammal.combined <- readRDS(file=paste(saveext,"/Integrate_with_RVT_stroma_withModuleScores.rds",sep=""))

Idents(D_s, cells=intersect(colnames(mammal.combined)[which(mammal.combined$prolifS1>mammal.combined$dS1)],ours)) <- "Prolf Strom"

mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo4,D_s), dims = 1:30, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_epithel.pdf",sep=""),width = 10, height = 8,p)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_3.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,2)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_1_2.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)


saveRDS(mammal.combined,file=paste(saveext,"/Alignwith_RVT_allEndometrial.rds",sep=""))


Markers <- c(
  "SOX9",
  "LGR5",
  "PGR",
  "ESR1",
  "MMP7",
  "CPM",
  "MKI67",
  "HMGB2",
  "PLAU",
  "IL32",
  "TNF",
  "WNT7A",
  "FOXJ1",
  "PIFO",
  "ABCG1",
  "SCGB2A2",
  "C2CD4A",
  "SLC18A2",
  "PAEP",
  "CXCL14",
  "SPP1",
  "PAX2",
  "VTCN1",
  "SLC26A7",
  "MSLN",
  "ACTA2",
 "PCOLCE",
  "MMP11",
  "ECM1",
  "FOXO1")


Markers2 <- c(
  "ACTA2",
  "C7",
  "IGF1",
  "PCOLCE",
  "MMP11",
  "ECM1",
  "FOXO1",
  "IL15",
  "CFD",
  "CEBPB",
  "PDGFA")



uID <- as.character(mammal.combined$ID3)
uID[which(uID%in%c("C1","C2","C3","C4","C6","C7","C8"))] <- "Ours"

mammal.combined$ID0 <- Idents(mammal.combined)

Idents(mammal.combined) <- paste(uID,mammal.combined$ID0)

 
Idents(mammal.combined) <- as.character(Idents(mammal.combined))
   

catorder <- c("proliferative Fibroblast C7","early-secretory Fibroblast C7","late-secretory Fibroblast C7",
              "proliferative dS","early-secretory dS","late-secretory dS","Ours Stroma",
              "proliferative eS","early-secretory eS","late-secretory eS","Ours Prolf Strom",
              "proliferative SOX9","early-secretory SOX9","late-secretory SOX9","Ours SOX9","Ours LRG5",
              "proliferative Lumenal","early-secretory Lumenal","late-secretory Lumenal","Ours Lumenal",
              "proliferative Glandular","early-secretory Glandular","late-secretory Glandular","Ours Glandular",
              "proliferative Ciliated","early-secretory Ciliated","late-secretory Ciliated","Ours Ciliated")


Idents(mammal.combined) <- factor(Idents(mammal.combined), levels = catorder)


list1 <- c("proliferative Ciliated","early-secretory Ciliated","late-secretory Ciliated","proliferative Glandular","Ours Ciliated",
           "early-secretory Glandular","late-secretory Glandular","Ours Glandular","proliferative Lumenal","early-secretory Lumenal","late-secretory Lumenal",
           "Ours Lumenal","proliferative SOX9","early-secretory SOX9","late-secretory SOX9","Ours SOX9","Ours LRG5","proliferative eS",
           "early-secretory eS","late-secretory eS","Ours Prolf Strom","proliferative dS","early-secretory dS","late-secretory dS","Ours Stroma")

DefaultAssay(mammal.combined) <- "RNA"
p<-DotPlot(mammal.combined, features = Markers, cols = c("lightgrey", "red"), assay = "RNA",  dot.scale =4,scale =FALSE )
ggsave(filename=paste(saveext,"/DimRed/EpitheliaDotPlot.pdf",sep=""),width = 10, height = 8,p, limitsize = FALSE)

mammal.combinedsub <- subset(mammal.combined,idents =  c("Ours Ciliated","Ours Glandular","Ours Lumenal","Ours SOX9","Ours LRG5","Ours Prolf Strom","Ours Stroma") )


p<-DotPlot(mammal.combinedsub, features = Markers, cols = c("lightgrey", "red"), idents = c("Ours Ciliated","Ours Glandular","Ours Lumenal","Ours SOX9","Ours LRG5","Ours Prolf Strom","Ours Stroma") , assay = "RNA",  dot.scale =4, scale =FALSE)
ggsave(filename=paste(saveext,"/DimRed/EpitheliaDotPlotSubset.pdf",sep=""),width = 10, height = 8,p, limitsize = FALSE)



p<-DotPlot(mammal.combinedsub, features = Markers2, cols = c("lightgrey", "red"), idents = c("Ours Ciliated","Ours Glandular","Ours Lumenal","Ours SOX9","Ours LRG5","Ours Prolf Strom","Ours Stroma") , assay = "RNA" ,  dot.scale =4, scale=FALSE)
ggsave(filename=paste(saveext,"/DimRed/StromaDotPlotSubset.pdf",sep=""),width = 10, height = 8,p, limitsize = FALSE)
p<-DotPlot(mammal.combined, features = Markers, cols = c("lightgrey", "red"), assay = "RNA",  dot.scale =4, scale=FALSE)
ggsave(filename=paste(saveext,"/DimRed/StromaDotPlot.pdf",sep=""),width = 10, height = 8,p, limitsize = FALSE)



set1 <- c("proliferative dS","early-secretory dS","late-secretory dS","proliferative eS","early-secretory eS","late-secretory eS")

set2 <- c("proliferative SOX9","early-secretory SOX9","late-secretory SOX9",
              "proliferative Lumenal","early-secretory Lumenal","late-secretory Lumenal",
              "proliferative Glandular","early-secretory Glandular","late-secretory Glandular",
              "proliferative Ciliated","early-secretory Ciliated","late-secretory Ciliated")

set3 <- c("Ours Stroma","Ours Prolf Strom")
set4 <- c("Ours SOX9","Ours LRG5","Ours Lumenal","Ours Glandular","Ours Ciliated")


mammal.combined2 <- mammal.combined
Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents=set1)) <- "Stroma"
Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents=set2)) <- "Epithelial"
Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents=set3)) <- "Stroma_Ours"
Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents=set4)) <- "Epithelial_Ours"




List1 <- FindMarkers(mammal.combined2,ident.1 = "Stroma_Ours",ident.2 = c("Epithelial_Ours"), test.use = "MAST")
List2 <- FindMarkers(mammal.combined2,ident.1 = "Stroma", ident.2 = "Epithelial", test.use = "MAST")


List3 <- FindMarkers(mammal.combined2,ident.1 = "Stroma_Ours",ident.2 = c("Epithelial_Ours"), test.use = "MAST",only.pos=TRUE )
List4 <- FindMarkers(mammal.combined2,ident.1 = "Stroma", ident.2 = "Epithelial", test.use = "MAST",only.pos=TRUE)

List5 <- FindMarkers(mammal.combined2,ident.2 = "Stroma_Ours",ident.1 = c("Epithelial_Ours"), test.use = "MAST" ,only.pos=TRUE)
List6 <- FindMarkers(mammal.combined2,ident.2 = c("Stroma"), ident.1 = "Epithelial", test.use = "MAST",only.pos=TRUE)



DefaultAssay(mammal.combined2) <- "RNA"
Av <- AverageExpression(mammal.combined2)
Cl1 <- List1
Ae5 <- as.data.frame(Av$RNA)
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_log2FC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log2(0.01))
Ae5$ColInd <- 0

Ae5[ which( (Ae5$FC1)>log(1.2) & Ae5$Pval1>-log2(0.01)), c("ColInd")] <- 1
Ae5[ which( (Ae5$FC1)< -log(1.2) & Ae5$Pval1>-log2(0.01)), c("ColInd")] <- 2
Ae5$ColInd <- as.factor(Ae5$ColInd)
Ae5$AvExp <- log(0.5*(Ae5$'Epithelial_Ours' + Ae5$'Stroma_Ours') +1)
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
  geom_hline(aes(yintercept = 0))

genes.to.label <- genes.to.label <-unique(c(intersect(rownames(List3),rownames(List4)),intersect(rownames(List5),rownames(List6))))

p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_Stroma_vs_epithelial.pdf",sep=""),width = 20, height = 10, plot = p1, useDingbats=FALSE)
dev.off()








set1 <- c("proliferative dS","early-secretory dS","late-secretory dS")
set2 <- c("proliferative eS","early-secretory eS","late-secretory eS")


mammal.combined2 <- mammal.combined
Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents=set1)) <- "dStroma"
Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents=set2)) <- "pStroma"
#Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents=set3)) <- "Stroma_Ours"
#Idents(mammal.combined2,cells=WhichCells(mammal.combined2,idents=set4)) <- "Epithelial_Ours"
#"Ours Stroma","Ours Prolf Strom"

List1 <- FindMarkers(mammal.combined2,ident.1 = "Ours Stroma",ident.2 = c("Ours Prolf Strom"), test.use = "MAST")
List2 <- FindMarkers(mammal.combined2,ident.1 = "dStroma", ident.2 = "pStroma", test.use = "MAST")


List3 <- FindMarkers(mammal.combined2,ident.1 = "Ours Stroma",ident.2 = c("Ours Prolf Strom"), test.use = "MAST",only.pos=TRUE )
List4 <- FindMarkers(mammal.combined2,ident.1 = "dStroma", ident.2 = "pStroma", test.use = "MAST",only.pos=TRUE)

List5 <- FindMarkers(mammal.combined2,ident.2 = "Ours Stroma",ident.1 = c("Ours Prolf Strom"), test.use = "MAST" ,only.pos=TRUE)
List6 <- FindMarkers(mammal.combined2,ident.2 = "dStroma", ident.1 = "pStroma", test.use = "MAST",only.pos=TRUE)



DefaultAssay(mammal.combined2) <- "RNA"
Av <- AverageExpression(mammal.combined2)
Cl1 <- List1
Ae5 <- as.data.frame(Av$RNA)
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_log2FC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log2(0.01))
Ae5$ColInd <- 0

Ae5[ which( (Ae5$FC1)>log(1.2) & Ae5$Pval1>-log2(0.01)), c("ColInd")] <- 1
Ae5[ which( (Ae5$FC1)< -log(1.2) & Ae5$Pval1>-log2(0.01)), c("ColInd")] <- 2
Ae5$ColInd <- as.factor(Ae5$ColInd)
Ae5$AvExp <- log(0.5*(Ae5$'Ours Prolf Strom' + Ae5$'Ours Stroma') +1)
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
  geom_hline(aes(yintercept = 0))

genes.to.label <- genes.to.label <-unique(c(intersect(rownames(List3),rownames(List4)),intersect(rownames(List5),rownames(List6))))

p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_Stroma_dec_vs_prolif.pdf",sep=""),width = 20, height = 10, plot = p1, useDingbats=FALSE)
dev.off()

Idents(mammal.combined) <- mammal.combined$ID0
mammal.combined1 <- subset(mammal.combined,idents=c("eS","dS","Fibroblast C7","Prolf Strom","Stroma"))

uID <- as.character(mammal.combined1$ID3)
uID[which(uID%in%c("C1","C2","C4","C3","C6","C7","C8"))] <- "C1"
mammal.combined1$IDX <- uID
data <- data.frame(x=as.factor(mammal.combined1$IDX),y=as.factor(mammal.combined1$ID0) )
p0 <- ggplot(data, aes(fill=y, x=x, y=y)) + geom_bar(position="fill", stat="identity") + theme_bw() + coord_cartesian(ylim = c(0, 1))
ggsave(filename=paste(saveext,"/BarPlot_Stroma_Prolf.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p0)


