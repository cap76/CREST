library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library("Matrix")
library(harmony)

set.seed(1)

saveext = "./EndometriumCPDB/"
dir.create(saveext)
dir.create(paste(saveext,"DimRed/",sep=""))
dir.create(paste(saveext,"Marker/",sep=""))

EndoD <- read.table('Endometrium/GSE111976_summary_10x_day_donor_ctype.csv',sep=",", header = T, row.names=1)
Endo_data <- readRDS(file = "Endometrium/GSE111976_ct_endo_10x.rds")
Endom <- CreateSeuratObject(counts = Endo_data, assay = "RNA",min.cells = 0, min.features = 0)
Endom$ID1 <- "Endometrium"
Endom$ID2 <- EndoD[,3] #Type
Endom$ID3 <- EndoD[,2] #Patient
Endom$Time <- EndoD[,1] #Patient Time

Endom$ID4 <- EndoD[,2] #Patient
Endom$ID5 <- EndoD[,2] #Patient
Idents(Endom) <- EndoD[,3]
Endom <- NormalizeData(Endom, verbose = FALSE)
Endom <- FindVariableFeatures(Endom, selection.method = "vst", nfeatures = 20000)

Endom <- ScaleData(Endom, verbose = FALSE)
Endom <- RunPCA(Endom, npcs = 20, verbose = FALSE)
Endom <- FindNeighbors(Endom, reduction = "pca", dims = 1:20)
Endom <- RunUMAP(Endom, reduction = "pca", dims = 1:20, n.neighbors = 47)
Endom <- FindNeighbors(Endom, reduction = "pca", dims = 1:20)


DimPlot(Endom,  reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE,pt.size = 2)
ggsave(filename=paste(saveext,"/All_UMAP_split_mergs",".pdf",sep=""),width = 100, height = 10,limitsize = FALSE, useDingbats=FALSE)

DimPlot(Endom,  reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, pt.size = 2)
ggsave(filename=paste(saveext,"/All_PCA_split_mergs",".pdf",sep=""),width = 100, height = 10,limitsize = FALSE, useDingbats=FALSE)


Endom <- Endom %>% RunHarmony("ID3", plot_convergence = TRUE)

Endom <- RunUMAP(Endom, reduction = "harmony", dims = 1:20)

DimPlot(Endom,  reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE,pt.size = 2)
ggsave(filename=paste(saveext,"/All_UMAP_split_harmony_mergs",".pdf",sep=""),width = 100, height = 10,limitsize = FALSE, useDingbats=FALSE)





saddsadsdads

#14 19 20 29 39 41 57 58 60 63


#Placenta <- readRDS(Placenta,file=paste(saveext,"/All_mergs_Seurat3",".pdf",sep=""))

#sdadsa



Idents(Endom) <- Endom$ID3
Data1 <- subset(Endom, idents="14") #No
Data2 <- subset(Endom, idents="19") #No
Data3 <- subset(Endom, idents="20") #Yes
Data4 <- subset(Endom, idents="29") #Yes
Data5 <- subset(Endom, idents="39") #Yes
Data6 <- subset(Endom, idents="41") #Yes
Data7 <- subset(Endom, idents="57") #Yes
Data8 <- subset(Endom, idents="58") #Yes
Data9 <- subset(Endom, idents="60") #Yes
Data10 <- subset(Endom, idents="63") #Yes


Idents(Data1) <- Data1$ID2
Idents(Data2) <- Data2$ID2
Idents(Data3) <- Data3$ID2
Idents(Data4) <- Data4$ID2
Idents(Data5) <- Data5$ID2
Idents(Data6) <- Data6$ID2
Idents(Data7) <- Data6$ID2
Idents(Data8) <- Data6$ID2
Idents(Data9) <- Data6$ID2
Idents(Data10) <- Data6$ID2


dir.create(paste(saveext,"/Embryo1/",sep=""))
writeMM(Data1@assays$RNA@data, file = paste(saveext,'/Embryo1/matrix.mtx',sep="") )
write(x = rownames(Data1@assays$RNA@data), file = paste(saveext,"/Embryo1/features.tsv",sep="") )
write(x = colnames(Data1@assays$RNA@data), file = paste(saveext,"/Embryo1/barcodes.tsv",sep="") )
Data1$cell_type <-  paste("Cl", Data1$ID2, sep="")
Data1@meta.data$Cell = rownames(Data1@meta.data)
df = Data1@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'Embryo1_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)

dir.create(paste(saveext,"/Embryo2/",sep=""))
writeMM(Data2@assays$RNA@data, file = paste(saveext,'/Embryo2/matrix.mtx',sep="") )
write(x = rownames(Data2@assays$RNA@data), file = paste(saveext,"/Embryo2/features.tsv",sep="") )
write(x = colnames(Data2@assays$RNA@data), file = paste(saveext,"/Embryo2/barcodes.tsv",sep="") )
Data2$cell_type <-  paste("Cl", Data2$ID2, sep="")
Data2@meta.data$Cell = rownames(Data2@meta.data)
df = Data2@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'Embryo2_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)

dir.create(paste(saveext,"/Embryo3/",sep=""))
writeMM(Data3@assays$RNA@data, file = paste(saveext,'/Embryo3/matrix.mtx',sep="") )
write(x = rownames(Data3@assays$RNA@data), file = paste(saveext,"/Embryo3/features.tsv",sep="") )
write(x = colnames(Data3@assays$RNA@data), file = paste(saveext,"/Embryo3/barcodes.tsv",sep="") )
Data3$cell_type <-  paste("Cl", Data3$ID2, sep="")
Data3@meta.data$Cell = rownames(Data3@meta.data)
df = Data3@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'Embryo3_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)

dir.create(paste(saveext,"/Embryo4/",sep=""))
writeMM(Data4@assays$RNA@data, file = paste(saveext,'/Embryo4/matrix.mtx',sep="") )
write(x = rownames(Data4@assays$RNA@data), file = paste(saveext,"/Embryo4/features.tsv",sep="") )
write(x = colnames(Data4@assays$RNA@data), file = paste(saveext,"/Embryo4/barcodes.tsv",sep="") )
Data4$cell_type <-  paste("Cl", Data4$ID2, sep="")
Data4@meta.data$Cell = rownames(Data4@meta.data)
df = Data4@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'Embryo4_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)

dir.create(paste(saveext,"/Embryo5/",sep=""))
writeMM(Data5@assays$RNA@data, file = paste(saveext,'/Embryo5/matrix.mtx',sep="") )
write(x = rownames(Data5@assays$RNA@data), file = paste(saveext,"/Embryo5/features.tsv",sep="") )
write(x = colnames(Data5@assays$RNA@data), file = paste(saveext,"/Embryo5/barcodes.tsv",sep="") )
Data5$cell_type <-  paste("Cl", Data5$ID2, sep="")
Data5@meta.data$Cell = rownames(Data5@meta.data)
df = Data5@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'Embryo5_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)

dir.create(paste(saveext,"/Embryo6/",sep=""))
writeMM(Data6@assays$RNA@data, file = paste(saveext,'/Embryo6/matrix.mtx',sep="") )
write(x = rownames(Data6@assays$RNA@data), file = paste(saveext,"/Embryo6/features.tsv",sep="") )
write(x = colnames(Data6@assays$RNA@data), file = paste(saveext,"/Embryo6/barcodes.tsv",sep="") )
Data6$cell_type <-  paste("Cl", Data6$ID2, sep="")
Data6@meta.data$Cell = rownames(Data6@meta.data)
df = Data6@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'Embryo6_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)

dir.create(paste(saveext,"/Embryo7/",sep=""))
writeMM(Data7@assays$RNA@data, file = paste(saveext,'/Embryo7/matrix.mtx',sep="") )
write(x = rownames(Data7@assays$RNA@data), file = paste(saveext,"/Embryo7/features.tsv",sep="") )
write(x = colnames(Data7@assays$RNA@data), file = paste(saveext,"/Embryo7/barcodes.tsv",sep="") )
Data7$cell_type <-  paste("Cl", Data7$ID2, sep="")
Data7@meta.data$Cell = rownames(Data7@meta.data)
df = Data7@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'Embryo7_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)

dir.create(paste(saveext,"/Embryo8/",sep=""))
writeMM(Data8@assays$RNA@data, file = paste(saveext,'/Embryo8/matrix.mtx',sep="") )
write(x = rownames(Data8@assays$RNA@data), file = paste(saveext,"/Embryo8/features.tsv",sep="") )
write(x = colnames(Data8@assays$RNA@data), file = paste(saveext,"/Embryo8/barcodes.tsv",sep="") )
Data8$cell_type <-  paste("Cl", Data8$ID2, sep="")
Data8@meta.data$Cell = rownames(Data8@meta.data)
df = Data8@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'Embryo8_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)

dir.create(paste(saveext,"/Embryo9/",sep=""))
writeMM(Data9@assays$RNA@data, file = paste(saveext,'/Embryo9/matrix.mtx',sep="") )
write(x = rownames(Data9@assays$RNA@data), file = paste(saveext,"/Embryo9/features.tsv",sep="") )
write(x = colnames(Data9@assays$RNA@data), file = paste(saveext,"/Embryo9/barcodes.tsv",sep="") )
Data9$cell_type <-  paste("Cl", Data9$ID2, sep="")
Data9@meta.data$Cell = rownames(Data9@meta.data)
df = Data9@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'Embryo9_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)

dir.create(paste(saveext,"/Embryo10/",sep=""))
writeMM(Data10@assays$RNA@data, file = paste(saveext,'/Embryo10/matrix.mtx',sep="") )
write(x = rownames(Data10@assays$RNA@data), file = paste(saveext,"/Embryo10/features.tsv",sep="") )
write(x = colnames(Data10@assays$RNA@data), file = paste(saveext,"/Embryo10/barcodes.tsv",sep="") )
Data10$cell_type <-  paste("Cl", Data10$ID2, sep="")
Data10@meta.data$Cell = rownames(Data10@meta.data)
df = Data10@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'Embryo10_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)



Endom <- Endom %>% RunHarmony("ID3", plot_convergence = TRUE)

Endom <- RunUMAP(Endom, reduction = "harmony", dims = 1:20)




DimPlot(Endom,  reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE,pt.size = 2)
ggsave(filename=paste(saveext,"/All_UMAP_split_harmony_mergs",".pdf",sep=""),width = 100, height = 10,limitsize = FALSE, useDingbats=FALSE)


