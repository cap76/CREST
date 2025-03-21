library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library("Matrix")

set.seed(1)

#Same folder as before
saveext = "./Align_with_new/"
dir.create(saveext)
dir.create(paste(saveext,"DimRed/",sep=""))
dir.create(paste(saveext,"Markers/",sep=""))

#Cell anotation and colours
cols <- c("#0c9cf5","#0979E3","#0556D1","#0233BF","#F04C04",
"#E04704","#CF4104","#BF3C04","#921FE6","#8717D1",
"#7C10BB","#7108A6","#877BD6","#6355B5",
"#3E2E94","#1A0873","BF0489","#AA8600",
"#E68600","#d17600","#BF9600","#BF7104",
"black","#921FE6")

cType <- c("EmDisc_d9","EmDisc_d11","EmDisc_d12","EmDisc_d14",
"Hyp_d9","Hyp_d11","Hyp_d12","Hyp_d14","CTB_d9","CTB_d11",
"CTB_d12","CTB_d14","STB_d9","STB_d11","STB_d12","STB_d14","EVT_d14",
"Ciliated","Unciliated epithelia 1","Unciliated epithelia 2","Stromal fibroblasts","Unciliated epithelia 3","Unciliated epithelia 4",
"EVT")

#Save the clustered data
mammal.combined <- readRDS(file = paste(saveext,"Seurat_combined2.rds",sep=""))

newID <- as.character(mammal.combined$Cells)
newID[which(mammal.combined$Cl05%in%c(0,10,12,16,18,19) & mammal.combined$species1=="New")] <- "Stromal fibroblasts"
newID[which(mammal.combined$Cl05%in%c(8) & mammal.combined$species1=="New")] <- "Ciliated"
newID[which(mammal.combined$Cl05%in%c(3,4,5,7) & mammal.combined$species1=="New")] <- "Unciliated epithelia 2"
newID[which(mammal.combined$Cl05%in%c(20) & mammal.combined$species1=="New")] <- "Unciliated epithelia 4"

newID[which(mammal.combined$Cl05%in%c(1) & mammal.combined$species1=="New")] <- "Unciliated epithelia 1"
newID[which(mammal.combined$Cl05%in%c(13) )] <- "Unciliated epithelia 3"

newID[which(mammal.combined$Cl05%in%c(17) & mammal.combined$species1=="New")] <- "Hyp_d14"
newID[which(mammal.combined$Cl05%in%c(19) & mammal.combined$species1=="New")] <- "EmDisc_d14"

newID[which(mammal.combined$Cl05%in%c(2,15) & mammal.combined$species1=="New")] <- "CTB_d14"
newID[which(mammal.combined$Cl05%in%c(6,9,11,14) & mammal.combined$species1=="New")] <- "STB_d14"
newID[which(mammal.combined$Cl05%in%c(21) & mammal.combined$species1=="New")] <- "EVT"
Idents(mammal.combined) <- newID

saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined_ano1.rds",sep=""))


colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- cols[colind]

DimPlot(mammal.combined, cols = coluse, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat_splitIICol",".pdf",sep=""),width = 60, height = 10, limitsize = FALSE, useDingbats=FALSE)
#DimPlot(mammal.combined, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"UMAP_SeuratIICol",".pdf",sep=""),width = 60, height = 10, limitsize = FALSE, useDingbats=FALSE)


#sadsad
#DimPlot(mammal.combined, cols = coluse,  reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"UMAP_Seurat_splitIICol",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)
#DimPlot(mammal.combined, cols = coluse, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"UMAP_SeuratIICol",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE)

DimPlot(mammal.combined,cols = coluse,  reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat_splitIICol",".pdf",sep=""),width = 60, height = 10, limitsize = FALSE, useDingbats=FALSE)
#DimPlot(mammal.combined, cols = coluse, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"PCA_SeuratIICol",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE)

#saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined_ano1.rds",sep=""))

