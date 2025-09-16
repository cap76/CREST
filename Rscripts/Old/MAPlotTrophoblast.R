library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(pheatmap)
library("readxl")
library(stringr)

set.seed(1) 

saveext = "~/Desktop/Data/Endometrial/InVitro/Matteo/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

D1 <- readRDS(file=paste(saveext,"MatteoD1.rds",sep=""))
D2 <- readRDS(file=paste(saveext,"MatteoD2.rds",sep=""))
D3 <- readRDS(file=paste(saveext,"MatteoD3.rds",sep=""))
D6 <- readRDS(file=paste(saveext,"MatteoD6.rds",sep=""))
D7 <- readRDS(file=paste(saveext,"MatteoD7.rds",sep=""))
D8 <- readRDS(file=paste(saveext,"MatteoD8.rds",sep=""))

DefaultAssay(D1) <- "RNA"
DefaultAssay(D2) <- "RNA"
DefaultAssay(D3) <- "RNA"
DefaultAssay(D6) <- "RNA"
DefaultAssay(D7) <- "RNA"
DefaultAssay(D8) <- "RNA"

Idents(D1) <- D1$Cells
Idents(D2) <- D2$Cells
Idents(D3) <- D3$Cells
Idents(D6) <- D6$Cells
Idents(D7) <- D7$Cells
Idents(D8) <- D8$Cells

#[1] Lumenal             SOX9P               Hyp_d14             Prolif              putExMes            ExMes_d14           Stromal fibroblasts STB_d14             putSTB             
#[10] SOX9LRG5            Am/EmDisc_d14       EmDisc_d14          Am_d14              Hyp/Am              Ciliated            CTB_d14             Glandular     

#D1 <- subset(D1,idents=c("Hyp_d14","putExMes","ExMes_d14","STB_d14","putSTB","Am/EmDisc_d14","EmDisc_d14","Hyp/Am","Am_d1","CTB_d14")) 


#Now load mamroset
#Now load in other datasets
Key20307 <- read.table("~/Desktop/Data/Embryonic/Marmoset/CRUK_SLX-20307/allQC20307.txt",sep="\t",header = T, row.names=1)
raw_counts20307 <- read.table("~/Desktop/Data/Embryonic/Marmoset/CRUK_SLX-20307/featurecountsAll_extended_SLX-20307.csv",sep=",",header = T, row.names=1)

Key20308 <- read.table("~/Desktop/Data/Embryonic/Marmoset/CRUK_SLX-20308/allQC20308.txt",sep="\t",header = T, row.names=1)
raw_counts20308 <- read.table("~/Desktop/Data/Embryonic/Marmoset/CRUK_SLX-20308/featurecountsAll_extended_SLX-20308.csv",sep=",",header = T, row.names=1)

KeyCS5 <- read.table("~/Desktop/Data/Embryonic/Marmoset/CS5/CS5Key_210122.csv",sep=",",header = T, row.names=1)
raw_countsCS5 <- read.table("~/Desktop/Data/Embryonic/Marmoset/CS5/featurecountsCS5.csv",sep=",",header = T, row.names=1)

KeyCS6 <- read.table("~/Desktop/Data/Embryonic/Marmoset/CS6/CS6Key_210402.csv",sep=",",header = T, row.names=1)
raw_countsCS6 <- read.table("~/Desktop/Data/Embryonic/Marmoset/CS6/featurecountsCS6.csv",sep=",",header = T, row.names=1)

KeyCS7 <- read.table("~/Desktop/Data/Embryonic/Marmoset/CS7/CS7Key_210217.csv",sep=",",header = T, row.names=1) 
raw_countsCS7 <- read.table("~/Desktop/Data/Embryonic/Marmoset/CS7/featurecountsCS7.csv",sep=",",header = T, row.names=1)

KeyCRUK <- read.table("~/Desktop/Data/Embryonic/Marmoset/CRUK_SLX-20304/CRUKKey2.csv",sep=",",header = T, row.names=1) 
raw_countsCRUK <- read.table("~/Desktop/Data/Embryonic/Marmoset/CRUK_SLX-20304/featurecounts-CRUK.csv",sep=",",header = T, row.names=1)

KeyCRUK2 <- read.table("~/Desktop/Data/Embryonic/Marmoset/CRUK_210301_A00489_0804_AH3FMGDRXY/allQC.csv",sep=",",header = T, row.names=1) 
raw_countsCRUK2 <- read.table("~/Desktop/Data/Embryonic/Marmoset/CRUK_210301_A00489_0804_AH3FMGDRXY/featurecounts.csv",sep=",",header = T, row.names=1)

KeyCRUK3 <- read.table("~/Desktop/Data/Embryonic/Marmoset/CRUK_210303_A00489_0808_BH3NYWDRXY/allQC.csv",sep=",",header = T, row.names=1) 
raw_countsCRUK3 <- read.table("~/Desktop/Data/Embryonic/Marmoset/CRUK_210303_A00489_0808_BH3NYWDRXY/featurecounts.csv",sep=",",header = T, row.names=1)

KeyPre <- read.table("~/Desktop/Data/Embryonic/Marmoset/CS1-3/CS1-3Key.csv",sep=",",header = T, row.names=1) 
raw_countsPre <- read.table("~/Desktop/Data/Embryonic/Marmoset/CS1-3/featurecountsCS1-3.csv",sep=",",header = T, row.names=1)


#Create Seurat arrays
marmoset_data_20307 <- CreateSeuratObject(counts = raw_counts20307[,which(Key20307$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_20307) <- Key20307$Primary.Annotation[which(Key20307$QC>0)]
marmoset_data_20307$LOC <- Key20307$Location[which(Key20307$QC>0)]
marmoset_data_20307$Stage <- "CS6"
marmoset_data_20307$Batch <- "New"
marmoset_data_20307 <- subset(marmoset_data_20307, idents = c("Tb_CS6","EmDisc_Gast_CS6","ExMes_CS6","Am_CS6","SYS_CS6","PGC_CS6","VE_CS6","EmDisc_CS6","Am_CS6_EmDisc","Stalk_CS6")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
marmoset_data_20307 <- NormalizeData(marmoset_data_20307, verbose = FALSE)
marmoset_data_20307$Dataset <- "InVivo"

marmoset_data_20308 <- CreateSeuratObject(counts = raw_counts20308[,which(Key20308$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_20308) <- Key20308$Primary.Annotation[which(Key20308$QC>0)]
marmoset_data_20308$LOC <- Key20308$Location[which(Key20308$QC>0)]
marmoset_data_20308$Stage <- "CS7"
marmoset_data_20308$Batch <- "New"
marmoset_data_20308 <- subset(marmoset_data_20308, idents = c("Tb_CS7","EmDisc_Gast_CS7","VE_CS7","EmDisc_CS7","Am_CS7","ExMes_stalk_CS7","PGC_CS7","SYS_CS7","ExMes_CS7","Am_CS7_EmDisc","ReStroma_CS7")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
marmoset_data_20308 <- NormalizeData(marmoset_data_20308, verbose = FALSE)
marmoset_data_20308$Dataset <- "InVivo"

marmoset_data_Pre <- CreateSeuratObject(counts = raw_countsPre[,which(KeyPre$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_Pre) <- KeyPre$Primary.lineage[which(KeyPre$QC>0)]
marmoset_data_Pre$LOC <- KeyPre$Location[which(KeyPre$QC>0)]
marmoset_data_Pre <- subset(marmoset_data_Pre, idents = c("Tb_CS3","Epi_CS3","ICM_CS3","cMor_CS3","Hyp_CS3","Zy_CS1","8-cell_CS2","4-cell_CS2"))
marmoset_data_Pre <- NormalizeData(marmoset_data_Pre, verbose = FALSE)
marmoset_data_Pre$Stage <- "Pre"
marmoset_data_Pre$Batch <- "Old"
marmoset_data_Pre$Dataset <- "InVivo"

marmoset_data_CRUK <- CreateSeuratObject(counts = raw_countsCRUK[,which(KeyCRUK$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_CRUK) <- KeyCRUK$Primary.Annotation[which(KeyCRUK$QC>0)]
marmoset_data_CRUK$LOC <- KeyCRUK$Loc[which(KeyCRUK$QC>0)]
marmoset_data_CRUK$Stage <- KeyCRUK$Stage[which(KeyCRUK$QC>0)]
marmoset_data_CRUK$Batch <- "New"
marmoset_data_CRUK <- subset(marmoset_data_CRUK, idents = c("Tb_CS5","EmDisc_CS6","Am_CS7","VE_CS5","EmDisc_CS7","ExMes_CS6","ExMes_CS7","Am_CS6","SYS_CS7","EmDisc_CS5","EmDisc_CS7_Am","VE_CS6","ExMes_CS5","Am_CS5_PGC","Am_CS5","SYS_CS5","Stalk_CS7","Stalk_CS6")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
marmoset_data_CRUK <- NormalizeData(marmoset_data_CRUK, verbose = FALSE)
marmoset_data_CRUK$Dataset <- "InVivo"

marmoset_data_CRUK2 <- CreateSeuratObject(counts = raw_countsCRUK2[,which(KeyCRUK2$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_CRUK2) <- KeyCRUK2$Primary.Annotation[which(KeyCRUK2$QC>0)]
marmoset_data_CRUK2$LOC <- KeyCRUK2$Loc[which(KeyCRUK2$QC>0)]
marmoset_data_CRUK2$Stage <- KeyCRUK2$Stage[which(KeyCRUK2$QC>0)]
marmoset_data_CRUK2$Batch <- "New"
marmoset_data_CRUK2 <- subset(marmoset_data_CRUK2, idents = c("EmDisc_CS5","EmDisc_CS6","Am_CS6","EmDisc_CS5","Am_CS6_EmDisc","EmDisc_CS6_Gast","PGC_CS6","EmDisc_CS6_Am","ExMes_CS6","SYS_CS6","VE_CS6","Tb_CS6","Stalk_CS6")) #,"Hyp_CS3","ICM_CS3","Tb_CS3","Epi_CS3"))
marmoset_data_CRUK2 <- NormalizeData(marmoset_data_CRUK2, verbose = FALSE)
marmoset_data_CRUK2Dataset <- "InVivo"

marmoset_data_CRUK3 <- CreateSeuratObject(counts = raw_countsCRUK3[,which(KeyCRUK3$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_CRUK3) <- KeyCRUK3$Primary.Annotation[which(KeyCRUK3$QC>0)]
marmoset_data_CRUK3$LOC <- KeyCRUK3$Loc[which(KeyCRUK3$QC>0)]
marmoset_data_CRUK3$Stage <- KeyCRUK3$Stage[which(KeyCRUK3$QC>0)]
marmoset_data_CRUK3$Batch <- "New"
marmoset_data_CRUK3 <- subset(marmoset_data_CRUK3, idents = c("ExMes_CS6","Am_CS6","Am_CS7","EmDisc_CS7","EmDisc_CS6","EmDisc_CS7_Gast","PGC_CS7","VE_CS7","VE_CS6","SYS_CS7","SYS_CS6","Stroma_CS6","Stalk_CS6")) #,"Stroma_CS6")) 
marmoset_data_CRUK3 <- NormalizeData(marmoset_data_CRUK3, verbose = FALSE)
marmoset_data_CRUK3$Dataset <- "InVivo"

marmoset_data_CS5 <- CreateSeuratObject(counts = raw_countsCS5[,which(KeyCS5$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_CS5) <- KeyCS5$Primary.Annotation[which(KeyCS5$QC>0)]
marmoset_data_CS5$LOC <- KeyCS5$Location[which(KeyCS5$QC>0)]
marmoset_data_CS5$Stage <- "CS5"
marmoset_data_CS5$Batch <- "Old"
marmoset_data_CS5 <- subset(marmoset_data_CS5, idents = c("Tb_CS5","ExMes_CS5","EmDisc_CS5","SYS_CS5","Am_CS5","VE_CS5","Am_CS5_PGC","EmDisc_CS5_Am","ReGland_CS5","Gland_CS5","ReStroma_CS5","Stroma_CS5")) #","Gland_CS5","ReGland_CS5","ReStroma_CS5")) 
marmoset_data_CS5 <- NormalizeData(marmoset_data_CS5, verbose = FALSE)
marmoset_data_CS5$Dataset <- "InVivo"

marmoset_data_CS6 <- CreateSeuratObject(counts = raw_countsCS6[,which(KeyCS6$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_CS6) <- KeyCS6$Primary.Anotation[which(KeyCS6$QC>0)]
marmoset_data_CS6$LOC <- KeyCS6$Location[which(KeyCS6$QC>0)]
marmoset_data_CS6$Stage <- "CS6"
marmoset_data_CS6$Batch <- "Old"
marmoset_data_CS6 <- subset(marmoset_data_CS6, idents = c("Tb_CS6","ExMes_CS6","Am_CS6_EmDisc","Am_CS6_PGC","EmDisc_CS6","SYS_CS6","Am_CS6","VE_CS6","PGC_CS6","EmDisc_CS6_Am","EmDisc_CS6_PGC","Stalk_CS6","EmDisc_gast_CS6","Gland_CS6","ReStroma_CS6","Stalk_CS6")) #,"Gland_CS6","Stroma_CS6","ReStroma_CS6"))
marmoset_data_CS6 <- NormalizeData(marmoset_data_CS6, verbose = FALSE)
marmoset_data_CS6$Dataset <- "InVivo"

marmoset_data_CS7 <- CreateSeuratObject(counts = raw_countsCS7[,which(KeyCS7$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(marmoset_data_CS7) <- KeyCS7$Primary.Annotation[which(KeyCS7$QC>0)]
marmoset_data_CS7$Stage <- "CS7"
marmoset_data_CS7$Batch <- "Old"
marmoset_data_CS7$LOC <- KeyCS7$Location[which(KeyCS7$QC>0)]
marmoset_data_CS7 <- subset(marmoset_data_CS7, idents = c("Tb_abembryonal_CS7","Tb_CS7","ExMes_CS7","EmDisc_CS7","SYS_CS7","Am_CS7","Stalk_CS7_PGC","EmDisc_CS7_PGC","ExMes_stalk_CS7","Stalk_CS7","Gland_CS7","ReStroma_CS7","Stroma_CS7","ReGland_CS7","Myo_CS7")) #"Gland_CS7","ReStroma_CS7","ReGland_CS7","Myo_CS7")) 
marmoset_data_CS7 <- NormalizeData(marmoset_data_CS7, verbose = FALSE)
marmoset_data_CS7$Dataset <- "InVivo"

marmoset_dataInVivo <- merge(marmoset_data_CS5, y = c(marmoset_data_Pre,marmoset_data_CS6,marmoset_data_CS7,marmoset_data_CRUK,marmoset_data_CRUK2,marmoset_data_CRUK3,marmoset_data_20307,marmoset_data_20308), project = "merged")
marmoset_dataInVivo$Dataset <- "2) Marmoset in vivo"

marmoset_dataInVivo2 <- marmoset_dataInVivo
marmoset_dataInVivo2 <- FindVariableFeatures(marmoset_dataInVivo2, selection.method = "vst", nfeatures = 20000)
marmoset_dataInVivo2 <- ScaleData(marmoset_dataInVivo2, verbose = FALSE)
marmoset_dataInVivo2 <- RunPCA(marmoset_dataInVivo2, npcs = 20, verbose = FALSE)
marmoset_dataInVivo2 <- RunUMAP(marmoset_dataInVivo2, reduction = "pca", dims = 1:20)
marmoset_dataInVivo2 <- RunTSNE(marmoset_dataInVivo2, reduction = "pca", dims = 1:20)
marmoset_dataInVivo2 <- FindNeighbors(marmoset_dataInVivo2, reduction = "pca", dims = 1:20)

marmoset_dataInVivo2 <- subset(marmoset_dataInVivo2, idents = c("Am_CS6_EmDisc","Am_CS7_EmDisc", "EmDisc_CS5_Am", "EmDisc_CS6_Am", "EmDisc_CS7_Am") , invert = TRUE)

newID <- as.character(Idents(marmoset_dataInVivo2))
newID[which(newID=="Tb_abembryonal_CS7")] <- "Tb_CS7"
newID[which(newID=="Am_CS5_PGC")] <- "PGC_CS5"
newID[which(newID=="EmDisc_CS5_Am")] <- "EmDisc_CS5"
newID[which(newID=="EmDisc_gast_CS6 ")] <- "EmDisc_CS6"
newID[which(newID=="EmDisc_Gast_CS6")] <- "EmDisc_CS6"
newID[which(newID=="EmDisc_stalk_CS6")] <- "Stalk_CS6"
newID[which(newID=="EmDisc_stalk_CS7")] <- "Stalk_CS7"
newID[which(newID=="EmDisc_CS6_PGC")] <- "PGC_CS6"
newID[which(newID=="Am_CS6_PGC")] <- "PGC_CS6"
newID[which(newID=="Stalk_CS7_PGC")] <- "PGC_CS7"
newID[which(newID=="EmDisc_CS7_PGC")] <- "PGC_CS7"
newID[which(newID=="ExMes_stalk_CS7")] <- "Stalk_CS7"
newID[which(newID=="EmDisc_CS7_Am")] <- "EmDisc_CS7"
newID[which(newID=="EmDisc_CS7_Am")] <- "EmDisc_CS7"
newID[which(newID=="Am_CS6_EmDisc")] <- "Am_CS6"
newID[which(newID=="EmDisc_CS6_Gast")] <- "EmDisc_CS6"
newID[which(newID=="EmDisc_gast_CS6")] <- "EmDisc_CS6"
newID[which(newID=="EmDisc_CS6_Am")] <- "EmDisc_CS6"
newID[which(newID=="EmDisc_Stalk_CS6")] <- "Stalk_CS6"
newID[which(newID=="EmDisc_Gast_CS7")] <- "EmDisc_CS7"
newID[which(newID=="EmDisc_CS7_Gast")] <- "EmDisc_CS7"
newID[which(newID=="EmDisc_CS6_Gast")] <- "EmDisc_CS6"
newID[which(newID=="Am_CS7_EmDisc")] <- "Am_CS7"
newID[which(newID=="cMor_CS3")] <- "cMor_CS2"

newID[which(newID=="ReStroma_CS5")] <- "Remodelled"
newID[which(newID=="ReStroma_CS6")] <- "Remodelled"
newID[which(newID=="ReStroma_CS7")] <- "Remodelled"
newID[which(newID=="Stroma_CS5")] <- "Stroma"
newID[which(newID=="Stroma_CS6")] <- "Stroma"
newID[which(newID=="Stroma_CS7")] <- "Stroma"
newID[which(newID=="Myo_CS7")] <- "Myo"


newID[which(newID=="ReGland_CS5")] <- "Remodelled"
newID[which(newID=="ReGland_CS6")] <- "Remodelled"
newID[which(newID=="ReGland_CS7")] <- "Remodelled"
newID[which(newID=="Gland_CS5")] <- "Gland"
newID[which(newID=="Gland_CS6")] <- "Gland"
newID[which(newID=="Gland_CS7")] <- "Gland"
newID[which(newID=="ReGland_CS5")] <- "Remodelled"
newID[which(newID=="ReGland_CS6")] <- "Remodelled"
newID[which(newID=="ReGland_CS7")] <- "Remodelled"

Idents(marmoset_dataInVivo2) <- newID

Ano <- read_excel("/Users/christopherpenfold/Desktop/Temp/Matteo_Anotations.xlsx", sheet = "Combined")

saveext = "~/Desktop/Data/Endometrial/InVitro/Matteo/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

#Get lists of markers
PrimedNaiveForm <- read_excel("/Users/christopherpenfold/Downloads/mmc3.xls")
PrimedNaive <- read_excel("/Users/christopherpenfold/Downloads/mmc2.xls")
AmnionMarkers <- read_excel("/Users/christopherpenfold/Downloads/AmnionMarkers.xlsx")
TbMarkers <- read_excel("/Users/christopherpenfold/Downloads/AnotatedAeDEupdated2.xlsx",sheet = 2) 


#Get the dataset
D <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/CorrectedFinalAno2.rds")
Idents(D) <- D$CorrectLabel
D$Cells <- Idents(D)

#
cl2 <- c("TCTCCGAGTAGGCTCC-1_2","CACTTCGCAGCTGCCA-1_2","AATCGACTCCAGCAAT-1_2","GCCCAGATCGTTGCCT-1_2","GTGAGCCAGCACGTCC-1_2")
Idents(D,cells =cl2) <- "putSTB"
cl1 <- c("AACTTCTGTAGCCCTG-1_3","ACTATTCGTCTGTGAT-1_3","ACTTAGGCAGTAACGG-1_3","AGATAGATCTCCCTAG-1_3","AGGATCTAGCGTATAA-1_3","AGGTCATCAGAAGTTA-1_3","ATACTTCTCATGCCAA-1_3","ATCCACCTCATTATCC-1_3","ATTCATCGTGTCTAAC-1_3","CAAGGGACATCTGGGC-1_3","CCACCATAGTGAGTTA-1_3","CGAGAAGCATTGACAC-1_3","CTGGTCTGTATGAGCG-1_3","GATTGGTGTCCAGGTC-1_3","GCCAGGTAGCTTGTGT-1_3","GCCTGTTTCATGGCCG-1_3","GCTACAACAGTAGAAT-1_3","GGAAGTGCATGACAGG-1_3","GGAATCTGTCCTGTCT-1_3","TACTGCCTCGACTCCT-1_3","TGGAACTCAGGTGTTT-1_3","TTGGGATCATGGATCT-1_3","TTTGACTAGGGTCTTT-1_3")
cl3 <- c("AACACACAGGGCAATC-1_3","AAGTTCGTCCACGTGG-1_3","ACACAGTTCGATACTG-1_3","ACACAGTTCTGCTTAT-1_3","ACACTGATCCATACTT-1_3","ACCAACATCATCACTT-1_3","ACCCTTGAGGTTTACC-1_3","ACGATGTCACAGCCAC-1_3","ACGTACACAGTTGCGC-1_3","ACGTTCCCAATCTGCA-1_3","ACTATCTTCGTTGTTT-1_3","ACTATTCAGTCTACCA-1_3","ACTTATCAGGTAAGGA-1_3","ACTTCGCGTCTAGATC-1_3","AGACCCGAGAGCAACC-1_3","AGCTACATCTCGTGAA-1_3","AGCTTCCTCGCAGTTA-1_3","AGTACCACAACGGCTC-1_3","ATAGGCTCACACGGAA-1_3","ATATCCTAGACTCAAA-1_3","ATCCTATAGACAGTCG-1_3","ATCGGATAGAGGGTAA-1_3","ATCGTGACAGACCAGA-1_3","ATCTCTAGTACCTATG-1_3","ATCTTCACACTTGAGT-1_3","ATGCATGCAGGGTCTC-1_3","ATTACTCCAATACCTG-1_3","ATTCTTGTCAGACCCG-1_3","ATTGTTCCAGCGTGCT-1_3","ATTTCACTCACACCGG-1_3","CACAACAAGTAAACAC-1_3","CATACTTAGGTCGTGA-1_3","CATGCCTCAGGGATAC-1_3","CATGCGGGTGTCTTCC-1_3","CATTGCCTCATCAGTG-1_3","CCCTGATAGCTCTGTA-1_3","CCGGTGAAGCAATTAG-1_3","CCGTAGGGTAAGTTAG-1_3","CCTACGTAGGACTTCT-1_3","CCTATCGAGCACTTTG-1_3","CCTATCGAGTTGCCCG-1_3","CCTCAACTCTTGCAAG-1_3","CCTCTCCTCGTTGTGA-1_3","CCTTTGGAGGTTGTTC-1_3","CGAAGGATCACCACAA-1_3","CGAAGTTAGTGCTCAT-1_3","CGAGGCTAGGACAACC-1_3","CGATGCGCATCATTGG-1_3","CGCAGGTCAGAAACCG-1_3","CGCAGGTTCGGCCTTT-1_3","CGGGTGTAGTTGTAGA-1_3","CGTAGTACAATACCTG-1_3","CTAACCCAGATAGTGT-1_3","CTAACCCAGGCACCAA-1_3","CTCAGGGCACGATTCA-1_3","CTCATCGCATCAACCA-1_3","CTCATGCCACACACGC-1_3","CTCCAACGTTCTCCAC-1_3","CTCCCTCCAGCACAAG-1_3","CTTCAATTCCTCTCTT-1_3","GAAATGAGTCTTGCGG-1_3","GACCCTTTCTAACACG-1_3","GACGCTGTCGAGCCAC-1_3","GAGACCCCACTGCGAC-1_3","GAGACTTCATCAGCGC-1_3","GAGTTACAGCAGCCCT-1_3","GATCACAGTCACTTCC-1_3","GATCATGAGTATCTGC-1_3","GCACATAGTCATCCCT-1_3","GCACGGTCAAAGCGTG-1_3","GCCCGAAGTTGTTTGG-1_3","GGACGTCCAGCCCACA-1_3","GGAGGTATCTAGACCA-1_3","GGCGTCACAATGTTGC-1_3","GTAACCACATTCTGTT-1_3","GTAGCTACAACAGAGC-1_3","GTAGGTTGTCCTCATC-1_3","GTCTACCCAACACGTT-1_3","GTGATGTCACCCTCTA-1_3","GTGGTTATCCGGACTG-1_3","GTGTAACCAAGCAATA-1_3","GTGTGATGTCATATGC-1_3","TAATTCCCACACCAGC-1_3","TACCTGCCATGCCGCA-1_3","TACCTGCGTGTTTCTT-1_3","TACGGTATCATGCAGT-1_3","TACTTACGTACCGGCT-1_3","TACTTCACATGAATCC-1_3","TAGACCAAGGAGACCT-1_3","TCAAGACCAATCTCTT-1_3","TCACAAGTCACGAGGA-1_3","TCACACCCAGGTTACT-1_3","TCACGCTTCAGGAAAT-1_3","TCATACTAGTATAACG-1_3","TCCTCGATCTGCGAGC-1_3","TCGAACACAAGAGCTG-1_3","TCGATTTCAATTTCCT-1_3","TCGCACTCAAATAGCA-1_3","TCGGGTGAGATTCGCT-1_3","TCGGTCTTCCACCTCA-1_3","TCGTGGGCAACGGCCT-1_3","TCTGGCTGTACTCAAC-1_3","TGAACGTAGTTTGAGA-1_3","TGACTCCTCAGTGCGC-1_3","TGAGTCATCATCTACT-1_3","TGATCTTCACACCTGG-1_3","TGATCTTCACACCTGG-1_3","TGCAGGCAGGAGGCAG-1_3","TGCCGAGCATCATTTC-1_3","TGGAGAGAGTATAACG-1_3","TGGTGATCAATCTCGA-1_3","TGTACAGGTTAAGTCC-1_3","TGTCCCAAGTGAATAC-1_3","TTCATGTCACCTTCGT-1_3","TTCTAGTGTAGGAAAG-1_3","TTCTTGACACGAAGAC-1_3","TTGAACGTCTGGCCAG-1_3","TTGACCCAGCTGTCCG-1_3","TTGACCCTCCTTGGAA-1_3","TTTACCATCATCGGGC-1_3","TTTATGCTCCGGACGT-1_3")
Idents(D,cells =cl1) <- "putSTB"
Idents(D,cells =cl3) <- "putSTB"
cl1 <- c("ACATCGACAGAGGGTT-1_7","ATATCCTAGAGAGCAA-1_7","ATCCCTGAGTGCACCC-1_7","CCTAAGATCCTTCGAC-1_7","GAGGGTACACTCATAG-1_7","GTAGTACCAGCGTTTA-1_7","TGCAGGCGTAGCGAGT-1_7","TTAGGGTGTACCATAC-1_7","AGGGTGATCCGAGATT-1_7")
Idents(D,cells =cl1) <- "putSTB"


newcl1 <- c("AGATAGATCTCCCTAG-1_3","ACTATTCGTCTGTGAT-1_3","GGAATCTGTCCTGTCT-1_3","AACTTCTGTAGCCCTG-1_3","TGGAACTCAGGTGTTT-1_3","CCACCATAGTGAGTTA-1_3","ATCCACCTCATTATCC-1_3","AGGTCATCAGAAGTTA-1_3","TACTGCCTCGACTCCT-1_3","TTGGGATCATGGATCT-1_3","GATTGGTGTCCAGGTC-1_3","TTTGACTAGGGTCTTT-1_3","ACTTAGGCAGTAACGG-1_3","GCCTGTTTCATGGCCG-1_3","CTGGTCTGTATGAGCG-1_3","GCTACAACAGTAGAAT-1_3","CAAGGGACATCTGGGC-1_3","GGAAGTGCATGACAGG-1_3","AGGATCTAGCGTATAA-1_3","CATTGCCGTAGCTGCC-1_3","GCCAGGTAGCTTGTGT-1_3","ATTCATCGTGTCTAAC-1_3")
newcl1 <- str_replace(newcl1,"-","-")

newcl2 <- c("TGTACAGGTTAAGTCC-1_3","GTGTGATGTCATATGC-1_3","TTCTAGTGTAGGAAAG-1_3","CGTAGTACAATACCTG-1_3","TGTCCCAAGTGAATAC-1_3","ATATCCTAGACTCAAA-1_3","ACTTATCAGGTAAGGA-1_3","ATTGTTCCAGCGTGCT-1_3","CCTATCGAGTTGCCCG-1_3","ATGCATGCAGGGTCTC-1_3","TGCAGGCAGGAGGCAG-1_3","ACACAGTTCTGCTTAT-1_3","CCCTGATAGCTCTGTA-1_3","ACGATGTCACAGCCAC-1_3","ACGTACACAGTTGCGC-1_3","GTGTAACCAAGCAATA-1_3","GTAACCACATTCTGTT-1_3","AGTACCACAACGGCTC-1_3","CATTGCCTCATCAGTG-1_3","ATTTCACTCACACCGG-1_3","ATCGGATAGAGGGTAA-1_3","GTGATGTCACCCTCTA-1_3","GAGACTTCATCAGCGC-1_3","CTTCAATTCCTCTCTT-1_3","ACCAACATCATCACTT-1_3","GACCCTTTCTAACACG-1_3","CCTCAACTCTTGCAAG-1_3","GCACGGTCAAAGCGTG-1_3","ATTACTCCAATACCTG-1_3","TCTGGCTGTACTCAAC-1_3","TCAAGACCAATCTCTT-1_3","CGAAGTTAGTGCTCAT-1_3","TGGAGAGAGTATAACG-1_3","GATCATGAGTATCTGC-1_3","TTGACCCTCCTTGGAA-1_3","TCGGTCTTCCACCTCA-1_3","GCCCGAAGTTGTTTGG-1_3","ATCCTATAGACAGTCG-1_3","GACGCTGTCGAGCCAC-1_3","GATCACAGTCACTTCC-1_3","CTCATGCCACACACGC-1_3","TACCTGCCATGCCGCA-1_3","CATACTTAGGTCGTGA-1_3","AGACCCGAGAGCAACC-1_3","TGAACGTAGTTTGAGA-1_3","ATCGTGACAGACCAGA-1_3","CATGCGGGTGTCTTCC-1_3","TCGCACTCAAATAGCA-1_3","GGAGGTATCTAGACCA-1_3","TCGATTTCAATTTCCT-1_3","CTCCAACGTTCTCCAC-1_3","CGCAGGTCAGAAACCG-1_3","TTGAACGTCTGGCCAG-1_3","CCTACGTAGGACTTCT-1_3","TCACACCCAGGTTACT-1_3","AGCTACATCTCGTGAA-1_3","GTAGCTACAACAGAGC-1_3","CGGGTGTAGTTGTAGA-1_3","ACTTCGCGTCTAGATC-1_3","TGGTGATCAATCTCGA-1_3","CGAGGCTAGGACAACC-1_3","GGCGTCACAATGTTGC-1_3","CCGGTGAAGCAATTAG-1_3","CCGTAGGGTAAGTTAG-1_3","AAGTTCGTCCACGTGG-1_3","GTCTACCCAACACGTT-1_3","ATAGGCTCACACGGAA-1_3","GCACATAGTCATCCCT-1_3","AGCTTCCTCGCAGTTA-1_3","CTCCCTCCAGCACAAG-1_3","ACGTTCCCAATCTGCA-1_3","ACCCTTGAGGTTTACC-1_3","TCGGGTGAGATTCGCT-1_3","TACTTACGTACCGGCT-1_3","CTCATCGCATCAACCA-1_3","TCCTCGATCTGCGAGC-1_3","TACCTGCGTGTTTCTT-1_3","AACACACAGGGCAATC-1_3","TCGAACACAAGAGCTG-1_3","CTAACCCAGATAGTGT-1_3","GAAATGAGTCTTGCGG-1_3","TACGGTATCATGCAGT-1_3","GTGGTTATCCGGACTG-1_3","TAGACCAAGGAGACCT-1_3")
newcl2 <- str_replace(newcl2,"-","-")
Idents(D,cells =newcl1) <- "putSTB"
Idents(D,cells =newcl2) <- "putSTB"


newcl_1 <- c("AGATAGATCTCCCTAG-1_3","ACTATTCGTCTGTGAT-1_3","GGAATCTGTCCTGTCT-1_3","AACTTCTGTAGCCCTG-1_3","TGGAACTCAGGTGTTT-1_3","CCACCATAGTGAGTTA-1_3","ATCCACCTCATTATCC-1_3","AGGTCATCAGAAGTTA-1_3","TACTGCCTCGACTCCT-1_3","TTGGGATCATGGATCT-1_3","GATTGGTGTCCAGGTC-1_3","TTTGACTAGGGTCTTT-1_3","ACTTAGGCAGTAACGG-1_3","GCCTGTTTCATGGCCG-1_3","CTGGTCTGTATGAGCG-1_3","GCTACAACAGTAGAAT-1_3","CAAGGGACATCTGGGC-1_3","GGAAGTGCATGACAGG-1_3","AGGATCTAGCGTATAA-1_3","CATTGCCGTAGCTGCC-1_3","GCCAGGTAGCTTGTGT-1_3","ATTCATCGTGTCTAAC-1_3")
newcl_2 <- c("GTGTGATGTCATATGC-1_3","TTCTAGTGTAGGAAAG-1_3","CGTAGTACAATACCTG-1_3","TGTCCCAAGTGAATAC-1_3","ATATCCTAGACTCAAA-1_3","ACTTATCAGGTAAGGA-1_3","ATTGTTCCAGCGTGCT-1_3","CCTATCGAGTTGCCCG-1_3","ATGCATGCAGGGTCTC-1_3","TGCAGGCAGGAGGCAG-1_3","ACACAGTTCTGCTTAT-1_3","CCCTGATAGCTCTGTA-1_3","ACGATGTCACAGCCAC-1_3","ACGTACACAGTTGCGC-1_3","GTGTAACCAAGCAATA-1_3","GTAACCACATTCTGTT-1_3","AGTACCACAACGGCTC-1_3","CATTGCCTCATCAGTG-1_3","ATTTCACTCACACCGG-1_3","ATCGGATAGAGGGTAA-1_3","GTGATGTCACCCTCTA-1_3","GAGACTTCATCAGCGC-1_3","CTTCAATTCCTCTCTT-1_3","ACCAACATCATCACTT-1_3","GACCCTTTCTAACACG-1_3","CCTCAACTCTTGCAAG-1_3","GCACGGTCAAAGCGTG-1_3","ATTACTCCAATACCTG-1_3","TCTGGCTGTACTCAAC-1_3","TCAAGACCAATCTCTT-1_3","CGAAGTTAGTGCTCAT-1_3","TGGAGAGAGTATAACG-1_3","GATCATGAGTATCTGC-1_3","TTGACCCTCCTTGGAA-1_3","TCGGTCTTCCACCTCA-1_3","GCCCGAAGTTGTTTGG-1_3","ATCCTATAGACAGTCG-1_3","GACGCTGTCGAGCCAC-1_3","GATCACAGTCACTTCC-1_3","CTCATGCCACACACGC-1_3","TACCTGCCATGCCGCA-1_3","CATACTTAGGTCGTGA-1_3","AGACCCGAGAGCAACC-1_3","TGAACGTAGTTTGAGA-1_3","ATCGTGACAGACCAGA-1_3","CATGCGGGTGTCTTCC-1_3","TCGCACTCAAATAGCA-1_3","GGAGGTATCTAGACCA-1_3","TCGATTTCAATTTCCT-1_3","CTCCAACGTTCTCCAC-1_3","CGCAGGTCAGAAACCG-1_3","TTGAACGTCTGGCCAG-1_3","CCTACGTAGGACTTCT-1_3","TCACACCCAGGTTACT-1_3","AGCTACATCTCGTGAA-1_3","GTAGCTACAACAGAGC-1_3","CGGGTGTAGTTGTAGA-1_3","ACTTCGCGTCTAGATC-1_3","TGGTGATCAATCTCGA-1_3","CGAGGCTAGGACAACC-1_3","GGCGTCACAATGTTGC-1_3","CCGGTGAAGCAATTAG-1_3","CCGTAGGGTAAGTTAG-1_3","AAGTTCGTCCACGTGG-1_3","GTCTACCCAACACGTT-1_3","ATAGGCTCACACGGAA-1_3","GCACATAGTCATCCCT-1_3","AGCTTCCTCGCAGTTA-1_3","CTCCCTCCAGCACAAG-1_3","ACGTTCCCAATCTGCA-1_3","ACCCTTGAGGTTTACC-1_3","TCGGGTGAGATTCGCT-1_3","TACTTACGTACCGGCT-1_3","CTCATCGCATCAACCA-1_3","TCCTCGATCTGCGAGC-1_3","TACCTGCGTGTTTCTT-1_3","AACACACAGGGCAATC-1_3","TCGAACACAAGAGCTG-1_3","CTAACCCAGATAGTGT-1_3","GAAATGAGTCTTGCGG-1_3","TACGGTATCATGCAGT-1_3","TACGGTATCATGCAGT-1_3","GTGGTTATCCGGACTG-1_3","TGTACAGGTTAAGTCC-1_3")
newcl_1 <- str_replace(newcl_1,"-","-")
newcl_2 <- str_replace(newcl_2,"-","-")
Idents(D,cells =newcl_1) <- "putSTB"
Idents(D,cells =newcl_2) <- "putSTB"


D$CorrectLabel <- Idents(D)


#Get the refined gene anotations
NewG1 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D1_emb.tsv", header=TRUE)
NewG2 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D2_emb.tsv", header=TRUE)
NewG3B <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D3B_emb.tsv", header=TRUE)
NewG6 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D6_emb.tsv", header=TRUE)
NewG8 <- read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/Genotypes/D8_emb.tsv", header=TRUE)

D$Cells <- D$CorrectLabel

#Subset on embryonic only and anotate
Dsubset1 <- subset(D,idents=c("ExMes_d14","STB_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","CTB_d14","putSTB","EVT_d14","Hyp_d14","Am/EmDisc","Hyp/Am","putExMes"))

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
Idents(Dsubset1, cells= paste(NewG6$barcode[which(NewG6$assignment=="0/1")],6,sep="_") ) <- "NC6_1"
Idents(Dsubset1, cells= paste(NewG6$barcode[which(NewG6$assignment=="1/0")],6,sep="_") ) <- "NC6_2"

Idents(Dsubset1, cells= paste(NewG8$barcode[which(NewG8$assignment=="0")],8,sep="_") ) <- "NC8_1"
Idents(Dsubset1, cells= paste(NewG8$barcode[which(NewG8$assignment=="1")],8,sep="_") ) <- "NC8_2"
Idents(Dsubset1, cells= paste(NewG8$barcode[which(NewG8$assignment=="0/1")],8,sep="_") ) <- "NC8_1"
Idents(Dsubset1, cells= paste(NewG8$barcode[which(NewG8$assignment=="1/0")],8,sep="_") ) <- "NC8_2"

Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0")],5,sep="_") ) <- "NC3_1"
Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/2")],5,sep="_") ) <- "NC3_1"
Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/1")],5,sep="_") ) <- "NC3_1"
Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1")],5,sep="_") ) <- "NC3_2"
Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1/2")],5,sep="_") ) <- "NC3_2"
Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1/0")],5,sep="_") ) <- "NC3_2"
Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="2")],5,sep="_") ) <- "NC3_3"
Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="2/0")],5,sep="_") ) <- "NC3_3"
Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="0/1")],5,sep="_") ) <- "NC2_1"
Idents(Dsubset1, cells= paste(NewG3B$barcode[which(NewG3B$assignment=="1/0")],5,sep="_") ) <- "NC2_2"

Dsubset1$Genotype2 <- Idents(Dsubset1)
Dsubset1$Cells <- Dsubset1$CorrectLabel
Idents(Dsubset1) <- Dsubset1$Cells

#Plot the UMAP/PCAs
p<-DimPlot(Dsubset1,  pt.size = 4, reduction = "umap", split.by = "Genotype2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Embryo_split.pdf",sep=""),width = 80, height = 8,p,limitsize = FALSE)
p<-DimPlot(Dsubset1,  pt.size = 4, reduction = "pca", split.by = "Genotype2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Embryo_split.pdf",sep=""),width = 80, height = 8,p,limitsize = FALSE)


Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="STB")]) <- "STB"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="Am")]) <- "Am"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="CTB")]) <- "CTB"
#Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="EmDisc/PGC")]) <- "EmDisc/PGC"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="ExMes")]) <- "ExMes"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="Hypoblast")]) <- "Hypoblast"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="STB1")]) <- "STB1"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="SYS")]) <- "SYS"
#Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="Am/PGC")]) <- "Am/PGC"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="EmDisc")]) <- "EmDisc"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="Am/EmDisc")]) <- "Am/EmDisc"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="EVT")]) <- "EVT"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="ExMes/SYS")]) <- "ExMes/SYS"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="pPGC")]) <- "PGC"







#Dsubset2 <- subset(mammal.combined,idents=c("Am_CS5","Am_CS6","Am_CS7","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am","Am/EmDisc","EmDisc","Hypoblast","VE_CS5","VE_CS6","VE_CS7"))
DefaultAssay(Dsubset1) <- "RNA"
DE1 <- FindMarkers(marmoset_dataInVivo2,ident.1 = c("EmDisc_CS5","Am_CS5","Am_CS7","EmDisc_CS7","Am_CS6","EmDisc_CS6"),ident.2 = c("Tb_CS5","Tb_CS6","Tb_CS7"), test.use = "MAST")
DE2 <- FindMarkers(Dsubset1,ident.1 = c("Am","EmDisc","Am/EmDisc"),ident.2 = c("STB","STB1"), test.use = "MAST", logfc.threshold = log(0))

Idents(Dsubset1,cells=WhichCells(Dsubset1,idents=c("EmDisc","Am","Am/EmDisc"))) <- "Am"
Idents(Dsubset1,cells=WhichCells(Dsubset1,idents=c("STB","STB1"))) <- "STB"

AvE <- AverageExpression(Dsubset1)
AvExp <- as.data.frame(AvE$RNA)
AvExp$'AvExp' <- NA
AvExp[,'AvExp'] <- (AvExp$STB+AvExp$Am)/2
AvExp$'log2FC' <- NA
AvExp[rownames(DE2),'log2FC'] <- DE2$avg_log2FC
AvExp[rownames(DE2),'pval'] <- DE2$p_val_adj
AvExp$'IsIn' <- 0
AvExp[which(AvExp$pval<0.05),'IsIn'] <- 1
AvExp$'IsIn' <- as.factor(AvExp$'IsIn')
genes.to.label1 = intersect(rownames(DE1),rownames(DE2))
p1 <- ggplot(AvExp, aes(log2(AvExp),log2FC)) + geom_point(aes(color=IsIn)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Average expression", y = "log2 FC")
ggsave(filename=paste(saveext,"AmEmDiscvsTb_marmosethuman.pdf",sep=""),width = 13, height = 13, plot = p1)

genes.to.label1 = c("GATA2","GATA3","CGA","TFAP2C")
p1 <- ggplot(AvExp, aes(log2(AvExp),log2FC)) + geom_point(aes(color=IsIn)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Average expression", y = "log2 FC")
ggsave(filename=paste(saveext,"AmEmDiscvsTb_marmosethuman2.pdf",sep=""),width = 13, height = 13, plot = p1)




write.table(as.data.frame(DE1),file=paste("DEMarm_Tb_vs_AmEmDisc.csv",sep=""))
write.table(as.data.frame(DE2),file=paste("DEOurs_Tb_vs_AmEmDisc.csv",sep=""))
#write.table(as.data.frame(DE2B),file=paste("DETyser_SYS_vs_VE.csv",sep=""))





#Dsubset2 <- subset(mammal.combined,idents=c("Am_CS5","Am_CS6","Am_CS7","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am","Am/EmDisc","EmDisc","Hypoblast","VE_CS5","VE_CS6","VE_CS7"))
DefaultAssay(Dsubset1) <- "RNA"
DE1 <- FindMarkers(marmoset_dataInVivo2,ident.1 = c("SYS_CS5","SYS_CS5","SYS_CS7"),ident.2 = c("Tb_CS5","Tb_CS6","Tb_CS7"), test.use = "MAST")
DE2 <- FindMarkers(Dsubset1,ident.1 = c("SYS"),ident.2 = c("STB","STB1"), test.use = "MAST", logfc.threshold = log(0))


#Idents(Dsubset1,cells=WhichCells(Dsubset1,idents=c("EmDisc","Am","Am/EmDisc"))) <- "Am"
Idents(Dsubset1,cells=WhichCells(Dsubset1,idents=c("STB","STB1"))) <- "STB"

AvE <- AverageExpression(Dsubset1)
AvExp <- as.data.frame(AvE$RNA)
AvExp$'AvExp' <- NA
AvExp[,'AvExp'] <- (AvExp$STB+AvExp$SYS)/2
AvExp$'log2FC' <- NA
AvExp[rownames(DE2),'log2FC'] <- DE2$avg_log2FC
AvExp[rownames(DE2),'pval'] <- DE2$p_val_adj
AvExp$'IsIn' <- 0
AvExp[which(AvExp$pval<0.05),'IsIn'] <- 1
AvExp$'IsIn' <- as.factor(AvExp$'IsIn')
genes.to.label1 = intersect(rownames(DE1),rownames(DE2))
p1 <- ggplot(AvExp, aes(log2(AvExp),log2FC)) + geom_point(aes(color=IsIn)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Average expression", y = "log2 FC")
ggsave(filename=paste(saveext,"SYSvsTb_marmosethuman.pdf",sep=""),width = 13, height = 13, plot = p1)

genes.to.label1 = c("GATA2","GATA3","CGA","TFAP2C")
p1 <- ggplot(AvExp, aes(log2(AvExp),log2FC)) + geom_point(aes(color=IsIn)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Average expression", y = "log2 FC")
ggsave(filename=paste(saveext,"SYSvsTb_marmosethuman2.pdf",sep=""),width = 13, height = 13, plot = p1)


write.table(as.data.frame(DE1),file=paste("DEMarm_Tb_vs_SYS.csv",sep=""))
write.table(as.data.frame(DE2),file=paste("DEOurs_Tb_vs_SYS.csv",sep=""))
