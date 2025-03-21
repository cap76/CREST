library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(pheatmap)
library(stringr)
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

TomMarkers <- c("PLA2G2A",
"DIO2",
"SCARA5",
"CXCL8",
"CXCL14",
"TIMP3",
"IL15",
"ECM1",
"AXL",
"WNT7A",
"PTGS1",
"HEY1",
"DNAI1",
"FOXJ1",
"SPP1",
"PAEP",
"DPP4",
"LIF",
"AREG",
"EREG",
"NEAT1",
"KCNQ1OT1",
"WNT5A")

#Load in the dataset
#D <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/CorrectedFinalAno.rds")


#Get the dataset
D <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/CorrectedFinalAno2.rds")
Idents(D) <- D$CorrectLabel


D$LongID <- paste(D$ID3,D$Genotype,sep="_")
p<-FeaturePlot(D,  pt.size = 1, features = "GATA3", reduction = "pca" , split.by = "LongID") 
ggsave(filename=paste(saveext,"/DimRed/PCA_GATA3_longID.pdf",sep=""),width = 200, height = 8,p,limitsize = FALSE)
p<-FeaturePlot(D,  pt.size = 1, features = "PADI1", reduction = "umap") 
ggsave(filename=paste(saveext,"/DimRed/UMAP_AllDatasets_nosplit_PADI1.pdf",sep=""),width = 10, height = 8,p)


p<-FeaturePlot(D,  pt.size = 1, features = "POU5F1", reduction = "pca" , split.by = "LongID") 
ggsave(filename=paste(saveext,"/DimRed/PCA_OCT4_longID.pdf",sep=""),width = 200, height = 8,p,limitsize = FALSE)

Idents(D) <- D$Cl05
p<-DimPlot(D,  pt.size = 4, reduction = "pca", split.by = "LongID", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Cl_longID.pdf",sep=""),width = 200, height = 8,p,limitsize = FALSE)


cl2 <- c("TCTCCGAGTAGGCTCC-1_2","CACTTCGCAGCTGCCA-1_2","AATCGACTCCAGCAAT-1_2","GCCCAGATCGTTGCCT-1_2","GTGAGCCAGCACGTCC-1_2")
Idents(D,cells =cl2) <- "putSTB"
cl1 <- c("AACTTCTGTAGCCCTG-1_3","ACTATTCGTCTGTGAT-1_3","ACTTAGGCAGTAACGG-1_3","AGATAGATCTCCCTAG-1_3","AGGATCTAGCGTATAA-1_3","AGGTCATCAGAAGTTA-1_3","ATACTTCTCATGCCAA-1_3","ATCCACCTCATTATCC-1_3","ATTCATCGTGTCTAAC-1_3","CAAGGGACATCTGGGC-1_3","CCACCATAGTGAGTTA-1_3","CGAGAAGCATTGACAC-1_3","CTGGTCTGTATGAGCG-1_3","GATTGGTGTCCAGGTC-1_3","GCCAGGTAGCTTGTGT-1_3","GCCTGTTTCATGGCCG-1_3","GCTACAACAGTAGAAT-1_3","GGAAGTGCATGACAGG-1_3","GGAATCTGTCCTGTCT-1_3","TACTGCCTCGACTCCT-1_3","TGGAACTCAGGTGTTT-1_3","TTGGGATCATGGATCT-1_3","TTTGACTAGGGTCTTT-1_3")
cl3 <- c("AACACACAGGGCAATC-1_3","AAGTTCGTCCACGTGG-1_3","ACACAGTTCGATACTG-1_3","ACACAGTTCTGCTTAT-1_3","ACACTGATCCATACTT-1_3","ACCAACATCATCACTT-1_3","ACCCTTGAGGTTTACC-1_3","ACGATGTCACAGCCAC-1_3","ACGTACACAGTTGCGC-1_3","ACGTTCCCAATCTGCA-1_3","ACTATCTTCGTTGTTT-1_3","ACTATTCAGTCTACCA-1_3","ACTTATCAGGTAAGGA-1_3","ACTTCGCGTCTAGATC-1_3","AGACCCGAGAGCAACC-1_3","AGCTACATCTCGTGAA-1_3","AGCTTCCTCGCAGTTA-1_3","AGTACCACAACGGCTC-1_3","ATAGGCTCACACGGAA-1_3","ATATCCTAGACTCAAA-1_3","ATCCTATAGACAGTCG-1_3","ATCGGATAGAGGGTAA-1_3","ATCGTGACAGACCAGA-1_3","ATCTCTAGTACCTATG-1_3","ATCTTCACACTTGAGT-1_3","ATGCATGCAGGGTCTC-1_3","ATTACTCCAATACCTG-1_3","ATTCTTGTCAGACCCG-1_3","ATTGTTCCAGCGTGCT-1_3","ATTTCACTCACACCGG-1_3","CACAACAAGTAAACAC-1_3","CATACTTAGGTCGTGA-1_3","CATGCCTCAGGGATAC-1_3","CATGCGGGTGTCTTCC-1_3","CATTGCCTCATCAGTG-1_3","CCCTGATAGCTCTGTA-1_3","CCGGTGAAGCAATTAG-1_3","CCGTAGGGTAAGTTAG-1_3","CCTACGTAGGACTTCT-1_3","CCTATCGAGCACTTTG-1_3","CCTATCGAGTTGCCCG-1_3","CCTCAACTCTTGCAAG-1_3","CCTCTCCTCGTTGTGA-1_3","CCTTTGGAGGTTGTTC-1_3","CGAAGGATCACCACAA-1_3","CGAAGTTAGTGCTCAT-1_3","CGAGGCTAGGACAACC-1_3","CGATGCGCATCATTGG-1_3","CGCAGGTCAGAAACCG-1_3","CGCAGGTTCGGCCTTT-1_3","CGGGTGTAGTTGTAGA-1_3","CGTAGTACAATACCTG-1_3","CTAACCCAGATAGTGT-1_3","CTAACCCAGGCACCAA-1_3","CTCAGGGCACGATTCA-1_3","CTCATCGCATCAACCA-1_3","CTCATGCCACACACGC-1_3","CTCCAACGTTCTCCAC-1_3","CTCCCTCCAGCACAAG-1_3","CTTCAATTCCTCTCTT-1_3","GAAATGAGTCTTGCGG-1_3","GACCCTTTCTAACACG-1_3","GACGCTGTCGAGCCAC-1_3","GAGACCCCACTGCGAC-1_3","GAGACTTCATCAGCGC-1_3","GAGTTACAGCAGCCCT-1_3","GATCACAGTCACTTCC-1_3","GATCATGAGTATCTGC-1_3","GCACATAGTCATCCCT-1_3","GCACGGTCAAAGCGTG-1_3","GCCCGAAGTTGTTTGG-1_3","GGACGTCCAGCCCACA-1_3","GGAGGTATCTAGACCA-1_3","GGCGTCACAATGTTGC-1_3","GTAACCACATTCTGTT-1_3","GTAGCTACAACAGAGC-1_3","GTAGGTTGTCCTCATC-1_3","GTCTACCCAACACGTT-1_3","GTGATGTCACCCTCTA-1_3","GTGGTTATCCGGACTG-1_3","GTGTAACCAAGCAATA-1_3","GTGTGATGTCATATGC-1_3","TAATTCCCACACCAGC-1_3","TACCTGCCATGCCGCA-1_3","TACCTGCGTGTTTCTT-1_3","TACGGTATCATGCAGT-1_3","TACTTACGTACCGGCT-1_3","TACTTCACATGAATCC-1_3","TAGACCAAGGAGACCT-1_3","TCAAGACCAATCTCTT-1_3","TCACAAGTCACGAGGA-1_3","TCACACCCAGGTTACT-1_3","TCACGCTTCAGGAAAT-1_3","TCATACTAGTATAACG-1_3","TCCTCGATCTGCGAGC-1_3","TCGAACACAAGAGCTG-1_3","TCGATTTCAATTTCCT-1_3","TCGCACTCAAATAGCA-1_3","TCGGGTGAGATTCGCT-1_3","TCGGTCTTCCACCTCA-1_3","TCGTGGGCAACGGCCT-1_3","TCTGGCTGTACTCAAC-1_3","TGAACGTAGTTTGAGA-1_3","TGACTCCTCAGTGCGC-1_3","TGAGTCATCATCTACT-1_3","TGATCTTCACACCTGG-1_3","TGATCTTCACACCTGG-1_3","TGCAGGCAGGAGGCAG-1_3","TGCCGAGCATCATTTC-1_3","TGGAGAGAGTATAACG-1_3","TGGTGATCAATCTCGA-1_3","TGTACAGGTTAAGTCC-1_3","TGTCCCAAGTGAATAC-1_3","TTCATGTCACCTTCGT-1_3","TTCTAGTGTAGGAAAG-1_3","TTCTTGACACGAAGAC-1_3","TTGAACGTCTGGCCAG-1_3","TTGACCCAGCTGTCCG-1_3","TTGACCCTCCTTGGAA-1_3","TTTACCATCATCGGGC-1_3","TTTATGCTCCGGACGT-1_3")
Idents(D,cells =cl1) <- "putSTB"
Idents(D,cells =cl3) <- "putSTB"
cl1 <- c("ACATCGACAGAGGGTT-1_7","ATATCCTAGAGAGCAA-1_7","ATCCCTGAGTGCACCC-1_7","CCTAAGATCCTTCGAC-1_7","GAGGGTACACTCATAG-1_7","GTAGTACCAGCGTTTA-1_7","TGCAGGCGTAGCGAGT-1_7","TTAGGGTGTACCATAC-1_7","AGGGTGATCCGAGATT-1_7")
Idents(D,cells =cl1) <- "putSTB"
newcl1 <- c("AGATAGATCTCCCTAG−1_3","ACTATTCGTCTGTGAT−1_3","GGAATCTGTCCTGTCT−1_3","AACTTCTGTAGCCCTG−1_3","TGGAACTCAGGTGTTT−1_3","CCACCATAGTGAGTTA−1_3","ATCCACCTCATTATCC−1_3","AGGTCATCAGAAGTTA−1_3","TACTGCCTCGACTCCT−1_3","TTGGGATCATGGATCT−1_3","GATTGGTGTCCAGGTC−1_3","TTTGACTAGGGTCTTT−1_3","ACTTAGGCAGTAACGG−1_3","GCCTGTTTCATGGCCG−1_3","CTGGTCTGTATGAGCG−1_3","GCTACAACAGTAGAAT−1_3","CAAGGGACATCTGGGC−1_3","GGAAGTGCATGACAGG−1_3","AGGATCTAGCGTATAA−1_3","CATTGCCGTAGCTGCC−1_3","GCCAGGTAGCTTGTGT−1_3","ATTCATCGTGTCTAAC−1_3")
newcl1 <- str_replace(newcl1,"−","-")
newcl2 <- c("TGTACAGGTTAAGTCC−1_3","GTGTGATGTCATATGC−1_3","TTCTAGTGTAGGAAAG−1_3","CGTAGTACAATACCTG−1_3","TGTCCCAAGTGAATAC−1_3","ATATCCTAGACTCAAA−1_3","ACTTATCAGGTAAGGA−1_3","ATTGTTCCAGCGTGCT−1_3","CCTATCGAGTTGCCCG−1_3","ATGCATGCAGGGTCTC−1_3","TGCAGGCAGGAGGCAG−1_3","ACACAGTTCTGCTTAT−1_3","CCCTGATAGCTCTGTA−1_3","ACGATGTCACAGCCAC−1_3","ACGTACACAGTTGCGC−1_3","GTGTAACCAAGCAATA−1_3","GTAACCACATTCTGTT−1_3","AGTACCACAACGGCTC−1_3","CATTGCCTCATCAGTG−1_3","ATTTCACTCACACCGG−1_3","ATCGGATAGAGGGTAA−1_3","GTGATGTCACCCTCTA−1_3","GAGACTTCATCAGCGC−1_3","CTTCAATTCCTCTCTT−1_3","ACCAACATCATCACTT−1_3","GACCCTTTCTAACACG−1_3","CCTCAACTCTTGCAAG−1_3","GCACGGTCAAAGCGTG−1_3","ATTACTCCAATACCTG−1_3","TCTGGCTGTACTCAAC−1_3","TCAAGACCAATCTCTT−1_3","CGAAGTTAGTGCTCAT−1_3","TGGAGAGAGTATAACG−1_3","GATCATGAGTATCTGC−1_3","TTGACCCTCCTTGGAA−1_3","TCGGTCTTCCACCTCA−1_3","GCCCGAAGTTGTTTGG−1_3","ATCCTATAGACAGTCG−1_3","GACGCTGTCGAGCCAC−1_3","GATCACAGTCACTTCC−1_3","CTCATGCCACACACGC−1_3","TACCTGCCATGCCGCA−1_3","CATACTTAGGTCGTGA−1_3","AGACCCGAGAGCAACC−1_3","TGAACGTAGTTTGAGA−1_3","ATCGTGACAGACCAGA−1_3","CATGCGGGTGTCTTCC−1_3","TCGCACTCAAATAGCA−1_3","GGAGGTATCTAGACCA−1_3","TCGATTTCAATTTCCT−1_3","CTCCAACGTTCTCCAC−1_3","CGCAGGTCAGAAACCG−1_3","TTGAACGTCTGGCCAG−1_3","CCTACGTAGGACTTCT−1_3","TCACACCCAGGTTACT−1_3","AGCTACATCTCGTGAA−1_3","GTAGCTACAACAGAGC−1_3","CGGGTGTAGTTGTAGA−1_3","ACTTCGCGTCTAGATC−1_3","TGGTGATCAATCTCGA−1_3","CGAGGCTAGGACAACC−1_3","GGCGTCACAATGTTGC−1_3","CCGGTGAAGCAATTAG−1_3","CCGTAGGGTAAGTTAG−1_3","AAGTTCGTCCACGTGG−1_3","GTCTACCCAACACGTT−1_3","ATAGGCTCACACGGAA−1_3","GCACATAGTCATCCCT−1_3","AGCTTCCTCGCAGTTA−1_3","CTCCCTCCAGCACAAG−1_3","ACGTTCCCAATCTGCA−1_3","ACCCTTGAGGTTTACC−1_3","TCGGGTGAGATTCGCT−1_3","TACTTACGTACCGGCT−1_3","CTCATCGCATCAACCA−1_3","TCCTCGATCTGCGAGC−1_3","TACCTGCGTGTTTCTT−1_3","AACACACAGGGCAATC−1_3","TCGAACACAAGAGCTG−1_3","CTAACCCAGATAGTGT−1_3","GAAATGAGTCTTGCGG−1_3","TACGGTATCATGCAGT−1_3","GTGGTTATCCGGACTG−1_3","TAGACCAAGGAGACCT−1_3")
newcl2 <- str_replace(newcl2,"−","-")
Idents(D,cells =newcl1) <- "putSTB"
Idents(D,cells =newcl2) <- "putSTB"
newcl_1 <- c("AGATAGATCTCCCTAG−1_3","ACTATTCGTCTGTGAT−1_3","GGAATCTGTCCTGTCT−1_3","AACTTCTGTAGCCCTG−1_3","TGGAACTCAGGTGTTT−1_3","CCACCATAGTGAGTTA−1_3","ATCCACCTCATTATCC−1_3","AGGTCATCAGAAGTTA−1_3","TACTGCCTCGACTCCT−1_3","TTGGGATCATGGATCT−1_3","GATTGGTGTCCAGGTC−1_3","TTTGACTAGGGTCTTT−1_3","ACTTAGGCAGTAACGG−1_3","GCCTGTTTCATGGCCG−1_3","CTGGTCTGTATGAGCG−1_3","GCTACAACAGTAGAAT−1_3","CAAGGGACATCTGGGC−1_3","GGAAGTGCATGACAGG−1_3","AGGATCTAGCGTATAA−1_3","CATTGCCGTAGCTGCC−1_3","GCCAGGTAGCTTGTGT−1_3","ATTCATCGTGTCTAAC−1_3")
newcl_2 <- c("GTGTGATGTCATATGC−1_3","TTCTAGTGTAGGAAAG−1_3","CGTAGTACAATACCTG−1_3","TGTCCCAAGTGAATAC−1_3","ATATCCTAGACTCAAA−1_3","ACTTATCAGGTAAGGA−1_3","ATTGTTCCAGCGTGCT−1_3","CCTATCGAGTTGCCCG−1_3","ATGCATGCAGGGTCTC−1_3","TGCAGGCAGGAGGCAG−1_3","ACACAGTTCTGCTTAT−1_3","CCCTGATAGCTCTGTA−1_3","ACGATGTCACAGCCAC−1_3","ACGTACACAGTTGCGC−1_3","GTGTAACCAAGCAATA−1_3","GTAACCACATTCTGTT−1_3","AGTACCACAACGGCTC−1_3","CATTGCCTCATCAGTG−1_3","ATTTCACTCACACCGG−1_3","ATCGGATAGAGGGTAA−1_3","GTGATGTCACCCTCTA−1_3","GAGACTTCATCAGCGC−1_3","CTTCAATTCCTCTCTT−1_3","ACCAACATCATCACTT−1_3","GACCCTTTCTAACACG−1_3","CCTCAACTCTTGCAAG−1_3","GCACGGTCAAAGCGTG−1_3","ATTACTCCAATACCTG−1_3","TCTGGCTGTACTCAAC−1_3","TCAAGACCAATCTCTT−1_3","CGAAGTTAGTGCTCAT−1_3","TGGAGAGAGTATAACG−1_3","GATCATGAGTATCTGC−1_3","TTGACCCTCCTTGGAA−1_3","TCGGTCTTCCACCTCA−1_3","GCCCGAAGTTGTTTGG−1_3","ATCCTATAGACAGTCG−1_3","GACGCTGTCGAGCCAC−1_3","GATCACAGTCACTTCC−1_3","CTCATGCCACACACGC−1_3","TACCTGCCATGCCGCA−1_3","CATACTTAGGTCGTGA−1_3","AGACCCGAGAGCAACC−1_3","TGAACGTAGTTTGAGA−1_3","ATCGTGACAGACCAGA−1_3","CATGCGGGTGTCTTCC−1_3","TCGCACTCAAATAGCA−1_3","GGAGGTATCTAGACCA−1_3","TCGATTTCAATTTCCT−1_3","CTCCAACGTTCTCCAC−1_3","CGCAGGTCAGAAACCG−1_3","TTGAACGTCTGGCCAG−1_3","CCTACGTAGGACTTCT−1_3","TCACACCCAGGTTACT−1_3","AGCTACATCTCGTGAA−1_3","GTAGCTACAACAGAGC−1_3","CGGGTGTAGTTGTAGA−1_3","ACTTCGCGTCTAGATC−1_3","TGGTGATCAATCTCGA−1_3","CGAGGCTAGGACAACC−1_3","GGCGTCACAATGTTGC−1_3","CCGGTGAAGCAATTAG−1_3","CCGTAGGGTAAGTTAG−1_3","AAGTTCGTCCACGTGG−1_3","GTCTACCCAACACGTT−1_3","ATAGGCTCACACGGAA−1_3","GCACATAGTCATCCCT−1_3","AGCTTCCTCGCAGTTA−1_3","CTCCCTCCAGCACAAG−1_3","ACGTTCCCAATCTGCA−1_3","ACCCTTGAGGTTTACC−1_3","TCGGGTGAGATTCGCT−1_3","TACTTACGTACCGGCT−1_3","CTCATCGCATCAACCA−1_3","TCCTCGATCTGCGAGC−1_3","TACCTGCGTGTTTCTT−1_3","AACACACAGGGCAATC−1_3","TCGAACACAAGAGCTG−1_3","CTAACCCAGATAGTGT−1_3","GAAATGAGTCTTGCGG−1_3","TACGGTATCATGCAGT−1_3","TACGGTATCATGCAGT−1_3","GTGGTTATCCGGACTG−1_3","TGTACAGGTTAAGTCC−1_3")
newcl_1 <- str_replace(newcl_1,"−","-")
newcl_2 <- str_replace(newcl_2,"−","-")
Idents(D,cells =newcl_1) <- "putSTB"
Idents(D,cells =newcl_2) <- "putSTB"
putPGC <- c("CAAGAGGTCCACCCTA-1_2","TTATTGCAGCACAAAT-1_2","TCGCACTAGCCAAGGT-1_2","GTAGATCAGCTGGCCT-1_2","ATCGATGGTAGCGTAG-1_2","GTGTTCCAGGGCAGGA-1_2","CAAGACTGTATTGCCA-1_2","TTCTTCCAGTCATCGT-1_2","CCTGCATCAAGTATCC-1_3","ACACCAACAACCGCTG-1_3","CCTCAACTCTTGCAAG-1_3","TAACTTCGTACAAAGT-1_3","AGCGCCACACGGCACT-1_3","TCATGTTGTCGACGCT-1_3","CAGTGCGAGTAGCCAG-1_3","ATGAGTCAGTCAATCC-1_3","TTACAGGGTCTTCTAT-1_3","CATCGTCTCGTTGTAG-1_6","AAAGAACAGTAAACGT-1_6","ACACGCGAGGAGGGTG-1_6","GCATCTCTCGGAATGG-1_6","GGAAGTGAGGCCTTGC-1_6")
Idents(D,cells=putPGC) <- "putPGC"

D$Cells <- Idents(D)
p<-DimPlot(D,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_AllDatasets.pdf",sep=""),width = 40, height = 8,p)
p<-DimPlot(D,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_AllDatasets.pdf",sep=""),width = 40, height = 8,p)

p<-DimPlot(D,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_AllDatasets_nosplit.pdf",sep=""),width = 10, height = 8,p)
p<-DimPlot(D,  pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_AllDatasets_nosplit.pdf",sep=""),width = 10, height = 8,p)

DefaultAssay(D) <- "RNA"

p<-FeaturePlot(D,  pt.size = 1, features = "RYBP", reduction = "pca" ) 
ggsave(filename=paste(saveext,"/DimRed/PCA_AllDatasets_nosplit_RYBP.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(D,  pt.size = 1, features = "PADI1", reduction = "umap") 
ggsave(filename=paste(saveext,"/DimRed/UMAP_AllDatasets_nosplit_PADI1.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(D,  pt.size = 1, features = "PADI2", reduction = "pca") 
ggsave(filename=paste(saveext,"/DimRed/PCA_AllDatasets_nosplit_PADI2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(D,  pt.size = 1, features = "PADI2", reduction = "umap") 
ggsave(filename=paste(saveext,"/DimRed/UMAP_AllDatasets_nosplit_PADI2.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(D,  pt.size = 1, features = "PADI3", reduction = "pca") 
ggsave(filename=paste(saveext,"/DimRed/PCA_AllDatasets_nosplit_PADI3.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(D,  pt.size = 1, features = "PADI3", reduction = "umap") 
ggsave(filename=paste(saveext,"/DimRed/UMAP_AllDatasets_nosplit_PADI3.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(D,  pt.size = 1, features = "PADI4", reduction = "pca") 
ggsave(filename=paste(saveext,"/DimRed/PCA_AllDatasets_nosplit_PADI4.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(D,  pt.size = 1, features = "PADI4", reduction = "umap") 
ggsave(filename=paste(saveext,"/DimRed/UMAP_AllDatasets_nosplit_PADI4.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(D,  pt.size = 1, features = "PADI6", reduction = "pca") 
ggsave(filename=paste(saveext,"/DimRed/PCA_AllDatasets_nosplit_PADI6.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(D,  pt.size = 1, features = "PADI6", reduction = "umap") 
ggsave(filename=paste(saveext,"/DimRed/UMAP_AllDatasets_nosplit_PADI6.pdf",sep=""),width = 10, height = 8,p)


Idents(D) <- D$Genotype
p<-DimPlot(D,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_AllDatasetsGT.pdf",sep=""),width = 40, height = 8,p)
p<-DimPlot(D,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_AllDatasetsGT.pdf",sep=""),width = 40, height = 8,p)

Dtest <- subset(D,idents=c("Embryonic_G1","Embryonic_G2","Embryonic_G4","Embryonic_G5"))
p<-DimPlot(Dtest,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_AllDatasetsGTsub.pdf",sep=""),width = 40, height = 8,p)
p<-DimPlot(Dtest,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_AllDatasetsGTsub.pdf",sep=""),width = 40, height = 8,p)
Idents(Dtest) <- Dtest$Cells
p<-DimPlot(Dtest,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_AllDatasetssub.pdf",sep=""),width = 40, height = 8,p)

p<-DimPlot(Dtest,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_AllDatasetssub.pdf",sep=""),width = 40, height = 8,p)

#correctedLabe <- as.character(D$Cells)
#correctedLabe[which(D$ID3=="C1" & D$Genotype%in%c("Embryonic_G1","Embryonic_G2","Embryonic_G4","Embryonic_G5") & D$Cells=="Cont")] <- "putExMes"
#correctedLabe[which(D$ID3=="C1" & D$Genotype%in%c("Embryonic_G1","Embryonic_G2","Embryonic_G4","Embryonic_G5") & D$Cells=="SOX9LRG5")] <- "Hyp/Am"
#correctedLabe[which(D$ID3=="C1" & D$Genotype%in%c("Embryonic_G1","Embryonic_G2","Embryonic_G4","Embryonic_G5") & D$Cells=="Prolif")] <- "Hyp/Am"
#correctedLabe[which(D$ID3=="C1" & D$Genotype%in%c("Embryonic_G1","Embryonic_G2","Embryonic_G4","Embryonic_G5") & D$Cells=="SOX9P")] <- "EmDisc_d14"
#correctedLabe[which(D$ID3=="C2" & D$Genotype%in%c("Embryonic_G1","Embryonic_G2","Embryonic_G4","Embryonic_G5") & D$Cells=="Cont")] <- "putExMes"
#correctedLabe[which(D$ID3=="C2" & D$Genotype%in%c("Embryonic_G1","Embryonic_G2","Embryonic_G4","Embryonic_G5") & D$Cells=="Lumenal")] <- "putExMes"
#correctedLabe[which(D$ID3=="C2" & D$Genotype%in%c("Embryonic_G1","Embryonic_G2","Embryonic_G4","Embryonic_G5") & D$Cells=="Prolif")] <- "Am/EmDisc"
#correctedLabe[which(D$ID3=="C2" & D$Genotype%in%c("Embryonic_G1","Embryonic_G2","Embryonic_G4","Embryonic_G5") & D$Cells=="SOX9P")] <- "EmDisc_d14"
#correctedLabe[which(D$ID3=="C4" & D$Genotype%in%c("Embryonic_G1","Embryonic_G2","Embryonic_G4","Embryonic_G5") & D$Cells=="SOX9LRG5")] <- "Am/EmDisc"
#correctedLabe[which(D$ID3=="C4" & D$Genotype%in%c("Embryonic_G1","Embryonic_G2","Embryonic_G4","Embryonic_G5") & D$Cells=="SOX9P")] <- "Am/EmDisc"
#correctedLabe[which(D$ID3=="C6" & D$Genotype%in%c("Embryonic_G1","Embryonic_G2","Embryonic_G4","Embryonic_G5") & D$Cells=="Cont")] <- "putExMes"
#correctedLabe[which(D$ID3=="C6" & D$Genotype%in%c("Embryonic_G1","Embryonic_G2","Embryonic_G4","Embryonic_G5") & D$Cells=="SOX9LRG5")] <- "Am/EmDisc"
#correctedLabe[which(D$ID3=="C6" & D$Genotype%in%c("Embryonic_G1","Embryonic_G2","Embryonic_G4","Embryonic_G5") & D$Cells=="Prolif")] <- "Am/EmDisc"
#correctedLabe[which(D$ID3=="C6" & D$Genotype%in%c("Embryonic_G1","Embryonic_G2","Embryonic_G4","Embryonic_G5") & D$Cells=="SOX9P")] <- "Am/EmDisc"
#correctedLabe[which(D$ID3=="C7" & D$Genotype%in%c("Embryonic_G1","Embryonic_G2","Embryonic_G4","Embryonic_G5") & D$Cells=="SOX9P")] <- "Am/EmDisc"
#correctedLabe[which(D$ID3=="C2" & D$Genotype%in% c("Epithelial_G4","Epithelial_G3","Stromal_G2","Epithelial_G2","Stromal_G4","Stromal_G6","Epithelial_G1","Epithelial_G6","Stromal_G3","Stromal_G1") & D$Cells=="EmDisc_d14")] <- "Prolif"
#correctedLabe[which(D$ID3=="C6" & D$Genotype%in% c("Epithelial_G4","Epithelial_G3","Stromal_G2","Epithelial_G2","Stromal_G4","Stromal_G6","Epithelial_G1","Epithelial_G6","Stromal_G3","Stromal_G1") & D$Cells=="EmDisc_d14")] <- "Stromal fibroblasts"

D$CorrectLabel <- D$Cells
Dtest2 <- subset(D,idents=c("Epithelial_G4","Epithelial_G3","Stromal_G2","Epithelial_G2","Stromal_G4","Stromal_G6","Epithelial_G1","Epithelial_G6","Stromal_G3","Stromal_G1"))
p<-DimPlot(Dtest2,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_AllDatasetsGTsub2.pdf",sep=""),width = 40, height = 8,p)

p<-DimPlot(Dtest2,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_AllDatasetsGTsub2.pdf",sep=""),width = 40, height = 8,p)
Idents(Dtest2) <- Dtest2$Cells
p<-DimPlot(Dtest2,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_AllDatasetssub2.pdf",sep=""),width = 40, height = 8,p)

p<-DimPlot(Dtest2,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_AllDatasetssub2.pdf",sep=""),width = 40, height = 8,p)






#Split by batch
Idents(D) <- paste(D$ID3,D$Cells,sep="_")
D1 <- subset(D,idents=c("C1_Lumenal","C1_SOX9P","C1_Prolif","C1_SOX9LRG5","C1_Ciliated","C1_Glandular"))
D2 <- subset(D,idents=c("C2_Lumenal","C2_SOX9P","C2_Prolif","C2_SOX9LRG5","C2_Ciliated","C2_Glandular"))
D3 <- subset(D,idents=c("C3_Lumenal","C3_SOX9P","C3_Prolif","C3_Ciliated","C3_Glandular"))
D4 <- subset(D,idents=c("C4_Lumenal","C4_SOX9P","C4_Prolif","C4_SOX9LRG5","C4_Ciliated","C4_Glandular"))
D6 <- subset(D,idents=c("C6_Lumenal","C6_SOX9P","C6_Prolif","C6_SOX9LRG5","C6_Ciliated","C6_Glandular"))
D7 <- subset(D,idents=c("C7_Lumenal","C7_SOX9P","C7_Prolif","C7_Ciliated","C7_Glandular"))
D8 <- subset(D,idents=c("C8_Lumenal","C8_SOX9P","C8_Prolif","C8_SOX9LRG5","C8_Ciliated","C8_Glandular"))
#rm(D)


#Also the assembloids data
#Load monkey
raw_counts1<-read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Assembloids/GSE168405_assembloid_counts.csv",sep=",", row.names=1, header=T)
BS<-read.table("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Assembloids/GSE168405_assembloid_metadata.csv",sep=",",header = T, row.names=1)

assemb  <- CreateSeuratObject(counts = raw_counts1, assay = "RNA",min.cells = 0, min.features = 0)
Idents(assemb) <- BS$Cluster
assemb <- NormalizeData(assemb, verbose = FALSE)
assemb <- FindVariableFeatures(assemb, selection.method = "vst", nfeatures = 20000)
assemb$ID1 <- "Assemb"
assemb$ID2 <- BS$Cluster
assemb$ID3 <- BS$SampleName
assemb$ID4 <- BS$Patient
assemb$ID5 <- BS$Treatment
Idents(assemb) <- BS$Cluster

Endo0 <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Placenta/Placentaharm.rds")
Idents(Endo0) <- Endo0$ID2
Endo0 <- subset(Endo0,idents=c("dS1","dS2","dS3","Epi1","Epi2"))

#Load in reference and split by stage
Endo <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Endometrium/other_endometrium_data.rds")
Idents(Endo) <- Endo$ID3
Endo <- subset(Endo,idents="early-secretory")
Idents(Endo) <- Endo$ID5
Endo <- subset(Endo,idents=c("Glandular","SOX9","Lumenal","Ciliated") ) #"dS","Fibroblast C7","eS"
List1 <- FindMarkers(Endo,ident.1 = "Lumenal",ident.2 = c("Glandular","SOX9","Ciliated"), test.use = "MAST", only.pos = TRUE )
List2 <- FindMarkers(Endo,ident.1 = "Ciliated",ident.2 = c("Lumenal","Glandular","SOX9"), test.use = "MAST", only.pos = TRUE )
List3 <- FindMarkers(Endo,ident.1 = "Glandular",ident.2 = c("Lumenal","Ciliated","SOX9"), test.use = "MAST", only.pos = TRUE )
List4 <- FindMarkers(Endo,ident.1 = "SOX9",ident.2 = c("Lumenal","Ciliated","Glandular"), test.use = "MAST", only.pos = TRUE )

Endo2 <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Endometrium/other_endometrium_data.rds")
Idents(Endo2) <- Endo2$ID3
Endo2 <- subset(Endo2,idents="proliferative")
Idents(Endo2) <- Endo2$ID5
Endo2 <- subset(Endo2,idents=c("Glandular","SOX9","Lumenal","Ciliated") ) #"dS","Fibroblast C7","eS"
List5 <- FindMarkers(Endo2,ident.1 = "Lumenal",ident.2 = c("Glandular","SOX9","Ciliated"), test.use = "MAST", only.pos = TRUE )
List6 <- FindMarkers(Endo2,ident.1 = "Ciliated",ident.2 = c("Lumenal","Glandular","SOX9"), test.use = "MAST", only.pos = TRUE )
List7 <- FindMarkers(Endo2,ident.1 = "Glandular",ident.2 = c("Lumenal","Ciliated","SOX9"), test.use = "MAST", only.pos = TRUE )
List8 <- FindMarkers(Endo2,ident.1 = "SOX9",ident.2 = c("Lumenal","Ciliated","Glandular"), test.use = "MAST", only.pos = TRUE )


Endo3 <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Endometrium/other_endometrium_data.rds")
Idents(Endo3) <- Endo3$ID3
Endo3 <- subset(Endo3,idents=c("early-secretory","proliferative","late-secretory") )
#Idents(Endo3) <- Endo3$ID5
List9 <- FindMarkers(Endo3,ident.1 = "early-secretory",ident.2 = c("proliferative"), test.use = "MAST",  only.pos = TRUE )
List10 <- FindMarkers(Endo3,ident.1 = "proliferative",ident.2 = c("early-secretory"), test.use = "MAST",  only.pos = TRUE)

#Gland markers up regulated in early secretory phase
List11 <- intersect(intersect(rownames(List7)[which(List7$p_val_adj<0.1 & List7$avg_log2FC > 0)], rownames(List9)[which(List9$p_val_adj<0.1)] ), rownames(D))




#Now do the same but look only at our datasets ...
D <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/CorrectedFinalAno2.rds")
Idents(D) <- D$CorrectLabel
cl2 <- c("TCTCCGAGTAGGCTCC-1_2","CACTTCGCAGCTGCCA-1_2","AATCGACTCCAGCAAT-1_2","GCCCAGATCGTTGCCT-1_2","GTGAGCCAGCACGTCC-1_2")
Idents(D,cells =cl2) <- "putSTB"
cl1 <- c("AACTTCTGTAGCCCTG-1_3","ACTATTCGTCTGTGAT-1_3","ACTTAGGCAGTAACGG-1_3","AGATAGATCTCCCTAG-1_3","AGGATCTAGCGTATAA-1_3","AGGTCATCAGAAGTTA-1_3","ATACTTCTCATGCCAA-1_3","ATCCACCTCATTATCC-1_3","ATTCATCGTGTCTAAC-1_3","CAAGGGACATCTGGGC-1_3","CCACCATAGTGAGTTA-1_3","CGAGAAGCATTGACAC-1_3","CTGGTCTGTATGAGCG-1_3","GATTGGTGTCCAGGTC-1_3","GCCAGGTAGCTTGTGT-1_3","GCCTGTTTCATGGCCG-1_3","GCTACAACAGTAGAAT-1_3","GGAAGTGCATGACAGG-1_3","GGAATCTGTCCTGTCT-1_3","TACTGCCTCGACTCCT-1_3","TGGAACTCAGGTGTTT-1_3","TTGGGATCATGGATCT-1_3","TTTGACTAGGGTCTTT-1_3")
cl3 <- c("AACACACAGGGCAATC-1_3","AAGTTCGTCCACGTGG-1_3","ACACAGTTCGATACTG-1_3","ACACAGTTCTGCTTAT-1_3","ACACTGATCCATACTT-1_3","ACCAACATCATCACTT-1_3","ACCCTTGAGGTTTACC-1_3","ACGATGTCACAGCCAC-1_3","ACGTACACAGTTGCGC-1_3","ACGTTCCCAATCTGCA-1_3","ACTATCTTCGTTGTTT-1_3","ACTATTCAGTCTACCA-1_3","ACTTATCAGGTAAGGA-1_3","ACTTCGCGTCTAGATC-1_3","AGACCCGAGAGCAACC-1_3","AGCTACATCTCGTGAA-1_3","AGCTTCCTCGCAGTTA-1_3","AGTACCACAACGGCTC-1_3","ATAGGCTCACACGGAA-1_3","ATATCCTAGACTCAAA-1_3","ATCCTATAGACAGTCG-1_3","ATCGGATAGAGGGTAA-1_3","ATCGTGACAGACCAGA-1_3","ATCTCTAGTACCTATG-1_3","ATCTTCACACTTGAGT-1_3","ATGCATGCAGGGTCTC-1_3","ATTACTCCAATACCTG-1_3","ATTCTTGTCAGACCCG-1_3","ATTGTTCCAGCGTGCT-1_3","ATTTCACTCACACCGG-1_3","CACAACAAGTAAACAC-1_3","CATACTTAGGTCGTGA-1_3","CATGCCTCAGGGATAC-1_3","CATGCGGGTGTCTTCC-1_3","CATTGCCTCATCAGTG-1_3","CCCTGATAGCTCTGTA-1_3","CCGGTGAAGCAATTAG-1_3","CCGTAGGGTAAGTTAG-1_3","CCTACGTAGGACTTCT-1_3","CCTATCGAGCACTTTG-1_3","CCTATCGAGTTGCCCG-1_3","CCTCAACTCTTGCAAG-1_3","CCTCTCCTCGTTGTGA-1_3","CCTTTGGAGGTTGTTC-1_3","CGAAGGATCACCACAA-1_3","CGAAGTTAGTGCTCAT-1_3","CGAGGCTAGGACAACC-1_3","CGATGCGCATCATTGG-1_3","CGCAGGTCAGAAACCG-1_3","CGCAGGTTCGGCCTTT-1_3","CGGGTGTAGTTGTAGA-1_3","CGTAGTACAATACCTG-1_3","CTAACCCAGATAGTGT-1_3","CTAACCCAGGCACCAA-1_3","CTCAGGGCACGATTCA-1_3","CTCATCGCATCAACCA-1_3","CTCATGCCACACACGC-1_3","CTCCAACGTTCTCCAC-1_3","CTCCCTCCAGCACAAG-1_3","CTTCAATTCCTCTCTT-1_3","GAAATGAGTCTTGCGG-1_3","GACCCTTTCTAACACG-1_3","GACGCTGTCGAGCCAC-1_3","GAGACCCCACTGCGAC-1_3","GAGACTTCATCAGCGC-1_3","GAGTTACAGCAGCCCT-1_3","GATCACAGTCACTTCC-1_3","GATCATGAGTATCTGC-1_3","GCACATAGTCATCCCT-1_3","GCACGGTCAAAGCGTG-1_3","GCCCGAAGTTGTTTGG-1_3","GGACGTCCAGCCCACA-1_3","GGAGGTATCTAGACCA-1_3","GGCGTCACAATGTTGC-1_3","GTAACCACATTCTGTT-1_3","GTAGCTACAACAGAGC-1_3","GTAGGTTGTCCTCATC-1_3","GTCTACCCAACACGTT-1_3","GTGATGTCACCCTCTA-1_3","GTGGTTATCCGGACTG-1_3","GTGTAACCAAGCAATA-1_3","GTGTGATGTCATATGC-1_3","TAATTCCCACACCAGC-1_3","TACCTGCCATGCCGCA-1_3","TACCTGCGTGTTTCTT-1_3","TACGGTATCATGCAGT-1_3","TACTTACGTACCGGCT-1_3","TACTTCACATGAATCC-1_3","TAGACCAAGGAGACCT-1_3","TCAAGACCAATCTCTT-1_3","TCACAAGTCACGAGGA-1_3","TCACACCCAGGTTACT-1_3","TCACGCTTCAGGAAAT-1_3","TCATACTAGTATAACG-1_3","TCCTCGATCTGCGAGC-1_3","TCGAACACAAGAGCTG-1_3","TCGATTTCAATTTCCT-1_3","TCGCACTCAAATAGCA-1_3","TCGGGTGAGATTCGCT-1_3","TCGGTCTTCCACCTCA-1_3","TCGTGGGCAACGGCCT-1_3","TCTGGCTGTACTCAAC-1_3","TGAACGTAGTTTGAGA-1_3","TGACTCCTCAGTGCGC-1_3","TGAGTCATCATCTACT-1_3","TGATCTTCACACCTGG-1_3","TGATCTTCACACCTGG-1_3","TGCAGGCAGGAGGCAG-1_3","TGCCGAGCATCATTTC-1_3","TGGAGAGAGTATAACG-1_3","TGGTGATCAATCTCGA-1_3","TGTACAGGTTAAGTCC-1_3","TGTCCCAAGTGAATAC-1_3","TTCATGTCACCTTCGT-1_3","TTCTAGTGTAGGAAAG-1_3","TTCTTGACACGAAGAC-1_3","TTGAACGTCTGGCCAG-1_3","TTGACCCAGCTGTCCG-1_3","TTGACCCTCCTTGGAA-1_3","TTTACCATCATCGGGC-1_3","TTTATGCTCCGGACGT-1_3")
Idents(D,cells =cl1) <- "putSTB"
Idents(D,cells =cl3) <- "putSTB"
cl1 <- c("ACATCGACAGAGGGTT-1_7","ATATCCTAGAGAGCAA-1_7","ATCCCTGAGTGCACCC-1_7","CCTAAGATCCTTCGAC-1_7","GAGGGTACACTCATAG-1_7","GTAGTACCAGCGTTTA-1_7","TGCAGGCGTAGCGAGT-1_7","TTAGGGTGTACCATAC-1_7","AGGGTGATCCGAGATT-1_7")
Idents(D,cells =cl1) <- "putSTB"
newcl1 <- c("AGATAGATCTCCCTAG−1_3","ACTATTCGTCTGTGAT−1_3","GGAATCTGTCCTGTCT−1_3","AACTTCTGTAGCCCTG−1_3","TGGAACTCAGGTGTTT−1_3","CCACCATAGTGAGTTA−1_3","ATCCACCTCATTATCC−1_3","AGGTCATCAGAAGTTA−1_3","TACTGCCTCGACTCCT−1_3","TTGGGATCATGGATCT−1_3","GATTGGTGTCCAGGTC−1_3","TTTGACTAGGGTCTTT−1_3","ACTTAGGCAGTAACGG−1_3","GCCTGTTTCATGGCCG−1_3","CTGGTCTGTATGAGCG−1_3","GCTACAACAGTAGAAT−1_3","CAAGGGACATCTGGGC−1_3","GGAAGTGCATGACAGG−1_3","AGGATCTAGCGTATAA−1_3","CATTGCCGTAGCTGCC−1_3","GCCAGGTAGCTTGTGT−1_3","ATTCATCGTGTCTAAC−1_3")
newcl1 <- str_replace(newcl1,"−","-")
newcl2 <- c("TGTACAGGTTAAGTCC−1_3","GTGTGATGTCATATGC−1_3","TTCTAGTGTAGGAAAG−1_3","CGTAGTACAATACCTG−1_3","TGTCCCAAGTGAATAC−1_3","ATATCCTAGACTCAAA−1_3","ACTTATCAGGTAAGGA−1_3","ATTGTTCCAGCGTGCT−1_3","CCTATCGAGTTGCCCG−1_3","ATGCATGCAGGGTCTC−1_3","TGCAGGCAGGAGGCAG−1_3","ACACAGTTCTGCTTAT−1_3","CCCTGATAGCTCTGTA−1_3","ACGATGTCACAGCCAC−1_3","ACGTACACAGTTGCGC−1_3","GTGTAACCAAGCAATA−1_3","GTAACCACATTCTGTT−1_3","AGTACCACAACGGCTC−1_3","CATTGCCTCATCAGTG−1_3","ATTTCACTCACACCGG−1_3","ATCGGATAGAGGGTAA−1_3","GTGATGTCACCCTCTA−1_3","GAGACTTCATCAGCGC−1_3","CTTCAATTCCTCTCTT−1_3","ACCAACATCATCACTT−1_3","GACCCTTTCTAACACG−1_3","CCTCAACTCTTGCAAG−1_3","GCACGGTCAAAGCGTG−1_3","ATTACTCCAATACCTG−1_3","TCTGGCTGTACTCAAC−1_3","TCAAGACCAATCTCTT−1_3","CGAAGTTAGTGCTCAT−1_3","TGGAGAGAGTATAACG−1_3","GATCATGAGTATCTGC−1_3","TTGACCCTCCTTGGAA−1_3","TCGGTCTTCCACCTCA−1_3","GCCCGAAGTTGTTTGG−1_3","ATCCTATAGACAGTCG−1_3","GACGCTGTCGAGCCAC−1_3","GATCACAGTCACTTCC−1_3","CTCATGCCACACACGC−1_3","TACCTGCCATGCCGCA−1_3","CATACTTAGGTCGTGA−1_3","AGACCCGAGAGCAACC−1_3","TGAACGTAGTTTGAGA−1_3","ATCGTGACAGACCAGA−1_3","CATGCGGGTGTCTTCC−1_3","TCGCACTCAAATAGCA−1_3","GGAGGTATCTAGACCA−1_3","TCGATTTCAATTTCCT−1_3","CTCCAACGTTCTCCAC−1_3","CGCAGGTCAGAAACCG−1_3","TTGAACGTCTGGCCAG−1_3","CCTACGTAGGACTTCT−1_3","TCACACCCAGGTTACT−1_3","AGCTACATCTCGTGAA−1_3","GTAGCTACAACAGAGC−1_3","CGGGTGTAGTTGTAGA−1_3","ACTTCGCGTCTAGATC−1_3","TGGTGATCAATCTCGA−1_3","CGAGGCTAGGACAACC−1_3","GGCGTCACAATGTTGC−1_3","CCGGTGAAGCAATTAG−1_3","CCGTAGGGTAAGTTAG−1_3","AAGTTCGTCCACGTGG−1_3","GTCTACCCAACACGTT−1_3","ATAGGCTCACACGGAA−1_3","GCACATAGTCATCCCT−1_3","AGCTTCCTCGCAGTTA−1_3","CTCCCTCCAGCACAAG−1_3","ACGTTCCCAATCTGCA−1_3","ACCCTTGAGGTTTACC−1_3","TCGGGTGAGATTCGCT−1_3","TACTTACGTACCGGCT−1_3","CTCATCGCATCAACCA−1_3","TCCTCGATCTGCGAGC−1_3","TACCTGCGTGTTTCTT−1_3","AACACACAGGGCAATC−1_3","TCGAACACAAGAGCTG−1_3","CTAACCCAGATAGTGT−1_3","GAAATGAGTCTTGCGG−1_3","TACGGTATCATGCAGT−1_3","GTGGTTATCCGGACTG−1_3","TAGACCAAGGAGACCT−1_3")
newcl2 <- str_replace(newcl2,"−","-")
Idents(D,cells =newcl1) <- "putSTB"
Idents(D,cells =newcl2) <- "putSTB"
newcl_1 <- c("AGATAGATCTCCCTAG−1_3","ACTATTCGTCTGTGAT−1_3","GGAATCTGTCCTGTCT−1_3","AACTTCTGTAGCCCTG−1_3","TGGAACTCAGGTGTTT−1_3","CCACCATAGTGAGTTA−1_3","ATCCACCTCATTATCC−1_3","AGGTCATCAGAAGTTA−1_3","TACTGCCTCGACTCCT−1_3","TTGGGATCATGGATCT−1_3","GATTGGTGTCCAGGTC−1_3","TTTGACTAGGGTCTTT−1_3","ACTTAGGCAGTAACGG−1_3","GCCTGTTTCATGGCCG−1_3","CTGGTCTGTATGAGCG−1_3","GCTACAACAGTAGAAT−1_3","CAAGGGACATCTGGGC−1_3","GGAAGTGCATGACAGG−1_3","AGGATCTAGCGTATAA−1_3","CATTGCCGTAGCTGCC−1_3","GCCAGGTAGCTTGTGT−1_3","ATTCATCGTGTCTAAC−1_3")
newcl_2 <- c("GTGTGATGTCATATGC−1_3","TTCTAGTGTAGGAAAG−1_3","CGTAGTACAATACCTG−1_3","TGTCCCAAGTGAATAC−1_3","ATATCCTAGACTCAAA−1_3","ACTTATCAGGTAAGGA−1_3","ATTGTTCCAGCGTGCT−1_3","CCTATCGAGTTGCCCG−1_3","ATGCATGCAGGGTCTC−1_3","TGCAGGCAGGAGGCAG−1_3","ACACAGTTCTGCTTAT−1_3","CCCTGATAGCTCTGTA−1_3","ACGATGTCACAGCCAC−1_3","ACGTACACAGTTGCGC−1_3","GTGTAACCAAGCAATA−1_3","GTAACCACATTCTGTT−1_3","AGTACCACAACGGCTC−1_3","CATTGCCTCATCAGTG−1_3","ATTTCACTCACACCGG−1_3","ATCGGATAGAGGGTAA−1_3","GTGATGTCACCCTCTA−1_3","GAGACTTCATCAGCGC−1_3","CTTCAATTCCTCTCTT−1_3","ACCAACATCATCACTT−1_3","GACCCTTTCTAACACG−1_3","CCTCAACTCTTGCAAG−1_3","GCACGGTCAAAGCGTG−1_3","ATTACTCCAATACCTG−1_3","TCTGGCTGTACTCAAC−1_3","TCAAGACCAATCTCTT−1_3","CGAAGTTAGTGCTCAT−1_3","TGGAGAGAGTATAACG−1_3","GATCATGAGTATCTGC−1_3","TTGACCCTCCTTGGAA−1_3","TCGGTCTTCCACCTCA−1_3","GCCCGAAGTTGTTTGG−1_3","ATCCTATAGACAGTCG−1_3","GACGCTGTCGAGCCAC−1_3","GATCACAGTCACTTCC−1_3","CTCATGCCACACACGC−1_3","TACCTGCCATGCCGCA−1_3","CATACTTAGGTCGTGA−1_3","AGACCCGAGAGCAACC−1_3","TGAACGTAGTTTGAGA−1_3","ATCGTGACAGACCAGA−1_3","CATGCGGGTGTCTTCC−1_3","TCGCACTCAAATAGCA−1_3","GGAGGTATCTAGACCA−1_3","TCGATTTCAATTTCCT−1_3","CTCCAACGTTCTCCAC−1_3","CGCAGGTCAGAAACCG−1_3","TTGAACGTCTGGCCAG−1_3","CCTACGTAGGACTTCT−1_3","TCACACCCAGGTTACT−1_3","AGCTACATCTCGTGAA−1_3","GTAGCTACAACAGAGC−1_3","CGGGTGTAGTTGTAGA−1_3","ACTTCGCGTCTAGATC−1_3","TGGTGATCAATCTCGA−1_3","CGAGGCTAGGACAACC−1_3","GGCGTCACAATGTTGC−1_3","CCGGTGAAGCAATTAG−1_3","CCGTAGGGTAAGTTAG−1_3","AAGTTCGTCCACGTGG−1_3","GTCTACCCAACACGTT−1_3","ATAGGCTCACACGGAA−1_3","GCACATAGTCATCCCT−1_3","AGCTTCCTCGCAGTTA−1_3","CTCCCTCCAGCACAAG−1_3","ACGTTCCCAATCTGCA−1_3","ACCCTTGAGGTTTACC−1_3","TCGGGTGAGATTCGCT−1_3","TACTTACGTACCGGCT−1_3","CTCATCGCATCAACCA−1_3","TCCTCGATCTGCGAGC−1_3","TACCTGCGTGTTTCTT−1_3","AACACACAGGGCAATC−1_3","TCGAACACAAGAGCTG−1_3","CTAACCCAGATAGTGT−1_3","GAAATGAGTCTTGCGG−1_3","TACGGTATCATGCAGT−1_3","TACGGTATCATGCAGT−1_3","GTGGTTATCCGGACTG−1_3","TGTACAGGTTAAGTCC−1_3")
newcl_1 <- str_replace(newcl_1,"−","-")
newcl_2 <- str_replace(newcl_2,"−","-")
Idents(D,cells =newcl_1) <- "putSTB"
Idents(D,cells =newcl_2) <- "putSTB"
putPGC <- c("CAAGAGGTCCACCCTA-1_2","TTATTGCAGCACAAAT-1_2","TCGCACTAGCCAAGGT-1_2","GTAGATCAGCTGGCCT-1_2","ATCGATGGTAGCGTAG-1_2","GTGTTCCAGGGCAGGA-1_2","CAAGACTGTATTGCCA-1_2","TTCTTCCAGTCATCGT-1_2","CCTGCATCAAGTATCC-1_3","ACACCAACAACCGCTG-1_3","CCTCAACTCTTGCAAG-1_3","TAACTTCGTACAAAGT-1_3","AGCGCCACACGGCACT-1_3","TCATGTTGTCGACGCT-1_3","CAGTGCGAGTAGCCAG-1_3","ATGAGTCAGTCAATCC-1_3","TTACAGGGTCTTCTAT-1_3","CATCGTCTCGTTGTAG-1_6","AAAGAACAGTAAACGT-1_6","ACACGCGAGGAGGGTG-1_6","GCATCTCTCGGAATGG-1_6","GGAAGTGAGGCCTTGC-1_6")
Idents(D,cells=putPGC) <- "putPGC"
D$Cells <- Idents(D)

D <- subset(D,idents=c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular","Stromal fibroblasts"))
Idents(D) <- D$Genotype
#D <- subset(D,idents=c("Stromal_G4","Stromal_G6","Stromal_G3","Stromal_G1","Stromal_G2","Embryonic_G1","Embryonic_G4","Embryonic_G2","Embryoic_G3"),invert=TRUE)
#Idents(D) <- D$Genotype
#D <- subset(D,idents=c("Stromal_G2"),invert=TRUE)
#Idents(D) <- D$Cells
Idents(D) <- D$Cells

DELISTA <- FindMarkers(D,ident.1=c("Prolif"), ident.2 = c("SOX9P","Lumenal","SOX9LRG5","Ciliated","Glandular"), test.use = "MAST", logfc.threshold = log(0))
DELISTB <- FindMarkers(D,ident.1=c("Lumenal"), ident.2 = c("SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular"), test.use = "MAST", logfc.threshold = log(0))
DELISTC <- FindMarkers(D,ident.1=c("Ciliated"), ident.2 = c("Lumenal","SOX9P","Prolif","SOX9LRG5","Glandular"), test.use = "MAST", logfc.threshold = log(0))
DELISTD <- FindMarkers(D,ident.1=c("Glandular"), ident.2 = c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated"), test.use = "MAST", logfc.threshold = log(0))
DELISTE <- FindMarkers(D,ident.1=c("SOX9LRG5","SOX9P"), ident.2 = c("Lumenal","Prolif","Ciliated","Glandular"), test.use = "MAST", logfc.threshold = log(0))

DELISTGeneral <- FindMarkers(D,ident.2=c("Stromal fibroblasts"), ident.1 = c("Prolif","SOX9P","Lumenal","SOX9LRG5","Ciliated","Glandular"), test.use = "MAST", only.pos = TRUE)




###
#1,2,3,4 - lumenal cilliated glandular sox9
#5,6,7,8

AvExp <- AverageExpression(D)

UberList <- c(
  intersect(rownames(List1)[which(List1$avg_log2FC > log2(1.2) & List1$p_val_adj<0.01)],
  rownames(List5)[which(List5$avg_log2FC > log2(1.2) & List5$p_val_adj<0.01)] ),
  intersect(rownames(List2)[which(List2$avg_log2FC > log2(1.5) & List2$p_val_adj<0.01)],
            rownames(List6)[which(List6$avg_log2FC > log2(1.5) & List6$p_val_adj<0.01)] ),
  intersect(rownames(List3)[which(List3$avg_log2FC > log2(1.5) & List3$p_val_adj<0.01)],
            rownames(List7)[which(List7$avg_log2FC > log2(1.5) & List7$p_val_adj<0.01)] ),
  intersect(rownames(List4)[which(List4$avg_log2FC > log2(1.2) & List4$p_val_adj<0.01)],
            rownames(List8)[which(List8$avg_log2FC > log2(1.2) & List8$p_val_adj<0.01)] )
)
Expmax <- apply(AvExp$RNA, 1, mean)  
Expmax <- Expmax[UberList]
#UberList <- intersect(UberList,rownames(DELISTGeneral)[which(DELISTGeneral$avg_log2FC > log2(1.2) & DELISTGeneral$p_val_adj<0.01)] )

UberList <- UberList[which(Expmax>.8)]

X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C1") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch1_EndometrialUberList",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch2_EndometrialUberList",".pdf",sep="") ,width=55,height=38)



X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C3","C4") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch3_EndometrialUberList",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C6") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch6_EndometrialUberList",".pdf",sep="") ,width=55,height=38)



X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C7") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch7_EndometrialUberList",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C8") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch8_EndometrialUberList",".pdf",sep="") ,width=55,height=38)





Exp <- AverageExpression(D)

UberList <- c(intersect(rownames(List1)[which(List1$avg_log2FC > log2(1.2) & List1$p_val_adj<0.01)],
            rownames(List5)[which(List5$avg_log2FC > log2(1.2) & List5$p_val_adj<0.01)] ))
Expmax <- apply(Exp$RNA, 1, mean)  
Expmax <- Expmax[UberList]
UberList <- UberList[which(Expmax>.8)]
X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C1") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch1_EndometrialLumenalList",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch2_EndometrialLumenalList",".pdf",sep="") ,width=55,height=38)



X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C3","C4") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch3_EndometrialLumenalList",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C6") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch6_EndometrialLumenalList",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C7") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch7_EndometrialLumenalList",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C8") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch8_EndometrialLumenalList",".pdf",sep="") ,width=55,height=38)







UberList <- c(intersect(rownames(List2)[which(List2$avg_log2FC > log2(1.2) & List2$p_val_adj<0.01)],
                        rownames(List6)[which(List6$avg_log2FC > log2(1.2) & List6$p_val_adj<0.01)] ))
Expmax <- apply(Exp$RNA, 1, mean)  
Expmax <- Expmax[UberList]
UberList <- UberList[which(Expmax>.8)]
X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C1") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch1_EndometrialCiliatedList",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch2_EndometrialCiliatedList",".pdf",sep="") ,width=55,height=38)



X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C3","C4") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch3_EndometrialCiliatedList",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C6") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch6_EndometrialCiliatedList",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C7") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch7_EndometrialCiliatedList",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C8") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch8_EndometrialCiliatedList",".pdf",sep="") ,width=55,height=38)







UberList <- c(intersect(rownames(List3)[which(List3$avg_log2FC > log2(1.2) & List3$p_val_adj<0.01)],
                        rownames(List7)[which(List7$avg_log2FC > log2(1.2) & List7$p_val_adj<0.01)] ))
Expmax <- apply(Exp$RNA, 1, mean)  
Expmax <- Expmax[UberList]
UberList <- UberList[which(Expmax>.8)]
X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C1") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch1_EndometrialGlandularList",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch2_EndometrialGlandularList",".pdf",sep="") ,width=55,height=38)



X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C3","C4") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch3_EndometrialGlandularList",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C6") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch6_EndometrialGlandulardList",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C7") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch7_EndometrialGlandularList",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C8") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch8_EndometrialGlandularList",".pdf",sep="") ,width=55,height=38)








UberList <- c(intersect(rownames(List4)[which(List4$avg_log2FC > log2(1.2) & List4$p_val_adj<0.01)],
                        rownames(List8)[which(List8$avg_log2FC > log2(1.2) & List8$p_val_adj<0.01)] ))
Expmax <- apply(Exp$RNA, 1, mean)  
Expmax <- Expmax[UberList]
UberList <- UberList[which(Expmax>.8)]
X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C1") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch1_EndometrialSOX9List",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch2_EndometrialSOX9List",".pdf",sep="") ,width=55,height=38)



X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C3","C4") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch3_EndometrialSOX9List",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C6") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch6_EndometrialSOX9dList",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C7") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch7_EndometrialSOX9List",".pdf",sep="") ,width=55,height=38)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C8") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Batch8_EndometrialSOX9List",".pdf",sep="") ,width=55,height=38)



















X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C1") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( c(rownames(DELISTA)[which(DELISTA$p_val_adj<0.1)], rownames(DELISTB)[which(DELISTB$p_val_adj<0.1)]  )  ,    intersect(rownames(List1)[which(List1$p_val_adj<0.1)], rownames(D) ) )   ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_All_Lumenal_Batch1Epithelia",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C1") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTC)[which(DELISTC$p_val_adj<0.1)] )  ,  intersect(rownames(List2)[which(List2$p_val_adj<0.1)]  , rownames(D) ) )   ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Ciliated_Batch1Epithelia",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C1") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTD)[which(DELISTD$p_val_adj<0.1)]  )  ,  intersect(  rownames(List3)[which(List3$p_val_adj<0.1)] , rownames(D) ) )  ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Glandular_Batch1Epithelia",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C1") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect(  rownames(List4)[which(List4$p_val_adj<0.1)] , rownames(D) )   ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_SOX9_Batch1Epithelia",".pdf",sep="") ,width=10,height=10)

#Batch 2
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( c(rownames(DELISTA)[which(DELISTA$p_val_adj<0.1)], rownames(DELISTB)[which(DELISTB$p_val_adj<0.1)]  )  ,    intersect(rownames(List1)[which(List1$p_val_adj<0.1)], rownames(D) ) )   ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Lumenal_Batch2Epithelia",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTC)[which(DELISTC$p_val_adj<0.1)] )  ,  intersect(rownames(List2)[which(List2$p_val_adj<0.1)]  , rownames(D) ) )   ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Ciliated_Batch2Epithelia",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTD)[which(DELISTD$p_val_adj<0.1)]  )  ,  intersect(  rownames(List3)[which(List3$p_val_adj<0.1)] , rownames(D) ) )  ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Glandular_Batch2Epithelia",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect(  rownames(List4)[which(List4$p_val_adj<0.1)] , rownames(D) )   ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_SOX9_Batch2Epithelia",".pdf",sep="") ,width=10,height=10)

#Batch 3
subs <- which(D$ID3%in%c("C3","C4") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( c(rownames(DELISTA)[which(DELISTA$p_val_adj<0.1)], rownames(DELISTB)[which(DELISTB$p_val_adj<0.1)]  )  ,    intersect(rownames(List1)[which(List1$p_val_adj<0.1)], rownames(D) ) )   ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Lumenal_Batch3Epithelia",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTC)[which(DELISTC$p_val_adj<0.1)] )  ,  intersect(rownames(List2)[which(List2$p_val_adj<0.1)]  , rownames(D) ) )   ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Ciliated_Batch3Epithelia",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTD)[which(DELISTD$p_val_adj<0.1)]  )  ,  intersect(  rownames(List3)[which(List3$p_val_adj<0.1)] , rownames(D) ) )  ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Glandular_Batch3Epithelia",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTE)[which(DELISTE$p_val_adj<0.1)]  )  ,  intersect(  rownames(List4)[which(List4$p_val_adj<0.1)] , rownames(D) ) )  ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_SOX9_Batch3Epithelia",".pdf",sep="") ,width=10,height=10)


#Batch 7
subs <- which(D$ID3%in%c("C7") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( c(rownames(DELISTA)[which(DELISTA$p_val_adj<0.1)], rownames(DELISTB)[which(DELISTB$p_val_adj<0.1)]  )  ,    intersect(rownames(List1)[which(List1$p_val_adj<0.1)], rownames(D) ) )   ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Lumenal_Batch7Epithelia",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTC)[which(DELISTC$p_val_adj<0.1)] )  ,  intersect(rownames(List2)[which(List2$p_val_adj<0.1)]  , rownames(D) ) )   ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Ciliated_Batch7Epithelia",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTD)[which(DELISTD$p_val_adj<0.1)]  )  ,  intersect(  rownames(List3)[which(List3$p_val_adj<0.1)] , rownames(D) ) )  ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Glandular_Batch7Epithelia",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTE)[which(DELISTE$p_val_adj<0.1)]  )  ,  intersect(  rownames(List4)[which(List4$p_val_adj<0.1)] , rownames(D) ) )  ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_SOX9_Batch7Epithelia",".pdf",sep="") ,width=10,height=10)

#Batch 8
subs <- which(D$ID3%in%c("C8") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( c(rownames(DELISTA)[which(DELISTA$p_val_adj<0.1)], rownames(DELISTB)[which(DELISTB$p_val_adj<0.1)]  )  ,    intersect(rownames(List1)[which(List1$p_val_adj<0.1)], rownames(D) ) )   ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Lumenal_Batch8Epithelia",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTC)[which(DELISTC$p_val_adj<0.1)] )  ,  intersect(rownames(List2)[which(List2$p_val_adj<0.1)]  , rownames(D) ) )   ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Ciliated_Batch8Epithelia",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTD)[which(DELISTD$p_val_adj<0.1)]  )  ,  intersect(  rownames(List3)[which(List3$p_val_adj<0.1)] , rownames(D) ) )  ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Glandular_Batch8Epithelia",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(D,assay = "RNA")) 
subs <- which(D$ID3%in%c("C2") & D$Cells%in%c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular") )
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTE)[which(DELISTE$p_val_adj<0.1)]  )  ,  intersect(  rownames(List4)[which(List4$p_val_adj<0.1)] , rownames(D) ) )  ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(D$Cells[subs])), Batch =  factor((D$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_SOX9_Batch8Epithelia",".pdf",sep="") ,width=10,height=10)


###All
X <- (GetAssayData(D,assay = "RNA")) 
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTE)[which(DELISTE$p_val_adj<0.1)]  )  ,  intersect(  rownames(List4)[which(List4$p_val_adj<0.1)] , rownames(D) ) )  ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(D$Cells)), Batch =  factor((D$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_SOX9_Epithelia",".pdf",sep="") ,width=40,height=30)

Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( c(rownames(DELISTA)[which(DELISTA$p_val_adj<0.1)], rownames(DELISTB)[which(DELISTB$p_val_adj<0.1)]  )  ,    intersect(rownames(List1)[which(List1$p_val_adj<0.1)], rownames(D) ) )   ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(D$Cells)), Batch =  factor((D$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Lumenal_Epithelia",".pdf",sep="") ,width=10,height=10)

Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTC)[which(DELISTC$p_val_adj<0.1)] )  ,  intersect(rownames(List2)[which(List2$p_val_adj<0.1)]  , rownames(D) ) )   ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(D$Cells)), Batch =  factor((D$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Ciliated_Epithelia",".pdf",sep="") ,width=10,height=10)


Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTD)[which(DELISTD$p_val_adj<0.1)]  )  ,  intersect(  rownames(List3)[which(List3$p_val_adj<0.1)] , rownames(D) ) )  ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(D$Cells)), Batch =  factor((D$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_Glandular_Epithelia",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(D,assay = "RNA")) 
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[smallglandlist ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(D$Cells)), Batch =  factor((D$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_GlandularSmallList_Epithelia",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(Endo3,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(rownames(DELISTD)[which(DELISTD$p_val_adj<0.1)]  )  ,  List11 )  ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(Endo3$ID5)), Batch =  factor((Endo3$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_GlandularSecRef_Epithelia",".pdf",sep="") ,width=100,height=15)


smallglandlist <- c("VMP1","NEAT1","SOD2","KLF6","ERRFI1","TLE1","PNPLA8","SCGB1D2","GABRP","SCGB2A1","CLDN4","TM4SF1","UBE2B","FTH1","FXYD3","DEFB1","PHLDA2","ADM","PHLDA2","ADM","NDUFA1","SEC61G","NOP10","AMD1","MIF","RASD1")
Xh <- t(scale(t(X)))
Xp <- Xh[smallglandlist,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(Endo3$ID5)), Batch =  factor((Endo3$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_GlandularSecRef2_Epithelia",".pdf",sep="") ,width=100,height=15)




#Align all data to reference...


write( D$ID3 , file=paste(saveext,"/EpithelialIdentsGenomeoldano.csv",sep="") )
write( as.character(Idents(D)) , file=paste(saveext,"/EpithelialIdentsoldano.csv",sep="") )
write( colnames(D) , file=paste(saveext,"/EpithelialIdentsIDoldano.csv",sep="") )


Dnocil <- subset(D,idents = c("Ciliated","Stromal fibroblasts"), invert = TRUE)
Idents(Endo3) <- Endo3$ID5
Endo3 <- subset(Endo3,idents=c("SOX9","Lumenal","Glandular"))
mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo3,Dnocil), dims = 1:30, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_nocil.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_nocil.pdf",sep=""),width = 10, height = 8,p)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_VentoTormo_split_nocil.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_VentoTormo.pdf",sep=""),width = 10, height = 8,p)
saveRDS(mammal.combined,file=paste(saveext,"/Integrate_with_RVT_stroma_eithelia_nocil_nocil.rds",sep=""))
#saveRDS(mammal.combined,file=paste(saveext,"/Integrate_with_RVT_stroma_eithelia_nocil.rds",sep=""))


mammal.combined <- readRDS(file=paste(saveext,"/Integrate_with_RVT_stroma_eithelia_nocil.rds",sep=""))

mammal.combine <- mammal.combined #readRDS(file=paste(saveext,"/Integrate_with_RVT_stroma_eithelia_nocil.rds",sep=""))
cilcells = WhichCells(mammal.combined,idents="Ciliated")

DefaultAssay(mammal.combine) <- "RNA"
Idents(mammal.combine) <- mammal.combine$ID3
ours <- WhichCells(mammal.combine,idents=c("C1","C2","C3","C4","C6","C7","C8"))
mammal.combined <- subset(mammal.combine,idents=c("C1","C2","C3","C4","C6","C7","C8"))


newIDs <- as.character(Idents(mammal.combined))
newIDs[which(mammal.combined$ID3 %in% c("C1","C2","C3","C4","C6","C7","C8") )] <- "Lumenal"

SOX <- WhichCells(mammal.combine,expression=SOX9> log(1+1) )
LRG <- WhichCells(mammal.combine,expression=LGR5> log(1+1) )
SPP <- WhichCells(mammal.combine,expression=SPP1> log(1+1) )
SCGB <- WhichCells(mammal.combine,expression=SCGB2A2> log(1+1) )
mammal.combined2 <- mammal.combined
Idents(mammal.combined) <- newIDs
Idents(mammal.combined, cells= intersect(SPP,ours)) <- "Glandular"
Idents(mammal.combined, cells=intersect(SCGB,ours)) <- "Glandular"
Idents(mammal.combined, cells=intersect(SOX,ours)) <- "SOX9"
Idents(mammal.combined, cells=intersect(LRG,ours)) <- "LRG5"
Idents(mammal.combined, cells=cilcells) <- "Ciliated"


p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_nocil_newlab.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_nocil_newlab.pdf",sep=""),width = 10, height = 8,p)


mammal.combinedd <- mammal.combined
mammal.combinedd$Cells <- Idents(mammal.combinedd)
Idents(mammal.combinedd) <- mammal.combinedd$ID3
mammal.combinedd <- subset(mammal.combinedd,idents=c("C1","C2","C3","C4","C6","C7","C8"))
Idents(mammal.combinedd) <- mammal.combinedd$Cells

write( mammal.combined$ID3 , file=paste(saveext,"/EpithelialIdentsGenome.csv",sep="") )
write( as.character(Idents(mammal.combined)) , file=paste(saveext,"/EpithelialIdents.csv",sep="") )
write( colnames(mammal.combined) , file=paste(saveext,"/EpithelialIdentsID.csv",sep="") )


newIDs <- as.character(Idents(mammal.combined))
newIDs[which(mammal.combined$ID3 %in% c("C1","C2","C3","C4","C6","C7","C8") )] <- "Lumenal"
cilcells = WhichCells(mammal.combined,idents="Ciliated")

SOX <- WhichCells(mammal.combine,expression=SOX9> log(5) )
LRG <- WhichCells(mammal.combine,expression=LGR5>log(2)  )
SPP <- WhichCells(mammal.combine,expression=SPP1>log(5)  )
SCGB <- WhichCells(mammal.combine,expression=SCGB2A2> log(5)  )

Idents(mammal.combined) <- newIDs
Idents(mammal.combined, cells= intersect(SPP,ours) ) <- "Glandular"
Idents(mammal.combined, cells=intersect(SCGB,ours)) <- "Glandular"
Idents(mammal.combined, cells=intersect(SOX,ours)) <- "SOX9"
Idents(mammal.combined, cells=intersect(LRG,ours)) <- "LRG5"
table(Idents(mammal.combined))

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_nocil_newlab2.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_nocil_newlab2.pdf",sep=""),width = 10, height = 8,p)

write( mammal.combinedd$ID3 , file=paste(saveext,"/EpithelialIdentsGenome.csv",sep="") )
write( as.character(Idents(mammal.combinedd)) , file=paste(saveext,"/EpithelialIdents.csv",sep="") )
write( colnames(mammal.combinedd) , file=paste(saveext,"/EpithelialIdentsID.csv",sep="") )

mammal.combinedd <- mammal.combined
mammal.combinedd$Cells <- Idents(mammal.combinedd)
Idents(mammal.combinedd) <- mammal.combinedd$ID3
mammal.combinedd <- subset(mammal.combinedd,idents=c("C1","C2","C3","C4","C6","C7","C8"))
Idents(mammal.combinedd) <- mammal.combinedd$Cells

write( mammal.combined$ID3 , file=paste(saveext,"/EpithelialIdentsGenome.csv",sep="") )
write( as.character(Idents(mammal.combined)) , file=paste(saveext,"/EpithelialIdents.csv",sep="") )
write( colnames(mammal.combined) , file=paste(saveext,"/EpithelialIdentsID.csv",sep="") )


List1 <- FindMarkers(Endo3,ident.1 = "SOX9",ident.2 = c("Lumenal","Glandular"), test.use = "MAST", only.pos = TRUE )
List2 <- FindMarkers(Endo3,ident.1 = "Lumenal",ident.2 = c("SOX9","Glandular"), test.use = "MAST", only.pos = TRUE )
List3 <- FindMarkers(Endo3,ident.1 = "Glandular",ident.2 = c("Lumenal","SOX9"), test.use = "MAST", only.pos = TRUE )


newdata1 <- List1[order(List1$avg_log2FC,decreasing = TRUE),]
newdata2 <- List2[order(List3$avg_log2FC,decreasing = TRUE),]
newdata3 <- List3[order(List3$avg_log2FC,decreasing = TRUE),]

mammal.combine <- mammal.combined
mammal.combine$ID6 <- Idents(mammal.combine)
DefaultAssay(mammal.combine) <- "RNA"
Idents(mammal.combine) <- mammal.combine$ID3
ours <- WhichCells(mammal.combine,idents=c("C1","C2","C3","C4","C6","C7","C8"))
mammal.combine <- subset(mammal.combine,idents=c("C1","C2","C3","C4","C6","C7","C8"))
Idents(mammal.combine) <- mammal.combine$ID6

UberList <- c( rownames(newdata1)[1:100], rownames(newdata2)[1:100], rownames(newdata3)[1:100]  )
X <- (GetAssayData(mammal.combine,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(mammal.combine) )    ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combine$ID6)), Batch =  factor((mammal.combine$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/EndometrialMarkers",".pdf",sep="") ,width=55,height=38)

Idents(Endo3) <- paste(Endo3$ID3,Endo3$ID5,sep="_")
List1 <- FindMarkers(Endo3,ident.1 = "late-secretory_Glandular",ident.2 = c("early-secretory_Lumenal"), test.use = "MAST", only.pos = TRUE )
newdata1 <- List1[order(List1$avg_log2FC,decreasing = TRUE),]




UberList <- c( rownames(newdata1)[1:300]  )
mammal.combine2 <- subset(mammal.combine,idents="Glandular")
X <- (GetAssayData(mammal.combine2,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(mammal.combine2) )    ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combine2$ID6)), Batch =  factor((mammal.combine2$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/EndometrialMarkersGland",".pdf",sep="") ,width=55,height=38)




#Align all data to reference...
Idents(Endo3) <- Endo3$ID5
Endo3 <- subset(Endo3,idents=c("eS","dS","SOX9","Ciliated","Lumenal","Fibroblast C7","Glandular"))
mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo3,D), dims = 1:30, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo.pdf",sep=""),width = 10, height = 8,p)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_VentoTormo_split.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_VentoTormo.pdf",sep=""),width = 10, height = 8,p)


saveRDS(mammal.combined,file=paste(saveext,"/Integrate_with_RVT_stroma_eithelia.rds",sep=""))
mammal.combined<-readRDS(file=paste(saveext,"/Integrate_with_RVT_stroma_eithelia.rds",sep=""))
ours <- colnames(mammal.combined)[which(mammal.combined$Dataset=="10X Ours")]
cilcells = WhichCells(mammal.combined,idents="Ciliated")

SOX <- WhichCells(mammal.combined,expression=SOX9> log(1+1) )
LRG <- WhichCells(mammal.combined,expression=LGR5>log(1+1)  )
SPP <- WhichCells(mammal.combined,expression=SPP1>log(1+1)  )
SCGB <- WhichCells(mammal.combined,expression=SCGB2A2> log(1+1)  )
Idents(mammal.combined) <- newIDs
Idents(mammal.combined, cells= intersect(SPP,ours) ) <- "Glandular"
Idents(mammal.combined, cells=intersect(SCGB,ours)) <- "Glandular"
Idents(mammal.combined, cells=intersect(SOX,ours)) <- "SOX9"
Idents(mammal.combined, cells=intersect(LRG,ours)) <- "LRG5"

mammal.combined$ID0 <- Idents(mammal.combined)
Idents(mammal.combined) <- mammal.combined$Dataset
onlyOurs <- subset(mammal.combined,idents = "10X Ours")
Idents(onlyOurs) <- onlyOurs$ID0

mammal.combined<-readRDS(file=paste(saveext,"/Integrate_with_RVT_stroma_eithelia.rds",sep=""))
DefaultAssay(mammal.combined) <- "RNA"
newlab <- as.character(mammal.combined$ID3)
newlab[which(newlab%in%c("C1","C2","C3","C4","C6","C7","C8"))] <- "Ours"
mammal.combined$ID0 <- newlab

DefaultAssay(mammal.combined) <- "RNA"
p1<-FeaturePlot(mammal.combined, features = "RYBP", pt.size = 4, reduction = "pca", split.by = "ID0", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Endometrial_RYBP.pdf",sep=""),width = 30, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(mammal.combined, features = "PADI1", pt.size = 4, reduction = "pca", split.by = "ID0", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/ndometrial_PADI1.pdf",sep=""),width = 30, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(mammal.combined, features = "PADI2", pt.size = 4, reduction = "pca", split.by = "ID0", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/ndometrial_PADI2.pdf",sep=""),width = 30, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(mammal.combined, features = "PADI3", pt.size = 4, reduction = "pca", split.by = "ID0", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/ndometrial_PADI3.pdf",sep=""),width = 30, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(mammal.combined, features = "PADI4", pt.size = 4, reduction = "pca", split.by = "ID0", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/ndometrial_PADI4.pdf",sep=""),width = 30, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(mammal.combined, features = "PADI6", pt.size = 4, reduction = "pca", split.by = "ID0", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/ndometrial_PADI6.pdf",sep=""),width = 30, height = 8,p1,limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", label = TRUE, split.by = "ID0", repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Epithelial_stroma_all.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)

Idents(Endo3) <- Endo3$ID5
Endo4 <- subset(Endo3,idents=c("SOX9","Ciliated","Lumenal","Glandular","eS","dS","Fibroblast C7"))
D_s <- subset(D,idents=c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular","Stromal fibroblasts"))

strom <- WhichCells(D_s,idents="Stromal fibroblasts")
cil <- WhichCells(D_s,idents="Ciliated")

ours <- colnames(D_s)

newIDs <- Idents(D_s)
Idents(D_s) <- newIDs
Idents(D_s, cells= intersect(SPP,ours)) <- "Glandular"
Idents(D_s, cells=intersect(SCGB,ours)) <- "Glandular"
Idents(D_s, cells=intersect(SOX,ours)) <- "SOX9"
Idents(D_s, cells=intersect(LRG,ours)) <- "LRG5"

Idents(D_s, cells=strom) <- "Stroma"
Idents(D_s, cells=cil) <- "Ciliated"


newIDs <- Idents(mammal.combined)
Idents(mammal.combined) <- newIDs
Idents(mammal.combined, cells= intersect(SPP,ours)) <- "Glandular"
Idents(mammal.combined, cells=intersect(SCGB,ours)) <- "Glandular"
Idents(mammal.combined, cells=intersect(SOX,ours)) <- "SOX9"
Idents(mammal.combined, cells=intersect(LRG,ours)) <- "LRG5"

Idents(mammal.combined, cells=strom) <- "Stroma"
Idents(mammal.combined, cells=cil) <- "Ciliated"

mammal.combined$nn <- Idents(mammal.combined)
Idents(mammal.combined) <- mammal.combined$ID3
mammal.combinedours <- subset(mammal.combined,idents=c("C1","C2","C3","C4","C6","C7","C8"))
Idents(mammal.combinedours) <- mammal.combinedours$nn
p<-DimPlot(mammal.combinedours,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_epithelial_stroma_ours.pdf",sep=""),width = 10, height = 8,p)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", label = TRUE, split.by = "ID3", repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_epithelial_stroma_all.pdf",sep=""),width = 70, height = 8,p,limitsize = FALSE)


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


p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_VentoTormo_split_epithel.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_VentoTormo_epithel.pdf",sep=""),width = 10, height = 8,p)

saveRDS(mammal.combined,file=paste(saveext,"/Integrate_with_RVT_eithelia.rds",sep=""))

mammal.combined <- readRDS(file=paste(saveext,"/Integrate_with_RVT_eithelia.rds",sep=""))

Idents(Endo4) <- paste(Endo4$ID3,Idents(Endo4),sep="")



p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_test.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)


List1 <- FindMarkers(Endo4,ident.2 = "proliferativeLumenal",ident.1 = c("late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )
List2 <- FindMarkers(Endo4,ident.1 = "proliferativeLumenal",ident.2 = c("late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )

List3 <- FindMarkers(Endo4,ident.1 = "proliferativeSOX9",ident.2 = c("proliferativeCiliated","proliferativeGlandular","proliferativeLumenal"), test.use = "MAST", only.pos = TRUE )
List4 <- FindMarkers(Endo4,ident.1 = "proliferativeCiliated",ident.2 = c("proliferativeLumenal"), test.use = "MAST", only.pos = TRUE )

List5 <- FindMarkers(Endo4,ident.1 = "proliferativeGlandular",ident.2 = c("proliferativeLumenal"), test.use = "MAST", only.pos = TRUE )
List6 <- FindMarkers(Endo4,ident.2 = "proliferativeGlandular",ident.1 = c("proliferativeLumenal"), test.use = "MAST", only.pos = TRUE )


newdata1 <- List1[order(List1$avg_log2FC,decreasing = TRUE),]
newdata2 <- List2[order(List2$avg_log2FC,decreasing = TRUE),]
newdata3 <- List3[order(List3$avg_log2FC,decreasing = TRUE),]
newdata4 <- List4[order(List4$avg_log2FC,decreasing = TRUE),]
newdata5 <- List5[order(List5$avg_log2FC,decreasing = TRUE),]
newdata6 <- List6[order(List6$avg_log2FC,decreasing = TRUE),]

DefaultAssay(mammal.combined) <- "RNA"

mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata1)[1:100] ), name = "Gland")
p1<-FeaturePlot(mammal.combined, features = "Gland1", pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Gland_VE.pdf",sep=""),width = 70, height = 8,p1,limitsize = FALSE)

mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata2)[1:100] ), name = "Luminal")
p1<-FeaturePlot(mammal.combined, features = "Luminal1", pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Lum_VE.pdf",sep=""),width = 70, height = 8,p1,limitsize = FALSE)

mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata3)[1:100] ), name = "SOX9P_")
p1<-FeaturePlot(mammal.combined, features = "SOX9P_1", pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SOX9_VE.pdf",sep=""),width = 70, height = 8,p1,limitsize = FALSE)

mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata4)[1:100] ), name = "Cil")
p1<-FeaturePlot(mammal.combined, features = "Cil1", pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Cil_VE.pdf",sep=""),width = 70, height = 8,p1,limitsize = FALSE)



mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata5)[1:100] ), name = "eGland")
p1<-FeaturePlot(mammal.combined, features = "eGland1", pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/earlyGland_VE.pdf",sep=""),width = 70, height = 8,p1,limitsize = FALSE)


mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata6)[1:100] ), name = "eLuminal")
p1<-FeaturePlot(mammal.combined, features = "eLuminal1", pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/earlyLumenal_VE.pdf",sep=""),width = 70, height = 8,p1,limitsize = FALSE)


ID6 <- as.character(mammal.combined$ID3)
ID6[which(ID6%in%c("C1","C2","C3","C4","C5","C6","C7","C8"))] <- "Ours"
mammal.combined$ID6 <- ID6

Idents(mammal.combined) <- mammal.combined$ID6
mammal.combinedours <- subset(mammal.combined,idents="Ours")
p1<-FeaturePlot(mammal.combinedours, features = "Gland1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Gland_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(mammal.combinedours, features = "Luminal1", pt.size = 4, reduction = "pca",label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Lum_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(mammal.combinedours, features = "SOX9P_1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SOX9_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(mammal.combinedours, features = "Cil1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Cil_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(mammal.combinedours, features = "eGland1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/earlyGland_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(mammal.combinedours, features = "eGland1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/earlyLumenal_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)



Idents(mammal.combined) <- mammal.combined$ID3
onlyOurs <- subset(mammal.combined,idents=c("C1","C2","C3","C4","C6","C7","C8"))
notOurs1 <- subset(mammal.combined,idents=c("proliferative","late-secretory"))
Idents(notOurs1) <- paste(notOurs1$ID3,notOurs1$ID5,sep="_")
p <- FeatureScatter(onlyOurs, feature1 = "Gland1", feature2 = "Luminal1", )
ggsave(filename=paste(saveext,"/DimRed/Ours_Gland_vs_Luminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p <- FeatureScatter(notOurs1, feature1 = "Gland1", feature2 = "Luminal1", )
ggsave(filename=paste(saveext,"/DimRed/NotOurs1_Gland_vs_Luminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


p <- FeatureScatter(onlyOurs, feature1 = "Cil1", feature2 = "Luminal1", )
ggsave(filename=paste(saveext,"/DimRed/Ours_Cil_vs_Luminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p <- FeatureScatter(notOurs1, feature1 = "Cil1", feature2 = "Luminal1", )
ggsave(filename=paste(saveext,"/DimRed/NotOurs1_Cil_vs_Luminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)



p <- FeatureScatter(onlyOurs, feature1 = "SOX9P_1", feature2 = "Luminal1", )
ggsave(filename=paste(saveext,"/DimRed/Ours_SOX9_vs_Luminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p <- FeatureScatter(notOurs1, feature1 = "SOX9P_1", feature2 = "Luminal1", )
ggsave(filename=paste(saveext,"/DimRed/NotOurs1_SOX9_vs_Luminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p <- FeatureScatter(onlyOurs, feature1 = "eGland1", feature2 = "Luminal1", )
ggsave(filename=paste(saveext,"/DimRed/Ours_earlyland_vs_Luminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p <- FeatureScatter(notOurs1, feature1 = "eGland1", feature2 = "Luminal1", )
ggsave(filename=paste(saveext,"/DimRed/NotOurs1_earlyGland_vs_Luminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p <- FeatureScatter(onlyOurs, feature1 = "eGland1", feature2 = "eLuminal1", )
ggsave(filename=paste(saveext,"/DimRed/Ours_earlyland_vs_eLuminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p <- FeatureScatter(notOurs1, feature1 = "eGland1", feature2 = "eLuminal1", )
ggsave(filename=paste(saveext,"/DimRed/NotOurs1_earlyGland_vs_eLuminal.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


Idents(onlyOurs) <- onlyOurs$ID6

Idents(onlyOurs, cells= intersect(SPP,ours)) <- "Glandular"
Idents(onlyOurs, cells=intersect(SCGB,ours)) <- "Glandular"
Idents(onlyOurs, cells=intersect(SOX,ours)) <- "SOX9"
Idents(onlyOurs, cells=intersect(LRG,ours)) <- "LRG5"

Idents(onlyOurs, cells=strom) <- "Stroma"
Idents(onlyOurs, cells=cil) <- "Ciliated"
Idents(onlyOurs, WhichCells(onlyOurs,idents="Ours")) <- "Lumenal"

p <- VlnPlot(onlyOurs,features = c("Luminal1","eLuminal1","eGland1","Gland1","Cil1"))
ggsave(filename=paste(saveext,"/DimRed/OursSplitVlnbyModuleScore.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)






write.table(onlyOurs$eGland1,file=paste(saveext,'eGlandular.csv',sep=""),sep=",")
write.table(onlyOurs$Gland1,file=paste(saveext,'Glandular.csv',sep=""),sep=",")
write.table(onlyOurs$SOX9P_1,file=paste(saveext,'SOX9P_1.csv',sep=""),sep=",")
write.table(onlyOurs$Cil1,file=paste(saveext,'Cil1.csv',sep=""),sep=",")
write.table(onlyOurs$Luminal1,file=paste(saveext,'Luminal.csv',sep=""),sep=",")
write.table(onlyOurs$eLuminal1,file=paste(saveext,'eLuminal.csv',sep=""),sep=",")

write.table(onlyOurs$ID3,file=paste(saveext,'ID.csv',sep=""),sep=",")


D1 <- GetAssayData(onlyOurs, assay = "RNA")



write.table(t(D1[ c(rownames(newdata1)[1:100] , rownames(newdata2)[1:100], rownames(newdata3)[1:100], rownames(newdata4)[1:100], rownames(newdata5)[1:100]   )   ,]),file=paste(saveext,'Markers.csv',sep=""),sep=",")


write.table(List1,file=paste(saveext,'Glandular.csv',sep=""),sep=",")
write.table(List2,file=paste(saveext,'Luminal.csv',sep=""),sep=",")
write.table(List3,file=paste(saveext,'Ciliated.csv',sep=""),sep=",")
write.table(List4,file=paste(saveext,'SOX9.csv',sep=""),sep=",")

#PAEP
#GPX3
#CXCL14
#NUPR1
#CRYAB
#RIMKLB
#IGFBP7
#VCAN
#S100A1
#SLC18A2

DefaultAssay(onlyOurs) <- "RNA"
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "PAEP", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_PAEP.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "GPX3",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_GPX3.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "CXCL14",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_CXCL14.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "NUPR1",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_NUPR1.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "CRYAB",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_CRYAB.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "RIMKLB",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RIMKLB.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "IGFBP7",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_IGFBP7.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "VCAN",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_VCAN.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "S100A1",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_S100A1.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "SLC18A2",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_SLC18A2.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)

DefaultAssay(notOurs1) <- "RNA"
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "PAEP", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_PAEP_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "GPX3",   split.by = "ID3",cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_GPX3_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "CXCL14",   split.by = "ID3",cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_CXCL14_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "NUPR1",   split.by = "ID3",cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_NUPR1_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "CRYAB",   split.by = "ID3",cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_CRYAB_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "RIMKLB",   split.by = "ID3",cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RIMKLB_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "IGFBP7",  split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_IGFBP7_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "VCAN",   split.by = "ID3",cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_VCAN_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "S100A1",  split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_S100A1_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "SLC18A2",  split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_SLC18A2_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)

p<-FeaturePlot(notOurs1,  reduction = "pca", features = "PAEP",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_PAEP_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "GPX3",   cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_GPX3_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "CXCL14",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_CXCL14_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "NUPR1",   cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_NUPR1_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "CRYAB",   cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_CRYAB_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "RIMKLB",   cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RIMKLB_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "IGFBP7",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_IGFBP7_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "VCAN",   cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_VCAN_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "S100A1",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_S100A1_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "SLC18A2",   cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_SLC18A2_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)


#TPPP3
#C20orf85
#IGFBP7
#CAPS
#C1orf194
#C9orf24
#RSPH1
#FAM183A
#MORN2
#AGR3
#C11orf88
#PIFO

p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "TPPP3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_TPPP3.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "C20orf85", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_C20orf85.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "CAPS", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_CAPS.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "C1orf194", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_C1orf194.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "C9orf24", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_C9orf24.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "RSPH1", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RSPH1.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "FAM183A", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_FAM183A.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "MORN2", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_MORN2.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "AGR3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_AGR3.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "C11orf88", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_C11orf88.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "PIFO",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_PIFO.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)


p<-FeaturePlot(notOurs1,  reduction = "pca", features = "TPPP3",  split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_TPPP3_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "C20orf85",  split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_C20orf85_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "CAPS",  split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_CAPS_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "C1orf194",  split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_C1orf194_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "C9orf24",  split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_C9orf24_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "RSPH1",  split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RSPH1_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "FAM183A",  split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_FAM183A_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "MORN2",  split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_MORN2_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "AGR3",  split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_AGR3_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "C11orf88",  split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_C11orf88_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "PIFO",  split.by = "ID3",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_PIFO_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)

p<-FeaturePlot(notOurs1,  reduction = "pca", features = "TPPP3",   cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_TPPP3_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "C20orf85",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_C20orf85_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "CAPS",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_CAPS_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "C1orf194", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_C1orf194_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "C9orf24",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_C9orf24_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "RSPH1",   cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RSPH1_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "FAM183A",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_FAM183A_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "MORN2",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_MORN2_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "AGR3",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_AGR3_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "C11orf88",   cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_C11orf88_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "PIFO",   cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_PIFO_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)


#MMP7
#VIM
#WFDC2
#RPS20
#ANXA1
#RPS17
#RPL13A
#PLAU
#RPL37A
#RPL23
#RPL31
#RPL7

p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "MMP7", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_MMP7.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "VIM", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_VIM.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "WFDC2",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_WFDC2.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "RPS20",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RPS20.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "ANXA1",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_ANXA1.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "RPS17",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RPS17.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "RPL13A",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RPL13A.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "PLAU",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_PLAU.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "RPL37A",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RPL37A.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "RPL23",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RPL23.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)


p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "SOX9",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_SOX9.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)

p<-FeaturePlot(notOurs1,  reduction = "pca", features = "MMP7",  split.by = "ID3",cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_MMP7_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "VIM",  split.by = "ID3",cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_VIM_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "WFDC2",   split.by = "ID3",cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_WFDC2_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "RPS20",   split.by = "ID3",cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RPS20_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "ANXA1",   split.by = "ID3",cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_ANXA1_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "RPS17",   split.by = "ID3",cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RPS17_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "RPL13A",  split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RPL13A_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "PLAU",   split.by = "ID3",cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_PLAU_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "RPL37A",  split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RPL37A_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "RPL23",   split.by = "ID3",cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RPL23_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)

p<-FeaturePlot(notOurs1,  reduction = "pca", features = "MMP7",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_MMP7_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "VIM", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_VIM_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "WFDC2",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_WFDC2_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "RPS20",   cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RPS20_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "ANXA1",   cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_ANXA1_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "RPS17",   cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RPS17_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "RPL13A",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RPL13A_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "PLAU",   cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_PLAU_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "RPL37A",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RPL37A_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs1,  reduction = "pca", features = "RPL23",   cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_PCA_RPL23_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)

p1<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Epithelial_all.pdf",sep=""),width = 70, height = 8,p1,limitsize = FALSE)


newIDs <- as.character(Idents(mammal.combined))
newIDs[which(mammal.combined$ID3 %in% c("C1","C2","C3","C4","C6","C7","C8") )] <- "Lumenal"
SOX <- WhichCells(mammal.combine,expression=SOX9> log(1) )
LRG <- WhichCells(mammal.combine,expression=LGR5>log(1)  )
SPP <- WhichCells(mammal.combine,expression=SPP1>log(1)  )
SCGB <- WhichCells(mammal.combine,expression=SCGB2A2> log(1)  )
Idents(mammal.combined) <- newIDs
Idents(mammal.combined,cells=SOX) <- "SOX"
Idents(mammal.combined,cells=LRG) <- "LRG"
Idents(mammal.combined,cells=SPP) <- "Glandular"
Idents(mammal.combined,cells=SCGB) <- "Glandular"
Idents(mammal.combined,cells=cil) <- "Ciliated"
mammal.combined$nn <- Idents(mammal.combined)
Idents(mammal.combined) <- mammal.combined$ID3
mammal.combinedours <- subset(mammal.combined,idents=c("C1","C2","C3","C4","C6","C7","C8"))
Idents(mammal.combinedours) <- mammal.combinedours$nn
p1<-DimPlot(mammal.combinedours,  pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Epithelial_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

DefaultAssay(mammal.combined) <- "RNA"

mammal.combined$nn <- Idents(mammal.combined)
Idents(mammal.combined) <- mammal.combined$ID3
onlyOurs <- subset(mammal.combined,idents=c("C1","C2","C3","C4","C7","C6","C8"))
Idents(onlyOurs) <- onlyOurs$nn
DefaultAssay(onlyOurs) <- "RNA"


onlyOurs$ID4 <- NULL
p1<-VlnPlot(onlyOurs,  idents = c("Ciliated","Lumenal","Glandular"),features = rownames(newdata1)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_set1.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(onlyOurs,   idents = c("Ciliated","Lumenal","Glandular"),features = rownames(newdata2)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_set2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(onlyOurs,   idents = c("Ciliated","Lumenal","Glandular"),features = rownames(newdata3)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_set3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(onlyOurs,   idents = c("Ciliated","Lumenal","Glandular"),features = rownames(newdata4)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_set4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

onlyOurs$ID1 <- NULL

p1<-VlnPlot(onlyOurs,  idents = c("Ciliated","Lumenal","Glandular"),features = rownames(newdata1)[41:80], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_set1_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(onlyOurs,   idents = c("Ciliated","Lumenal","Glandular"),features = rownames(newdata2)[41:80],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_set2_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(onlyOurs,   idents = c("Ciliated","Lumenal","Glandular"),features = rownames(newdata3)[41:80],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_set3_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(onlyOurs,   idents = c("Ciliated","Lumenal","Glandular"),features = rownames(newdata4)[41:80],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_set4_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)




DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- FindClusters(mammal.combined, resolution = 0.1)
p1<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Epithelial_Cl.pdf",sep=""),width = 70, height = 8,p1,limitsize = FALSE)
mammal.combined$newCl <- Idents(mammal.combined)

Idents(mammal.combined) <- mammal.combined$ID3
clusteredmammal.combinednotours1 <- subset(mammal.combined,idents=c("proliferative"))
clusteredmammal.combinednotours2 <- subset(mammal.combined,idents=c("early-secretory"))
clusteredmammal.combinednotours3 <- subset(mammal.combined,idents=c("late-secretory"))
clusteredmammal.combinednotours4 <- subset(mammal.combined,idents=c("C1","C2","C3","C4","C6","C7","C8"))

data <- data.frame(x=as.factor(clusteredmammal.combinednotours1$ID5),y=as.factor(clusteredmammal.combinednotours1$newCl) )
p0 <- ggplot(data, aes(fill=y, x=x, y=y)) + geom_bar(position="fill", stat="identity") + theme_bw() + coord_cartesian(ylim = c(0, 1))
ggsave(filename=paste(saveext,"/BarPlot_EpithelialClbyType_Prolf.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p0)

data <- data.frame(x=as.factor(clusteredmammal.combinednotours2$ID5),y=as.factor(clusteredmammal.combinednotours2$newCl) )
p0 <- ggplot(data, aes(fill=y, x=x, y=y)) + geom_bar(position="fill", stat="identity") + theme_bw() + coord_cartesian(ylim = c(0, 1))
ggsave(filename=paste(saveext,"/BarPlot_EpithelialClbyType_eSec.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p0)

data <- data.frame(x=as.factor(clusteredmammal.combinednotours3$ID5),y=as.factor(clusteredmammal.combinednotours3$newCl) )
p0 <- ggplot(data, aes(fill=y, x=x, y=y)) + geom_bar(position="fill", stat="identity") + theme_bw() + coord_cartesian(ylim = c(0, 1))
ggsave(filename=paste(saveext,"/BarPlot_EpithelialClbyType_lSec.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p0)

data <- data.frame(x=as.factor(clusteredmammal.combinednotours4$ID5),y=as.factor(clusteredmammal.combinednotours4$newCl), z=as.factor(clusteredmammal.combinednotours4$ID3) )
p0 <- ggplot(data, aes(fill=y, x=z, y=y)) + geom_bar(position="fill", stat="identity") + theme_bw() + coord_cartesian(ylim = c(0, 1))
ggsave(filename=paste(saveext,"/BarPlot_EpithelialClbyType_ours.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p0)

clusteredmammal.combinednotours5 <- subset(mammal.combined,idents=c("proliferative","early-secretory","late-secretory"))

data <- data.frame(x=as.factor(paste(clusteredmammal.combinednotours5$ID3,clusteredmammal.combinednotours5$ID5,sep="")),y=as.factor(clusteredmammal.combinednotours5$newCl) )
p0 <- ggplot(data, aes(fill=x, x=y, y=x)) + geom_bar(position="fill", stat="identity") + theme_bw() + coord_cartesian(ylim = c(0, 1))
ggsave(filename=paste(saveext,"/BarPlot_EpithelialClbyType_ByCluster.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p0)

data <- data.frame(x=as.factor(clusteredmammal.combinednotours4$ID3),y=as.factor(clusteredmammal.combinednotours4$newCl), z=as.factor(clusteredmammal.combinednotours4$ID3) )
p0 <- ggplot(data, aes(fill=y, x=z, y=y)) + geom_bar(position="fill", stat="identity") + theme_bw() + coord_cartesian(ylim = c(0, 1))
ggsave(filename=paste(saveext,"/BarPlot_EpithelialClbyType_ByClusterOurs.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p0)


uID <- as.character(mammal.combined$newCl) #paste(mammal.combined$newCl,mammal.combined$newCl,sep="_")
uID[which(mammal.combined$ID3%in%c("proliferative","early-secretory","late-secretory"))] <- paste(mammal.combined$ID3[which(mammal.combined$ID3%in%c("proliferative","early-secretory","late-secretory"))],mammal.combined$ID5[which(mammal.combined$ID3%in%c("proliferative","early-secretory","late-secretory"))], sep ="_")
Idents(mammal.combined) <- uID

AvExp <- AverageExpression(mammal.combined)
C1 <- cor(AvExp$integrated,method = "pearson")
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(as.data.frame(C1),color =  redblue1(50), fontsize = 12,  border_color = NA,cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/ClusterCorrelations",".pdf",sep="") ,width=55,height=38)


ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T) +
  ggtitle("Studying 4 species..") +
  theme_ipsum() +
  xlab("")


write.table(Idents(mammal.combined),file=paste(saveext,"/EpithelialClusters.csv",sep=""))


AVE2 <- data.frame(ID=mammal.combined$ID3,Cell=mammal.combined$ID5,eGl=mammal.combined$eGland1,Gl=mammal.combined$Gland1,Cil=mammal.combined$Cil1,SO=mammal.combined$SOX9P_1)

#colnames(AVE2) <- c("Prolif","Gland","Lumenal")
#AVE2$Sum <- log2(AvExp3[rownames(AVE2),"Sum"]+1)
#AVE2$Delta <- AvExp3[rownames(AVE2),"Delta"]
#AVE2$MAX <- AvExp3[rownames(AVE2),"MAX"]
#AVE2$Sum[which(is.na(AVE2$Sum)==1)] <- 0
#AVE2$Delta[which(is.na(AVE2$Delta)==1)] <- 0
#AVE2$Delta <- 5*(AVE2$Delta/max(AVE2$Delta))
#AVE2$Sum <- 5*(AVE2$Sum/max(AVE2$Sum))
#AVE2<- AVE2[which(AVE2$MAX>1),]
#AVE2 <- AVE2[intersect(list2,rownames(AVE2)),]
#library(ggtern)
#p1 <- ggtern(data=AVE2, aes(x =Prolif,z=Gland,y=Lumenal)) + geom_point(size=AVE2$MAX, color = "black", position= position_jitter_tern(x=0.01, y=0.01, z=0.01))
#p1 <- p1+ theme_showgrid()
#p1 <- p1 + annotate(geom  = 'text', x = (AVE2$Prolif)/(AVE2$Prolif+AVE2$Gland+AVE2$Lumenal),
#                    z= (AVE2$Gland)/(AVE2$Prolif+AVE2$Gland+AVE2$Lumenal),
#                    y= (AVE2$Lumenal)/(AVE2$Prolif+AVE2$Gland+AVE2$Lumenal),label = rownames(AVE2),color = c("black"))
#ggsave(filename=paste(saveext,"GGTern_Markers_ProlifGlandLumenal.pdf",sep=""),width = 20, height = 20, plot = p1)


Endo3 <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Endometrium/other_endometrium_data.rds")
Idents(Endo3) <- Endo3$ID3
Endo3 <- subset(Endo3,idents=c("early-secretory","proliferative","late-secretory") )
Idents(Endo3) <- Endo3$ID5
Endo4 <- subset(Endo3,idents=c("Fibroblast C7","dS","eS"))
D_s <- subset(D,idents=c("Stromal fibroblasts"))

mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo4,D_s), dims = 1:30, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_stroma.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_stroma.pdf",sep=""),width = 10, height = 8,p)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_stroma_1_3.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)


p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_VentoTormo_stroma.pdf",sep=""),width = 10, height = 8,p)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_VentoTormo_split_stroma.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)

saveRDS(mammal.combined,file=paste(saveext,"/Integrate_with_RVT_stroma.rds",sep=""))
mammal.combined <- readRDS(file=paste(saveext,"/Integrate_with_RVT_stroma.rds",sep=""))


tempID <- mammal.combined$ID5
tempID[which(is.na(tempID))] <- "Stromal"
Idents(mammal.combined) <- tempID
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_stroma_1_3.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)


p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_stromatest.pdf",sep=""),width = 10, height = 8,p)

Idents(Endo4) <- paste(Endo4$ID3,Idents(Endo4),sep="")

List1 <- FindMarkers(Endo4,ident.2 = "proliferativeeS",ident.1 = c("late-secretorydS"), test.use = "MAST", only.pos = TRUE )
List2 <- FindMarkers(Endo4,ident.1 = "proliferativeeS",ident.2 = c("late-secretorydS"), test.use = "MAST", only.pos = TRUE )
List3 <- FindMarkers(Endo4,ident.1 = "proliferativeFibroblast C7",ident.2 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List4 <- FindMarkers(Endo4,ident.1 = "late-secretoryFibroblast C7",ident.2 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List5 <- FindMarkers(Endo4,ident.2 = "proliferativeFibroblast C7",ident.1 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List6 <- FindMarkers(Endo4,ident.2 = "late-secretoryFibroblast C7",ident.1 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )

Idents(Endo3) <- paste(Endo3$ID3,Idents(Endo3),sep="")

Iunique(Idents(Endo3))


saveRDS(List1,file="~/Desktop/ProlifStroma_vs_SecdStroma.rds")
saveRDS(List2,file="~/Desktop/SecdStroma_vs_ProlifStroma.rds")


saveRDS(List3,file="~/Desktop/ProlfiFib_vs_ProlifStroma.rds")
saveRDS(List4,file="~/Desktop/ProlifStroma_vs_ProlfiFib.rds")

saveRDS(List5,file="~/Desktop/SecFib_vs_SecStroma.rds")
saveRDS(List6,file="~/Desktop/SecStroma_vs_SecFib.rds")


#There were 49.36% and 34.1% stromal cells in the control and endometriosis groups, respectively. The proportion of endothelial and immune cells was 8.43% and 15.43% in the control group, respectively, and 6.4% and 26.3% in the endometriosis group, respectively.


newdata1 <- List1[order(List1$avg_log2FC,decreasing = TRUE),]
newdata2 <- List2[order(List3$avg_log2FC,decreasing = TRUE),]
newdata3 <- List3[order(List3$avg_log2FC,decreasing = TRUE),]
newdata4 <- List4[order(List4$avg_log2FC,decreasing = TRUE),]
newdata5 <- List5[order(List5$avg_log2FC,decreasing = TRUE),]
newdata6 <- List6[order(List6$avg_log2FC,decreasing = TRUE),]


Idents(mammal.combined) <- mammal.combined$ID3
onlyOurs <- subset(mammal.combined,idents=c("C1","C2","C3","C4","C7","C6","C8"))
notOurs <- subset(mammal.combined,idents=c("proliferative","late-secretory"))


DefaultAssay(onlyOurs) <- "RNA"


p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "PAEP",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_PAEP.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "SOD2",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_SOD2.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "CFD", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_CFD.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "TIMP3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_TIMP3.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "JUND",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_JUND.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "HSPA1A",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_HSPA1A.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "CD81", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_CD81.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "C11orf96",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_C11orf96.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "S100A6",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_S100A6.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "SERPING1",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_SERPING1.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)


p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "BEX2",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_BEX2.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "ADAM12", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_ADAM12.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "AGPAT2", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_AGPAT21.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "ANKRD29",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_ANKRD29.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "ATP5L", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_ATP5L.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "ANXA5",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_AANXA5.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(onlyOurs,  reduction = "pca", features = "C14orf2", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_C14orf2.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)


DefaultAssay(notOurs) <- "RNA"
p<-FeaturePlot(notOurs,  reduction = "pca", features = "PAEP", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_PAEP_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "SOD2", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_SOD2_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "CFD", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_CFD_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "TIMP3", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_TIMP3_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "JUND", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_JUND_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "HSPA1A", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_HSPA1A_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "CD81", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_CD81_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "C11orf96", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_C11orf96_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "S100A6", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_S100A6_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "SERPING1", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_SERPING1_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "BEX2", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_BEX2_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "ADAM12", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_ADAM12_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "AGPAT2", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_AGPAT21_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "ANKRD29", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_ANKRD29_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "ATP5L", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_ATP5L_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "ANXA5", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_AANXA5_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "C14orf2", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_C14orf2_notours.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)


p<-FeaturePlot(notOurs,  reduction = "pca", features = "PAEP",cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_PAEP_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "SOD2", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_SOD2_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "CFD", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_CFD_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "TIMP3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_TIMP3_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "JUND", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_JUND_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "HSPA1A", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_HSPA1A_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "CD81",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_CD81_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "C11orf96",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_C11orf96_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "S100A6", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_S100A6_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "SERPING1",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_SERPING1_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "BEX2",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_BEX2_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "ADAM12",cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_ADAM12_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "AGPAT2", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_AGPAT21_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "ANKRD29",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_ANKRD29_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "ATP5L",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_ATP5L_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "ANXA5",  cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_AANXA5_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(notOurs,  reduction = "pca", features = "C14orf2", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_C14orf2_notourscb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)



DefaultAssay(mammal.combined) <- "RNA"
Idents(mammal.combined) <- mammal.combined$ID3
mammal.combined3 <- subset(mammal.combined,idents=c("C1","C2","C3","C4","C6","C7","C8"))
mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata1)[1:100] ), name = "dS")
p1<-FeaturePlot(mammal.combined3, features = "dS1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Stroma_dS_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata2)[1:100] ), name = "prolifS")
p1<-FeaturePlot(mammal.combined3, features = "prolifS1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Stroma_prolifS_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata3)[1:100] ), name = "prolifC7")
p1<-FeaturePlot(mammal.combined3, features = "prolifC71", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Stroma_prolifFibC7_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata4)[1:100] ), name = "lateC7")
p1<-FeaturePlot(mammal.combined3, features = "lateC71", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Stroma_lateFibC7_ours.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)


p <- FeatureScatter(mammal.combined3, feature1 = "dS1", feature2 = "prolifS1", )
ggsave(filename=paste(saveext,"/DimRed/Ours_Stroma.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

saveRDS(mammal.combined,file=paste(saveext,"/Integrate_with_RVT_stroma_withModuleScores.rds",sep=""))


mammal.combined <- readRDS(file=paste(saveext,"/Integrate_with_RVT_stroma_withModuleScores.rds",sep=""))



mammal.combined4 <- subset(mammal.combined,idents=c("proliferative","early-secretory","late-secretory"))
mammal.combined4 <- AddModuleScore(mammal.combined4, features = list( rownames(newdata1)[1:100] ), name = "dS")
mammal.combined4 <- AddModuleScore(mammal.combined4, features = list( rownames(newdata2)[1:100] ), name = "prolifS")
mammal.combined4 <- AddModuleScore(mammal.combined4, features = list( rownames(newdata3)[1:100] ), name = "prolifC7")
mammal.combined4 <- AddModuleScore(mammal.combined4, features = list( rownames(newdata4)[1:100] ), name = "lateC7")


#mammal.combined3 <- subset(mammal.combined,idents=c("proliferative","early-secretory","late-secretory"))
mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata1)[1:100] ), name = "dS")
mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata2)[1:100] ), name = "prolifS")
mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata3)[1:100] ), name = "prolifC7")
mammal.combined3 <- AddModuleScore(mammal.combined3, features = list( rownames(newdata4)[1:100] ), name = "lateC7")

#mammal.combined3 <- subset(mammal.combined,idents=c("proliferative","early-secretory","late-secretory"))
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata1)[1:100] ), name = "dS")
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata2)[1:100] ), name = "prolifS")
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata3)[1:100] ), name = "prolifC7")
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata4)[1:100] ), name = "lateC7")


Idents(mammal.combined) <- mammal.combined3$ID3
onlyOurs <- subset(mammal.combined,idents=c("C1","C2","C3","C4","C7","C6","C8"))
uID <- as.character( Idents(onlyOurs) )
uID[1:length(uID)] <- "dS"
uID[which(onlyOurs$prolifS1 > onlyOurs$dS1)] <- "pS"
Idents(onlyOurs) <- uID
DefaultAssay(onlyOurs) <- "RNA"

uID <- as.character( Idents(mammal.combined) )
uID[1:length(uID)] <- "dS"
uID[which(mammal.combined$prolifS1 > mammal.combined$dS1)] <- "pS"
uID[which(mammal.combined$ID3%in%c("proliferative","early-secretory","late-secretory"))] <- "Reference"
Idents(mammal.combined) <- uID

DefaultAssay(mammal.combined) <- "RNA"

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_stroma_samescale.pdf",sep=""),width = 10, height = 8,p, limitsize = FALSE)

mammal.combined$Delta <- mammal.combined$dS1 - mammal.combined$prolifS1
twosplit <- as.character(mammal.combined$ID3)
twosplit[which(twosplit%in%c("C1","C2","C3","C4","C6","C7","C8"))] <- "ours"
twosplit[which(twosplit%in%c("proliferative","early-secretory","late-secretory"))] <- "Ref"
mammal.combined$twosplit <- twosplit
Idents(mammal.combined) <- mammal.combined$Delta

p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "Delta", split.by = "twosplit", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_Delta_SvdS.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "Delta", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_Delta_SvdScb.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)

Idents(mammal.combined) <- mammal.combined$twosplit
mammal.combinedsss <- subset(mammal.combined,idents=c("ours"))
p<-FeaturePlot(mammal.combinedsss,  reduction = "pca", features = "Delta", split.by = "twosplit", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_Delta_SvdS2.pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combinedsss,  reduction = "pca", features = "Delta", cols =  c("blue", "red"), pt.size = 2 ) # + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_Delta_SvdScb2.pdf",sep=""),width = 10, height = 8, limitsize = FALSE,p)


uID <- as.character( Idents(mammal.combinedsss) )
uID[1:length(uID)] <- "dS"
uID[which(mammal.combinedsss$prolifS1 > mammal.combinedsss$dS1)] <- "pS"
Idents(mammal.combinedsss) <- uID
DefaultAssay(mammal.combinedsss) <- "RNA"

mammal.combinedsss$ID4 <- NULL
p<-VlnPlot(mammal.combinedsss,  features = rownames(newdata1)[1:40], cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_decidualised.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p)


p<-VlnPlot(mammal.combinedsss,  features = rownames(newdata2)[1:40], cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_notdecidualised.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p)


p <- FeatureScatter(mammal.combined4, feature1 = "dS1", feature2 = "prolifS1", )
ggsave(filename=paste(saveext,"/DimRed/Ref_Stroma.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

write.table(mammal.combined3$dS1,file=paste(saveext,"/Ours_Stroma_Dec.csv",sep=""))
write.table(mammal.combined3$prolifS1,file=paste(saveext,"/Ours_Stroma_Prolif.csv",sep=""))
write.table(colnames(mammal.combined3),file=paste(saveext,"/StromaID.csv",sep=""))
write.table(mammal.combined3$ID3,file=paste(saveext,"/StromaGenome.csv",sep=""))


UberList <- c(rownames(List1)[1:100],rownames(List2)[1:100])
X <- (GetAssayData(D_s,assay = "RNA")) 
Idents(D) <- D$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[ intersect( UberList, rownames(D) )    ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(D$Cells)), Batch =  factor((D$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Fullheamtap_Decidualisaiton",".pdf",sep="") ,width=55,height=38)

mammal.combined <- readRDS(file=paste(saveext,"/Integrate_with_RVT_eithelia.rds",sep=""))
write.table(List1,file=paste(saveext,'DecidualisedStroma.csv',sep=""),sep=",")
write.table(List2,file=paste(saveext,'UnDecidualisedStroma.csv',sep=""),sep=",")

DefaultAssay(mammal.combined) <- "RNA"
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "PAEP", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_PAEP.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "SOD2", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_SOD2.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "CFD", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_CFD.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "TIMP3", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_TIMP3.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "JUND", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_JUND.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "HSPA1A", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_HSPA1A.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "CD81", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_CD81.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "C11orf96", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_C11orf96.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "S100A6", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_S100A6.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)

PAEP
SOD2
C11orf96
CFD
TIMP3
JUND
HSPA1A
SERPING1
CD81
S100A6



p<-DimPlot(onlyOurs,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_stroma_onlyours.pdf",sep=""),width = 10, height = 8,p, limitsize = FALSE)


p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "SFRP4", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_SFRP4.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "MMP11", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_MMP11.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "CRABP2", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_CRABP2.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "ECM1", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_ECM1.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "PGRMC1", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_PGRMC1.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "TGFBI", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_TGFBI.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)


p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "RBP7", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_RBP7.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "SFRP1", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_SFRP1.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "HIST1H4C", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_HIST1H4C.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)

p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "PAMR1", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_PCA_PAMR1.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)

SFRP4
MMP11
CRABP2
ECM1
PGRMC1
TGFBI
RBP7
SFRP1
HIST1H4C
PAMR1

DefaultAssay(mammal.combined) <- "RNA"
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "PIFO", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_PIFO.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)

p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "FOXJ1", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_FOXJ1.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)

p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "TPPP3", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_TPPP3.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)

p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "TP73", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_TP73.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)


p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "MUC12", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_MUC12.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)

p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "HES6", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_HES6.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)

p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "LAMP3", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_CD208.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)


p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "CCNO", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_CCNO.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)

p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "MMP7", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_MMP7.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)

p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "ESR1", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_ESR1.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)


p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "WNT7A", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_WNT7A.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)


p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "ESR1", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_ESR1.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)


p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "MKI67", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_MKI67.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)


p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "KRT17", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_KRT17.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)





p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "ESR1", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_ESR1.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)



p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "IHH", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_IHH.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)



p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "PAEP", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_PAEP.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)



p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "PTGS1", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_PTGS1.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)


p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "KRT5", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_KRT5.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)


p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "CXCL8", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_CXCL8.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)



p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "HES1", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_HES1.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)


p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "HEY1", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_HEY1.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)


p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "CSRNP1", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_CSRNP1.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)



p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "FOXO1", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_FOXO1.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)




p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "SOX9", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_SOX9.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)



p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "LGR5", split.by = "ID3", cols =  c("lightgrey", "red"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_LGR5.pdf",sep=""),width = 70, height = 8, limitsize = FALSE,p)

DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- FindClusters(mammal.combined, resolution = .05)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_Cl.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)

Idents(mammal.combined) <- mammal.combined$Cl
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_VentoTormo_split_epithel_Cl.pdf",sep=""),width = 70, height = 8,p, limitsize = FALSE)



#We are here ...
DefaultAssay(mammal.combined) <- "RNA"

mammal.combined$Cl <- Idents(mammal.combined)
Idents(mammal.combined) <- mammal.combined$ID3

mammal.combined1 <- subset(mammal.combined,idents="C1")
Idents(mammal.combined1) <- mammal.combined1$Cl
p <- VlnPlot(mammal.combined1, features = c("PIFO","FOXJ1","TPPP3","TP73","MUC12","HES6","LAMP3","CCNO","MMP7","WNT7A","MKI67","KRT17","ESR1","IHH","PAEP","PTGS1","KRT5","CXCL8","HES1", "HEY1","CSRNP1","SCGB2A2","FOXO1","SOX9","LGR5") )
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_VlnByCl_Batch1.pdf",sep=""),width = 30, height = 30,p, limitsize = FALSE)

mammal.combined1 <- subset(mammal.combined,idents="C2")
Idents(mammal.combined1) <- mammal.combined1$Cl
p <- VlnPlot(mammal.combined1, features = c("PIFO","FOXJ1","TPPP3","TP73","MUC12","HES6","LAMP3","CCNO","MMP7","WNT7A","MKI67","KRT17","ESR1","IHH","PAEP","PTGS1","KRT5","CXCL8","HES1", "HEY1","CSRNP1","SCGB2A2","FOXO1","SOX9","LGR5") )
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_VlnByCl_Batch2.pdf",sep=""),width = 30, height = 30,p, limitsize = FALSE)


mammal.combined1 <- subset(mammal.combined,idents=c("C3","C4") )
Idents(mammal.combined1) <- mammal.combined1$Cl
p <- VlnPlot(mammal.combined1, features = c("PIFO","FOXJ1","TPPP3","TP73","MUC12","HES6","LAMP3","CCNO","MMP7","WNT7A","MKI67","KRT17","ESR1","IHH","PAEP","PTGS1","KRT5","CXCL8","HES1", "HEY1","CSRNP1","SCGB2A2","FOXO1","SOX9","LGR5") )
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_VlnByCl_Batch3.pdf",sep=""),width = 30, height = 30,p, limitsize = FALSE)


mammal.combined1 <- subset(mammal.combined,idents="C6")
Idents(mammal.combined1) <- mammal.combined1$Cl
p <- VlnPlot(mammal.combined1, features = c("PIFO","FOXJ1","TPPP3","TP73","MUC12","HES6","LAMP3","CCNO","MMP7","WNT7A","MKI67","KRT17","ESR1","IHH","PAEP","PTGS1","KRT5","CXCL8","HES1", "HEY1","CSRNP1","SCGB2A2","FOXO1","SOX9","LGR5") )
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_VlnByCl_Batch6.pdf",sep=""),width = 30, height = 30,p, limitsize = FALSE)

mammal.combined1 <- subset(mammal.combined,idents="C7")
Idents(mammal.combined1) <- mammal.combined1$Cl
p <- VlnPlot(mammal.combined1, features = c("PIFO","FOXJ1","TPPP3","TP73","MUC12","HES6","LAMP3","CCNO","MMP7","WNT7A","MKI67","KRT17","ESR1","IHH","PAEP","PTGS1","KRT5","CXCL8","HES1", "HEY1","CSRNP1","SCGB2A2","FOXO1","SOX9","LGR5") )
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_VlnByCl_Batch7.pdf",sep=""),width = 30, height = 30,p, limitsize = FALSE)

mammal.combined1 <- subset(mammal.combined,idents="C8")
Idents(mammal.combined1) <- mammal.combined1$Cl
p <- VlnPlot(mammal.combined1, features = c("PIFO","FOXJ1","TPPP3","TP73","MUC12","HES6","LAMP3","CCNO","MMP7","WNT7A","MKI67","KRT17","ESR1","IHH","PAEP","PTGS1","KRT5","CXCL8","HES1", "HEY1","CSRNP1","SCGB2A2","FOXO1","SOX9","LGR5") )
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_VlnByCl_Batch8.pdf",sep=""),width = 30, height = 30,p, limitsize = FALSE)



mammal.combined1 <- subset(mammal.combined,idents="early-secretory")
Idents(mammal.combined1) <- mammal.combined1$Cl
p <- VlnPlot(mammal.combined1, features = c("PIFO","FOXJ1","TPPP3","TP73","MUC12","HES6","LAMP3","CCNO","MMP7","WNT7A","MKI67","KRT17","ESR1","IHH","PAEP","PTGS1","KRT5","CXCL8","HES1", "HEY1","CSRNP1","SCGB2A2","FOXO1","SOX9","LGR5") )
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_VlnByCl_EarlySec.pdf",sep=""),width = 30, height = 30,p, limitsize = FALSE)


mammal.combined1 <- subset(mammal.combined,idents="proliferative")
Idents(mammal.combined1) <- mammal.combined1$Cl
p <- VlnPlot(mammal.combined1, features = c("PIFO","FOXJ1","TPPP3","TP73","MUC12","HES6","LAMP3","CCNO","MMP7","WNT7A","MKI67","KRT17","ESR1","IHH","PAEP","PTGS1","KRT5","CXCL8","HES1", "HEY1","CSRNP1","SCGB2A2","FOXO1","SOX9","LGR5") )
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_split_epithel_VlnByCl_Prolif.pdf",sep=""),width = 30, height = 30,p, limitsize = FALSE)


mammal.combined <- readRDS(file=paste(saveext,"/Integrate_with_RVT_eithelia.rds",sep=""))
mammal.combined <- subset(mammal.combined,idents = "Ciliated",invert= TRUE)

mammal.combined$Cells <- Idents(mammal.combined)
Idents(mammal.combined) <- mammal.combined$ID3
mammal.combined <- subset(mammal.combined,idents = "C2",invert= TRUE)
Idents(mammal.combined) <- mammal.combined$Cells

C2outlier1 <- c("AAAGTGACAATACGCT-1_3", "AAGACTCAGGTTCACT-1_3", "AATGCCAAGCGACCCT-1_3", "ACCAAACTCGAGCCAC-1_3", "ACCTGAATCTAACGCA-1_3", "AGACCCGTCTTACGGA-1_3", "AGAGAATTCTCTGCTG-1_3",
"AGCCAGCCACCTCAGG-1_3", "AGCTCAAGTTCAAACC-1_3", "AGGATCTTCCAAATGC-1_3", "AGGTAGGGTGACACAG-1_3", "AGTCAACCACTGATTG-1_3", "ATAGGCTGTTGGATCT-1_3", "ATCCACCCACTTGTCC-1_3",
"ATCGGCGGTAAGAACT-1_3", "ATGTCTTTCACGATAC-1_3", "ATTGGGTCACGTGAGA-1_3", "ATTTCTGCACTAAACC-1_3", "CAACGGCCAGATCCAT-1_3", "CACACAAGTGGTAATA-1_3", "CACATGACAGCTAACT-1_3",
"CACATGATCTACGGTA-1_3", "CAGATACCAATCCAGT-1_3", "CAGCAATGTTGCATGT-1_3", "CCGGGTATCCAATCCC-1_3", "CCGGTGACAGTCAACT-1_3", "CCTCACAGTCGTTGGC-1_3", "CCTGTTGAGCCTCAAT-1_3",
"CGATGCGTCACCTTGC-1_3", "CGGAACCTCACCGCTT-1_3", "CTACATTTCCGTATGA-1_3", "CTCAACCAGTCGTTAC-1_3", "CTCATGCCACGTAACT-1_3", "CTCCACACAGCTACAT-1_3", "CTCGAGGTCATAGAGA-1_3",
"CTCTCAGGTCCCGGTA-1_3", "CTCTGGTAGACAACAT-1_3", "CTGCCATCAGGCACAA-1_3", "CTGCCATCAGTCTGGC-1_3", "CTGTGAACAGCTTCCT-1_3", "GAACTGTCAACCGGAA-1_3", "GACCTTCTCCTACCAC-1_3",
"GACGTTAAGCGTGAAC-1_3", "GACTATGGTTACGGAG-1_3", "GAGTGAGGTTGTCTAG-1_3", "GATGACTGTCTCCCTA-1_3", "GATGCTACATAGATCC-1_3", "GCAACCGCAACACGTT-1_3", "GCACATAAGTGAGGCT-1_3",
"GCACGGTGTCCAGCAC-1_3", "GCAGCTGGTCAGTCTA-1_3", "GCATCTCAGACTCCGC-1_3", "GCCCAGAAGTCACTAC-1_3", "GCTGAATTCGCTAAAC-1_3", "GCTGGGTGTGGAACAC-1_3", "GGAATGGAGGCCCACT-1_3",
"GGTTGTATCAAAGGAT-1_3", "GTAACACCAAGTCGTT-1_3", "GTAGTACAGCGGTAAC-1_3", "GTCATTTTCCATTTGT-1_3", "GTGATGTTCCTTATAC-1_3", "GTGGTTAGTGTTCGTA-1_3", "GTTGCTCTCCATCTCG-1_3",
"GTTTACTGTAGTTACC-1_3", "TAACTTCTCAAACCTG-1_3", "TAAGCCATCCGTTGGG-1_3", "TACCTCGTCAACCGAT-1_3", "TATATCCGTCAGATTC-1_3", "TATATCCTCGAGCTGC-1_3", "TCAAGACTCCATTGTT-1_3",
"TCAAGCAAGGATACAT-1_3", "TCACGGGCAGGTGACA-1_3", "TCCACGTGTGAATGAT-1_3", "TCCATCGCACCGCTAG-1_3", "TCCATGCTCCGTACGG-1_3", "TCCTCGACACGTACTA-1_3", "TCCTTCTGTGCCGTAC-1_3",
"TCGACCTAGCATCAAA-1_3", "TCGCTCACAGGTCCGT-1_3", "TCTTTGACACTGTCCT-1_3", "TGAACGTCACTTGAGT-1_3", "TGCTGAAGTACTTCCC-1_3", "TGTCAGAGTATCGCGC-1_3", "TGTGATGAGTAGGCCA-1_3",
"TTCGCTGCACCCTCTA-1_3", "TTGGTTTCACAGTACT-1_3", "TTTATGCAGATGGTAT-1_3")
  
C2outlier2 <- c("AAACCCAGTAGCGCTC-1_3", "AAAGAACCAACAACAA-1_3", "AAATGGACACTTTAGG-1_3", "AACAAAGGTGCAGATG-1_3", "AACACACTCATAAGGA-1_3", "AACCAACAGCAACTTC-1_3", "AAGACTCGTAAGAACT-1_3",
"AAGCATCTCGTCGATA-1_3", "AAGGTAACAAGCACCC-1_3", "AATCGTGAGTAACCTC-1_3", "AATTCCTTCACGATCA-1_3", "ACAAGCTCAGCCTTCT-1_3", "ACAGAAACACCCTTAC-1_3", "ACAGAAAGTGAATAAC-1_3",
"ACCCAAAGTCATAACC-1_3", "ACGTACAGTGCGGTAA-1_3", "ACGTCCTAGAGAAGGT-1_3", "AGACACTCATAACTCG-1_3", "AGATGCTCAGCACCCA-1_3", "AGCGCCATCTGCGGGT-1_3", "AGCGTCGGTTGAGGAC-1_3",
"AGGACGATCAGGAAAT-1_3", "AGTAGTCAGGTCTACT-1_3", "AGTCATGCATCGAACT-1_3", "AGTCATGGTGACTCGC-1_3", "AGTCTCCCATCGAACT-1_3", "AGTGACTAGTCATCGT-1_3", "ATACCGAGTTAGGAGC-1_3",
"ATCACAGAGCTGACTT-1_3", "ATCATTCAGATTGACA-1_3", "ATCGCCTTCAGACAAA-1_3", "ATGGATCCAACAACAA-1_3", "ATGGATCCATGTGGCC-1_3", "ATTTCTGGTAGTCTGT-1_3", "CAACGGCCAGAACATA-1_3",
"CACACAATCTAGCAAC-1_3", "CACTTCGCAGACGGAT-1_3", "CACTTCGGTTCTTGCC-1_3", "CAGATACAGATACGAT-1_3", "CAGATCAGTCAGTCCG-1_3", "CAGCAATTCGAACCAT-1_3", "CATAAGCTCTATACGG-1_3",
"CATACAGAGCCGCACT-1_3", "CATGCTCGTATGCTAC-1_3", "CCACACTCACCGTACG-1_3", "CCCTCAAAGAGTCTGG-1_3", "CCGGACAGTACCGTGC-1_3", "CCTTTGGCAGGCACAA-1_3", "CGGACACAGACTCATC-1_3",
"CGGGACTCAGCAGGAT-1_3", "CGGGACTTCAATCCAG-1_3", "CGGGTGTGTGTCATGT-1_3", "CGTGCTTCATCCGTGG-1_3", "CGTTAGACATGAAGGC-1_3", "CTAACCCTCTAACACG-1_3", "CTACCCAAGCGCATCC-1_3",
"CTATAGGTCATTGAGC-1_3", "CTCATCGCAGGCATGA-1_3", "CTCCATGTCTGATGGT-1_3", "CTCGAGGCATACAGAA-1_3", "CTCGAGGGTTTCGACA-1_3", "CTGCCTAAGACTCCGC-1_3", "CTGCTCAAGTACCGGA-1_3",
"CTGTAGACAAGTATCC-1_3", "CTTAGGAGTGGGCTTC-1_3", "GAAGAATTCGCGTTTC-1_3", "GAATCACCAGTGTGCC-1_3", "GACCTTCTCAGTGTCA-1_3", "GAGACTTTCATGTCTT-1_3", "GAGCCTGGTAGACGGT-1_3",
"GAGGCCTAGACATACA-1_3", "GAGTCTATCTTACTGT-1_3", "GAGTGAGTCCATACTT-1_3", "GATTCTTAGCATCCCG-1_3", "GATTGGTTCGTAGGGA-1_3", "GCACTAAAGTGGTGAC-1_3", "GCAGCTGAGCACGTCC-1_3",
"GCCAGCACAAATCGTC-1_3", "GCCATTCCACCACATA-1_3", "GCCGTGACAGACCATT-1_3", "GCCTGTTAGCACCGTC-1_3", "GCGAGAACAAGGATGC-1_3", "GCTTTCGCACTACCCT-1_3", "GGAGATGCAGGGAATC-1_3",
"GGCAGTCAGCAGCAGT-1_3", "GGCTTTCAGTACCATC-1_3", "GGGTAGATCTACAGGT-1_3", "GGTCACGTCTGGCCAG-1_3", "GGTTGTAAGTCTCCTC-1_3", "GTAGGAGTCCCGAGTG-1_3", "GTGAGCCTCCAAGGGA-1_3",
"GTGAGTTCAAAGGGTC-1_3", "GTGCAGCAGTGGTCAG-1_3", "GTGCAGCGTTGCGGAA-1_3", "GTGCTTCAGTGCGACA-1_3", "GTGCTTCTCTGCGGCA-1_3", "GTGGAAGCATGTGTCA-1_3", "GTGTTAGGTCAATCTG-1_3",
"GTTACCCCAGGCAATG-1_3", "GTTAGTGGTCAACATC-1_3", "GTTGCTCCATCAGCTA-1_3", "GTTTACTCATTCCTCG-1_3", "TACATTCTCTGATGGT-1_3", "TACCCACGTGAAAGTT-1_3", "TAGACTGAGTCAGCGA-1_3",
"TAGGTACCACAGAGAC-1_3", "TAGGTTGAGCATCAAA-1_3", "TATCCTAAGGCCTGAA-1_3", "TCAAGCAAGAGCCCAA-1_3", "TCACGCTCAAACGAGC-1_3", "TCACGGGGTCGTCTCT-1_3", "TCAGCCTCAGCTGTGC-1_3",
"TCATCCGGTGACGTCC-1_3", "TCATTGTAGCCTCTCT-1_3", "TCATTTGGTGCCTGCA-1_3", "TCCACCACAGACATCT-1_3", "TCGTAGACACGAAGAC-1_3", "TCTACATTCGTCTAAG-1_3", "TCTGTCGCAACTCCCT-1_3",
"TGAGCATTCCGCAACG-1_3", "TGAGGGAAGCCTAGGA-1_3", "TGATCAGAGATCCCGC-1_3", "TGATGCACAGTAGTTC-1_3", "TGCGATAGTCTTCCGT-1_3", "TGGAACTGTGTCTCCT-1_3", "TGGTGATCAACGCATT-1_3",
"TGGTTAGGTCAAGCGA-1_3", "TGTCCCAGTTCGGCTG-1_3", "TGTTTGTCATGTGCTA-1_3", "TTACTGTTCGTAATGC-1_3", "TTCATGTGTTGCAACT-1_3", "TTCCGGTCAAGTGATA-1_3", "TTCGGTCAGCCGTTGC-1_3",
"TTGAGTGTCTGTCGCT-1_3", "TTGGGATGTGTTGAGG-1_3", "TTGTTGTTCCAGTACA-1_3", "TTTCACATCTAGAACC-1_3")

library(destiny)
mammal.combined <- subset(mammal.combined,cells=C2outlier1,invert=TRUE)
mammal.combined <- subset(mammal.combined,cells=C2outlier2,invert=TRUE)
DN1 <- as.data.frame(mammal.combined[["pca"]]@cell.embeddings)[,1:20]
dm <- DiffusionMap(DN1)
DM <- dm@eigenvectors*100
rownames(DM) <- colnames(mammal.combined)
mammal.combined[["dm"]] <- CreateDimReducObject(embeddings = DM, key = "DM_", assay = DefaultAssay(mammal.combined))
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial",".pdf",sep=""),width = 80, height = 8, limitsize = FALSE,p)




DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
mammal.combined <- FindClusters(mammal.combined, resolution = 1)
mammal.combined$Cl <- Idents(mammal.combined)
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial",".pdf",sep=""),width = 80, height = 8, limitsize = FALSE,p)

Idents(mammal.combined) <- mammal.combined$ID3
mammal.combined <- subset(mammal.combined,idents = "C2",invert= FALSE)
#Idents(mammal.combined) <- mammal.combined$Cells
Idents(mammal.combined) <- mammal.combined$Cl #colnames(mammal.combined)

p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 12, label.size = 12, label = TRUE, repel = TRUE) + NoLegend()
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial",".pdf",sep=""),width = 30, height = 30, limitsize = FALSE,p)

#p<-FeaturePlot(DsubsetTb,  reduction = "dm", features = "HLA-G", cols =  c("lightgrey", "black"), pt.size = 4)
#ggsave(filename=paste(saveext,"/DimRed/Tbsubset3_PCA_DM.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)



#Plots
TomList <- c("PLA2G2A",
"DIO2",
"SCARA5",
"CXCL8",
"CXCL14",
"TIMP3",
"IL15",
"ECM1",
"AXL",
"WNT7A",
"PTGS1",
"HEY1",
"DNAI1",
"FOXJ1",
"SPP1",
"PAEP",
"DPP4",
"LIF",
"AREG",
"EREG",
"NEAT1",
"KCNQ1OT1",
"WNT5A")







#Align each batch with the reference
mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D1), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC1_to_VentoTormo.pdf",sep=""),width = 20, height = 8,p)

mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D2), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC2_to_VentoTormo.pdf",sep=""),width = 20, height = 8,p)

mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D3), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC3_to_VentoTormo.pdf",sep=""),width = 20, height = 8,p)

mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D4), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC4_to_VentoTormo.pdf",sep=""),width = 20, height = 8,p)

mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D6), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC6_to_VentoTormo.pdf",sep=""),width = 20, height = 8,p)

mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D7), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC7_to_VentoTormo.pdf",sep=""),width = 20, height = 8,p)

mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D8), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC8_to_VentoTormo.pdf",sep=""),width = 20, height = 8,p)




#Allign all data to the reference
D <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/CorrectedFinalAno.rds")
D$Cells <- Idents(D)
#Idents(D) <- paste(D$ID3,D$Cells,sep="_")
D1 <- subset(D,idents=c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular"))
rm(D)

mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D1), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo.pdf",sep=""),width = 90, height = 8,p,limitsize = FALSE)
Dsubset1 <- FindClusters(mammal.combined, resolution = 0.3)
p<-DimPlot(Dsubset1, pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_Cl.pdf",sep=""),width = 62, height = 10, useDingbats = FALSE, limitsize = FALSE,p)




#Do it again but for proliferative stage
D <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/CorrectedFinalAno.rds")
D$Cells <- Idents(D)

Idents(D) <- paste(D$ID3,D$Cells,sep="_")
D1 <- subset(D,idents=c("C1_Lumenal","C1_SOX9P","C1_Prolif","C1_SOX9LRG5","C1_Ciliated","C1_Glandular"))
D2 <- subset(D,idents=c("C2_Lumenal","C2_SOX9P","C2_Prolif","C2_SOX9LRG5","C2_Ciliated","C2_Glandular"))
D3 <- subset(D,idents=c("C3_Lumenal","C3_SOX9P","C3_Prolif","C3_Ciliated","C3_Glandular"))
D4 <- subset(D,idents=c("C4_Lumenal","C4_SOX9P","C4_Prolif","C4_SOX9LRG5","C4_Ciliated","C4_Glandular"))
D6 <- subset(D,idents=c("C6_Lumenal","C6_SOX9P","C6_Prolif","C6_SOX9LRG5","C6_Ciliated","C6_Glandular"))
D7 <- subset(D,idents=c("C7_Lumenal","C7_SOX9P","C7_Prolif","C7_Ciliated","C7_Glandular"))
D8 <- subset(D,idents=c("C8_Lumenal","C8_SOX9P","C8_Prolif","C8_SOX9LRG5","C8_Ciliated","C8_Glandular"))
rm(D)

Endo <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Endometrium/other_endometrium_data.rds")
Idents(Endo) <- Endo$ID3
Endo <- subset(Endo,idents="proliferative")
Idents(Endo) <- Endo$ID5
Endo <- subset(Endo,idents=c("Glandular","SOX9","Lumenal","Ciliated") ) #"dS","Fibroblast C7","eS"

mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D1), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC1_to_VentoTormo_proliferative.pdf",sep=""),width = 20, height = 8,p)

mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D2), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC2_to_VentoTormo_proliferative.pdf",sep=""),width = 20, height = 8,p)


mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D3), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC3_to_VentoTormo_proliferative.pdf",sep=""),width = 20, height = 8,p)



mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D4), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC4_to_VentoTormo_proliferative.pdf",sep=""),width = 20, height = 8,p)


mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D6), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC6_to_VentoTormo_proliferative.pdf",sep=""),width = 20, height = 8,p)



mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D7), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC7_to_VentoTormo_proliferative.pdf",sep=""),width = 20, height = 8,p)

mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D8), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC8_to_VentoTormo_proliferative.pdf",sep=""),width = 20, height = 8,p)


#
D <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/CorrectedFinalAno.rds")
D$Cells <- Idents(D)
#Idents(D) <- paste(D$ID3,D$Cells,sep="_")
D1 <- subset(D,idents=c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular"))
rm(D)

mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D1), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_proliferative.pdf",sep=""),width = 90, height = 8,p,limitsize = FALSE)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_proliferative_dim_1_3.pdf",sep=""),width = 90, height = 8,p,limitsize = FALSE)


Dsubset1 <- FindClusters(mammal.combined, resolution = 0.3)
p<-DimPlot(Dsubset1, pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_VentoTormo_Cl_proliferative.pdf",sep=""),width = 62, height = 10, useDingbats = FALSE, limitsize = FALSE,p)





Endo <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Endometrium/other_endometrium_data.rds")
#Idents(Endo) <- Endo$ID3
#Endo <- subset(Endo,idents="proliferative")
Idents(Endo) <- Endo$ID5
Endo <- subset(Endo,idents=c("Glandular","SOX9","Lumenal","Ciliated") ) #"dS","Fibroblast C7","eS"

#
D <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/CorrectedFinalAno.rds")
D$Cells <- Idents(D)
#Idents(D) <- paste(D$ID3,D$Cells,sep="_")
D1 <- subset(D,idents=c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular"))
rm(D)


mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D1), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_AllVentoTormo.pdf",sep=""),width = 90, height = 8,p,limitsize = FALSE)


mammal.combined$Cells2 <- Idents(mammal.combined)

OrderID <- as.character(mammal.combined$ID3)
OrderID[which(OrderID=="proliferative")] <- "1) proliferative"
OrderID[which(OrderID=="early-secretory")] <- "2) early-secretory"
OrderID[which(OrderID=="early-mid-secretory")] <- "3) early-mid-secretory"
OrderID[which(OrderID=="mid-secretory")] <- "4) mid-secretory"
OrderID[which(OrderID=="late-secretory")] <- "5) late-secretory"

#OrderID[which(OrderID=="C2")] <- "5) C2"



mammal.combined$OrderID <- OrderID

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "OrderID", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_AllVentoTormo.pdf",sep=""),width = 90, height = 8,p,limitsize = FALSE)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "OrderID", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_AllVentoTormo_1_3.pdf",sep=""),width = 90, height = 8,p,limitsize = FALSE)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "OrderID", label = TRUE, repel = TRUE, dims = c(2,3)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_AllVentoTormo_2_3.pdf",sep=""),width = 90, height = 8,p,limitsize = FALSE)


p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "OrderID", label = TRUE, repel = TRUE, dims = c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_AllVentoTormo_1_4.pdf",sep=""),width = 90, height = 8,p,limitsize = FALSE)


ExtendedID <- paste(mammal.combined$ID3,mammal.combined$Cells2,sep="_")
Idents(mammal.combined) <- ExtendedID


AvExp <- AverageExpression(mammal.combined)


C1 <- cor( as.matrix( log(AvExp$RNA+1) ), as.matrix( log(AvExp$RNA +1) ), method = "pearson")
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1, color =  redblue1(20),  display_numbers = round(C1,2), border_color = NA, cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/Markers/RawExpressionCorr",".pdf",sep="") ,width=35,height=35)


C1 <- cor( as.matrix( log(AvExp$integrated+1) ), as.matrix( log(AvExp$integrated +1) ), method = "pearson")
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1, color =  redblue1(20),  display_numbers = round(C1,2), border_color = NA, cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/Markers/NormalisedExpressionCorr",".pdf",sep="") ,width=35,height=35)

list1 <- c("proliferative_SOX9",  #2370
"early-secretory_SOX9", #1246
"early-mid-secretory_SOX9", #24
#"mid-secretory_SOX9", #3
#"late-secretory_SOX9", #3
"proliferative_Lumenal", #999
"early-secretory_Lumenal", #909
"early-mid-secretory_Lumenal", #12091
"mid-secretory_Lumenal", #240
#"late-secretory_Lumenal", #13
"proliferative_Glandular", #588
"early-secretory_Glandular", #1575
"early-mid-secretory_Glandular", #6543  
"mid-secretory_Glandular", #3284
"late-secretory_Glandular", #1110
"proliferative_Ciliated",   #505
"early-secretory_Ciliated", #107 
"early-mid-secretory_Ciliated", #2333 
"mid-secretory_Ciliated", #57
"late-secretory_Ciliated") #187


list2 <-c("C1_Prolif", #560 
"C2_Prolif", #1904
"C3_Prolif", #89 
"C4_Prolif", #492
"C6_Prolif", #936 
"C7_Prolif", #257
"C8_Prolif", #1595
"C1_SOX9P", #304
"C2_SOX9P", #1933 
"C3_SOX9P", #41
"C4_SOX9P", # 420 
"C6_SOX9P", #783 
"C7_SOX9P", #185 
"C8_SOX9P", #1403
#"C1_SOX9LRG5", #12
#"C2_SOX9LRG5", #8 
#"C4_SOX9LRG5", #2
#"C6_SOX9LRG5", #3 
#"C8_SOX9LRG5", #3
"C1_Lumenal", #361
"C2_Lumenal", #797
"C3_Lumenal", #113
"C4_Lumenal", #531
"C6_Lumenal", #1973
"C7_Lumenal", #305
"C8_Lumenal", #1027
#"C1_Glandular", #1
#"C2_Glandular", #1
#"C3_Glandular", #4 
"C4_Glandular", #72 
"C6_Glandular", #15
#"C7_Glandular", #3
#"C8_Glandular", #3
#"C1_Ciliated", #4 
"C2_Ciliated", #104
#"C3_Ciliated", #3 
"C4_Ciliated", #31 
"C6_Ciliated", #42 
#"C7_Ciliated", #2
"C8_Ciliated") #14 




C1 <- cor( as.matrix( log(AvExp$RNA[,list1]+1) ), as.matrix( log(AvExp$RNA[,list2] +1) ), method = "pearson")
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1, color =  redblue1(20),  display_numbers = round(C1,2), border_color = NA, cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/Markers/RawExpressionCrossCorr",".pdf",sep="") ,width=35,height=35)

C1 <- cor( as.matrix( log(AvExp$integrated[,list1]+1) ), as.matrix( log(AvExp$integrated[,list2] +1) ), method = "pearson")
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1, color =  redblue1(20),  display_numbers = round(C1,2), border_color = NA, cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/Markers/NormExpressionCrossCorr",".pdf",sep="") ,width=35,height=35)


C1 <- cor( as.matrix( log(AvExp$RNA[,list1]+1) ), as.matrix( log(AvExp$RNA[,list1] +1) ), method = "pearson")
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1, color =  redblue1(20),  display_numbers = round(C1,2), border_color = NA, cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/Markers/RAWExpressionRefCorr",".pdf",sep="") ,width=35,height=35)



C1 <- cor( as.matrix( log(AvExp$RNA[,list2]+1) ), as.matrix( log(AvExp$RNA[,list2] +1) ), method = "pearson")
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1, color =  redblue1(20),  display_numbers = round(C1,2), border_color = NA, cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/Markers/RAWExpressionMatteoCorr",".pdf",sep="") ,width=35,height=35)


Idents(mammal.combined) <- mammal.combined$ID3
AvExp <- AverageExpression(mammal.combined)



C1 <- cor( as.matrix( log(AvExp$RNA+1) ), as.matrix( log(AvExp$RNA +1) ), method = "pearson")
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1, color =  redblue1(20),  display_numbers = round(C1,2), border_color = NA, cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/Markers/RawExpressionBulkedCorr",".pdf",sep="") ,width=10,height=10)



Endo <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Endometrium/other_endometrium_data.rds")
#Idents(Endo) <- Endo$ID3
#Endo <- subset(Endo,idents="proliferative")
Idents(Endo) <- Endo$ID5
Endo <- subset(Endo,idents=c("Glandular","SOX9","Lumenal","Ciliated") ) #"dS","Fibroblast C7","eS"


Endo <- ScaleData(Endo, verbose = FALSE)
Endo <- RunPCA(Endo, npcs = 20, verbose = FALSE)
Endo <- RunUMAP(Endo, reduction = "pca", dims = 1:20)
Endo <- FindNeighbors(Endo, reduction = "pca", dims = 1:20)#

p<-DimPlot(Endo,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_VentoTormo_All.pdf",sep=""),width = 20, height = 8,p)




Endo <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Endometrium/other_endometrium_data.rds")
#Idents(Endo) <- Endo$ID3
#Endo <- subset(Endo,idents="proliferative")
Idents(Endo) <- Endo$ID5
Endo <- subset(Endo,idents=c("Glandular","SOX9","Lumenal","Ciliated") ) #"dS","Fibroblast C7","eS"


mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D1), dims = 1:30, anchor.features = 5000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_All_VentoTormo_proliferative.pdf",sep=""),width = 90, height = 8,p,limitsize = FALSE)


p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_All_VentoTormo_proliferative.pdf",sep=""),width = 90, height = 8,p,limitsize = FALSE)

Dsubset1 <- FindClusters(mammal.combined, resolution = 0.3)
p<-DimPlot(Dsubset1, pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_All_VentoTormo_Cl_proliferative.pdf",sep=""),width = 62, height = 10, useDingbats = FALSE, limitsize = FALSE,p)



DefaultAssay(mammal.combined) <- "RNA"


p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "SOX9", split.by = "ID3", cols =  c("lightgrey", "black"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_SOX9.pdf",sep=""),width = 62, height = 10, limitsize = FALSE,p)


p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "LGR5", split.by = "ID3", cols =  c("lightgrey", "black"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_LGR5.pdf",sep=""),width = 62, height = 10, limitsize = FALSE,p)


p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "PAEP", split.by = "ID3", cols =  c("lightgrey", "black"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_PAEP.pdf",sep=""),width = 62, height = 10, limitsize = FALSE,p)



p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "CXCL8", split.by = "ID3", cols =  c("lightgrey", "black"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_CXCL8.pdf",sep=""),width = 62, height = 10, limitsize = FALSE,p)


p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "KLF6", split.by = "ID3", cols =  c("lightgrey", "black"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/PCA_KKLF6.pdf",sep=""),width = 62, height = 10, limitsize = FALSE,p)





p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_All_VentoTormo_proliferative.pdf",sep=""),width = 90, height = 8,p,limitsize = FALSE)






mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo,D1), dims = 1:30, anchor.features = 3000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareAllC_to_All_VentoTormo_3k.pdf",sep=""),width = 90, height = 8,p,limitsize = FALSE)


p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_CompareAllC_to_All_VentoTormo_3k.pdf",sep=""),width = 90, height = 8,p,limitsize = FALSE)





Idents(Endo3) <- Endo3$ID5
Endo3 <- subset(Endo3,idents = c("Glandular","SOX9","Ciliated","Lumenal"))
mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo3,D1), dims = 1:30, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC1_to_VentoTormo_proliferative.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)
DN1 <- as.data.frame(mammal.combined[["pca"]]@cell.embeddings)[,1:20]
dm <- DiffusionMap(DN1)
DM <- dm@eigenvectors*100
rownames(DM) <- colnames(mammal.combined)
mammal.combined[["dm"]] <- CreateDimReducObject(embeddings = DM, key = "DM_", assay = DefaultAssay(mammal.combined))
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial_C1",".pdf",sep=""),width = 30, height = 8, limitsize = FALSE,p)



Endo4 <- Endo3
Idents(Endo3) <- Endo3$ID5

Endo3 <- subset(Endo3,idents = c("Glandular","SOX9","Ciliated","Lumenal"))

D2p 
D2p <- subset(D2,cells=outlier,invert=TRUE)
mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo3,D2p), dims = 1:30, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC2_to_VentoTormo_proliferative.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)
DN1 <- as.data.frame(mammal.combined[["pca"]]@cell.embeddings)[,1:20]
dm <- DiffusionMap(DN1)
DM <- dm@eigenvectors*100
rownames(DM) <- colnames(mammal.combined)
mammal.combined[["dm"]] <- CreateDimReducObject(embeddings = DM, key = "DM_", assay = DefaultAssay(mammal.combined))
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial_C2",".pdf",sep=""),width = 30, height = 8, limitsize = FALSE,p)
mammal.combined <- FindClusters(mammal.combined, resolution = .1)
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial_C2_Cl",".pdf",sep=""),width = 30, height = 8, limitsize = FALSE,p)






Idents(Endo3) <- Endo3$ID5
Endo3 <- subset(Endo3,idents = c("Glandular","SOX9","Ciliated","Lumenal"))
D3p <- subset(D3,cells=outlier,invert=TRUE)
mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo3,D3p), dims = 1:30, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC3_to_VentoTormo_proliferative.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)
DN1 <- as.data.frame(mammal.combined[["pca"]]@cell.embeddings)[,1:20]
dm <- DiffusionMap(DN1)
DM <- dm@eigenvectors*100
rownames(DM) <- colnames(mammal.combined)
mammal.combined[["dm"]] <- CreateDimReducObject(embeddings = DM, key = "DM_", assay = DefaultAssay(mammal.combined))
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial_C3",".pdf",sep=""),width = 30, height = 8, limitsize = FALSE,p)
mammal.combined <- FindClusters(mammal.combined, resolution = .1)
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial_C3_Cl",".pdf",sep=""),width = 30, height = 8, limitsize = FALSE,p)




Idents(Endo3) <- Endo3$ID5
Endo3 <- subset(Endo3,idents = c("Glandular","SOX9","Ciliated","Lumenal"))
D6p <- subset(D6,cells=outlier,invert=TRUE)
mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo3,D6p), dims = 1:30, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC6_to_VentoTormo_proliferative.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)
DN1 <- as.data.frame(mammal.combined[["pca"]]@cell.embeddings)[,1:20]
dm <- DiffusionMap(DN1)
DM <- dm@eigenvectors*100
rownames(DM) <- colnames(mammal.combined)
mammal.combined[["dm"]] <- CreateDimReducObject(embeddings = DM, key = "DM_", assay = DefaultAssay(mammal.combined))
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial_C6",".pdf",sep=""),width = 30, height = 8, limitsize = FALSE,p)
mammal.combined <- FindClusters(mammal.combined, resolution = .1)
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial_C6_Cl",".pdf",sep=""),width = 30, height = 8, limitsize = FALSE,p)




Idents(Endo3) <- Endo3$ID5
Endo3 <- subset(Endo3,idents = c("Glandular","SOX9","Ciliated","Lumenal"))
D7p <- subset(D7,cells=outlier,invert=TRUE)
mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo3,D7p), dims = 1:30, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC7_to_VentoTormo_proliferative.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)
DN1 <- as.data.frame(mammal.combined[["pca"]]@cell.embeddings)[,1:20]
dm <- DiffusionMap(DN1)
DM <- dm@eigenvectors*100
rownames(DM) <- colnames(mammal.combined)
mammal.combined[["dm"]] <- CreateDimReducObject(embeddings = DM, key = "DM_", assay = DefaultAssay(mammal.combined))
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial_C7",".pdf",sep=""),width = 30, height = 8, limitsize = FALSE,p)
mammal.combined <- FindClusters(mammal.combined, resolution = .1)
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial_C7_Cl",".pdf",sep=""),width = 30, height = 8, limitsize = FALSE,p)




Idents(Endo3) <- Endo3$ID5
Endo3 <- subset(Endo3,idents = c("Glandular","SOX9","Ciliated","Lumenal"))
D8p <- subset(D8,cells=outlier,invert=TRUE)
mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo3,D8p), dims = 1:30, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC8_to_VentoTormo_proliferative.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)
DN1 <- as.data.frame(mammal.combined[["pca"]]@cell.embeddings)[,1:20]
dm <- DiffusionMap(DN1)
DM <- dm@eigenvectors*100
rownames(DM) <- colnames(mammal.combined)
mammal.combined[["dm"]] <- CreateDimReducObject(embeddings = DM, key = "DM_", assay = DefaultAssay(mammal.combined))
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial_C8",".pdf",sep=""),width = 30, height = 8, limitsize = FALSE,p)
mammal.combined <- FindClusters(mammal.combined, resolution = .1)
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial_C8_Cl",".pdf",sep=""),width = 30, height = 8, limitsize = FALSE,p)





#Do it again but for proliferative stage
D <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/CorrectedFinalAno.rds")
D$Cells <- Idents(D)
D <- subset(D,idents=c("Lumenal","SOX9P","Prolif","SOX9LRG5","Ciliated","Glandular"))
Idents(D) <- paste(D$ID3)
D <- subset(D,idents=c("C2"), invert=TRUE)
Idents(D) <- D$Cells
Endo4 <- Endo3
Idents(Endo3) <- Endo3$ID5
Endo3 <- subset(Endo3,idents = c("Glandular","SOX9","Ciliated","Lumenal"))

Endo3$IDX <- Endo3$ID3
D$IDX <- "Our dataset"
mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo3,D), dims = 1:30, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "IDX", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC2_to_VentoTormo_proliferative.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)
DN1 <- as.data.frame(mammal.combined[["pca"]]@cell.embeddings)[,1:20]
dm <- DiffusionMap(DN1)
DM <- dm@eigenvectors*100
rownames(DM) <- colnames(mammal.combined)
mammal.combined[["dm"]] <- CreateDimReducObject(embeddings = DM, key = "DM_", assay = DefaultAssay(mammal.combined))
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "IDX", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial_notC2",".pdf",sep=""),width = 30, height = 8, limitsize = FALSE,p)
mammal.combined <- FindClusters(mammal.combined, resolution = .1)
#p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
#ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial_notC2_Cl",".pdf",sep=""),width = 30, height = 8, limitsize = FALSE,p)


DefaultAssay(mammal.combined) <- "RNA" 
p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "PAEP", split.by = "IDX", cols =  c("lightgrey", "red"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/TSNE_PAEP_alignstages.pdf",sep=""),width = 30, height = 8,p)

p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "SOX9", split.by = "IDX", cols =  c("lightgrey", "red"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/TSNE_SOX9_alignstages.pdf",sep=""),width = 30, height = 8,p)

p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "SPP1", split.by = "IDX", cols =  c("lightgrey", "red"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/TSNE_SPP1_alignstages.pdf",sep=""),width = 30, height = 8,p)

p<-FeaturePlot(mammal.combined,  reduction = "pca", features = "SCGB2A2", split.by = "IDX", cols =  c("lightgrey", "red"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/TSNE_SCGB2A2_alignstages.pdf",sep=""),width = 30, height = 8,p)





#Do it again but for proliferative stage
D <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/CorrectedFinalAno.rds")
D$Cells <- Idents(D)
D <- subset(D,idents=c("Lumenal","SOX9P","Prolif","SOX9LRG5","Glandular"))
Idents(D) <- paste(D$ID3)
D <- subset(D,idents=c("C2"), invert=TRUE)
Idents(D) <- D$Cells



Endo3 <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Endometrium/other_endometrium_data.rds")
Idents(Endo3) <- Endo3$ID3
Endo3 <- subset(Endo3,idents=c("early-secretory","proliferative") )
#Idents(Endo3) <- Endo3$ID5
Idents(Endo3) <- Endo3$ID5
Endo3 <- subset(Endo3,idents = c("Glandular","SOX9","Lumenal"))
mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo3,D), dims = 1:30, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC2_to_VentoTormo_proliferative_2.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)
DN1 <- as.data.frame(mammal.combined[["pca"]]@cell.embeddings)[,1:20]
dm <- DiffusionMap(DN1)
DM <- dm@eigenvectors*100
rownames(DM) <- colnames(mammal.combined)
mammal.combined[["dm"]] <- CreateDimReducObject(embeddings = DM, key = "DM_", assay = DefaultAssay(mammal.combined))
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial_notC2_2",".pdf",sep=""),width = 30, height = 8, limitsize = FALSE,p)
mammal.combined <- FindClusters(mammal.combined, resolution = .1)



D <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/CorrectedFinalAno.rds")
D$Cells <- Idents(D)
D <- subset(D,idents=c("Lumenal","SOX9P","Prolif","SOX9LRG5","Glandular","Stromal fibroblasts","Ciliated"))
Idents(D) <- paste(D$ID3)
D <- subset(D,idents=c("C2"), invert=TRUE)
Idents(D) <- D$Cells



Endo3 <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Endometrium/other_endometrium_data.rds")
Idents(Endo3) <- Endo3$ID3
Endo3 <- subset(Endo3,idents=c("early-secretory","proliferative") )
#Idents(Endo3) <- Endo3$ID5
Endo4 <- Endo3
Idents(Endo3) <- Endo3$ID5
Endo3 <- subset(Endo3,idents = c("Glandular","SOX9","Lumenal","Fibroblast C7","Ciliated","eS","dS"))

mammal.anchors <- FindIntegrationAnchors(object.list = list(Endo3,D), dims = 1:30, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_CompareC2_to_VentoTormo_proliferative_3.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)
DN1 <- as.data.frame(mammal.combined[["pca"]]@cell.embeddings)[,1:20]
dm <- DiffusionMap(DN1)
DM <- dm@eigenvectors*100
rownames(DM) <- colnames(mammal.combined)
mammal.combined[["dm"]] <- CreateDimReducObject(embeddings = DM, key = "DM_", assay = DefaultAssay(mammal.combined))
p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 2, split.by = "ID3", label = TRUE, repel = TRUE) +xlim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/DM_Epithelial_notC2_3",".pdf",sep=""),width = 30, height = 8, limitsize = FALSE,p)
mammal.combined <- FindClusters(mammal.combined, resolution = .1)




#Now do the vln plots
Endo3 <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Endometrium/other_endometrium_data.rds")
Idents(Endo3) <- Endo3$ID3
Endo3 <- subset(Endo3,idents=c("early-secretory","late-secretory","proliferative") )
Idents(Endo3) <- Endo3$ID5
Endo4 <- Endo3
Endo4 <- subset(Endo4,idents = c("Glandular","SOX9","Lumenal","Fibroblast C7","Ciliated","eS","dS"))

mammal.combined <- readRDS(file=paste(saveext,"/Integrate_with_RVT_eithelia.rds",sep=""))

mammal.combined <- readRDS(file=paste(saveext,"/Integrate_with_RVT_stroma_eithelia.rds",sep=""))
cilcells = WhichCells(mammal.combined,idents="Ciliated")



Idents(mammal.combined) <- paste(mammal.combined$ID3,mammal.combined$ID5,sep="_")
mammal.combine <- mammal.combined
Idents(mammal.combined) <- mammal.combine$ID3
ours <- WhichCells(mammal.combine,idents=c("C1","C2","C3","C4","C6","C7","C8"))

onlyOurs <- subset(mammal.combined,idents=c("C1","C2","C3","C4","C7","C6","C8"))
Idents(onlyOurs) <- onlyOurs$Cells
onlyOurs <- subset(onlyOurs,idents = c("Stromal fibroblasts"),invert = TRUE)

SOX <- WhichCells(mammal.combined,expression=SOX9> log(1+1) )
LRG <- WhichCells(mammal.combined,expression=LGR5> log(1+1) )
SPP <- WhichCells(mammal.combined,expression=SPP1> log(1+1) )
SCGB <- WhichCells(mammal.combined,expression=SCGB2A2> log(1+1) )

Idents(onlyOurs, cells= ours) <- "Lumenal"
Idents(onlyOurs, cells= intersect(SPP,ours)) <- "Glandular"
Idents(onlyOurs, cells=intersect(SCGB,ours)) <- "Glandular"
Idents(onlyOurs, cells=intersect(SOX,ours)) <- "SOX9"
Idents(onlyOurs, cells=intersect(LRG,ours)) <- "LRG5"
Idents(onlyOurs, cells=cilcells) <- "Ciliated"

Idents(Endo4) <- paste(Endo4$ID3,Idents(Endo4),sep="")
List1 <- FindMarkers(Endo4,ident.1 = c("late-secretoryGlandular"),ident.2 = "proliferativeLumenal", test.use = "MAST", only.pos = TRUE )
List2 <- FindMarkers(Endo4,ident.1 = "proliferativeLumenal",ident.2 = c("late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )
List3 <- FindMarkers(Endo4,ident.1 = "proliferativeSOX9",ident.2 = c("proliferativeCiliated","proliferativeGlandular","proliferativeLumenal"), test.use = "MAST", only.pos = TRUE )
List4 <- FindMarkers(Endo4,ident.1 = "proliferativeCiliated",ident.2 = c("proliferativeLumenal"), test.use = "MAST", only.pos = TRUE )

List4o <- FindMarkers(Endo4,ident.2 = "proliferativeCiliated",ident.1 = c("proliferativeLumenal"), test.use = "MAST", only.pos = TRUE )


List5 <- FindMarkers(Endo4,ident.1 = "proliferativeGlandular",ident.2 = c("proliferativeLumenal"), test.use = "MAST", only.pos = TRUE )
List6 <- FindMarkers(Endo4,ident.1 = c("proliferativeLumenal"),ident.2 = "proliferativeGlandular", test.use = "MAST", only.pos = TRUE )


newdata1 <- List1[order(List1$avg_log2FC,decreasing = TRUE),]
newdata2 <- List2[order(List2$avg_log2FC,decreasing = TRUE),]
newdata3 <- List3[order(List3$avg_log2FC,decreasing = TRUE),]
newdata4 <- List4[order(List4$avg_log2FC,decreasing = TRUE),]
newdata4o <- List4o[order(List4o$avg_log2FC,decreasing = TRUE),]

newdata5 <- List5[order(List5$avg_log2FC,decreasing = TRUE),]
newdata6 <- List6[order(List6$avg_log2FC,decreasing = TRUE),]

DefaultAssay(mammal.combined)
mammal.combined$ID4 <- NULL
p1<-VlnPlot(mammal.combined,  idents = c("Ciliated","Lumenal","Glandular"),features = rownames(newdata1)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_set1.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(mammal.combined,   idents = c("Ciliated","Lumenal","Glandular"),features = rownames(newdata2)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_set2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(mammal.combined,   idents = c("Ciliated","Lumenal","Glandular"),features = rownames(newdata3)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_set3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(mammal.combined,   idents = c("Ciliated","Lumenal","Glandular"),features = rownames(newdata4)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_set4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

p1<-VlnPlot(Endo4,  idents = c("proliferativeCiliated","late-secretoryCiliated","proliferativeLumenal","proliferativeGlandular","late-secretoryGlandular"),features = rownames(newdata1)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_RVT_set1.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(Endo4,   idents = c("proliferativeCiliated","late-secretoryCiliated","proliferativeLumenal","proliferativeGlandular","late-secretoryGlandular"),features = rownames(newdata2)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_RVT_set2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(Endo4,   idents = c("proliferativeCiliated","late-secretoryCiliated","proliferativeLumenal","proliferativeGlandular","late-secretoryGlandular"),features = rownames(newdata3)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_RVT_set3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(Endo4,   idents = c("proliferativeCiliated","late-secretoryCiliated","proliferativeLumenal","proliferativeGlandular","late-secretoryGlandular"),features = rownames(newdata4)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_RVT_set4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

DefaultAssay(Endo0) <- "RNA"
p1<-VlnPlot(Endo0,  idents = c("Epi1","Epi2"),features = rownames(newdata1)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_RVT2_set1.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(Endo0,   idents = c("Epi1","Epi2"),features = rownames(newdata2)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_RVT2_set2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(Endo0,   idents = c("Epi1","Epi2"),features = rownames(newdata3)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_RVT2_set3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(Endo0,   idents = c("Epi1","Epi2"),features = rownames(newdata4)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_RVT2_set4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

DefaultAssay(assemb) <- "RNA"
p1<-VlnPlot(assemb,  idents = c("EpS1","EpS2","EpS3","EpS4","EpS5"),features = rownames(newdata1)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_set1.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(assemb,   idents = c("EpS1","EpS2","EpS3","EpS4","EpS5"),features = rownames(newdata2)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_set2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(assemb,   idents = c("EpS1","EpS2","EpS3","EpS4","EpS5"),features = rownames(newdata3)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_set3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(assemb,   idents = c("EpS1","EpS2","EpS3","EpS4","EpS5"),features = rownames(newdata4)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_set4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

mammal.subset <- subset(mammal.combined,idents=c("Lumenal","Glandular","Ciliated"))

DefaultAssay(onlyOurs) <- "RNA"
AllData3 <- merge(onlyOurs,y=c(Endo4,Endo0,assemb),project = "Merged")

Idents(AllData3,cells=WhichCells(AllData3,idents=c("Lumenal"))) <- "11) Lumenal"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("proliferativeLumenal"))) <- "12) proliferativeLumenal"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("early-secretoryLumenal"))) <- "13) early-secretoryLumenal"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("Glandular"))) <- "21) Glandular"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("proliferativeGlandular"))) <- "22) proliferativeGlandular"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("early-secretoryGlandular"))) <- "23) early-secretoryGlandular"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("late-secretoryGlandular"))) <- "24) late-secretoryGlandular"

Idents(AllData3,cells=WhichCells(AllData3,idents=c("Ciliated"))) <- "31) Ciliated"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("proliferativeCiliated"))) <- "32) proliferativeCiliated"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("early-secretoryCiliated"))) <- "33) early-secretoryCiliated"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("late-secretoryCiliated"))) <- "34) late-secretoryCiliated"

Idents(AllData3,cells=WhichCells(AllData3,idents=c("Epi1"))) <- "41) Epi1"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("Epi2"))) <- "42) Epi2"

Idents(AllData3,cells=WhichCells(AllData3,idents=c("EpS1"))) <- "51) EpS1"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("EpS2"))) <- "52) EpS2"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("EpS3"))) <- "53) EpS3"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("EpS4"))) <- "54) EpS4"
Idents(AllData3,cells=WhichCells(AllData3,idents=c("EpS5"))) <- "55) EpS5"

AllData4 <- subset(AllData3,idents=c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"))
Idents(AllData4) <- factor(x = Idents(AllData4), levels = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"))

AllData4$ID1 <- NULL

#remotes::install_version("SeuratObject", "4.1.4")

DefaultAssay(AllData4) <- "RNA"
p1<-VlnPlot(AllData4,  idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata1)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_set1.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata2)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_set2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata3)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_set3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata4)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_set4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)


p1<-VlnPlot(AllData4,  idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata1)[41:80], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_set1_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata2)[41:80],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_set2_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata3)[41:80],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_set3_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata4)[41:80],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_set4_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

DefaultAssay(AllData4) <- "RNA"
DefaultAssay(AllData5) <- "RNA"

Dats <- subset(AllData4,idents=c("21) Glandular","11) Lumenal","31) Ciliated"))
Dats$ID0 <- Idents(Dats)
Idents(Dats) <- Dats$Dataset
Dats <- subset(Dats,idents="10X Ours")
Idents(Dats) <- Dats$ID0

List1B <- FindMarkers(AllData4,ident.1 = "21) Glandular", ident.2 = "11) Lumenal", test.use = "MAST", only.pos = TRUE) #, min.cells.group = 1, 
                      #min.cells.feature = 1,
                      #min.pct = 0,
                      #logfc.threshold = 0)
List2Bs <- FindMarkers(AllData4,ident.1 = "11) Lumenal", ident.2 = c("21) Glandular"), test.use = "MAST", only.pos = TRUE )#, min.cells.group = 1, 
                      #min.cells.feature = 1,
                      #min.pct = 0,
                      #logfc.threshold = 0)
List4Bs <- FindMarkers(AllData4,ident.1 = "31) Ciliated",ident.2 = c("11) Lumenal"), test.use = "MAST", only.pos = TRUE) #, min.cells.group = 1, 
                      #min.cells.feature = 1,
                      #min.pct = 0,
                      #logfc.threshold = 0)

List4Bo <- FindMarkers(AllData4,ident.2 = "31) Ciliated",ident.1 = c("11) Lumenal"), test.use = "MAST", only.pos = TRUE)#, min.cells.group = 1, 
                       #min.cells.feature = 1,
                       #min.pct = 0,
                       #logfc.threshold = 0 )


#List1B <- FindMarkers(Dats,ident.1 = "21) Glandular", ident.2 = "11) Lumenal", test.use = "MAST", only.pos = TRUE )
#List2B <- FindMarkers(Dats,ident.1 = "11) Lumenal", ident.2 = c("21) Glandular"), test.use = "MAST", only.pos = TRUE )
#List4B <- FindMarkers(Dats,ident.1 = "31) Ciliated",ident.2 = c("11) Lumenal"), test.use = "MAST", only.pos = TRUE )

#List4Bo <- FindMarkers(Dats,ident.2 = "31) Ciliated",ident.1 = c("11) Lumenal"), test.use = "MAST", only.pos = TRUE )


newdata1B <- List1B[order(List1B$avg_log2FC,decreasing = TRUE),]
newdata2B <- List2Bs[order(List2Bs$avg_log2FC,decreasing = TRUE),]
newdata4Bs <- List4Bs[order(List4Bs$avg_log2FC,decreasing = TRUE),]

cl1 <- intersect(rownames(List1B),rownames(List1))
cl2 <- intersect(rownames(List2Bs),rownames(List2))
cl4 <- intersect(rownames(List4Bs),rownames(List4))
cl4o <- intersect(rownames(List4Bo),rownames(List4o))

cl5 <- intersect(rownames(List1B),rownames(List5))
cl6 <- intersect(rownames(List2Bs),rownames(List6))

newdata1t <- List1[cl1,]
newdata1t <- newdata1t[order(newdata1t$avg_log2FC,decreasing = TRUE),]
newdata2t <- List2[cl2,]
newdata2t <- newdata2t[order(newdata2t$avg_log2FC,decreasing = TRUE),]
#newdata3t <- List3[cl2,]
#newdata3t <- List3[order(List3$avg_log2FC,decreasing = TRUE),]
newdata4t <- List4[cl4,]
newdata4t <- newdata4t[order(newdata4t$avg_log2FC,decreasing = TRUE),]
newdata5t <- List5[cl5,]
newdata5t <- newdata5t[order(newdata5t$avg_log2FC,decreasing = TRUE),]
newdata6t <- List6[cl6,]
newdata6t <- newdata6t[order(newdata6t$avg_log2FC,decreasing = TRUE),]


newdata1t2 <- List1B[cl1,]
newdata1t2 <- newdata1t2[order(newdata1t2$avg_log2FC,decreasing = TRUE),]
newdata2t2 <- List2Bs[cl2,]
newdata2t2 <- newdata2t2[order(newdata2t2$avg_log2FC,decreasing = TRUE),]
newdata4t2 <- List4Bs[cl4,]
newdata4t2 <- newdata4t2[order(newdata4t2$avg_log2FC,decreasing = TRUE),]
newdata5t2 <- List1B[cl5,]
newdata5t2 <- newdata5t2[order(newdata5t2$avg_log2FC,decreasing = TRUE),]
newdata6t2 <- List2Bs[cl6,]
newdata6t2 <- newdata6t2[order(newdata6t2$avg_log2FC,decreasing = TRUE),]



DefaultAssay(AllData4) <- "RNA"
p1<-VlnPlot(AllData4,  idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata1t)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_commonset1.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata2t)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_commonset2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata4t)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_commonassemb_set4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata5t)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_commonassemb_set5.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata6t)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_commonassemb_set6.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

p1<-VlnPlot(AllData4,  idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata1t2)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_commonset1_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata2t2)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_assemb_commonset2_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata4t2)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_commonassemb_set4_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata5t2)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_commonassemb_set5_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("11) Lumenal","21) Glandular","31) Ciliated","12) proliferativeLumenal","22) proliferativeGlandular","32) proliferativeCiliated","13) early-secretoryLumenal","23) early-secretoryGlandular","33) early-secretoryCiliated","24) late-secretoryGlandular","34) late-secretoryCiliated","41) Epi1","42) Epi2","51) EpS1","52) EpS2","53) EpS3","54) EpS4","55) EpS5"),features = rownames(newdata6t2)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Epithelial_commonassemb_set6_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

Dats <- subset(AllData4,idents=c("21) Glandular","11) Lumenal","31) Ciliated"))
DefaultAssay(Dats) <- "RNA"

#List1Ca <- FindMarkers(Dats,ident.2 = "21) Glandular", ident.1 = "11) Lumenal", test.use = "MAST", only.pos = TRUE)
List1C <- FindMarkers(AllData4,ident.1 = "21) Glandular", ident.2 = "11) Lumenal", test.use = "MAST", only.pos = FALSE)#, min.cells.group = 1, 
  #                    min.cells.feature = 1,
  #                    min.pct = 0,
  #                    logfc.threshold = 0)
List4C <- FindMarkers(AllData4,ident.1 = c("31) Ciliated"),ident.2 = c("11) Lumenal","21) Glandular"), test.use = "MAST", only.pos = FALSE)#, min.cells.group = 1, 
   #                   min.cells.feature = 1,
  #                    min.pct = 0,
  #                    logfc.threshold = 0)



#AllData5 <- AllData4
#Idents(AllData5)
Av <- AverageExpression(AllData4)
Cl1 <- List1C
Ae5 <- as.data.frame(Av$RNA)
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5[rownames(Cl1),"Pval1"] <- -log2(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_log2FC
pospos1 <- which( abs(Ae5$FC1)>log(1.1) & Ae5$Pval1>-log2(0.01))
Ae5$ColInd <- 0

Ae5[ which( (Ae5$FC1)>log(1.1) & Ae5$Pval1> -log2(0.01)), c("ColInd")] <- 1
Ae5[ which( (Ae5$FC1)< -log(1.1) & Ae5$Pval1> -log2(0.01)), c("ColInd")] <- 2
Ae5$ColInd <- as.factor(Ae5$ColInd)
Ae5$AvExp <- log(0.5*(Ae5$'21) Glandular' + Ae5$'11) Lumenal') +1)
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
  geom_hline(aes(yintercept = 0))

genes.to.label <- genes.to.label <- unique(c(cl1,cl2,cl5,cl6) )

p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_Gland_vs_Lumenal_newnew.pdf",sep=""),width = 20, height = 10, plot = p1, useDingbats=FALSE)
dev.off()


AllData5 <- AllData4
Idents(AllData5,cells=WhichCells(AllData5,idents=c("21) Glandular"))) <- "11) Lumenal"



#AllData5 <- AllData4
#Idents(AllData5)
Av <- AverageExpression(AllData5)
Cl1 <- List4C
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
Ae5$AvExp <- log(0.5*(Ae5$'31) Ciliated' + Ae5$'11) Lumenal') +1)
p1 <- ggplot(Ae5[,], aes(AvExp,FC1)) + geom_point(aes(color=ColInd)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#cc1337','#652582')) +
  geom_hline(aes(yintercept = 0))

genes.to.label <- genes.to.label <- unique(c(cl4,cl4o) )

p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = 'black')
p1 <- p1 + labs(x = "log Average Exp", y = "log FC")
ggsave(filename=paste(saveext,"MA_Ciliated_vs_Other_newnew.pdf",sep=""),width = 20, height = 10, plot = p1, useDingbats=FALSE)
dev.off()


DefaultAssay(AllData4) <- "RNA"
p1<-FeatureScatter(AllData4,feature1 = "KRT17",feature2 = "SCGB2A1", cells = WhichCells(AllData4,idents = c("21) Glandular","11) Lumenal")))
ggsave(filename=paste(saveext,"MA_KRT17_SCGB2A1_CiliatedandGland_vs_Luminal.pdf",sep=""),width = 10, height = 10, plot = p1, useDingbats=FALSE)


Idents(AllData4,cells=WhichCells(AllData4,idents=c("42) Epi2"))) <- "41) Epi1"




#
List_imp1 <- FindMarkers(AllData4,ident.1 = c("41) Epi1"),ident.2 = c("24) late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )
List_imp2 <- FindMarkers(AllData4,ident.2 = c("41) Epi1"),ident.1 = c("24) late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )

List_imp1 <- List_imp1[order(List_imp1$avg_log2FC,decreasing = TRUE),]
List_imp2 <- List_imp2[order(List_imp2$avg_log2FC,decreasing = TRUE),]

List_imp3 <- FindMarkers(AllData4,ident.1 = c("21) Glandular"),ident.2 = c("24) late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )
List_imp4 <- FindMarkers(AllData4,ident.2 = c("21) Glandular"),ident.1 = c("24) late-secretoryGlandular"), test.use = "MAST", only.pos = TRUE )

List_imp3 <- List_imp3[order(List_imp3$avg_log2FC,decreasing = TRUE),]
List_imp4 <- List_imp4[order(List_imp4$avg_log2FC,decreasing = TRUE),]

cl1 <- intersect(rownames(List_imp1),rownames(List_imp3))
cl2 <- intersect(rownames(List_imp2),rownames(List_imp4))

newdata1t <- List_imp1[cl1,]
newdata1t <- newdata1t[order(newdata1t$avg_log2FC,decreasing = TRUE),]
newdata2t <- List_imp2[cl2,]
newdata2t <- newdata2t[order(newdata2t$avg_log2FC,decreasing = TRUE),]

newdata3t <- List_imp3[cl1,]
newdata3t <- newdata3t[order(newdata3t$avg_log2FC,decreasing = TRUE),]
newdata4t <- List_imp4[cl2,]
newdata4t <- newdata4t[order(newdata4t$avg_log2FC,decreasing = TRUE),]

write.table(as.data.frame(newdata1t),file=paste(saveext,"GlandImplantationDE1.rds",sep=""),sep=",",quote = FALSE)
write.table(as.data.frame(newdata2t),file=paste(saveext,"GlandImplantationDE2.rds",sep=""),sep=",",quote = FALSE)
write.table(as.data.frame(newdata3t),file=paste(saveext,"GlandImplantationDE3.rds",sep=""),sep=",",quote = FALSE)
write.table(as.data.frame(newdata4t),file=paste(saveext,"GlandImplantationDE4.rds",sep=""),sep=",",quote = FALSE)

AllData5 <- AllData4
AllData5 <- subset(AllData5,idents=c("11) Lumenal","21) Glandular","41) Epi1","12) proliferativeLumenal","22) proliferativeGlandular","24) late-secretoryGlandular"))
Idents(AllData5) <- factor(x = Idents(AllData5), levels = c("11) Lumenal","21) Glandular","41) Epi1","12) proliferativeLumenal","22) proliferativeGlandular","24) late-secretoryGlandular"))


DefaultAssay(AllData5) <- "RNA"
AllData5$ID4 <- NULL
AllData5$ID1 <- NULL
AllData5$ID2 <- NULL
AllData5$ID3 <- NULL

p1<-VlnPlot(AllData5,  idents = c("11) Lumenal","21) Glandular","41) Epi1","12) proliferativeLumenal","22) proliferativeGlandular","24) late-secretoryGlandular"),features = rownames(newdata1t)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/AllImpl_DownReg1.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

p1<-VlnPlot(AllData5,  idents = c("11) Lumenal","21) Glandular","41) Epi1","12) proliferativeLumenal","22) proliferativeGlandular","24) late-secretoryGlandular"),features = rownames(newdata2t)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/AllImpl_DownReg2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

p1<-VlnPlot(AllData5,  idents = c("11) Lumenal","21) Glandular","41) Epi1","12) proliferativeLumenal","22) proliferativeGlandular","24) late-secretoryGlandular"),features = rownames(newdata3t)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/AllImpl_DownReg3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

p1<-VlnPlot(AllData5,  idents = c("11) Lumenal","21) Glandular","41) Epi1","12) proliferativeLumenal","22) proliferativeGlandular","24) late-secretoryGlandular"),features = rownames(newdata4t)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/AllImpl_DownReg4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)


onlyOurs <- subset(mammal.combined,idents=c("C1","C2","C3","C4","C7","C6","C8"))
Idents(onlyOurs) <- onlyOurs$Cells
onlyOurs <- subset(onlyOurs,idents = c("Stromal fibroblasts"),invert = FALSE)

DefaultAssay(onlyOurs) <- "RNA"
AllData3 <- merge(onlyOurs,y=c(Endo4,Endo0,assemb),project = "Merged")

AllData4 <- subset(AllData3,idents=c("proliferativeeS","early-secretoryeS","late-secretoryeS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","Stromal fibroblasts"))
Idents(AllData4) <- factor(x = Idents(AllData4), levels = c("proliferativeeS","early-secretoryeS","late-secretoryeS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","Stromal fibroblasts"))




#List_imp1 <- FindMarkers(AllData3,ident.1 = c("dS1","dS2","dS3"),ident.2 = c("late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
#List_imp2 <- FindMarkers(AllData3,ident.2 = c("dS1","dS2","dS3"),ident.1 = c("late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
#List_imp1 <- List_imp1[order(List_imp1$avg_log2FC,decreasing = TRUE),]
#List_imp2 <- List_imp2[order(List_imp2$avg_log2FC,decreasing = TRUE),]

List_imp1 <- FindMarkers(AllData4,ident.1 = c("dS1","dS2","dS3"),ident.2 = c("late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List_imp2 <- FindMarkers(AllData4,ident.2 = c("dS1","dS2","dS3"),ident.1 = c("late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )

saveRDS(List_imp1,file=paste(saveext,"/Stroma_implantation_vs_Sectretory.rds",sep=""))
saveRDS(List_imp2,file=paste(saveext,"/Stroma_Sectretory_vs_implantation.rds",sep=""))


List_imp3 <- FindMarkers(AllData4,ident.1 = c("dS1","dS2","dS3"),ident.2 = c("proliferativedS","proliferativeeS"), test.use = "MAST", only.pos = TRUE )
List_imp4 <- FindMarkers(AllData4,ident.2 = c("dS1","dS2","dS3"),ident.1 = c("proliferativedS","proliferativeeS"), test.use = "MAST", only.pos = TRUE )
saveRDS(List_imp3,file=paste(saveext,"/Stroma_implantation_vs_prolif.rds",sep=""))
saveRDS(List_imp4,file=paste(saveext,"/Stroma_Sectretory_vs_prolif.rds",sep=""))

List_imp1 <- List_imp1[order(List_imp1$avg_log2FC,decreasing = TRUE),]
List_imp2 <- List_imp2[order(List_imp2$avg_log2FC,decreasing = TRUE),]

List_imp3 <- FindMarkers(AllData4,ident.1 = c("Stromal fibroblasts"),ident.2 = c("late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List_imp4 <- FindMarkers(AllData4,ident.2 = c("Stromal fibroblasts"),ident.1 = c("late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE)

List_imp3 <- List_imp3[order(List_imp3$avg_log2FC,decreasing = TRUE),]
List_imp4 <- List_imp4[order(List_imp4$avg_log2FC,decreasing = TRUE),]

cl1 <- intersect(rownames(List_imp1),rownames(List_imp3))
cl2 <- intersect(rownames(List_imp2),rownames(List_imp4))

newdata1t <- List_imp1[cl1,]
newdata1t <- newdata1t[order(newdata1t$avg_log2FC,decreasing = TRUE),]
newdata2t <- List_imp2[cl2,]
newdata2t <- newdata2t[order(newdata2t$avg_log2FC,decreasing = TRUE),]

newdata3t <- List_imp3[cl1,]
newdata3t <- newdata3t[order(newdata3t$avg_log2FC,decreasing = TRUE),]
newdata4t <- List_imp4[cl2,]
newdata4t <- newdata4t[order(newdata4t$avg_log2FC,decreasing = TRUE),]

AllData5 <- AllData4
Idents(AllData5,cells=WhichCells(AllData5,idents=c("dS2","dS3"))) <- "dS1"
Idents(AllData5,cells=WhichCells(AllData5,idents=c("proliferativeeS","proliferativedS"))) <- "pStrom"
Idents(AllData5,cells=WhichCells(AllData5,idents=c("late-secretorydS","late-secretoryeS"))) <- "secStrom"



write.table(as.data.frame(newdata1t),file=paste(saveext,"StromaImplantationDE1.rds",sep=""),sep=",",quote = FALSE)
write.table(as.data.frame(newdata2t),file=paste(saveext,"StromaImplantationDE2.rds",sep=""),sep=",",quote = FALSE)
write.table(as.data.frame(newdata3t),file=paste(saveext,"StromaImplantationDE3.rds",sep=""),sep=",",quote = FALSE)
write.table(as.data.frame(newdata4t),file=paste(saveext,"StromaImplantationDE4.rds",sep=""),sep=",",quote = FALSE)


DefaultAssay(AllData5) <- "RNA"
AllData5$ID4 <- NULL
AllData5$ID1 <- NULL
AllData5$ID2 <- NULL
AllData5$ID3 <- NULL


Idents(AllData5) <- factor(x = Idents(AllData5), levels = c("pStrom","secStrom","Stromal fibroblasts","dS1","early-secretorydS","early-secretoryeS"))


p1<-VlnPlot(AllData5,  idents = c("pStrom","secStrom","Stromal fibroblasts","dS1"),features = rownames(newdata1t)[1:37], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/AllImplStroma_DownReg1.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

p1<-VlnPlot(AllData5,  idents = c("pStrom","secStrom","Stromal fibroblasts","dS1"),features = rownames(newdata2t)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/AllImplStroma_DownReg2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

p1<-VlnPlot(AllData5,  idents = c("pStrom","secStrom","Stromal fibroblasts","dS1"),features = rownames(newdata3t)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/AllImplStroma_DownReg3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)

p1<-VlnPlot(AllData5,  idents = c("pStrom","secStrom","Stromal fibroblasts","dS1"),features = rownames(newdata4t)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/AllImplStroma_DownReg4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)





#Not do stroma
List7 <- FindMarkers(Endo4,ident.2 = "proliferativeeS",ident.1 = c("late-secretorydS"), test.use = "MAST", only.pos = TRUE )
List8 <- FindMarkers(Endo4,ident.1 = "proliferativeeS",ident.2 = c("late-secretorydS"), test.use = "MAST", only.pos = TRUE )
List9 <- FindMarkers(Endo4,ident.1 = "proliferativeFibroblast C7",ident.2 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List10 <- FindMarkers(Endo4,ident.1 = "late-secretoryFibroblast C7",ident.2 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List11 <- FindMarkers(Endo4,ident.2 = "proliferativeFibroblast C7",ident.1 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List12 <- FindMarkers(Endo4,ident.2 = "late-secretoryFibroblast C7",ident.1 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )


List13 <- FindMarkers(Endo4,ident.2 = c("proliferativeeS","early-secretoryeS","late-secretoryeS"),ident.1 = c("proliferativedS","early-secretorydS","late-secretorydS"), test.use = "MAST", only.pos = TRUE )
List14 <- FindMarkers(Endo4,ident.1 = c("proliferativeeS","early-secretoryeS","late-secretoryeS"),ident.2 = c("proliferativedS","early-secretorydS","late-secretorydS"), test.use = "MAST", only.pos = TRUE )


newdata7 <- List7[order(List7$avg_log2FC,decreasing = TRUE),]
newdata8 <- List8[order(List8$avg_log2FC,decreasing = TRUE),]
newdata9 <- List9[order(List9$avg_log2FC,decreasing = TRUE),]
newdata10 <- List10[order(List10$avg_log2FC,decreasing = TRUE),]
newdata11 <- List11[order(List11$avg_log2FC,decreasing = TRUE),]
newdata12 <- List12[order(List12$avg_log2FC,decreasing = TRUE),]

newdata13 <- List13[order(List13$avg_log2FC,decreasing = TRUE),]
newdata14 <- List14[order(List14$avg_log2FC,decreasing = TRUE),]


#Now load our stroma

mammal.combined <- readRDS(file=paste(saveext,"/Integrate_with_RVT_stroma_withModuleScores.rds",sep=""))
Idents(mammal.combined) <- mammal.combined$ID3
onlyOurs <- subset(mammal.combined,idents=c("C1","C2","C3","C4","C7","C6","C8"))
uID <- as.character( Idents(onlyOurs) )
uID[1:length(uID)] <- "dS"
uID[which(onlyOurs$prolifS1 > onlyOurs$dS1)] <- "pS"
Idents(onlyOurs) <- uID
DefaultAssay(onlyOurs) <- "RNA"

uID <- as.character( Idents(onlyOurs) )
uID[1:length(uID)] <- "dS"
uID[which(onlyOurs$prolifS1 > onlyOurs$dS1)] <- "pS"
Idents(onlyOurs) <- uID
DefaultAssay(onlyOurs) <- "RNA"

AllData3 <- merge(onlyOurs,y=c(Endo4,Endo0,assemb),project = "Merged")

AllData4 <- subset(AllData3,idents=c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"))
Idents(AllData4) <- factor(x = Idents(AllData4), levels = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"))


DefaultAssay(AllData4) <- "RNA"
AllData4$ID4 <- NULL
AllData4$ID1 <- NULL
AllData4$ID2 <- NULL
AllData4$ID3 <- NULL

p1<-VlnPlot(AllData4,  idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata7)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set1.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata8)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
#p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata3)[1:40],  pt.size = 2)
#ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
#p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata4)[1:40],  pt.size = 2)
#ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,  idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata13)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata14)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)



p1<-VlnPlot(AllData4,  idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata7)[41:80], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set1_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata8)[41:80],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set2_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
#p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata3)[1:40],  pt.size = 2)
#ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
#p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata4)[1:40],  pt.size = 2)
#ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,  idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata13)[41:80], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set3_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata14)[41:80],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set4_2.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)



AllData3 <- merge(onlyOurs,y=c(Endo4,Endo0,assemb),project = "Merged")

AllData4 <- subset(AllData3,idents=c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"))
Idents(AllData4) <- factor(x = Idents(AllData4), levels = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"))


Idents(AllData4,cells=WhichCells(AllData4,idents=c("proliferativeeS","early-secretoryeS","late-secretoryeS"))) <- "proliferativeeS"
Idents(AllData4,cells=WhichCells(AllData4,idents=c("proliferativedS","early-secretorydS","late-secretorydS"))) <- "late-secretorydS"
Idents(AllData4,cells=WhichCells(AllData4,idents=c("dS2","dS3"))) <- "dS2"

AllData4 <- subset(AllData4,idents=c("pS","proliferativeeS","dS","late-secretorydS","dS1","dS2","SS1","SS2","SS3","SS4","SS5"))



Idents(AllData4) <- factor(x = Idents(AllData4), levels = c("pS","dS","proliferativeeS","late-secretorydS","dS1","dS2","SS1","SS2","SS3","SS4","SS5"))

DefaultAssay(AllData4) <- "RNA"
AllData4$ID4 <- NULL
AllData4$ID1 <- NULL
AllData4$ID2 <- NULL
AllData4$ID3 <- NULL

p1<-VlnPlot(AllData4,  idents = c("pS","proliferativeeS","dS","late-secretorydS","dS1","dS2","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata7)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set1_merge.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","dS","late-secretorydS","dS1","dS2","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata8)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set2_merge.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
#p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata3)[1:40],  pt.size = 2)
#ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set3.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
#p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","early-secretoryeS","late-secretoryeS","dS","proliferativedS","early-secretorydS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata4)[1:40],  pt.size = 2)
#ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set4.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,  idents = c("pS","proliferativeeS","dS","late-secretorydS","dS1","dS2","dS3","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata13)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set3_merge.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","dS","late-secretorydS","dS1","dS2","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata14)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set4_merge.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)



#Not do stroma
List7 <- FindMarkers(Endo4,ident.2 = "proliferativeeS",ident.1 = c("late-secretorydS"), test.use = "MAST", only.pos = TRUE )
List8 <- FindMarkers(Endo4,ident.1 = "proliferativeeS",ident.2 = c("late-secretorydS"), test.use = "MAST", only.pos = TRUE )
List9 <- FindMarkers(Endo4,ident.1 = "proliferativeFibroblast C7",ident.2 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List10 <- FindMarkers(Endo4,ident.1 = "late-secretoryFibroblast C7",ident.2 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List11 <- FindMarkers(Endo4,ident.2 = "proliferativeFibroblast C7",ident.1 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )
List12 <- FindMarkers(Endo4,ident.2 = "late-secretoryFibroblast C7",ident.1 = c("proliferativeeS","proliferativedS","early-secretorydS","late-secretorydS","late-secretoryeS"), test.use = "MAST", only.pos = TRUE )


List13 <- FindMarkers(Endo4,ident.2 = c("proliferativeeS","early-secretoryeS","late-secretoryeS"),ident.1 = c("proliferativedS","early-secretorydS","late-secretorydS"), test.use = "MAST", only.pos = TRUE )
List14 <- FindMarkers(Endo4,ident.1 = c("proliferativeeS","early-secretoryeS","late-secretoryeS"),ident.2 = c("proliferativedS","early-secretorydS","late-secretorydS"), test.use = "MAST", only.pos = TRUE )

List15 <- FindMarkers(onlyOurs,ident.2 = c("pS"),ident.1 = c("dS"), test.use = "MAST", only.pos = TRUE )
List16 <- FindMarkers(onlyOurs,ident.1 = c("pS"),ident.2 = c("dS"), test.use = "MAST", only.pos = TRUE )

int1 <- intersect(rownames(List13),rownames(List15))
int2 <- intersect(rownames(List14),rownames(List16))

newdata15 <- List13[int1,]
  
newdata15 <- newdata15[order(newdata15$avg_log2FC,decreasing = TRUE),]

newdata16 <- List14[int2,]
newdata16 <- newdata16[order(newdata16$avg_log2FC,decreasing = TRUE),]




p1<-VlnPlot(AllData4,  idents = c("pS","proliferativeeS","dS","late-secretorydS","dS1","dS2","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata15)[1:40], pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set3_mergecommon.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)
p1<-VlnPlot(AllData4,   idents = c("pS","proliferativeeS","dS","late-secretorydS","dS1","dS2","SS1","SS2","SS3","SS4","SS5"),features = rownames(newdata16)[1:40],  pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Stroma_assemb_set4_mergecommon.pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p1)



onlyOurs <- subset(mammal.combined,cells=ours)
notOurs <- subset(mammal.combined,cells=ours,invert=TRUE)

DefaultAssay(onlyOurs) <- "RNA"
DefaultAssay(notOurs) <- "RNA"
List1 <- FindMarkers(onlyOurs,ident.2 = "Lumenal",ident.1 = c("Glandular"), test.use = "MAST")
List2 <- FindMarkers(onlyOurs,ident.2 = "Lumenal",ident.1 = c("Ciliated"), test.use = "MAST")
List3 <- FindMarkers(onlyOurs,ident.2 = "Glandular",ident.1 = c("Ciliated"), test.use = "MAST")
List4 <- FindMarkers(notOurs,ident.2 = c("proliferative_Lumenal","early-secretory_Lumenal","late-secretory_Lumenal"),ident.1 = c("proliferative_Glandular","early-secretory_Glandular","late-secretory_Glandular"), test.use = "MAST", only.pos = TRUE )
List5 <- FindMarkers(notOurs,ident.2 = c("proliferative_Lumenal","early-secretory_Lumenal","late-secretory_Lumenal"),ident.1 = c("Ciliated"), test.use = "MAST", only.pos = TRUE )
List6 <- FindMarkers(notOurs,ident.2 = c("proliferative_Glandular","early-secretory_Glandular","late-secretory_Glandular"),ident.1 = c("Ciliated"), test.use = "MAST", only.pos = TRUE )


newdata1 <- List1[order(List1$avg_log2FC,decreasing = TRUE),]
newdata2 <- List2[order(List2$avg_log2FC,decreasing = TRUE),]
newdata3 <- List3[order(List3$avg_log2FC,decreasing = TRUE),]
newdata4 <- List4[order(List4$avg_log2FC,decreasing = TRUE),]
newdata5 <- List5[order(List5$avg_log2FC,decreasing = TRUE),]
newdata6 <- List6[order(List6$avg_log2FC,decreasing = TRUE),]


AvE <- AverageExpression(onlyOurs)
AvExp <- as.data.frame(AvE$RNA)
AvExp$'AvExp' <- NA
AvExp[,'AvExp'] <- (AvExp$Lumenal+AvExp$Glandular)/2
AvExp$'log2FC' <- NA

AvExp[rownames(List1),'log2FC'] <- List1$avg_log2FC
AvExp[rownames(List1),'pval'] <- List1$p_val_adj
AvExp$'IsIn' <- 0
AvExp[which(AvExp$pval<0.05),'IsIn'] <- 1
AvExp$'IsIn' <- as.factor(AvExp$'IsIn')

genes.to.label1 = unique(c(intersect(rownames(List1),rownames(List4))))

p1 <- ggplot(AvExp, aes(log2(AvExp),log2FC)) + geom_point(aes(color=IsIn)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Average expression", y = "log2 FC")
ggsave(filename=paste(saveext,"LumenGland_marmosethuman.pdf",sep=""),width = 13, height = 13, plot = p1)




AvExp <- as.data.frame(AvE$RNA)
AvExp$'AvExp' <- NA
AvExp[,'AvExp'] <- (AvExp$Lumenal+AvExp$Ciliated)/2
AvExp$'log2FC' <- NA

AvExp[rownames(List2),'log2FC'] <- List2$avg_log2FC
AvExp[rownames(List2),'pval'] <- List2$p_val_adj
AvExp$'IsIn' <- 0
AvExp[which(AvExp$pval<0.05),'IsIn'] <- 1
AvExp$'IsIn' <- as.factor(AvExp$'IsIn')

genes.to.label1 = unique(c(intersect(rownames(List2),rownames(List5))))

p1 <- ggplot(AvExp, aes(log2(AvExp),log2FC)) + geom_point(aes(color=IsIn)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Average expression", y = "log2 FC")
ggsave(filename=paste(saveext,"LumenCiliated_marmosethuman.pdf",sep=""),width = 13, height = 13, plot = p1)


AvExp <- as.data.frame(AvE$RNA)
AvExp$'AvExp' <- NA
AvExp[,'AvExp'] <- (AvExp$Glandular+AvExp$Ciliated)/2
AvExp$'log2FC' <- NA

AvExp[rownames(List3),'log2FC'] <- List3$avg_log2FC
AvExp[rownames(List3),'pval'] <- List3$p_val_adj
AvExp$'IsIn' <- 0
AvExp[which(AvExp$pval<0.05),'IsIn'] <- 1
AvExp$'IsIn' <- as.factor(AvExp$'IsIn')

genes.to.label1 = unique(c(intersect(rownames(List3),rownames(List6))))

p1 <- ggplot(AvExp, aes(log2(AvExp),log2FC)) + geom_point(aes(color=IsIn)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Average expression", y = "log2 FC")
ggsave(filename=paste(saveext,"GlandularCiliated_marmosethuman.pdf",sep=""),width = 13, height = 13, plot = p1)





