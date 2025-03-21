library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(pheatmap)
library("readxl")
library(stringr)

set.seed(1) 

set.seed(1) 


Ano <- read_excel("/Users/christopherpenfold/Desktop/Matteo_Anotations.xlsx", sheet = "Combined")

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


#Plot the UMAP/PCAs
p<-DimPlot(Dsubset1,  pt.size = 4, reduction = "umap", split.by = "Genotype2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Embryo_split.pdf",sep=""),width = 80, height = 8,p,limitsize = FALSE)
p<-DimPlot(Dsubset1,  pt.size = 4, reduction = "pca", split.by = "Genotype2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Embryo_split.pdf",sep=""),width = 80, height = 8,p,limitsize = FALSE)


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


marmoset_dataInVivo2$Lab <- newID
Idents(marmoset_dataInVivo2) <- newID

annotationL <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS2","ICM_CS3","Epi_CS3",
                 "EmDisc_CS5","EmDisc_CS6","EmDisc_CS7",
                 "PGC_CS5","PGC_CS6","PGC_CS7",
                 "Am_CS5","Am_CS6","Am_CS7",
                 "Hyp_CS3",
                 "VE_CS5","VE_CS6","VE_CS7",
                 "SYS_CS5","SYS_CS6","SYS_CS7",
                 "ExMes_CS5","ExMes_CS6","ExMes_CS7","Stalk_CS6","Stalk_CS7",
                 "Tb_CS3",
                 "Tb_CS5","Tb_CS6","Tb_CS7","Stroma","Gland","Remodelled","Myo")

Idents(marmoset_dataInVivo2) <- factor(Idents(marmoset_dataInVivo2), levels = annotationL)


marmoset_dataInVivo2 <- subset(marmoset_dataInVivo2,idents=c("ExMes_CS5","ExMes_CS6","ExMes_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Stalk_CS6","Stalk_CS7"))

#"Tb_CS5","Tb_CS6","Tb_CS7"
#"Tb_CS3","Epi_CS3","Hyp_CS3"
marmoset_dataInVivo2$Dataset <- "Marmoset" 

Dsubset2 <- subset(Dsubset1,idents=c("STB","STB1","CTB","EVT","putSTB"),invert=TRUE)

mammal.anchors <- FindIntegrationAnchors(object.list = list(Dsubset1,marmoset_dataInVivo2), reduction = "cca", normalization.method = c("LogNormalize"), dims = 1:20, anchor.features = 5000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll.pdf",sep=""),width = 20, height = 8,p,limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_1_3.pdf",sep=""),width = 20, height = 8,p,limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_umap.pdf",sep=""),width = 20, height = 8,p,limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_umap_nosplit.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_nosplit.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_tsne_nosplit.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


mammal.combined <- readRDS(file=paste(saveext,"/BatchAll_marmosetintegration.rds",sep=""))




cType <- c("Am/EmDisc","Hypoblast","ExMes/SYS","Am","EmDisc","PGC","VE","SYS","ExMes","None","Other","Gland_CS5","Gland_CS6","Gland","Stroma_CS5","Stroma","Remodelled_CS5","Remodelled_CS6","Remodelled","Stalk_CS6","2307","2308","Other","Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS2","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","Stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","Am","Amnoid_bead","BMP_MEF","BMP_noMEF","EmD","EmDisc","ActA_MEF","ActA_noMEF","SB43_MEF","CHIR_MEF","FGF_noMEF","Am_CS5_PGC","Am_CS5_ExMesPGC","EmDisc_CS5_Am","PGC_CS6","Stalk_CS7_PGC","EmDisc_CS7_PGC","EmDisc_CS7_Am","Am_CS6_EmDisc","EmDisc_CS5_Gast","EmDisc_CS6_Am","EmDisc_CS6_Gast","EmDisc_CS7_Gast","Stalk_CS7_PGC","Gland_CS6_","ReStroma_CS5","ReStroma_CS6","ReStroma_CS7","Stroma_CS5","Stroma_CS6","Stroma_CS7","VE_CS6","VE_CS7","Am_CS7_EmDisc","EmDisc_Gast_CS6","EmDisc_Gast_CS7","Other")
BaseCol <- c("#8779e0","#fa530a","#e3c709","#877bd6","#0c9cf5","#E6E600","#F04C04", "#e08402","#e6c800","#B3B2B2","#f5f2d0","#B3B2B2","#A5A4A3","#969593","#DFDFDF","#D0D0D0","#878684","#797775","#6A6866","#754C24","lightgrey","lightgrey","lightgrey","#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#603813","#E6E600","#BFBF04","#999903","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#7b3294","#c2a5cf","#a6dba0","#008837","#ca0020","#f4a582","#92c5de","#0571b0","#d8b365","#5ab4ac","#4d4d4d","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#0233BF","#5F54C7","#0c9cf5","#0767DA","#0767DA","#0233BF","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0",'#D74404',"black","#5F54C7","#0767DA","#0233BF","#754C24,","#f5f2d0")


colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

p<-DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_All_ID1.pdf",sep=""),width = 20, height = 8,p)


p<-DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_All_ID1_1_3.pdf",sep=""),width = 20, height = 8,p)


p<-DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, dims = c(2,3)) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_All_ID1_2_3.pdf",sep=""),width = 20, height = 8,p)


p<-DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, dims = c(1,4)) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_All_ID1_1_4.pdf",sep=""),width = 20, height = 8,p)

DE1 <- FindMarkers(marmoset_dataInVivo,ident.1 = c("VE_CS5","VE_CS6","VE_CS7"), ident.2 = c("SYS_CS5","SYS_CS6","SYS_CS7"),test.use = "MAST", only.pos = TRUE)
DE2 <- FindMarkers(marmoset_dataInVivo,ident.2 = c("VE_CS5","VE_CS6","VE_CS7"), ident.1 = c("SYS_CS5","SYS_CS6","SYS_CS7"),test.use = "MAST", only.pos = TRUE)

newdata1 <- DE1[order(DE1$avg_log2FC,decreasing = TRUE),]
newdata2 <- DE2[order(DE2$avg_log2FC,decreasing = TRUE),]

DefaultAssay(mammal.combined) <- "RNA"
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata1)[1:100] ), name = "VEScore")
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata2)[1:100] ), name = "SYSScore")
DefaultAssay(mammal.combined) <- "integrated"


mammal.combined2 <- subset(mammal.combined,idents=c("VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","SYS","Hypoblast"))


library(destiny)
DN1 <- as.data.frame(mammal.combined2[["pca"]]@cell.embeddings)[,1:20]
dm <- DiffusionMap(DN1)
DM <- dm@eigenvectors*100
rownames(DM) <- colnames(mammal.combined2)
mammal.combined2[["dm"]] <- CreateDimReducObject(embeddings = DM, key = "DM_", assay = DefaultAssay(mammal.combined2))

p <- DimPlot(mammal.combined2, reduction = "dm", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/VE_SYS_DM",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE,p)

p <- DimPlot(mammal.combined2, reduction = "dm", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) +xlim(-5,5) +ylim(-7,7) 
ggsave(filename=paste(saveext,"/DimRed/VE_SYS_DM2",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE,p)


colind <- integer( length( levels(Idents(mammal.combined2)) )  )
for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
}
coluse <- BaseCol[colind]


cluster_letters <- as.factor(cluster_letters)
names(cluster_letters)=rownames(mammal.combined@meta.data)
names[which(mammal.combined=="")]
mammal.combined <- AddMetaData(object = mammal.combined, metadata = cluster_letters,col.name = 'cell.orig')

#DimPlot(mammal.combined,  reduction = "umap", pt.size = 2, shape.by = "cell.orig", split.by = "species", label = TRUE, repel = TRUE)  + NoLegend()


p <- DimPlot(mammal.combined2, cols = coluse, shape.by="species", reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/VE_SYS_PCA_rr",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE,p)


p<-FeaturePlot(mammal.combined2,  reduction = "pca", features = "VEScore1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/VE_SYS_PCA_VEScore.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined2,  reduction = "pca", features = "SYSScore1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/VE_SYS_PCA_SYSScore.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)


mammal.combined3 <- subset(mammal.combined2,idents=c("Hypoblast","SYS"))
p<-FeaturePlot(mammal.combined3,  reduction = "pca", features = "VEScore1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/VE_SYS_PCA_VEScoreours.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined3,  reduction = "pca", features = "SYSScore1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/VE_SYS_PCA_SYSScoreours.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)



DefaultAssay(mammal.combined3) <- "RNA"

DE3 <- FindMarkers(mammal.combined3,ident.1 = c("Hypoblast"), ident.2 = c("SYS"),test.use = "MAST", only.pos = TRUE)
DE4 <- FindMarkers(mammal.combined3,ident.2 = c("Hypoblast"), ident.1 = c("SYS"),test.use = "MAST", only.pos = TRUE)


newdata3 <- DE3[order(DE3$avg_log2FC,decreasing = TRUE),]
newdata4 <- DE4[order(DE4$avg_log2FC,decreasing = TRUE),]


uList1 <- intersect(rownames(newdata3),rownames(newdata1))
uList2 <- intersect(rownames(newdata4),rownames(newdata2))

newdata3 <- DE3[uList1,]
newdata4 <- DE4[uList2,]

newdata3 <- newdata3[order(newdata3$avg_log2FC,decreasing = TRUE),]
newdata4 <- newdata4[order(newdata4$avg_log2FC,decreasing = TRUE),]


newdata1 <- DE1[uList1,]
newdata2 <- DE2[uList2,]

newdata1 <- newdata1[order(newdata1$avg_log2FC,decreasing = TRUE),]
newdata2 <- newdata2[order(newdata2$avg_log2FC,decreasing = TRUE),]



#p1<-FeaturePlot(mammal.combined, features = "Gland1", pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/Gland_VE.pdf",sep=""),width = 70, height = 8,p1,limitsize = FALSE)

#"LEFTY2" "CER1"   "APOE"   "AUTS2"  "FZD5"   "APOC1"  "FN1"    "CDH2"   "ROR1"   "FLT1" 
p1<-FeaturePlot(mammal.combined2, features = "LEFTY2", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_LEFTY2.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "CER1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_CER1.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "APOE", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_APOE.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "AUTS2", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_AUTS2.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


p1<-FeaturePlot(mammal.combined2, features = "FZD5", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_FZD5.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "APOC1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_APOC1.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


p1<-FeaturePlot(mammal.combined2, features = "FN1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_FN1.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)



p1<-FeaturePlot(mammal.combined2, features = "LEFTY2", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_LEFTY2.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "CER1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_CER1.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "APOE", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_APOE.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "AUTS2", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_AUTS2.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


p1<-FeaturePlot(mammal.combined2, features = "FZD5", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_FZD5.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "APOC1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_APOC1.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


p1<-FeaturePlot(mammal.combined2, features = "FN1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_FN1.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)



DefaultAssay(mammal.combined2) <- "RNA"

p1<-FeaturePlot(mammal.combined2, features = "C15orf48", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_C15orf48.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "SLPI", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_SLPI.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "S100A9", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_S100A9.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "TAGLN", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_TAGLN.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "SAT1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_SAT1.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)




p1<-FeaturePlot(mammal.combined2, features = "G0S2", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_G0S2.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


p1<-FeaturePlot(mammal.combined2, features = "MT1X", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_MT1X.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


p1<-FeaturePlot(mammal.combined2, features = "C19orf33", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_C19orf33.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


p1<-FeaturePlot(mammal.combined2, features = "S100A6", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_S100A6.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)





p1<-FeaturePlot(mammal.combined2, features = "ANXA2", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_ANXA2.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


p1<-FeaturePlot(mammal.combined2, features = "B2M", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_B2M.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


p1<-FeaturePlot(mammal.combined2, features = "ANXA1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_ANAXA1.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "COTL1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_COTL1.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


p1<-FeaturePlot(mammal.combined2, features = "NEDD9", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_NEDD9.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)



p1<-VlnPlot(mammal.combined2, features = rownames(newdata4)[1:40])
ggsave(filename=paste(saveext,"/DimRed/VESYS_test.pdf",sep=""),width = 100, height = 100,p1,limitsize = FALSE)


p1<-VlnPlot(mammal.combined2, features = rownames(newdata4)[41:80])
ggsave(filename=paste(saveext,"/DimRed/VESYS_test2.pdf",sep=""),width = 100, height = 100,p1,limitsize = FALSE)





p1<-VlnPlot(mammal.combined2, features = rownames(newdata2)[1:40])
ggsave(filename=paste(saveext,"/DimRed/VESYS_test3.pdf",sep=""),width = 100, height = 100,p1,limitsize = FALSE)


p1<-VlnPlot(mammal.combined2, features = rownames(newdata2)[41:80])
ggsave(filename=paste(saveext,"/DimRed/VESYS_test4.pdf",sep=""),width = 100, height = 100,p1,limitsize = FALSE)


mammal.combined4 <- mammal.combined2
Idents(mammal.combined4,cells=WhichCells(mammal.combined4,idents=c("VE_CS5","VE_CS6","VE_CS7")) ) <- "VE_Marm"
Idents(mammal.combined4,cells=WhichCells(mammal.combined4,idents=c("SYS_CS5","SYS_CS6","SYS_CS7")) ) <- "SYS_Marm"  

AvExp <- AverageExpression(mammal.combined4)
  

Dplot <- AvExp$RNA[c(rownames(newdata1)[1:40], rownames(newdata2)[1:40]), ]

redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Dplot,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA, cluster_rows=TRUE, cluster_cols=TRUE, scale = "row",  filename = paste(saveext,"/DimRed/SYSMarkers",".pdf",sep="") ,width=10,height=10)



rawCS7 <- read.table("/Users/christopherpenfold/Desktop/Data/Embryonic/Human/Human_CS7/featurecountsCS7_hg19.csv",sep=",",header = T, row.names=1)
h2BS <- read.table("/Users/christopherpenfold/Desktop/Data/Embryonic/Human/Human_CS7/CS7Meta.csv",sep=",",header = T, row.names=1)
human_dataA2 <- CreateSeuratObject(counts = rawCS7, assay = "RNA",min.cells = 0, min.features = 0)
Idents(human_dataA2) <- h2BS$Cluster.assigned.lineage
human_dataA2$ID <- Idents(human_dataA2)
human_dataA2$ID2 <- Idents(human_dataA2)
human_dataA2$species <- "3) Human (in vitro)"
human_dataA2$divergence1 <- "3) Ape"
human_dataA2 <- NormalizeData(human_dataA2, verbose = FALSE)
human_dataA2 <- FindVariableFeatures(human_dataA2, selection.method = "vst", nfeatures = 20000)
#human_dataA2 <- subset(human_dataA2,idents=c("YS mesoderm","Ectoderm"))



DE5 <- FindMarkers(human_dataA2,ident.1 = c("DE (Caudal)","DE (Rostal)"), ident.2 = c("YS Endoderm"),test.use = "MAST", only.pos = TRUE)
DE6 <- FindMarkers(human_dataA2,ident.2 = c("DE (Caudal)","DE (Rostal)"), ident.1 = c("YS Endoderm"),test.use = "MAST", only.pos = TRUE)
#YS mesoderm

newdata5 <- DE5[order(DE1$avg_log2FC,decreasing = TRUE),]
newdata6 <- DE6[order(DE1$avg_log2FC,decreasing = TRUE),]


DefaultAssay(mammal.combined3) <- "RNA"

DE3 <- FindMarkers(mammal.combined3,ident.1 = c("Hypoblast"), ident.2 = c("SYS"),test.use = "MAST", only.pos = TRUE)
DE4 <- FindMarkers(mammal.combined3,ident.2 = c("Hypoblast"), ident.1 = c("SYS"),test.use = "MAST", only.pos = TRUE)


newdata3 <- DE3[order(DE3$avg_log2FC,decreasing = TRUE),]
newdata4 <- DE4[order(DE4$avg_log2FC,decreasing = TRUE),]

uList1 <- intersect(rownames(newdata3),rownames(newdata5))
uList2 <- intersect(rownames(newdata4),rownames(newdata6))


newdata5 <- DE5[uList1,]
newdata6 <- DE6[uList2,]

newdata5 <- newdata5[order(newdata5$avg_log2FC,decreasing = TRUE),]
newdata6 <- newdata6[order(newdata6$avg_log2FC,decreasing = TRUE),]


"CTSV"    "WFDC2"   "STMN1"   "CCKBR"   "KRT18"   "KRT8"    "CRABP1"  "TBCB"    "PAPSS1"  "FAM159B"



uList1 <- intersect(rownames(newdata3),rownames(newdata1))
uList2 <- intersect(rownames(newdata4),rownames(newdata2))

newdata3 <- DE3[uList1,]
newdata4 <- DE4[uList2,]

newdata3 <- newdata3[order(newdata3$avg_log2FC,decreasing = TRUE),]
newdata4 <- newdata4[order(newdata4$avg_log2FC,decreasing = TRUE),]


"CTSV"    "WFDC2"   "STMN1"   "CCKBR"   "KRT18"   "KRT8"    "CRABP1"  "TBCB"    "PAPSS1"  "FAM159B"


p1<-FeaturePlot(mammal.combined2, features = "CRABP1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_CRABP1.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "TBCB", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_TBCB.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "PAPSS1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_PAPSS1.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "FAM159B", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_FAM159B.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


p1<-FeaturePlot(mammal.combined2, features = "CTSV", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_CTSV.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "WFDC2", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_WFDC2.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "STMN1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_STMN1.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "CCKBR", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_CCKBR.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


p1<-FeaturePlot(mammal.combined2, features = "KRT18", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_KRT18.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined2, features = "KRT8", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/VESYS_KRT8.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)









DE5 <- FindMarkers(human_dataA2,ident.1 = c("DE (Caudal)","DE (Rostal)"), ident.2 = c("YS Endoderm"),test.use = "MAST", only.pos = TRUE)
DE6 <- FindMarkers(human_dataA2,ident.2 = c("DE (Caudal)","DE (Rostal)"), ident.1 = c("YS Endoderm"),test.use = "MAST", only.pos = TRUE)
#YS mesoderm

newdata5 <- DE5[order(DE5$avg_log2FC,decreasing = TRUE),]
newdata6 <- DE6[order(DE6$avg_log2FC,decreasing = TRUE),]



DefaultAssay(mammal.combined) <- "RNA"
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata5)[1:100] ), name = "VEScoreb")
mammal.combined <- AddModuleScore(mammal.combined, features = list( rownames(newdata6)[1:100] ), name = "SYSScoreb")
DefaultAssay(mammal.combined) <- "integrated"


mammal.combined2 <- subset(mammal.combined,idents=c("VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","SYS","Hypoblast"))


colind <- integer( length( levels(Idents(mammal.combined2)) )  )
for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
}
coluse <- BaseCol[colind]

p <- DimPlot(mammal.combined2, cols = coluse, reduction = "pca", pt.size = 2, label.size = 2, label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/VE_SYS_PCA",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE,p)

p<-FeaturePlot(mammal.combined2,  reduction = "pca", features = "VEScoreb1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/VE_SYS_PCA_VEScorBe.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined2,  reduction = "pca", features = "SYSScoreb1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/VE_SYS_PCA_SYSScoreB.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)


mammal.combined3 <- subset(mammal.combined2,idents=c("Hypoblast","SYS"))
p<-FeaturePlot(mammal.combined3,  reduction = "pca", features = "VEScoreb1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/VE_SYS_PCA_VEScoreoursb.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)
p<-FeaturePlot(mammal.combined3,  reduction = "pca", features = "SYSScoreb1", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/VE_SYS_PCA_SYSScoreoursb.pdf",sep=""),width = 12, height = 10, useDingbats = FALSE, limitsize = FALSE,p)



