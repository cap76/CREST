library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(pheatmap)
library("readxl")
library(stringr)

set.seed(1) 


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
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="STB1")]) <- "STB"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="putSTB")]) <- "STB"

Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="Am")]) <- "Am"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="Am_2")]) <- "Am"

Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="CTB")]) <- "CTB"
#Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="EmDisc/PGC")]) <- "EmDisc/PGC"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="ExMes")]) <- "ExMes"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="Hypoblast")]) <- "Hypoblast"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="STB1")]) <- "STB"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="SYS")]) <- "SYS"
#Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="Am/PGC")]) <- "Am/PGC"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="EmDisc")]) <- "EmDisc"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="EmDisc_2")]) <- "EmDisc"

Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="Am/EmDisc")]) <- "Am/EmDisc"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="EVT")]) <- "EVT"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="ExMes/SYS")]) <- "ExMes/SYS"
Idents(Dsubset1,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="pPGC")]) <- "PGC"

Idents(Dsubset1,WhichCells(Dsubset1,idents=c("putSTB"))) <- "STB"

#Plot the UMAP/PCAs
p<-DimPlot(Dsubset1,  pt.size = 4, reduction = "umap", split.by = "Genotype2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Embryo_split.pdf",sep=""),width = 80, height = 8,p,limitsize = FALSE)
p<-DimPlot(Dsubset1,  pt.size = 4, reduction = "pca", split.by = "Genotype2", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Embryo_split.pdf",sep=""),width = 80, height = 8,p,limitsize = FALSE)

p<-DimPlot(D,  pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Embryo_nosplit.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


#Plot the UMAP/PCAs
p<-DimPlot(Dsubset1,  pt.size = 4, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Embryo_nosplit.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p<-DimPlot(Dsubset1,  pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Embryo_nosplit.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


Dsubset2 <- subset(Dsubset1,idents=c("PGC","EmDisc","Am","Am/EmDisc"))

p<-FeatureScatter(Dsubset2, feature1 = "SOX17", feature2 = "TFAP2C",  pt.size = 4) 
ggsave(filename=paste(saveext,"/DimRed/PCA_SOX17_TFAP2C.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


#MiniList <- c("SOX17","TFAP2A","TFAP2C","PRDM1","POU5F1","NANOS3","KIT","DND1","SOX15","KLF4")



MiniList <- c("SOX17","TFAP2A","TFAP2C","PRDM1","POU5F1","KLF4")
X <- GetAssayData(Dsubset2)
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(Dsubset2) )), Batch =  factor(( Idents(Dsubset2)))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/PGC_markerplot2",".pdf",sep="") ,width=15,height=3)



MiniList1 <- c("JARID2","CASC15","LINC01194","MT1G","MT1H","NLGN4X","RIMS2","CD9","DPPA4","CLU","NAV3","CDH1","CLDN4","PRLR","DMD","CLDN6","EZR","LGR5", "LHX1","OPHN1","SLC3A2","SLC2A3","TANC2","TCIM","COL5A2","COL6A3","TNC","SPARC","COL1A2","COL1A1","COL3A1","DCN","LUM","COL6A1","COL6A2","DPP4","CST3","DNAJC15","DNMT3A","DPPA3","OTX2","POU5F1","PDGFA","PODXL","SFRP1","SOX15","GABRP","VTCN1","WNT6","BAMBI","TFAP2C","TFAP2A","GATA3","CGA","SOX17","GATA4","CER1","LEFTY1","LEFTY2","APOA1","APOB","LINC00261","HHEX","FOXF1","PDGFRA","HAND2","BST2","GATA6","HGF","PRDM1","SNAI1","VIM","FN1","SOX17","TFAP2A","TFAP2C","PRDM1","POU5F1","KLF4")
MiniList2 <- c("DNMT3A","DPPA3","OTX2","NOTUM","JAM3","TACSTD2","ENPEP","GATA2","GATA3","TFAP2C","TFAP2A","TEAD3","TP63","OVOL1","CCKBR","CGA","CGB5","CGB8","SDC1","HLA-C","HLA-G","ASCL2","SNAI1","VIM","FN1","ELF5","MKI67","CDH1","FSTL1","FSTL3","KRT19","KRT7","SPARC","LAMC1","CTNNAL1")



MiniList <- c("SOX17","TFAP2A","TFAP2C","PRDM1","POU5F1","KLF4")
X <- GetAssayData(Dsubset2)
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList1 ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(Dsubset2) )), Batch =  factor(( Idents(Dsubset2)))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Emb_markerplot2",".pdf",sep="") ,width=15,height=20)



AvE <- AverageExpression(Dsubset1)
Expression1 <- AvE$integrated[,c("CTB","STB","EVT","ExMes","Hypoblast","SYS","Am","Am/EmDisc","EmDisc","PGC")]
Expression2 <- AvE$RNA[,c("CTB","STB","EVT","ExMes","Hypoblast","SYS","Am","Am/EmDisc","EmDisc","PGC")]

MiniList2 <- c("DNMT3A","DPPA3","OTX2","NOTUM","JAM3","TACSTD2","ENPEP","GATA2","GATA3","TFAP2C","TFAP2A","TEAD3","TP63","OVOL1","CCKBR","CGA","CGB5","CGB8","SDC1","HLA-C","HLA-G","ASCL2","SNAI1","VIM","FN1","ELF5","MKI67","CDH1","FSTL1","FSTL3","KRT19","KRT7","SPARC","LAMC1","CTNNAL1")
MiniList1 <- c("JARID2","CASC15","LINC01194","MT1G","MT1H","NLGN4X","RIMS2","CD9","DPPA4","CLU","NAV3","CDH1","CLDN4","PRLR","DMD","CLDN6","EZR","LGR5", "LHX1","OPHN1","SLC3A2","SLC2A3","TANC2","TCIM","COL5A2","COL6A3","TNC","SPARC","COL1A2","COL1A1","COL3A1","DCN","LUM","COL6A1","COL6A2","DPP4","CST3","DNAJC15","DNMT3A","DPPA3","OTX2","POU5F1","PDGFA","PODXL","SFRP1","SOX15","GABRP","VTCN1","WNT6","BAMBI","TFAP2C","TFAP2A","GATA3","CGA","SOX17","GATA4","CER1","LEFTY1","LEFTY2","APOA1","APOB","LINC00261","HHEX","FOXF1","PDGFRA","HAND2","BST2","GATA6","HGF","PRDM1","SNAI1","VIM","FN1")


Markers1 <- c("OPHN1",
              "EZR",
              "CD34",
              "TEAD3",
              "JAM3",
              "GATA2",
              "GATA3",
              "CGB8",
              "CGB5",
              "CGA",
              "SDC1",
              "ERVW-1",
              "PRDM6",
              "TBX3",
              "DIO2",
              "NOTUM",
              "HLA-G",
              "ASCL2",
              "SNAI1",
              "HGF",
              "SNAI2",
              "HAND1",
              "HAND2",
              "PDGFRA",
              "BST2",
              "SPARC",
              "TNC",
              "COL1A1",
              "COL1A2",
              "COL5A2",              
              "COL6A1",
              "COL6A2",
              "COL6A3",
              "GATA6",
              "GATA4",
              "CER1",
              "NODAL",
              "LEFTY1",
              "LEFTY2",
              "APOA1",
              "FOXF1",
              "SOX17",
              "PRDM1",
              "POU5F1",
              "TFAP2A",
              "TFAP2C",
              "DMD",
              "WNT6",
              "VTCN1",
              "BAMBI",
              "AKAP12",
              "PODXL",
              "SFRP1",
              "SOX15",
              "CD9",
              "CDH1",
              "DPPA4",
              "PDGFA",
              "GRID2",
              "TERF1",
              "L1TD1",
              "GAL",
              "AC009654.1",
              "GREB1L",
              "CD93",
              "CLDN6",
              "CLDN4")



mat_breaks <- seq(-3, 3, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap(Expression2[Markers1,],color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_selectmarkers.pdf",sep=""),scale="row",width=10,height=32)

Markers2 <- c("OPHN1",
              "CD34",
              "TEAD3",
              "JAM3",
              "GATA2",
              "GATA3",
              "CGB8",
              "CGB5",
              "CGA",
              "SDC1",
              "ERVW-1",
              "PRDM6",
              "TBX3",
              "DIO2",
              "NOTUM",
              "HLA-G",
              "ASCL2",
              "SNAI1",
              "HGF",
              "SNAI2",
              "HAND1",
              "HAND2",
              "PDGFRA",
              "BST2",
              "SPARC",
              "TNC",
              "COL1A1",
              "COL1A2",
              "COL5A2",              
              "COL6A3",
              "GATA4",
              "NODAL",
              "LEFTY1",
              "APOA1",
              "FOXF1",
              "SOX17",
              "PRDM1",
              "POU5F1",
              "TFAP2A",
              "TFAP2C",
              "DMD",
              "WNT6",
              "VTCN1",
              "BAMBI",
              "AKAP12",
              "PODXL",
              "SFRP1",
              "SOX15",
              "DPPA4",
              "L1TD1",
              "GAL",
              "GREB1L",
              "CLDN6",
              "CLDN4")

mat_breaks <- seq(-3, 3, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap(Expression1[Markers2,],color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_selectmarkers_normexpression.pdf",sep=""),scale="row",width=10,height=32)


pheatmap(Expression2[setdiff(MiniList1 ,Markers1),],color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_selectmarkers2.pdf",sep=""),scale="row",width=10,height=32)
pheatmap(Expression2[setdiff(MiniList2, Markers1),],color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_selectmarkers3.pdf",sep=""),scale="row",width=10,height=32)



#X <- GetAssayData(mammal.combined, assay = "RNA")
Xh <- t(scale(t(AvE$RNA)))
Xp <- Xh[MiniList1,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(mammal.combined))), ID = factor(( mammal.combined$ID3)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,  cluster_rows=TRUE, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/EmbryoLineageLabelledExtendedList_1",".pdf",sep="") ,width=40,height=10)



#X <- GetAssayData(mammal.combined, assay = "RNA")
Xh <- t(scale(t(AvE$RNA)))
Xp <- Xh[MiniList1,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(mammal.combined))), ID = factor(( mammal.combined$ID3)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/EmbryoLineageLabelledExtendedList_1",".pdf",sep="") ,width=40,height=10)



Dsubset2 <- subset(Dsubset1,idents=c("ExMes","EmDisc","SYS","STB","Hypoblast","CTB","EmDisc","PGC","Am","EVT","Am/EmDisc"))
Idents(Dsubset2,cells=WhichCells(Dsubset2,idents = "CTB")) <- "1) CTB"
Idents(Dsubset2,cells=WhichCells(Dsubset2,idents = "STB")) <- "2) STB"
Idents(Dsubset2,cells=WhichCells(Dsubset2,idents = "EVT")) <- "3) EVT"

Idents(Dsubset2,cells=WhichCells(Dsubset2,idents = "ExMes")) <- "4) ExMes"
Idents(Dsubset2,cells=WhichCells(Dsubset2,idents = "Hypoblast")) <- "5) Hypo"
Idents(Dsubset2,cells=WhichCells(Dsubset2,idents = "SYS")) <- "6) SYS"

Idents(Dsubset2,cells=WhichCells(Dsubset2,idents = "EmDisc")) <- "7) EmDisc"
Idents(Dsubset2,cells=WhichCells(Dsubset2,idents = "Am/EmDisc")) <- "8) Am/EmDisc"
Idents(Dsubset2,cells=WhichCells(Dsubset2,idents = "Am")) <- "9) Am"
Idents(Dsubset2,cells=WhichCells(Dsubset2,idents = "PGC")) <- "99) PGC"

uID <- as.character(Idents(Dsubset2))
uID <- as.factor(uID)
Idents(Dsubset2) <- uID

DefaultAssay(Dsubset2) <- "RNA"
p <- VlnPlot(Dsubset2,features = Markers1)
ggsave(filename=paste(saveext,"/DimRed/AllVln",".pdf",sep=""),width = 100, height = 100, limitsize = FALSE,p)

