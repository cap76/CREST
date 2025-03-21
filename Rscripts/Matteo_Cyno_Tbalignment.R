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

MiniList2b <- c("TACSTD2","HLA-C","KRT19","KRT7",'DNMT3A',"SNAI1","JAM3","CTNNAL1","MKI67","ENPEP","OTX2","ELF5","FN1","FSTL1","SPARC","VIM","CDH1","LAMC1","TP63","CGA","CGB5","CGB8","DPPA3",'CCKBR',"GATA3","OVOL1","TFAP2C","GATA2","SDC1","TFAP2A","TEAD3","NOTUM","ASCL2","HLA-G","FSTL3")



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

#Idents(Dsubset1,WhichCells(Dsubset1,idents="STB_d14")) <- "STB"
#Idents(Dsubset1,WhichCells(Dsubset1,idents="STB1")) <- "STB"

Dsubset2 <- subset(Dsubset1,idents=c("STB","CTB","EVT"))
Dsubset2B <- subset(Dsubset1,idents=c("SYS","Hypoblast","ExMes","ExMes/SYS"))

cyBS2<-read.table("/Users/christopherpenfold/Desktop/Data/Embryonic/Cyno/cyinvitLab.csv",sep=",", header = T, row.names=1)
raw_counts2<-read.table("/Users/christopherpenfold/Desktop/Data/Embryonic/Cyno/cyInVitData.csv",sep=",",header = T, row.names=1)
cynomolgous_dataB <- CreateSeuratObject(counts = raw_counts2, assay = "RNA",min.cells = 0, min.features = 0)
cynomolgous_dataB$species <- "cynomolgous"
cynomolgous_dataB$divergence1 <- "cynomolgous"
Idents(cynomolgous_dataB) <- cyBS2$Type
cynomolgous_dataB <- subset(cynomolgous_dataB, subset = nFeature_RNA > 0)

cynomolgous_dataB <- subset(cynomolgous_dataB,idents=c("SYS_CS6","SYS_CS7","VE_CS5","ExMes_CS6","ExMes_CS5","ExMes_CS7"),invert=FALSE)
cynomolgous_dataB <- NormalizeData(cynomolgous_dataB, verbose = FALSE)
cynomolgous_dataB <- FindVariableFeatures(cynomolgous_dataB, selection.method = "vst", nfeatures = 20000)
cynomolgous_dataB$Dataset <- "Cyno"



test <- readRDS("/Users/christopherpenfold/Desktop/Data/Irene/Data/DylanDM2.rds")
test$ID0 <- Idents(test)
Idents(test) <- test$species

Dsubset3 <- subset(test,idents=c("Cyno (in vitro)"))
Dsubset4 <- subset(test,idents=c("newCyno"))
Dsubset5 <- subset(test,idents=c("HumanCS5-7"))

Idents(Dsubset3) <- Dsubset3$ID0
Idents(Dsubset4) <- Dsubset4$ID0
Idents(Dsubset5) <- Dsubset5$ID0

DefaultAssay(Dsubset2) <- "RNA"
DefaultAssay(Dsubset3) <- "RNA"
DefaultAssay(Dsubset4) <- "RNA"
DefaultAssay(Dsubset5) <- "RNA"


Dsubset2$Dataset <- "Ours"
Dsubset3$Dataset <- "Cyno (10X)"
Dsubset4$Dataset <- "Cyno (SS2)"
Dsubset5$Dataset <- "Human (SS2)"


mammal.anchors <- FindIntegrationAnchors(object.list = list(Dsubset2,Dsubset3,Dsubset4,Dsubset5), reduction = "cca",normalization.method = c("LogNormalize"), dims = 1:20, anchor.features = 3000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_cyno_align_PCA.pdf",sep=""),width = 20, height = 8,p)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_cyno_align_PCA_1_3.pdf",sep=""),width = 20, height = 8,p)
#saveRDS(mammal.combined,file=paste(saveext,"/Batch1_marmosetintegration.rds",sep=""))

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tb_allhuman_cyno_align_PCA.pdf",sep=""),width = 20, height = 8,p)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/Tb_allhuman_cyno_align_PCA_1_3.pdf",sep=""),width = 20, height = 8,p)
#saveRDS(mammal.combined,file=paste(saveext,"/Batch1_marmosetintegration.rds",sep=""))


mammal.anchors <- FindIntegrationAnchors(object.list = list(Dsubset2,Dsubset3,Dsubset4), reduction = "cca",normalization.method = c("LogNormalize"), dims = 1:20, anchor.features = 3000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tb_cyno_align_PCA.pdf",sep=""),width = 20, height = 8,p)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/Tb_cyno_align_PCA_1_3.pdf",sep=""),width = 20, height = 8,p)
#saveRDS(mammal.combined,file=paste(saveext,"/Batch1_marmosetintegration.rds",sep=""))

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tb_cyno_align_PCA.pdf",sep=""),width = 20, height = 8,p)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/Tb_cyno_align_PCA_1_3.pdf",sep=""),width = 20, height = 8,p)
#saveRDS(mammal.combined,file=paste(saveext,"/Batch1_marmosetintegration.rds",sep=""))



DefaultAssay(mammal.combined) <- "RNA"
p1<-FeaturePlot(mammal.combined, features = "RYBP", pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_all_RYBP.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined, features = "PADI1", pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_all_PADI1.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined, features = "PADI2", pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_all_PADI2.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined, features = "PADI3", pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_all_PADI3.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined, features = "PADI4", pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_all_PADI4.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

#p1<-FeaturePlot(mammal.combined, features = "PADI5", pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
#ggsave(filename=paste(saveext,"/DimRed/Tb_all_PADI5.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined, features = "PADI6", pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_all_PADI6.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


mammal.combined$ID0 <- Idents(mammal.combined)
Idents(mammal.combined) <- paste(mammal.combined$Dataset,mammal.combined$ID0,sep="_")

DefaultAssay(mammal.combined) <- "RNA" #DE1 ss DE4 is 10X
DE1 <- FindMarkers(mammal.combined,ident.1 = c("Cyno (10X)_CTB"),ident.2 = c("Cyno (10X)_EVT","Cyno (10X)_STB"),test.use="MAST", only.pos = TRUE)
DE2 <- FindMarkers(mammal.combined,ident.1 = c("Cyno (10X)_STB"),ident.2 = c("Cyno (10X)_EVT","Cyno (10X)_CTB"),test.use="MAST", only.pos = TRUE)
DE3 <- FindMarkers(mammal.combined,ident.1 = c("Cyno (10X)_EVT"),ident.2 = c("Cyno (10X)_STB","Cyno (10X)_CTB"),test.use="MAST", only.pos = TRUE)


DE4 <- FindMarkers(mammal.combined,ident.1 = c("Cyno (SS2)_CTB"),ident.2 = c("Cyno (SS2)_EVT","Cyno (SS2)_STB"),test.use="MAST", only.pos = TRUE)
DE5 <- FindMarkers(mammal.combined,ident.1 = c("Cyno (SS2)_STB"),ident.2 = c("Cyno (SS2)_EVT","Cyno (SS2)_CTB"),test.use="MAST", only.pos = TRUE)
DE6 <- FindMarkers(mammal.combined,ident.1 = c("Cyno (SS2)_EVT"),ident.2 = c("Cyno (SS2)_STB","Cyno (SS2)_CTB"),test.use="MAST", only.pos = TRUE)

DE1a <- DE1[-c(grep("RPL",rownames(DE1)),grep("RPS",rownames(DE1))),]
DE2a <- DE2[-c(grep("RPL",rownames(DE2)),grep("RPS",rownames(DE2))),]
DE3a <- DE3 #[-c(grep("RPL",rownames(DE3)),grep("RPS",rownames(DE3))),]
DE4a <- DE4[-c(grep("RPL",rownames(DE4)),grep("RPS",rownames(DE4))),]
DE5a <- DE5[-c(grep("RPL",rownames(DE5)),grep("RPS",rownames(DE5))),]
DE6a <- DE6[-c(grep("RPL",rownames(DE6)),grep("RPS",rownames(DE6))),]


newdata1 <- DE1[order(DE1$avg_log2FC,decreasing = TRUE),]
newdata2 <- DE2[order(DE2$avg_log2FC,decreasing = TRUE),]
newdata3 <- DE3[order(DE3$avg_log2FC,decreasing = TRUE),]
newdata4 <- DE4[order(DE4$avg_log2FC,decreasing = TRUE),]
newdata5 <- DE5[order(DE5$avg_log2FC,decreasing = TRUE),]
newdata6 <- DE6[order(DE6$avg_log2FC,decreasing = TRUE),]


newdata1a <- DE1a[order(DE1a$avg_log2FC,decreasing = TRUE),]
newdata2a <- DE2a[order(DE2a$avg_log2FC,decreasing = TRUE),]
newdata3a <- DE3a[order(DE3a$avg_log2FC,decreasing = TRUE),]
newdata4a <- DE4a[order(DE4a$avg_log2FC,decreasing = TRUE),]
newdata5a <- DE5a[order(DE5a$avg_log2FC,decreasing = TRUE),]
newdata6a <- DE6a[order(DE6a$avg_log2FC,decreasing = TRUE),]

mammal.combined <- AddModuleScore(mammal.combined, features =  list( rownames(newdata1)[1:100] ) , name = "CTB10X")
mammal.combined <- AddModuleScore(mammal.combined, features =  list( rownames(newdata2)[1:100] ) , name = "STB10X")
mammal.combined <- AddModuleScore(mammal.combined, features =  list( rownames(newdata3)[1:100] ) , name = "EVT10X")
mammal.combined <- AddModuleScore(mammal.combined, features =  list( rownames(newdata4)[1:100] ) , name = "CTBSS")
mammal.combined <- AddModuleScore(mammal.combined, features =  list( rownames(newdata5)[1:100] ) , name = "STBSS")
mammal.combined <- AddModuleScore(mammal.combined, features =  list( rownames(newdata6)[1:100] ) , name = "EVTSS")
DefaultAssay(mammal.combined) <- "integrated"

#library(destiny)
#DN1 <- as.data.frame(mammal.combined[["pca"]]@cell.embeddings)[,1:20]
#dm <- DiffusionMap(DN1)
#DM <- dm@eigenvectors*100
#rownames(DM) <- colnames(mammal.combined)
#mammal.combined[["dm"]] <- CreateDimReducObject(embeddings = DM, key = "DM_", assay = DefaultAssay(mammal.combined))
#p <- DimPlot(mammal.combined, reduction = "dm", pt.size = 2, label.size = 4, split.by = "Dataset", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/Tb_allPrimate_DM",".pdf",sep=""),width = 20, height = 8, limitsize = FALSE,p)


OnlyOurs <- subset(mammal.combined,idents=c("Ours_STB","Ours_CTB","Ours_EVT"))

p1<-FeaturePlot(OnlyOurs, features = "CTB10X1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_CTB10X_ModuleScore.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

p1<-FeaturePlot(OnlyOurs, features = "STB10X1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_STB10X_ModuleScore.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(OnlyOurs, features = "EVT10X1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_EVT10X_ModuleScore.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)


newID1 <- as.character(Idents(OnlyOurs))
newID2 <- as.character(Idents(OnlyOurs))

#newID1 is SS
newID1[which(OnlyOurs$CTB10X1>OnlyOurs$STB10X1 & OnlyOurs$CTB10X1>OnlyOurs$EVT10X1)] <- "CTB"
newID1[which(OnlyOurs$STB10X1>OnlyOurs$CTB10X1 & OnlyOurs$STB10X1>OnlyOurs$EVT10X1)] <- "STB"
newID1[which(OnlyOurs$EVT10X1>OnlyOurs$STB10X1 & OnlyOurs$EVT10X1>OnlyOurs$CTB10X1)] <- "EVT"

#New ID2 is 10X
newID2[which(OnlyOurs$STBSS1>OnlyOurs$CTBSS1 & OnlyOurs$STBSS1>OnlyOurs$EVTSS1)] <- "STB"
newID2[which(OnlyOurs$EVTSS1>OnlyOurs$CTBSS1 & OnlyOurs$EVTSS1>OnlyOurs$STBSS1)] <- "EVT"
newID2[which(OnlyOurs$CTBSS1>OnlyOurs$STBSS1 & OnlyOurs$CTBSS1>OnlyOurs$EVTSS1)] <- "CTB"
OnlyOurs$newID1 <- newID1
OnlyOurs$newID2 <- newID2
Idents(OnlyOurs) <- newID1
p1<-DimPlot(OnlyOurs, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tb_ModuleScoreID_10X.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

p1<-DimPlot(OnlyOurs, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/Tb_ModuleScoreID_10X_1_3.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)


p1<-DimPlot(OnlyOurs, pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tb_ModuleScoreID_split_10X.pdf",sep=""),width = 60, height = 8,p1,limitsize = FALSE)

Idents(OnlyOurs) <- newID2
p1<-DimPlot(OnlyOurs, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tb_ModuleScoreID_SS2.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

p1<-DimPlot(OnlyOurs, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/Tb_ModuleScoreID_SS2_1_3.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)


p1<-DimPlot(OnlyOurs, pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tb_ModuleScoreID_split_SS2.pdf",sep=""),width = 60, height = 8,p1,limitsize = FALSE)

Idents(OnlyOurs) <- paste(newID1,newID2,sep="_")
p1<-DimPlot(OnlyOurs, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tb_ModuleScoreID_SS210X.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

newID3 <- paste(newID1,newID2,sep="_")
Idents(OnlyOurs) <- paste(newID1,newID2,sep="_")
p1<-DimPlot(OnlyOurs, pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tb_ModuleScoreID_split_SS210X.pdf",sep=""),width = 60, height = 8,p1,limitsize = FALSE)

MegaList1 <- c(rownames(newdata1a)[1:100],rownames(newdata2a)[1:100],rownames(newdata3a)[1:100] )
MegaList2 <- c(rownames(newdata4a)[1:100],rownames(newdata5a)[1:100],rownames(newdata6a)[1:100] )

DefaultAssay(OnlyOurs) <- "RNA"
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm6<-pheatmap(GetAssayData(OnlyOurs)[MegaList1,], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=TRUE, cluster_cols=TRUE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1",".pdf",sep="") ,width=30,height=30)
pm7<-pheatmap(GetAssayData(OnlyOurs)[MegaList2,],breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=TRUE, cluster_cols=TRUE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb2",".pdf",sep="") ,width=30,height=30)


oL2 <- c("DNMT3A",
         "TP63",
         "ELF5",
         "CDH1",
         "VIM",
         "CTNNAL1",
         "LAMC1",
         "TACSTD2",
         "DPPA3",
         "CCKBR",
         "ENPEP",
         "CGB8",
         "CGB5",
         "SDC1",
         "CGA",
         "OVOL1",
         "KRT19",
         "MKI67",
         "OTX2",
         "FN1",
         "FSTL3",
         "KRT7",
         "JAM3",
         "FSTL1",
         "NOTUM",
         "ASCL2",
         "HLA-G")

D <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Placenta/Placentaharm.rds")
Idents(D) <- D$ID2
#subs1 <- which(D$ID2%in%c("VCT","SCT","EVT"))
AvExp2 <- AverageExpression(D)
Expr2 <- AvExp2$RNA[oL2,c("VCT","SCT","EVT")]
#X <- (GetAssayData(D,assay = "RNA")) 
#Xh <- t(scale(t(X)))
#Xp <- Xh[MiniList2,]
#Xp <- na.omit(Xp)[,subs1]
Expr2 <- t(scale(t(Expr2)))
#Expr2 <- na.omit(Expr2)

Idents(OnlyOurs,cells=WhichCells(OnlyOurs,idents=c("CTB_STB","CTB_EVT"))) <- "CTB_CTB"

AvExp <- AverageExpression(OnlyOurs)
Expr <- AvExp$RNA[oL2,]
Expr <- t(scale(t(Expr)))
#Expr <- na.omit(Expr)
Expr <- Expr[,c("CTB_CTB","STB_STB","EVT_EVT")]
Expr <- t(scale(t(Expr)))

Exps <- cbind(Expr,Expr2)

redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm6<-pheatmap(Expr2,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Tb_VT_MinListPB1",".pdf",sep="") ,width=4,height=4)

pm6<-pheatmap(Expr,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Tb_VT_MinListPB2",".pdf",sep="") ,width=4,height=4)

Idents(OnlyOurs) <- OnlyOurs$newID1
DefaultAssay(OnlyOurs) <- "RNA"

DE7 <- FindMarkers(OnlyOurs,ident.1 = c("CTB"),ident.2 = c("EVT","STB"),test.use="MAST", only.pos = TRUE)
DE8 <- FindMarkers(OnlyOurs,ident.1 = c("STB"),ident.2 = c("EVT","CTB"),test.use="MAST", only.pos = TRUE)
DE9 <- FindMarkers(OnlyOurs,ident.1 = c("EVT"),ident.2 = c("STB","CTB"),test.use="MAST", only.pos = TRUE)
Idents(OnlyOurs) <- OnlyOurs$newID2
DE10 <- FindMarkers(OnlyOurs,ident.1 = c("CTB"),ident.2 = c("EVT","STB"),test.use="MAST", only.pos = TRUE)
DE11 <- FindMarkers(OnlyOurs,ident.1 = c("STB"),ident.2 = c("EVT","CTB"),test.use="MAST", only.pos = TRUE)
DE12 <- FindMarkers(OnlyOurs,ident.1 = c("EVT"),ident.2 = c("STB","CTB"),test.use="MAST", only.pos = TRUE)


LList1 <- intersect(rownames(DE1a),rownames(DE7)) #
LList2 <- intersect(rownames(DE2a),rownames(DE8))
LList3 <- intersect(rownames(DE3a),rownames(DE9))

LList4 <- intersect(rownames(DE4a),rownames(DE7))
LList5 <- intersect(rownames(DE5a),rownames(DE8))
LList6 <- intersect(rownames(DE6a),rownames(DE9))

LList7 <- intersect(rownames(DE1a),rownames(DE10))
LList8 <- intersect(rownames(DE2a),rownames(DE11))
LList9 <- intersect(rownames(DE3a),rownames(DE12))

LList10 <- intersect(rownames(DE4a),rownames(DE10))
LList11 <- intersect(rownames(DE5a),rownames(DE11))
LList12 <- intersect(rownames(DE6a),rownames(DE12))

LLList1 <- intersect(intersect(LList1,LList4),intersect(LList7,LList10))
LLList2 <- intersect(intersect(LList2,LList5),intersect(LList8,LList11))
LLList3 <- intersect(intersect(LList3,LList6),intersect(LList9,LList12))


newdata1 <- DE7[LList1,]
newdata1 <- newdata1[order(newdata1$avg_log2FC,decreasing = TRUE),]
newdata2 <- DE8[LList2,]
newdata2 <- newdata2[order(newdata2$avg_log2FC,decreasing = TRUE),]
newdata3 <- DE9[LList3,]
newdata3 <- newdata3[order(newdata3$avg_log2FC,decreasing = TRUE),]

newdata4 <- DE7[LList4,]
newdata4 <- newdata4[order(newdata4$avg_log2FC,decreasing = TRUE),]
newdata5 <- DE8[LList5,]
newdata5 <- newdata5[order(newdata5$avg_log2FC,decreasing = TRUE),]
newdata6 <- DE9[LList6,]
newdata6 <- newdata6[order(newdata6$avg_log2FC,decreasing = TRUE),]


newdata7 <- DE10[LList7,]
newdata7 <- newdata7[order(newdata7$avg_log2FC,decreasing = TRUE),]
newdata8 <- DE11[LList8,]
newdata8 <- newdata8[order(newdata8$avg_log2FC,decreasing = TRUE),]
newdata9 <- DE12[LList9,]
newdata9 <- newdata9[order(newdata9$avg_log2FC,decreasing = TRUE),]

newdata10 <- DE10[LList10,]
newdata10 <- newdata10[order(newdata10$avg_log2FC,decreasing = TRUE),]
newdata11 <- DE11[LList11,]
newdata11 <- newdata11[order(newdata11$avg_log2FC,decreasing = TRUE),]
newdata12 <- DE12[LList12,]
newdata12 <- newdata12[order(newdata12$avg_log2FC,decreasing = TRUE),]



newdata13 <- DE7[LLList1,]
newdata13 <- newdata13[order(newdata13$avg_log2FC,decreasing = TRUE),]
newdata14 <- DE8[LLList2,]
newdata14 <- newdata14[order(newdata14$avg_log2FC,decreasing = TRUE),]
newdata15 <- DE9[LLList3,]
newdata15 <- newdata15[order(newdata15$avg_log2FC,decreasing = TRUE),]



newdata16 <- DE10[LLList1,]
newdata16 <- newdata16[order(newdata16$avg_log2FC,decreasing = TRUE),]
newdata17 <- DE11[LLList2,]
newdata17 <- newdata17[order(newdata17$avg_log2FC,decreasing = TRUE),]
newdata18 <- DE12[LLList3,]
newdata18 <- newdata18[order(newdata18$avg_log2FC,decreasing = TRUE),]


MegaList1 <- c(rownames(newdata1)[1:100],rownames(newdata2)[1:100],rownames(newdata3)[1:100] )
MegaList2 <- c(rownames(newdata4)[1:100],rownames(newdata5)[1:100],rownames(newdata6)[1:100] )
MegaList3 <- c(rownames(newdata7)[1:100],rownames(newdata7)[1:100],rownames(newdata9)[1:100] )
MegaList4 <- c(rownames(newdata10)[1:100],rownames(newdata11)[1:100],rownames(newdata12)[1:100] )


MegaList1a <- c(rownames(newdata1a)[1:500],rownames(newdata2a)[1:500],rownames(newdata3a)[1:500] )



MegaList1A <- c(rownames(newdata1)[1:160],rownames(newdata2)[1:160],rownames(newdata3)[1:160] )
MegaList2B <- c(rownames(newdata4)[1:160],rownames(newdata5)[1:160],rownames(newdata6)[1:160] )
MegaList3C <- c(rownames(newdata7)[1:160],rownames(newdata7)[1:160],rownames(newdata9)[1:160] )
MegaList4D <- c(rownames(newdata10)[1:160],rownames(newdata11)[1:160],rownames(newdata12)[1:160] )




MegaList5 <- c(rownames(newdata13)[1:50],rownames(newdata14)[1:50],rownames(newdata15)[1:50] )
MegaList6 <- c(rownames(newdata16)[1:50],rownames(newdata17)[1:50],rownames(newdata17)[1:50] )

DefaultAssay(OnlyOurs) <- "RNA"
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)


annotation_col = data.frame(Stage = factor(OnlyOurs$newID1))
rownames(annotation_col) <- colnames(OnlyOurs)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))

annotation_col2 = data.frame(Stage = factor(OnlyOurs$newID2))
rownames(annotation_col2) <- colnames(OnlyOurs)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))


pm6<-pheatmap(GetAssayData(OnlyOurs)[MegaList1,], breaks = my_breaks, annotation_col = annotation_col, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=TRUE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1",".pdf",sep="") ,width=30,height=30)
pm7<-pheatmap(GetAssayData(OnlyOurs)[MegaList2,],breaks = my_breaks, annotation_col = annotation_col,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=TRUE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb2",".pdf",sep="") ,width=30,height=30)
pm7<-pheatmap(GetAssayData(OnlyOurs)[MegaList3,],breaks = my_breaks, annotation_col = annotation_col2,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=TRUE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb3",".pdf",sep="") ,width=30,height=30)
pm7<-pheatmap(GetAssayData(OnlyOurs)[MegaList4,],breaks = my_breaks, annotation_col = annotation_col2, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=TRUE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb4",".pdf",sep="") ,width=30,height=30)

pm7<-pheatmap(GetAssayData(OnlyOurs)[MegaList5,],breaks = my_breaks, annotation_col = annotation_col2, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=TRUE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb5",".pdf",sep="") ,width=30,height=30)


Idents(OnlyOurs) <- OnlyOurs$newID1
Av1 <- AverageExpression(OnlyOurs)
#Av1$R
Idents(OnlyOurs) <- OnlyOurs$newID2
Av2 <- AverageExpression(OnlyOurs)


pm6<-pheatmap(Av1$RNA[MegaList1a,c("CTB","STB","EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1a_",".pdf",sep="") ,width=6,height=30)



aa <- rownames(newdata1a)[1:500]
bb <- rownames(newdata2a)[1:500]
cc <- rownames(newdata3a)[1:500]

MegaList1again <- c(aa[which(rownames(newdata1a)[1:100]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)],
                    bb[which(rownames(newdata2a)[1:100]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)],
                    cc[which(rownames(newdata3a)[1:100]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)])



#MegaList1again <- c(rownames(newdata1)[1:100][which(rownames(newdata1)[1:100]%in%rownames(Dmat)==TRUE)],
#                    rownames(newdata2)[1:100][which(rownames(newdata2)[1:100]%in%rownames(Dmat)==TRUE)],
#                    rownames(newdata3)[1:100][which(rownames(newdata3)[1:100]%in%rownames(Dmat)==TRUE)])
pm6<-pheatmap(Av1$RNA[MegaList1again,c("CTB","STB","EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1again_",".pdf",sep="") ,width=6,height=30)


aa <- rownames(newdata1a)[1:500]
bb <- rownames(newdata2a)[1:500]
cc <- rownames(newdata3a)[1:500]

MegaList1again <- c(aa[which(rownames(newdata1a)[1:500]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)],
                    bb[which(rownames(newdata2a)[1:500]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)],
                    cc[which(rownames(newdata3a)[1:500]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)])



#MegaList1again <- c(rownames(newdata1)[1:100][which(rownames(newdata1)[1:100]%in%rownames(Dmat)==TRUE)],
#                    rownames(newdata2)[1:100][which(rownames(newdata2)[1:100]%in%rownames(Dmat)==TRUE)],
#                    rownames(newdata3)[1:100][which(rownames(newdata3)[1:100]%in%rownames(Dmat)==TRUE)])
pm6<-pheatmap(Av1$RNA[MegaList1again,c("CTB","STB","EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1a_",".pdf",sep="") ,width=6,height=30)



pm6<-pheatmap(Av1$RNA[MegaList1,c("CTB","STB","EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1_",".pdf",sep="") ,width=6,height=30)
pm7<-pheatmap(Av1$RNA[MegaList2,c("CTB","STB","EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb2_",".pdf",sep="") ,width=6,height=30)
pm8<-pheatmap(Av1$RNA[MegaList3,c("CTB","STB","EVT")], breaks = my_breaks,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb3_",".pdf",sep="") ,width=6,height=30)
pm9<-pheatmap(Av1$RNA[MegaList4,c("CTB","STB","EVT")], breaks = my_breaks,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb4_",".pdf",sep="") ,width=6,height=30)


pm10<-pheatmap(Av2$RNA[MegaList1,c("CTB","STB","EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb5_",".pdf",sep="") ,width=6,height=30)
pm11<-pheatmap(Av2$RNA[MegaList2,c("CTB","STB","EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb6_",".pdf",sep="") ,width=6,height=30)
pm12<-pheatmap(Av2$RNA[MegaList3,c("CTB","STB","EVT")], breaks = my_breaks,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb7_",".pdf",sep="") ,width=6,height=30)
pm13<-pheatmap(Av2$RNA[MegaList4,c("CTB","STB","EVT")], breaks = my_breaks,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb8_",".pdf",sep="") ,width=6,height=30)

pm13<-pheatmap(Av2$RNA[MegaList5,c("CTB","STB","EVT")], breaks = my_breaks,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb9_ml",".pdf",sep="") ,width=6,height=30)
pm13<-pheatmap(Av2$RNA[MegaList6,c("CTB","STB","EVT")], breaks = my_breaks,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb10_ml",".pdf",sep="") ,width=6,height=30)



order1 <- c("HINT1",
            "GSTP1",
            "S100A10",
            "COMMD6",
            "MIF",
            "FAU",
            "LGALS1",
            "CHCHD2",
            "TXN",
            "HSPE1",
            "BTF3",
            "SNRPG",
            "LSM7",
            "TPI1",
            "NDUFA1",
            "NEDD8",
            "GUK1",
            "RAN",
            "RBX1",
            "UBA52",
            "TOMM7",
            "HMGB1",
            "APRT",
            "ATP6V1F",
            "SNRPE",
            "PRELID1",
            "MINOS1",
            "POLR2I",
            "RGS10",
            "SLC25A5",
            "SOD1",
            "MGST3",
            "GAPDH",
            "PCBD1",
            "SNRPD1",
            "NDUFB1",
            "NDUFC1",
            "SNRPD2",
            "SNRPF",
            "RHOC",
            "NSA2",
            "RANBP1",
            "PEBP1",
            "TIMM13",
            "PGK1",
            "PSMB3",
            "EEF1B2",
            "HMGA1",
            "YBX1",
            "NDUFB7",
            "SIVA1",
            "CKS2",
            "LSM3",
            "SLIRP",
            "TAGLN2",
            "GCSH",
            "EMC6",
            "C1QBP",
            "PPP1R14B",
            "NHP2",
            "SNRPB",
            "DYNLL1",
            "ENY2",
            "SVIP",
            "WDR83OS",
            "PFDN2",
            "CD81",
            "PPA1",
            "EEF1D",
            "LAMTOR5",
            "ACP1",
            "GTF2A2",
            "PRDX2",
            "TOMM5",
            "PDCD5",
            "FAM162A",
            "FABP5",
            "NDUFV3",
            "NDUFB8",
            "CKS1B",
            "HEBP2",
            "SPCS1",
            "STMN1",
            "LSM2",
            "PSMB1",
            "CYSTM1",
            "RWDD1",
            "PHB2",
            "ERH",
            "HMGN3",
            "NPM1",
            "LSM4",
            "NDUFS7",
            "UBE2S",
            "MDH2",
            "TMEM14C",
            "TXNL4A",
            "TMEM14B",
            "TIMM10",
            "LAMTOR2",
            "MAGI1",
            "CPM",
            "ZNF292",
            "COBLL1",
            "MAN1A2",
            "ADAMTS20",
            "CBLB",
            "LYN",
            "CGA",
            "MTUS1",
            "NHSL1",
            "MAML3",
            "ABTB2",
            "SLC43A2",
            "NEBL",
            "DIRC2",
            "AFF1",
            "GADD45G",
            "LEP",
            "PLEKHA6",
            "FAR2",
            "PGF",
            "GSE1",
            "ESRRG",
            "NUCB2",
            "TTC7B",
            "SDC1",
            "SHANK2",
            "TACC2",
            "RASA1",
            "PHLPP1",
            "PDIA3",
            "ATP2B4",
            "DAB2",
            "DOCK9",
            "DYSF",
            "GAB2",
            "CSGALNACT1",
            "RBM47",
            "FHDC1",
            "MAST4",
            "FDX1",
            "BASP1",
            "CDH5",
            "FAM135A",
            "HERPUD1",
            "RAB31",
            "KRT18",
            "PACSIN2",
            "ZMYND8",
            "GCNT1",
            "ERC1",
            "FNDC3A",
            "BCAS3",
            "FYN",
            "GSK3B",
            "MAOA",
            "PSG6",
            "ZFHX3",
            "VGLL3",
            "HSPB8",
            "AAK1",
            "TNS3",
            "LGR4",
            "CCR7",
            "GDF15",
            "FAM184A",
            "AEBP2",
            "TGIF1",
            "EPB41L3",
            "CLIC3",
            "TBX3",
            "RNF111",
            "CHSY1",
            "SGPP1",
            "SLC38A1",
            "GAN",
            "MAN1C1",
            "CEBPB",
            "ZFAT",
            "S100P",
            "HEXIM1",
            "SP6",
            "CCSAP",
            "NUFIP2",
            "KRT8",
            "MMP15",
            "DIP2B",
            "PABPC1",
            "SLC7A2",
            "SGK1",
            "PRLR",
            "MED13",
            "KIAA1217",
            "FAM43A",
            "HDAC5",
            "CAP2",
            "INSL4",
            "SIPA1L1",
            "AGAP1",
            "LAIR2",
            "HPGD",
            "TAGLN",
            "INHBA",
            "TIMP3",
            "ITGB1",
            "DIO2",
            "ACTB",
            "FN1",
            "TMSB10",
            "COL4A2",
            "ITGA2",
            "S100A6",
            "S100A4",
            "COL4A1",
            "CYR61",
            "PFN1",
            "PALLD",
            "RASGRF2",
            "HINT1",
            "NOTUM",
            "SH3BGRL3",
            "MYH9",
            "FSTL1",
            "TPM2",
            "SOX4",
            "WFDC2",
            "TUBB6",
            "TAGLN2",
            "LRRFIP1",
            "XYLT1",
            "IGFBP3",
            "RBPMS",
            "CTSL",
            "MYL9",
            "ANXA2",
            "PDLIM5",
            "LGALS3",
            "MYL12A",
            "FAM155A",
            "PXDN",
            "ACTR3",
            "HMGB1",
            "ARPC2",
            "TCF7L2",
            "CAPN2",
            "SEMA3C",
            "CST3",
            "IL1B",
            "CD151",
            "YBX1",
            "FLNB",
            "PON2",
            "ITGA5",
            "HSPG2",
            "TRIO",
            "MSN",
            "EHBP1",
            "CAP1",
            "TNFRSF12A",
            "ANXA3",
            "TLN1",
            "SPATS2L",
            "OCIAD2",
            "PICALM",
            "FER",
            "PDGFC",
            "SRI",
            "TUBA1C",
            "VCL",
            "LIMS1",
            "NRIP1",
            "LAMA4",
            "LAMC1",
            "WNT7A",
            "CTSB",
            "ITM2B",
            "VEGFA",
            "COL1A2",
            "P4HA2",
            "ST3GAL6",
            "HTRA1",
            "EFNB2",
            "MGLL",
            "MTDH",
            "MFAP5",
            "FAM3C",
            "PLK2",
            "GNAQ",
            "TSPAN14",
            "ID1",
            "SIPA1L3",
            "TEAD1",
            "PLXNB2",
            "GALNT1",
            "YWHAB",
            "DOCK5",
            "DDX5",
            "ANTXR2",
            "IER3")



EndoRef <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Placenta/Placentaharm.rds")
Idents(EndoRef) <- EndoRef$ID2
EndoRef <- subset(EndoRef,idents=c("VCT","EVT","SCT"),invert=FALSE)

Av3 <- AverageExpression(EndoRef)

#meh <- Av3$RNA[order1,c("VCT","SCT","EVT")]


MegaList1 <- c(rownames(newdata1)[1:100],rownames(newdata2)[1:100],rownames(newdata3)[1:100] )
MegaList2 <- c(rownames(newdata4)[1:100],rownames(newdata5)[1:100],rownames(newdata6)[1:100] )
MegaList3 <- c(rownames(newdata7)[1:100],rownames(newdata7)[1:100],rownames(newdata9)[1:100] )
MegaList4 <- c(rownames(newdata10)[1:100],rownames(newdata11)[1:100],rownames(newdata12)[1:100] )


MegaList1a <- c(rownames(newdata1a)[1:500],rownames(newdata2a)[1:500],rownames(newdata3a)[1:500] )


MegaList1again <- c(rownames(newdata1)[1:100][which(rownames(newdata1)[1:100]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)],
rownames(newdata2)[1:100][which(rownames(newdata3)[1:100]%in%rownames(D)==TRUE)],
rownames(newdata3)[1:100][which(rownames(newdata3)[1:100]%in%rownames(D)==TRUE)])

#MegaList1p <- c(rownames(newdata1)[1:100],rownames(newdata2)[1:100],rownames(newdata3)[1:100] )


aa <- rownames(newdata1)[1:100]
bb <- rownames(newdata2)[1:100]
cc <- rownames(newdata3)[1:100]

MegaList1again <- c(aa[which(rownames(newdata1)[1:100]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)],
                    bb[which(rownames(newdata2)[1:100]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)],
                    cc[which(rownames(newdata3)[1:100]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)])


pm6<-pheatmap(Av3$RNA[MegaList1again,c("VCT","SCT","EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1_inRVT",".pdf",sep="") ,width=6,height=30)




aa <- rownames(newdata1a)[1:500]
bb <- rownames(newdata2a)[1:500]
cc <- rownames(newdata3a)[1:500]

MegaList1again <- c(aa[which(rownames(newdata1a)[1:100]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)],
                    bb[which(rownames(newdata2a)[1:100]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)],
                    cc[which(rownames(newdata3a)[1:100]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)])

pm6<-pheatmap(Av3$RNA[MegaList1again,c("VCT","SCT","EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1a_inRVT",".pdf",sep="") ,width=6,height=30)



#pm6<-pheatmap(Av3$RNA[MegaList1again,c("VCT","SCT","EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1_inRVT",".pdf",sep="") ,width=6,height=30)





NotOurs1 <- subset(mammal.combined,idents=c("Cyno (10X)_STB","Cyno (10X)_EVT","Cyno (10X)_CTB"))
NotOurs2 <- subset(mammal.combined,idents=c("Cyno (SS2)_STB","Cyno (SS2)_CTB","Cyno (SS2)_EVT"))
DefaultAssay(NotOurs1) <- "RNA"
DefaultAssay(NotOurs2) <- "RNA"

Av4 <- AverageExpression(NotOurs1)
Av5 <- AverageExpression(NotOurs2)


MegaList1again <- c(rownames(newdata1)[1:100][which(rownames(newdata1)[1:100]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)],
                    rownames(newdata2)[1:100][which(rownames(newdata3)[1:100]%in%rownames(D)==TRUE)],
                    rownames(newdata3)[1:100][which(rownames(newdata3)[1:100]%in%rownames(D)==TRUE)])


pm7<-pheatmap(Av4$RNA[MegaList1again,c("Cyno (10X)_CTB","Cyno (10X)_STB","Cyno (10X)_EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1_CynoRef1",".pdf",sep="") ,width=6,height=30)
pm8<-pheatmap(Av5$RNA[MegaList1again,c("Cyno (SS2)_CTB","Cyno (SS2)_STB","Cyno (SS2)_EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1_CynoRef2",".pdf",sep="") ,width=6,height=30)




aa <- rownames(newdata1a)[1:500]
bb <- rownames(newdata2a)[1:500]
cc <- rownames(newdata3a)[1:500]

MegaList1again <- c(aa[which(rownames(newdata1a)[1:100]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)],
                    bb[which(rownames(newdata2a)[1:100]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)],
                    cc[which(rownames(newdata3a)[1:100]%in%intersect(rownames(D),rownames(Dmat) )==TRUE)])


pm7<-pheatmap(Av4$RNA[MegaList1again,c("Cyno (10X)_CTB","Cyno (10X)_STB","Cyno (10X)_EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1_CynoRef1a",".pdf",sep="") ,width=6,height=30)
pm8<-pheatmap(Av5$RNA[MegaList1again,c("Cyno (SS2)_CTB","Cyno (SS2)_STB","Cyno (SS2)_EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1_CynoRef2a",".pdf",sep="") ,width=6,height=30)


pm7<-pheatmap(Av4$RNA[order1,c("Cyno (10X)_CTB","Cyno (10X)_STB","Cyno (10X)_EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1_CynoRef1",".pdf",sep="") ,width=6,height=30)
pm8<-pheatmap(Av5$RNA[order1,c("Cyno (SS2)_CTB","Cyno (SS2)_STB","Cyno (SS2)_EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1_CynoRef2",".pdf",sep="") ,width=6,height=30)

DefaultAssay(Dsubset2) <- "RNA"
Av6 <- AverageExpression(Dsubset2)

pm8<-pheatmap(Av6$RNA[order1,c("CTB","STB","EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1_oldAno",".pdf",sep="") ,width=6,height=30)

pm7<-pheatmap(Av4$RNA[MegaList1a,c("Cyno (10X)_CTB","Cyno (10X)_STB","Cyno (10X)_EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1_CynoRef1a",".pdf",sep="") ,width=6,height=30)
pm8<-pheatmap(Av5$RNA[MegaList1a,c("Cyno (SS2)_CTB","Cyno (SS2)_STB","Cyno (SS2)_EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1_CynoRef2a",".pdf",sep="") ,width=6,height=30)

pm6<-pheatmap(Av3$RNA[MegaList1a,c("VCT","SCT","EVT")], breaks = my_breaks, color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", filename = paste(saveext,"/DimRed/FullListTb1a_inRVT",".pdf",sep="") ,width=6,height=30)

#Module score plots

p1<-FeaturePlot(OnlyOurs, features = "CTBSS1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_CTBSS_ModuleScore.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(OnlyOurs, features = "STBSS1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_STBSS_ModuleScore.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)
p1<-FeaturePlot(OnlyOurs, features = "EVTSS1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_EVTSS_ModuleScore.pdf",sep=""),width = 10, height = 8,p1,limitsize = FALSE)

write.table(OnlyOurs@meta.data,file="TBMetaData.csv",sep=",",quote = FALSE)

oL2 <- c("DNMT3A",
         "TP63",
         "ELF5",
         "CDH1",
         "VIM",
         "CTNNAL1",
         "LAMC1",
         "TACSTD2",
         "DPPA3",
         "CCKBR",
         "ENPEP",
         "CGB8",
         "CGB5",
         "SDC1",
         "CGA",
         "OVOL1",
         "KRT19",
         "MKI67",
         "OTX2",
         "FN1",
         "FSTL3",
         "KRT7",
         "JAM3",
         "FSTL1",
         "NOTUM",
         "ASCL2",
         "HLA-G")



Idents(mammal.combined) <- paste(mammal.combined$Dataset,Idents(mammal.combined),sep="_")







p1<-FeaturePlot(mammal.combined, features = "HLA-G", pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_cyno_align_PCAHLA-G.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


p1<-FeaturePlot(mammal.combined, features = "ASCL2", pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_cyno_align_PCA_ASCL2.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)



p1<-FeaturePlot(mammal.combined, features = "FSTL1", pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_cyno_align_PCA_FSTL1.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


p1<-FeaturePlot(mammal.combined, features = "CGA", pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_cyno_align_PCA_CGA.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


p1<-FeaturePlot(mammal.combined, features = "HOPX", pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_cyno_align_PCA_HOPX.pdf",sep=""),width = 20, height = 8,p1,limitsize = FALSE)


mammal.combined$ID0 <- Idents(mammal.combined)


DefaultAssay(mammal.combined) <= "integrated"
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindClusters(mammal.combined, graph.name = "integrated_snn", resolution = 0.1)


p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_cyno_align_Cluster.pdf",sep=""),width = 20, height = 8,p)

DefaultAssay(mammal.combined) <- "RNA"
DE <- FindMarkers(mammal.combined,ident.1 = "1",ident.2 = "3",test.use="MAST")

write.table(as.data.frame(DE),file = "~/Desktop/DE_CTB_EVT.csv")


list <- c("CGA",
"HOPX",
"NPPB",
"LOC102137201",
"LOC102137969",
"KRT23",
"DAB2",
"COL19A1",
"ADAMTS20",
"LOC102138598",
"TFRC",
"COBLL1",
"LINC02506",
"NPPA",
"SH3KBP1",
"ZNF292",
"PPP1R14A",
"MALAT1",
"LOC102116649",
"RBM24",
"MAML3",
"MAN1A2",
"CD164",
"MLLT1",
"SLC24A3",
"HSD3B1",
"LOC101925912",
"KIF21A",
"TBX3",
"RPS6KA5",
"MYCNUT",
"ASB2",
"CYP19A1",
"SBSN",
"CBLB",
"SLC2A1",
"RASA1",
"MTUS1",
"MT-ND3",
"AFF1",
"SGK1",
"LOC101865755",
"SLC40A1",
"NEAT1",
"MAGI1",
"LOC102120640",
"PGF",
"MT-ATP6",
"RBM47",
"TRIM58",
"LINC02484",
"VSIG1",
"ARHGAP24",
"ESRRG",
"ERVW-1",
"LEP",
"SAT1",
"SIAH1",
"LOC101867199",
"LINC01949")


p1<-FeaturePlot(mammal.combined, features = list, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_cyno_align_PCA_NPPA.pdf",sep=""),width = 20, height = 200,p1,limitsize = FALSE)


listb <-c("MMP12",
"ISM2",
"SERPINB6",
"PLAC8",
"NOTUM",
"SRI",
"C15orf48",
"SERPINE1",
"CXCL5",
"TPM1",
"S100A6",
"ESAM",
"TAGLN2",
"WFDC2",
"PPP1R14B",
"PALLD",
"STMN1",
"ARPC1B",
"CTSL",
"TAGLN",
"PFN1",
"TINAGL1",
"TUBB2B",
"ARPC2",
"TLN1",
"ARPC5",
"HSD17B2",
"ITM2B",
"LIMS3",
"TMSB10",
"TPM2",
"ECM1")


p1<-FeaturePlot(mammal.combined, features = listb, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_cyno_align_PCA_allevt.pdf",sep=""),width = 20, height = 120,p1,limitsize = FALSE)





listc<-c("PLA2G7",
"LAPTM4B",
"NCAM1",
"SEPTIN11",
"JAM2",
"MMP2",
"RPL22L1",
"PLK2",
"LIMA1",
"CCN2",
"MIF",
"ACTG1",
"EFHD2",
"TUBA1B",
"PLOD2",
"MGST1",
"SRM",
"TUBB6",
"PECAM1",
"BNIP3",
"H2AZ1",
"ACTB",
"MFAP2",
"FERMT2",
"SERPINE2",
"TBCA",
"CREM",
"H2AZ2",
"FST")



p1<-FeaturePlot(mammal.combined, features = listc, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_cyno_align_PCA_allevt2.pdf",sep=""),width = 20, height = 120,p1,limitsize = FALSE)



DE <- FindMarkers(mammal.combined,ident.1 = "0",ident.2 = c("3","1"),test.use="MAST")

write.table(as.data.frame(DE),file = "~/Desktop/DE_CTB_STBEVT.csv")




listd <- c("NCL",
"PAGE4",
"RANBP1",
"HMGB2",
"ANP32B",
"TOP2A",
"PTPMT1",
"PTTG1",
"HSP90AB1",
"PHGDH",
"CENPF",
"YDJC",
"SERBP1",
"SNRPD1",
"DNMT1",
"TFDP2",
"H1-5",
"APRT",
"SLC25A5",
"NPM3",
"MRPS26",
"LSM5",
"NASP",
"DUT",
"LY6E",
"MRPS34",
"S100A10",
"MRPL23",
"LSM6",
"PPA1",
"CDC20")


p1<-FeaturePlot(mammal.combined, features = listd, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/Tb_human_cyno_align_PCA_allctb.pdf",sep=""),width = 20, height = 120,p1,limitsize = FALSE)

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
#human_dataA2 <- subset(human_dataA2,idents=c("YS mesoderm","YS Endoderm"))
human_dataA2 <- subset(human_dataA2,idents=c("YS mesoderm","YS Endoderm","DE (Caudal)","DE (Rostal)"))

human_dataA2$Dataset <- "Human CS7"


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
marmoset_dataInVivo2 <- subset(marmoset_dataInVivo2,idents=c("ExMes_CS5","ExMes_CS6","ExMes_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7"))

Dsubset2B$Dataset <- "Ours"
marmoset_dataInVivo2$Dataset <- "Marmoset"
mammal.anchors <- FindIntegrationAnchors(object.list = list(Dsubset2B,cynomolgous_dataB,marmoset_dataInVivo2), reduction = "cca",normalization.method = c("LogNormalize"), dims = 1:20, anchor.features = 3000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYS_human_cyno_align_PCA.pdf",sep=""),width = 30, height = 8,p)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/SYS_human_cyno_align_PCA_1_3.pdf",sep=""),width = 30, height = 8,p)
#saveRDS(mammal.combined,file=paste(saveext,"/Batch1_marmosetintegration.rds",sep=""))

DE <- FindMarkers(mammal.combined,ident.1 = c("SYS_CS5","SYS_CS6","SYS_CS7"),ident.2 = c("ExMes_CS5","ExMes_CS6","ExMes_CS7"),test.use="MAST")

write.table(as.data.frame(DE),file = "~/Desktop/DE_SYS_ExMes.csv")



listaa <- c("APOA1",
"TTR",
"SPINK1",
"FABP1",
"LRP2",
"GDF3",
"FAM155A",
"EPCAM",
"SYNE2",
"DPPA3",
"AK4",
"GPC3",
"SERPINE2",
"S100A16",
"CDH2",
"SLC44A5",
"CCKBR",
"ADAMTS19",
"CYTL1",
"MSMB",
"APOB",
"LGR5",
"NRG1",
"ENPP1",
"S100A6",
"WFDC2",
"KYNU",
"MTTP",
"ANXA4",
"WWC2",
"FAT3",
"AUTS2",
"VTN",
"EPSTI1",
"HPGD",
"KLHL14",
"FGF12",
"ADM",
"CTSV",
"OTX2")

listbb <- c("HGF",
"COL1A1",
"POSTN",
"COL3A1",
"COL6A3",
"IGFBP3",
"RBMS3",
"COL1A2",
"HAPLN1",
"COL6A1",
"ZFPM2",
"HMGA2",
"LAMA2",
"TGFBI",
"LSAMP",
"ARHGDIB",
"ADAMTS20",
"HES1",
"COL6A2",
"SST",
"VIM",
"CNN3",
"COLEC10",
"CHSY3",
"UST",
"TEK",
"ITGA1",
"PTPRG",
"DIP2C",
"BAMBI",
"PMP22",
"DPYSL3",
"LIN28B",
"TBX20",
"RSPO2",
"ADAMTS6",
"PRICKLE1",
"FGF7",
"FAM46A",
"PDLIM3")

DefaultAssay(mammal.combined) <- "RNA"

p1<-FeaturePlot(mammal.combined, features = listaa, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/SYSEXmes_human_cyno_align_PCA_allsys.pdf",sep=""),width = 20, height = 120,p1,limitsize = FALSE)




p1<-FeaturePlot(mammal.combined, features = listbb, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/EXmesSYS_human_cyno_align_PCA_allsys.pdf",sep=""),width = 20, height = 120,p1,limitsize = FALSE)



Dsubset2B <- subset(Dsubset2B,idents = "Hypoblast", invert = TRUE)
cynomolgous_dataB <- subset(cynomolgous_dataB,idents = c("VE_CS5"), invert = TRUE)
marmoset_dataInVivo2 <- subset(marmoset_dataInVivo2,idents = c("VE_CS5","VE_CS6","VE_CS7"), invert = TRUE)




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
#human_dataA2 <- subset(human_dataA2,idents=c("YS mesoderm","YS Endoderm"))
human_dataA2 <- subset(human_dataA2,idents=c("YS mesoderm","YS Endoderm","DE (Caudal)","DE (Rostal)"))

human_dataA2$Dataset <- "Human CS7"


mammal.anchors <- FindIntegrationAnchors(object.list = list(Dsubset2B,cynomolgous_dataB,marmoset_dataInVivo2), reduction = "cca",normalization.method = c("LogNormalize"), dims = 1:20, anchor.features = 3000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYS_noVE_human_cyno_align_PCA.pdf",sep=""),width = 30, height = 8,p)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/SYS_noVE_human_cyno_align_PCA_1_3.pdf",sep=""),width = 30, height = 8,p)
#saveRDS(mammal.combined,file=paste(saveext,"/Batch1_marmosetintegration.rds",sep=""))




DefaultAssay(mammal.combined) <- "RNA"

p1<-FeaturePlot(mammal.combined, features = listaa, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/SYSEXmes_noVE_human_cyno_align_PCA_allsys.pdf",sep=""),width = 20, height = 120,p1,limitsize = FALSE)

p1<-FeaturePlot(mammal.combined, features = listbb, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, cols = c("lightgrey","red")) 
ggsave(filename=paste(saveext,"/DimRed/EXmesSYS_noVE_human_cyno_align_PCA_allsys.pdf",sep=""),width = 20, height = 120,p1,limitsize = FALSE)

Idents(mammal.combined) <- paste(mammal.combined$Dataset,Idents(mammal.combined),sep="_")

DE1 <- FindMarkers(mammal.combined,ident.1 = c("Marmoset_SYS_CS5","Marmoset_SYS_CS6","Marmoset_SYS_CS7"),ident.2 = c("Marmoset_ExMes_CS5","Marmoset_ExMes_CS6","Marmoset_ExMes_CS7","Marmoset_VE_CS6","Marmoset_VE_CS5","Marmoset_VE_CS7"),test.use="MAST", only.pos = TRUE)
DE2 <- FindMarkers(mammal.combined,ident.1 = c("Marmoset_ExMes_CS5","Marmoset_ExMes_CS6","Marmoset_ExMes_CS7"),ident.2 = c("Marmoset_VE_CS7","Marmoset_SYS_CS5","Marmoset_SYS_CS6","Marmoset_SYS_CS7","Marmoset_VE_CS6","Marmoset_VE_CS5","Marmoset_VE_CS7"),test.use="MAST", only.pos = TRUE)
DE3 <- FindMarkers(mammal.combined,ident.1 = c("Marmoset_VE_CS5","Marmoset_VE_CS6","Marmoset_VE_CS7"),ident.2 = c("Marmoset_SYS_CS5","Marmoset_SYS_CS6","Marmoset_SYS_CS7","Marmoset_ExMes_CS5","Marmoset_ExMes_CS6","Marmoset_ExMes_CS7"),test.use="MAST")

DE4 <- FindMarkers(mammal.combined,ident.1 = c("Cyno_SYS_CS6","Cyno_SYS_CS7"),ident.2 = c("Cyno_ExMes_CS5","Cyno_ExMes_CS6","Cyno_ExMes_CS7","Cyno_VE_CS5"),test.use="MAST", only.pos = TRUE)
DE5 <- FindMarkers(mammal.combined,ident.1 = c("Cyno_ExMes_CS5","Cyno_ExMes_CS6","Cyno_ExMes_CS7"),ident.2 = c("Cyno_VE_CS5","Cyno_SYS_CS6","Cyno_SYS_CS7","Cyno_VE_CS5"),test.use="MAST", only.pos = TRUE)
DE6 <- FindMarkers(mammal.combined,ident.1 = c("Cyno_VE_CS5"),ident.2 = c("Cyno_SYS_CS6","Cyno_SYS_CS7","Cyno_ExMes_CS5","Cyno_ExMes_CS6","Cyno_ExMes_CS7"),test.use="MAST")


DE7 <- FindMarkers(human_dataA2,ident.1 = c("YS Endoderm"),ident.2 = c("YS mesoderm","DE (Caudal)","DE (Rostal)"),test.use="MAST")
DE8 <- FindMarkers(human_dataA2,ident.1 = c("YS mesoderm"),ident.2 = c("YS Endoderm","DE (Caudal)","DE (Rostal)"),test.use="MAST")
DE9 <- FindMarkers(human_dataA2,ident.1 = c("DE (Caudal)","DE (Rostal)"),ident.2 = c("YS mesoderm","YS Endoderm"),test.use="MAST")



newdata1 <- DE1[order(DE1$avg_log2FC,decreasing = TRUE),]
newdata2 <- DE2[order(DE2$avg_log2FC,decreasing = TRUE),]
newdata3 <- DE3[order(DE3$avg_log2FC,decreasing = TRUE),]
newdata4 <- DE4[order(DE4$avg_log2FC,decreasing = TRUE),]
newdata5 <- DE5[order(DE5$avg_log2FC,decreasing = TRUE),]
newdata6 <- DE6[order(DE6$avg_log2FC,decreasing = TRUE),]
newdata7 <- DE7[order(DE7$avg_log2FC,decreasing = TRUE),]
newdata8 <- DE8[order(DE8$avg_log2FC,decreasing = TRUE),]
newdata9 <- DE9[order(DE9$avg_log2FC,decreasing = TRUE),]

DefaultAssay(mammal.combined) <- "RNA"
mammal.combined <- AddModuleScore(mammal.combined, features =  list( rownames(newdata1)[1:100] ) , name = "SYSMark")
mammal.combined <- AddModuleScore(mammal.combined, features =  list( rownames(newdata2)[1:100] ) , name = "ExMesMarm")
mammal.combined <- AddModuleScore(mammal.combined, features =  list( rownames(newdata3)[1:100] ) , name = "VEMarm")
mammal.combined <- AddModuleScore(mammal.combined, features =  list( rownames(newdata4)[1:100] ) , name = "SYSCyno")
mammal.combined <- AddModuleScore(mammal.combined, features =  list( rownames(newdata5)[1:100] ) , name = "ExMesCyno")
mammal.combined <- AddModuleScore(mammal.combined, features =  list( rownames(newdata6)[1:100] ) , name = "VECyno")
mammal.combined <- AddModuleScore(mammal.combined, features =  list( rownames(newdata7)[1:100] ) , name = "SYSHuman")
mammal.combined <- AddModuleScore(mammal.combined, features =  list( rownames(newdata8)[1:100] ) , name = "ExMesHuman")
mammal.combined <- AddModuleScore(mammal.combined, features =  list( rownames(newdata9)[1:100] ) , name = "VEHuman")

onlyOurs <- subset(mammal.combined,idents=c("Ours_SYS","Ours_ExMes","Ours_Hypoblast","Ours_ExMes/SYS"))

p1<-FeaturePlot(onlyOurs, features = "SYSMark1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_MarmOurs_SYSMarmList.pdf",sep=""),width = 10, height = 8,p1)
p1<-FeaturePlot(onlyOurs, features = "ExMesMarm1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_MarmOurs_ExMesMarmList.pdf",sep=""),width = 10, height = 8,p1)
p1<-FeaturePlot(onlyOurs, features = "VEMarm1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_MarmOurs_VEMarmList.pdf",sep=""),width = 10, height = 8,p1)

p1<-FeaturePlot(onlyOurs, features = "SYSCyno1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_CynoOurs_SYSCynoList.pdf",sep=""),width = 10, height = 8,p1)
p1<-FeaturePlot(onlyOurs, features = "ExMesCyno1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_CynoOurs_ExMesCynoList.pdf",sep=""),width = 10, height = 8,p1)
p1<-FeaturePlot(onlyOurs, features = "VECyno1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_CynoOurs_VECynoList.pdf",sep=""),width = 10, height = 8,p1)

p1<-FeaturePlot(onlyOurs, features = "SYSHuman1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_PrimateOurs_SYSPHumanList.pdf",sep=""),width = 10, height = 8,p1)
p1<-FeaturePlot(onlyOurs, features = "ExMesHuman1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_PrimateOurs_ExMesHumanList.pdf",sep=""),width = 10, height = 8,p1)
p1<-FeaturePlot(onlyOurs, features = "VEHuman1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_PrimateOurs_VEHumanList.pdf",sep=""),width = 10, height = 8,p1)


newID1 <- as.character(Idents(onlyOurs))
newID1[which(onlyOurs$SYSMark1>onlyOurs$ExMesMarm1 & onlyOurs$SYSMark1>onlyOurs$VEMarm1)] <- "SYS"
newID1[which(onlyOurs$ExMesMarm1>onlyOurs$SYSMark1 & onlyOurs$ExMesMarm1>onlyOurs$VEMarm1)] <- "ExMes"
newID1[which(onlyOurs$VEMarm1>onlyOurs$ExMesMarm1 & onlyOurs$VEMarm1>onlyOurs$SYSMark1)] <- "VE"

newID2 <- as.character(Idents(onlyOurs))
newID2[which(onlyOurs$SYSHuman1>onlyOurs$ExMesHuman1 & onlyOurs$SYSHuman1>onlyOurs$VEHuman1)] <- "SYS"
newID2[which(onlyOurs$ExMesHuman1>onlyOurs$SYSHuman1 & onlyOurs$ExMesHuman1>onlyOurs$VEHuman1)] <- "ExMes"
newID2[which(onlyOurs$VEHuman1>onlyOurs$ExMesHuman1 & onlyOurs$VEHuman1>onlyOurs$SYSHuman1)] <- "VE"

newID3 <- as.character(Idents(onlyOurs))
newID3[which(onlyOurs$SYSCyno1>onlyOurs$ExMesCyno1 & onlyOurs$SYSCyno1>onlyOurs$VECyno1)] <- "SYS"
newID3[which(onlyOurs$ExMesCyno1>onlyOurs$SYSCyno1 & onlyOurs$ExMesCyno1>onlyOurs$VECyno1)] <- "ExMes"
newID3[which(onlyOurs$VECyno1>onlyOurs$ExMesCyno1 & onlyOurs$VECyno1>onlyOurs$SYSCyno1)] <- "VE"

Idents(onlyOurs) <- newID1
p1<-DimPlot(onlyOurs, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_newAnotation_Marm.pdf",sep=""),width = 10, height = 8,p1)
Idents(onlyOurs) <- newID2
p1<-DimPlot(onlyOurs, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_newAnotation_Human.pdf",sep=""),width = 10, height = 8,p1)
Idents(onlyOurs) <- newID3
p1<-DimPlot(onlyOurs, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_newAnotation_Cyno.pdf",sep=""),width = 10, height = 8,p1)



EmbModel <- read_excel("/Users/christopherpenfold/Downloads/41586_2023_6604_MOESM3_ESM.xlsx", sheet = "Human_SEM_markers")

JHL1 <- EmbModel$gene[which(EmbModel$`cluster annotation`=="Yolk-Sac / Hypoblast")]
JHL2 <- EmbModel$gene[which(EmbModel$`cluster annotation`=="ExEM")]
JHL3 <- EmbModel$gene[which(EmbModel$`cluster annotation`=="Secondary Yolk-Sac")]

JHL4 <- EmbModel$gene[which(EmbModel$`cluster annotation`=="Epiblast")]
JHL5 <- EmbModel$gene[which(EmbModel$`cluster annotation`=="Posterior Epiblast")]
JHL6 <- EmbModel$gene[which(EmbModel$`cluster annotation`=="Epiblast committed")]
JHL7 <- EmbModel$gene[which(EmbModel$`cluster annotation`=="Amnion")]


onlyOurs <- AddModuleScore(onlyOurs, features =  list( JHL3[1:150] ) , name = "SYSJH")
onlyOurs <- AddModuleScore(onlyOurs, features =  list( JHL2[1:150] ) , name = "ExMesJH")
onlyOurs <- AddModuleScore(onlyOurs, features =  list( JHL1[1:150] ) , name = "VEJH")
onlyOurs <- AddModuleScore(onlyOurs, features =  list( JHL4[1:150] ) , name = "EpJH")
onlyOurs <- AddModuleScore(onlyOurs, features =  list( JHL5[1:150] ) , name = "PEpJH")
onlyOurs <- AddModuleScore(onlyOurs, features =  list( JHL6[1:150] ) , name = "cEpJH")
onlyOurs <- AddModuleScore(onlyOurs, features =  list( JHL7[1:150] ) , name = "AmJH")

p1<-FeaturePlot(onlyOurs, features = "SYSJH1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_PrimateOurs_SYSJHList.pdf",sep=""),width = 10, height = 8,p1)
p1<-FeaturePlot(onlyOurs, features = "ExMesJH1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_PrimateOurs_ExMesJHList.pdf",sep=""),width = 10, height = 8,p1)
p1<-FeaturePlot(onlyOurs, features = "VEJH1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_PrimateOurs_VEJHList.pdf",sep=""),width = 10, height = 8,p1)

p1<-FeaturePlot(onlyOurs, features = "EpJH1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_PrimateOurs_EpiJHList.pdf",sep=""),width = 10, height = 8,p1)
p1<-FeaturePlot(onlyOurs, features = "PEpJH1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_PrimateOurs_pEpiJHList.pdf",sep=""),width = 10, height = 8,p1)
p1<-FeaturePlot(onlyOurs, features = "cEpJH1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_PrimateOurs_cEpiJHList.pdf",sep=""),width = 10, height = 8,p1)
p1<-FeaturePlot(onlyOurs, features = "AmJH1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_PrimateOurs_AmJHList.pdf",sep=""),width = 10, height = 8,p1)


newID4 <- as.character(Idents(onlyOurs))
newID4[which(onlyOurs$SYSJH1>onlyOurs$ExMesJH1 & onlyOurs$SYSJH1>onlyOurs$VEJH1)] <- "SYS"
newID4[which(onlyOurs$ExMesJH1>onlyOurs$SYSJH1 & onlyOurs$ExMesJH1>onlyOurs$VEJH1)] <- "ExMes"
newID4[which(onlyOurs$VEJH1>onlyOurs$SYSJH1 & onlyOurs$VEJH1>onlyOurs$ExMesJH1)] <- "VE"

Idents(onlyOurs) <- newID4
p1<-DimPlot(onlyOurs, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYSVEExMes_newAnotation_JH.pdf",sep=""),width = 10, height = 8,p1)



Dsubset2B <- subset(Dsubset2B,idents=c("SYS","Hypoblast"))
marmoset_dataInVivo2 <- subset(marmoset_dataInVivo2,idents=c("VE_CS5","SYS_CS5","SYS_CS6","VE_CS6","SYS_CS7","VE_CS7"))
Dsubset2B$Dataset <- "Ours"
marmoset_dataInVivo2$Dataset <- "Marmoset"
mammal.anchors <- FindIntegrationAnchors(object.list = list(Dsubset2B,marmoset_dataInVivo2), reduction = "cca",normalization.method = c("LogNormalize"), dims = 1:20, anchor.features = 3000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/SYS_marm_align_PCA.pdf",sep=""),width = 20, height = 8,p)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE, dims = c(1,3)) 
ggsave(filename=paste(saveext,"/DimRed/SYS_marm_align_PCA_1_3.pdf",sep=""),width = 20, height = 8,p)
#saveRDS(mammal.combined,file=paste(saveext,"/Batch1_marmosetintegration.rds",sep=""))

