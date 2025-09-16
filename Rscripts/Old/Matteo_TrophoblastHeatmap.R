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


saveRDS(Dsubset1,"/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/CorrectedFinalAno3.rds")

newID <- paste(Idents(Dsubset2),Dsubset2$ID3,sep="_")
Idents(Dsubset2) <- newID
DefaultAssay(Dsubset2) <- "RNA"
AvExp <- AverageExpression(Dsubset2)
nC <- cor( log2(AvExp$RNA[,c("CTB_C1","STB_C1","STB1_C1")] +1), method = "pearson")
my_breaks <- seq(0.8, 1, length.out = 50)
pm6<-pheatmap(nC,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=TRUE, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/OurTBCoor_C1",".pdf",sep="") ,width=12,height=12)

nC <- cor( log2(AvExp$RNA[,c("CTB_C2","STB_C2","EVT_C2")] +1), method = "pearson")
my_breaks <- seq(0.8, 1, length.out = 50)
pm6<-pheatmap(nC,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=TRUE, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/OurTBCoor_C2",".pdf",sep="") ,width=12,height=12)


nC <- cor( log2(AvExp$RNA[,c("CTB_C2","STB_C2","EVT_C2")] +1), method = "pearson")
my_breaks <- seq(0.8, 1, length.out = 50)
pm6<-pheatmap(nC,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=TRUE, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/OurTBCoor_C2",".pdf",sep="") ,width=12,height=12)



nC <- cor( log2(AvExp$RNA[,c("CTB_C4","STB_C4","EVT_C4","STB_C3","EVT_C3")] +1), method = "pearson")
my_breaks <- seq(0.8, 1, length.out = 50)
pm6<-pheatmap(nC,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=TRUE, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/OurTBCoor_C3",".pdf",sep="") ,width=12,height=12)


nC <- cor( log2(AvExp$RNA[,c("CTB_C6","STB_C6","EVT_C6","STB1_C6")] +1), method = "pearson")
my_breaks <- seq(0.8, 1, length.out = 50)
pm6<-pheatmap(nC,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=TRUE, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/OurTBCoor_C6",".pdf",sep="") ,width=12,height=12)



nC <- cor( log2(AvExp$RNA[,c("CTB_C7","STB_C7","EVT_C7")] +1), method = "pearson")
my_breaks <- seq(0.8, 1, length.out = 50)
pm6<-pheatmap(nC,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=TRUE, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/OurTBCoor_C7",".pdf",sep="") ,width=12,height=12)


nC <- cor( log2(AvExp$RNA[,c("CTB_C8","STB_C8","EVT_C8")] +1), method = "pearson")
my_breaks <- seq(0.8, 1, length.out = 50)
pm6<-pheatmap(nC,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=TRUE, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/OurTBCoor_C8",".pdf",sep="") ,width=12,height=12)

oL <- c("DNMT3A",
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
"TFAP2C",
"GATA3",
"CGB8",
"CGB5",
"SDC1",
"TFAP2A",
"CGA",
"OVOL1",
"TEAD3",
"KRT19",
"GATA2",
"MKI67",
"OTX2",
"FN1",
"FSTL3",
"KRT7",
"JAM3",
"FSTL1",
"HLA-C",
"SNAI1",
"SPARC",
"NOTUM",
"ASCL2",
"HLA-G")


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


AvExp <- AverageExpression(Dsubset1)
Expr <- AvExp$RNA[oL2,]
Expr <- t(scale(t(Expr)))
#Expr <- na.omit(Expr)
Expr <- Expr[,c("CTB","STB","EVT")]
Expr <- t(scale(t(Expr)))

Exps <- cbind(Expr,Expr2)

redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm6<-pheatmap(Expr2,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Tb_VT_MinListPB1",".pdf",sep="") ,width=4,height=4)

pm6<-pheatmap(Expr,color =  redblue1(50), fontsize = 4,  border_color = NA,  cluster_rows=FALSE, cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Tb_VT_MinListPB2",".pdf",sep="") ,width=4,height=4)




Data <- GetAssayData(Dsubset2)
Data <- Data[oL2,]
Data <- t(scale(t(Data)))

Datas <- Data[,which(Dsubset2$ID3=="C1")]
annotation_col = data.frame(Stage = factor(droplevels(Idents(Dsubset2)[which(Dsubset2$ID3=="C1")])))
my_breaks <- seq(-2, 2, length.out = 50)
pm6<-pheatmap(Datas,color =  redblue1(50), breaks = my_breaks, fontsize = 4,  border_color = NA,  cluster_rows=FALSE, annotation_col = annotation_col, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/Tb_AllList_C1",".pdf",sep="") ,width=10,height=4)

Datas <- Data[,which(Dsubset2$ID3=="C2")]
annotation_col = data.frame(Stage = factor(droplevels(Idents(Dsubset2)[which(Dsubset2$ID3=="C2")])))
my_breaks <- seq(-2, 2, length.out = 50)
pm6<-pheatmap(Datas,color =  redblue1(50), breaks = my_breaks, fontsize = 4,  border_color = NA,  cluster_rows=FALSE, annotation_col = annotation_col, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/Tb_AllList_C2",".pdf",sep="") ,width=20,height=4)


Datas <- Data[,which(Dsubset2$ID3%in%c("C3","C4"))]
annotation_col = data.frame(Stage = factor(droplevels(Idents(Dsubset2)[which(Dsubset2$ID3%in%c("C3","C4"))])))
my_breaks <- seq(-2, 2, length.out = 50)
pm6<-pheatmap(Datas,color =  redblue1(50), breaks = my_breaks, fontsize = 4,  border_color = NA,  cluster_rows=FALSE, annotation_col = annotation_col, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/Tb_AllList_C3",".pdf",sep="") ,width=10,height=4)


Datas <- Data[,which(Dsubset2$ID3%in%c("C6"))]
annotation_col = data.frame(Stage = factor(droplevels(Idents(Dsubset2)[which(Dsubset2$ID3%in%c("C6"))])))
my_breaks <- seq(-2, 2, length.out = 50)
pm6<-pheatmap(Datas,color =  redblue1(50), breaks = my_breaks, fontsize = 4,  border_color = NA,  cluster_rows=FALSE, annotation_col = annotation_col, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/Tb_AllList_C6",".pdf",sep="") ,width=10,height=4)


Datas <- Data[,which(Dsubset2$ID3%in%c("C7"))]
annotation_col = data.frame(Stage = factor(droplevels(Idents(Dsubset2)[which(Dsubset2$ID3%in%c("C7"))])))
my_breaks <- seq(-2, 2, length.out = 50)
pm6<-pheatmap(Datas,color =  redblue1(50), breaks = my_breaks, fontsize = 4,  border_color = NA,  cluster_rows=FALSE, annotation_col = annotation_col, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/Tb_AllList_C7",".pdf",sep="") ,width=6,height=4)


Datas <- Data[,which(Dsubset2$ID3%in%c("C8"))]
annotation_col = data.frame(Stage = factor(droplevels(Idents(Dsubset2)[which(Dsubset2$ID3%in%c("C8"))])))
my_breaks <- seq(-2, 2, length.out = 50)
pm6<-pheatmap(Datas,color =  redblue1(50), breaks = my_breaks, fontsize = 4,  border_color = NA,  cluster_rows=FALSE, annotation_col = annotation_col, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/Tb_AllList_C8",".pdf",sep="") ,width=8,height=4)


#annotation_col = data.frame(Stage = factor(droplevels( Idents(Dsubset2) )), Batch =  factor(( Idents(Dsubset2)))  )
#redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
#my_breaks <- seq(-2, 2, length.out = 50)
##pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/PGC_markerplot2",".pdf",sep="") ,width=15,height=3)



Datas <- Data[,which(Dsubset2$Genotype2=="NC2_1")]
annotation_col = data.frame(Stage = factor(droplevels(Idents(Dsubset2)[which(Dsubset2$Genotype2=="NC2_1")])))
my_breaks <- seq(-2, 2, length.out = 50)
pm6<-pheatmap(Datas,color =  redblue1(50), breaks = my_breaks, fontsize = 4,  border_color = NA,  cluster_rows=FALSE, annotation_col = annotation_col, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/Tb_AllList_C2_G1",".pdf",sep="") ,width=20,height=4)


Datas <- Data[,which(Dsubset2$Genotype2=="NC2_2")]
annotation_col = data.frame(Stage = factor(droplevels(Idents(Dsubset2)[which(Dsubset2$Genotype2=="NC2_2")])))
my_breaks <- seq(-2, 2, length.out = 50)
pm6<-pheatmap(Datas,color =  redblue1(50), breaks = my_breaks, fontsize = 4,  border_color = NA,  cluster_rows=FALSE, annotation_col = annotation_col, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/Tb_AllList_C2_G2",".pdf",sep="") ,width=20,height=4)



Dsubset3 <- subset(Dsubset2,idents=c("STB","CTB"))
Dsubset4 <- subset(Dsubset2,idents=c("CTB","EVT"))
Dsubset5 <- subset(Dsubset2,idents=c("STB","EVT"))



Data <- GetAssayData(Dsubset3)
Data <- Data[oL2,]
Data <- t(scale(t(Data)))

Datas <- Data
annotation_col = data.frame(Stage = factor(droplevels(Idents(Dsubset3))))
my_breaks <- seq(-2, 2, length.out = 50)
pm6<-pheatmap(Datas,color =  redblue1(50), breaks = my_breaks, fontsize = 4,  border_color = NA,  cluster_rows=FALSE, annotation_col = annotation_col, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/Tb_AllList_STB_CTB",".pdf",sep="") ,width=60,height=4)


Data <- GetAssayData(Dsubset4)
Data <- Data[oL2,]
Data <- t(scale(t(Data)))
Datas <- Data
annotation_col = data.frame(Stage = factor(droplevels(Idents(Dsubset4))))
my_breaks <- seq(-2, 2, length.out = 50)
pm6<-pheatmap(Datas,color =  redblue1(50), breaks = my_breaks, fontsize = 4,  border_color = NA,  cluster_rows=FALSE, annotation_col = annotation_col, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/Tb_AllList_EVT_CTB",".pdf",sep="") ,width=30,height=4)


Data <- GetAssayData(Dsubset5)
Data <- Data[oL2,]
Data <- t(scale(t(Data)))
Datas <- Data
annotation_col = data.frame(Stage = factor(droplevels(Idents(Dsubset5))))
my_breaks <- seq(-2, 2, length.out = 50)
pm6<-pheatmap(Datas,color =  redblue1(50), breaks = my_breaks, fontsize = 4,  border_color = NA,  cluster_rows=FALSE, annotation_col = annotation_col, cluster_cols=TRUE, filename = paste(saveext,"/DimRed/Tb_AllList_EVT_STB",".pdf",sep="") ,width=30,height=4)


