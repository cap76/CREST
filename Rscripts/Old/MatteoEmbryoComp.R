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

#Subset lists based on those in our assay
TbList1 <- intersect(rownames(D),TbMarkers$Gene[which(as.numeric(TbMarkers$'HumanSS2_CTB_vs_STB_adjpval')<0.05)])
TbList2 <- intersect(rownames(D),TbMarkers$Gene[which(as.numeric(TbMarkers$'HumanSS2_CTB_vs_EVT_adjpval')<0.05)])
TbList3 <- intersect(rownames(D),TbMarkers$Gene[which(as.numeric(TbMarkers$'HumanSS2_STB_vs_EVT_adjpval')<0.05)])
PrimedNaiveFormList <- intersect(rownames(D),PrimedNaiveForm$Gene[which(PrimedNaiveForm$P.Value<0.0001)]  )                                                                                                                                      
PrimedNaiveList <- intersect(rownames(D),PrimedNaive$Gene[which(PrimedNaive$FDR<0.0001)] )
hsTEList <- intersect(rownames(D),AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsTE")] )
hsEAm_TEList1 <- intersect(rownames(D),AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E + hsTE")] )
hsEAm_TEList2 <- intersect(rownames(D),AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E (+ hsTE)")] )
hsEAm1 <- intersect(rownames(D),AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E (+ cyAME-L)")] )
hsEAm2 <- intersect(rownames(D),AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E + cyAME-L")] )
hsEAm3 <- intersect(rownames(D),AmnionMarkers$Gene[which(AmnionMarkers$marker=="cyAME-L")] )
hsLAm_TEList <- intersect(rownames(D),AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsTE + cyAME-L")] )

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

#Some misanotation - fix
#Idents(Dsubset1) <- Dsubset1$Genotype
#Dsubset1 <- subset(Dsubset1,idents=c("Epithelial_G3","Stromal_G4","Epithelial_G1"),invert=TRUE)
#Idents(Dsubset1) <- Dsubset1$Cells

#Load in the reference dataset
Dref <- readRDS("/Users/christopherpenfold/Desktop/Data/Embryonic/Human/MagdaDataset/MagdaHumanDataHarmonyhg38_labelled.rds")
Dref$Cells2 <- Idents(Dref)
Idents(Dref) <- Dref$ID3

Dref1 <- subset(Dref,idents="B1")
Dref2 <- subset(Dref,idents="B2")
Dref3 <- subset(Dref,idents="B3")
Dref4 <- subset(Dref,idents="B4")
Dref5 <- subset(Dref,idents="B5")
Dref6 <- subset(Dref,idents="B6")

Idents(Dref1) <- Dref1$Cells2
Idents(Dref2) <- Dref2$Cells2
Idents(Dref3) <- Dref3$Cells2
Idents(Dref4) <- Dref4$Cells2
Idents(Dref5) <- Dref5$Cells2
Idents(Dref6) <- Dref6$Cells2
Idents(Dsubset1) <- Dsubset1$ID3
D1 <- subset(Dsubset1,idents="C1")
D2 <- subset(Dsubset1,idents="C2")
D3 <- subset(Dsubset1,idents="C3")
D4 <- subset(Dsubset1,idents="C4")
D6 <- subset(Dsubset1,idents="C6")
D7 <- subset(Dsubset1,idents="C7")
D8 <- subset(Dsubset1,idents="C8")
Idents(D1) <- D1$Cells
Idents(D2) <- D2$Cells
Idents(D3) <- D3$Cells
Idents(D4) <- D4$Cells
Idents(D6) <- D6$Cells
Idents(D7) <- D7$Cells
Idents(D8) <- D8$Cells

DefaultAssay(Dref3) <-"RNA"
DefaultAssay(Dsubset1) <-"RNA"

#Merge the individual reference datasets based on RNA 
MatteoRef <- merge(Dref1,y=c(Dref2,Dref3,Dref4,Dref5,Dref6),project = "Mes")
MatteoRef <- FindVariableFeatures(MatteoRef, selection.method = "vst", nfeatures = 20000)
MatteoRef <- ScaleData(MatteoRef)
MatteoRef <- RunPCA(MatteoRef, npcs = 20, verbose = FALSE)
MatteoRef <- RunUMAP(MatteoRef, reduction = "pca", dims = 1:20)
MatteoRef <- FindNeighbors(MatteoRef, reduction = "pca", dims = 1:2)

#Plot Matteo's reference data only as I think there are batch issues
p<-DimPlot(MatteoRef,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_MatteoRefEmbryo.pdf",sep=""),width = 40, height = 8,p)
p<-DimPlot(MatteoRef,  pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_MatteoRefEmbryo.pdf",sep=""),width = 40, height = 8,p)
#1-4, 5-6 batch differences

#Split by batch and take a look
MatteoRef1 <- merge(Dref1,y=c(Dref2,Dref3,Dref4),project = "Mes")
MatteoRef2 <- merge(Dref5,y=c(Dref6),project = "Mes")

MatteoRef1 <- FindVariableFeatures(MatteoRef1, selection.method = "vst", nfeatures = 20000)
MatteoRef1 <- ScaleData(MatteoRef1)
MatteoRef1 <- RunPCA(MatteoRef1, npcs = 20, verbose = FALSE)
MatteoRef1 <- RunUMAP(MatteoRef1, reduction = "pca", dims = 1:20)
MatteoRef1 <- FindNeighbors(MatteoRef1, reduction = "pca", dims = 1:2)

MatteoRef2 <- FindVariableFeatures(MatteoRef2, selection.method = "vst", nfeatures = 20000)
MatteoRef2 <- ScaleData(MatteoRef2)
MatteoRef2 <- RunPCA(MatteoRef2, npcs = 20, verbose = FALSE)
MatteoRef2 <- RunUMAP(MatteoRef2, reduction = "pca", dims = 1:20)
MatteoRef2 <- FindNeighbors(MatteoRef2, reduction = "pca", dims = 1:2)

#Some colours for plots
cType	<-	c("Am_d14","Am/EmDisc_d14","EmDisc_d14","Am/EmDisc",
           "EmDisc_d9","EmDisc_d11","EmDisc_d12",
           "Hyp_d9","Hyp_d11","Hyp_d12","Hyp_d14","Hyp/Am",
           "CTB_d9","CTB_d11","CTB_d12","CTB_d14",
           "STB_d9","STB_d11","STB_d12","STB_d14","putSTB",
           "EVT_d14","ExMes_d14","putExMes")

BaseCol	<-	c("#0c9cf5","#0767DA","#0233BF","#0767DA",
             "#0c9cf5","#0767DA","#0233BF",
             "#F04C04","#D74404","#BF3C04","#BF3C04","#BF3C04",
             "#877bd6","#5F54C7","#1A0873","#1A0873",
             "#921FE6","#8017c2","#7108a6","#7108a6","#7108a6",
             "#BF0489","#e6c800","#e6c800")

Idents(Dsubset1) <- Dsubset1$Cells
mammal.anchors <- FindIntegrationAnchors(object.list = list(MatteoRef1,Dsubset1), dims = 1:30, anchor.features = 3000, k.filter = 30)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#


colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

p<-DimPlot(mammal.combined, cols = coluse,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Embryo_Compare_to_MatteoInVit1.pdf",sep=""),width = 40, height = 8,p)

p<-DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Embryo_Compare_to_MatteoInVit1.pdf",sep=""),width = 40, height = 8,p)

Dsubset1$Cells <- Dsubset1$CorrectLabel
Idents(Dsubset1) <- Dsubset1$Cells
mammal.anchors <- FindIntegrationAnchors(object.list = list(MatteoRef2,Dsubset1), dims = 1:30, anchor.features = 3000, k.filter = 30)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20 )
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

p<-DimPlot(mammal.combined, cols = coluse,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Embryo_Compare_to_MatteoInVit2.pdf",sep=""),width = 40, height = 8,p)
p<-DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Embryo_Compare_to_MatteoInVit2.pdf",sep=""),width = 40, height = 8,p)

mammal.combined <- ScaleData(Dsubset1, verbose = FALSE)
mammal.combined <- FindVariableFeatures(mammal.combined, selection.method = "vst", nfeatures = 20000)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)#

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]



p<-DimPlot(mammal.combined, cols = coluse,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Embryo.pdf",sep=""),width = 40, height = 8,p)
p<-DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Embryo.pdf",sep=""),width = 40, height = 8,p)

             
Idents(mammal.combined) <- mammal.combined$Genotype
p<-DimPlot(mammal.combined,   pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Embryo_Genotype.pdf",sep=""),width = 40, height = 8,p)

Idents(mammal.combined) <- mammal.combined$Genotype2
p<-DimPlot(mammal.combined,   pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Embryo_Genotype2.pdf",sep=""),width = 40, height = 8,p)

Idents(mammal.combined) <- mammal.combined$Cells
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

#Clean up the labels
ID4 <- as.character(mammal.combined$Genotype2)
ID4[which(ID4=="NC1_1")] <- "C1_G1"
ID4[which(ID4=="NC1_2")] <- "C1_G2"
ID4[which(ID4=="C1")] <- "C1_Unlabelled"
ID4[which(ID4=="NC2_1")] <- "C2_G1"
ID4[which(ID4=="NC2_2")] <- "C2_G2"
ID4[which(ID4=="NC3_1")] <- "C3_G1"
ID4[which(ID4=="NC3_2")] <- "C3_G2"
ID4[which(ID4=="NC3_3")] <- "C3_G3"
ID4[which(ID4=="C3")] <- "C3_Unlabelled"
ID4[which(ID4=="NC6_1")] <- "C6_G1"
ID4[which(ID4=="NC6_2")] <- "C6_G2"
ID4[which(ID4=="C7")] < "C7"
ID4[which(ID4=="NC8_1")] <- "C8_G1"
ID4[which(ID4=="NC8_2")] <- "C8_G2"
mammal.combined$ID4 <- ID4

p<-DimPlot(mammal.combined, cols = coluse,  pt.size = 4, reduction = "pca", split.by = "ID4", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Embryo_splitbyGT.pdf",sep=""),width = 60, height = 8,p,limitsize = FALSE)

Dhat <- GetAssayData(mammal.combined,assay = "RNA")
Dhat2 <- exp(GetAssayData(mammal.combined,assay = "RNA"))-1

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
  "DCN",
  "MSLN")

list1<-c(
  "PDGFRA",
  "SOX17",
  "POU5F1",
  "SOX15",
  "TFAP2A",
  "TFAP2C",
  "WNT6",
  "CGA",
  "DAB2",
  "DIO2",
  "GATA2",
  "GATA3",
  "HLA-G",
  "JAM3",
  "DPPA3")


list1 <- c("JAM3","GATA2","GATA3","TFAP2A","TFAP2C","CGB8","CGB5","CGA","ERVW-1","PRDM6","TBX3","DIO2","NOTUM","HLA-G","ASCL2","SNAI1",
              "HGF","SNAI2","HAND1","HAND2","PDGFRA","BST2","GATA6","GATA4","CER1","NODAL","LEFTY1","LEFTY2","APOA1","WNT6","VTCN1","BAMBI","AKAP12","PODXL","SFRP1","SOX15","PDGFA","POU5F1")

Idents(mammal.combined) <- mammal.combined$ID4
Dsubs1 <- subset(mammal.combined,idents=c("C1_G1"))
Dsubs2 <- subset(mammal.combined,idents=c("C1_G2"))
Dsubs3 <- subset(mammal.combined,idents=c("C2_G1"))
Dsubs4 <- subset(mammal.combined,idents=c("C2_G2"))
Dsubs5 <- subset(mammal.combined,idents=c("C3_G1"))
Dsubs6 <- subset(mammal.combined,idents=c("C3_G2"))
Dsubs7 <- subset(mammal.combined,idents=c("C3_G3"))
Dsubs8 <- subset(mammal.combined,idents=c("C6_G1"))
Dsubs9 <- subset(mammal.combined,idents=c("C6_G2"))
Dsubs10 <- subset(mammal.combined,idents=c("C7"))
Dsubs11 <- subset(mammal.combined,idents=c("C8_G1"))
Dsubs12 <- subset(mammal.combined,idents=c("C8_G2"))

RN <- rownames(GetAssayData(mammal.combined,assay = "RNA"))

UberList <- c( intersect(PrimedNaiveForm$Gene[1:100], RN),
                      intersect(PrimedNaive$Gene[1:100],  RN),
                      intersect(PrimedNaiveForm$Gene[1:100],  RN),
                      intersect(AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsTE")[1:100]], RN),
                      intersect(AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E + hsTE")[1:100]], RN),
                      intersect(AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E (+ hsTE)")[1:100]], RN),
                      intersect(AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E (+ cyAME-L)")[1:100]], RN),
                      intersect(AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E + cyAME-L")[1:100]], RN),
                      intersect(AmnionMarkers$Gene[which(AmnionMarkers$marker=="cyAME-L")[1:100]], RN) )

#Now get on and do some heatmaps
Idents(mammal.combined) <- mammal.combined$CorrectLabel
DELIST1 <- FindMarkers(mammal.combined, ident.1 = c("EmDisc_d14","Am_d14","Am/EmDisc_d14","Am/EmDisc","Hyp/Am"), ident.2 = c("CTB_d14","STB_d14","EVT_d14","putSTB"), test.use = "MAST" )

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells)), Batch =  factor((mammal.combined$ID3))  )
#rownames(annotation_col) <- colnames( names(mammal.combined$Cells) )
#mycolors <- (c("#BF3C04","#e6c800","#7108a6","#0767DA","#0233FB","#0c9cf5","#1A0873","#BF0489"))
#names(mycolors) <- c("Hyp_d14","ExMes_d14","STB_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","CTB_d14", "EVT_d14")
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Idents(mammal.combined) <- mammal.combined$Cells
#DELIST1 <- FindMarkers(mammal.combined, ident.1 = c("EmDisc_d14","Am_d14","Am/EmDisc_d14"), ident.2 = c("CTB_d14","STB_d14","EVT_d14"), test.use = "MAST" )
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells)), Batch =  factor((mammal.combined$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells)), Batch =  factor((mammal.combined$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells)) , Batch =  factor((mammal.combined$ID3)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells)), Batch =  factor((mammal.combined$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells)), Batch =  factor((mammal.combined$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells)), Batch =  factor((mammal.combined$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList",".pdf",sep="") ,width=10,height=10)



X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells)), Batch =  factor((mammal.combined$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells)) , Batch =  factor((mammal.combined$ID3)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2",".pdf",sep="") ,width=10,height=10)



X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells)), Batch =  factor((mammal.combined$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3",".pdf",sep="") ,width=10,height=10)



mammal.combined1 <- subset(mammal.combined, idents = c("EmDisc_d14","Am_d14","Am/EmDisc_d14","Am/EmDisc","Hyp/Am"))
mammal.combined2 <- subset(mammal.combined, idents = c("STB_d14","CTB_d14","EVT_d14","putSTB"))





X <- (GetAssayData(mammal.combined1,assay = "RNA")) 
Idents(mammal.combined1) <- mammal.combined1$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined1$Cells)) , Batch =  factor((mammal.combined1$ID3)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_EmDiscOnly",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined1,assay = "RNA")) 
Idents(mammal.combined1) <- mammal.combined1$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined1$Cells)) , Batch =  factor((mammal.combined1$ID3)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_EmDiscOnly",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined1,assay = "RNA")) 
Idents(mammal.combined1) <- mammal.combined1$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined1$Cells)) , Batch =  factor((mammal.combined1$ID3)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_EmDiscOnly",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined1,assay = "RNA")) 
Idents(mammal.combined1) <- mammal.combined1$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined1$Cells)), Batch =  factor((mammal.combined1$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_EmDiscOnly",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined1,assay = "RNA")) 
Idents(mammal.combined1) <- mammal.combined1$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined1$Cells)) , Batch =  factor((mammal.combined1$ID3)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_EmDiscOnly",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined1,assay = "RNA")) 
Idents(mammal.combined1) <- mammal.combined1$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined1$Cells)) , Batch =  factor((mammal.combined1$ID3)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_EmDiscOnly",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined1,assay = "RNA")) 
Idents(mammal.combined1) <- mammal.combined1$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined1$Cells)), Batch =  factor((mammal.combined1$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_EmDiscOnly",".pdf",sep="") ,width=10,height=10)









X <- (GetAssayData(mammal.combined2,assay = "RNA")) 
Idents(mammal.combined2) <- mammal.combined2$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined2$Cells)), Batch =  factor((mammal.combined2$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_TbOnly",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined2,assay = "RNA")) 
Idents(mammal.combined2) <- mammal.combined2$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined2$Cells)) , Batch =  factor((mammal.combined2$ID3)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_ETbOnly",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined2,assay = "RNA")) 
Idents(mammal.combined2) <- mammal.combined2$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined2$Cells)) , Batch =  factor((mammal.combined2$ID3)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_TbOnly",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined2,assay = "RNA")) 
Idents(mammal.combined2) <- mammal.combined2$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined2$Cells)) , Batch =  factor((mammal.combined2$ID3)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_TbOnly",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined2,assay = "RNA")) 
Idents(mammal.combined2) <- mammal.combined2$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined2$Cells)) , Batch =  factor((mammal.combined2$ID3)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_TbOnly",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined2,assay = "RNA")) 
Idents(mammal.combined2) <- mammal.combined2$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined2$Cells)) , Batch =  factor((mammal.combined2$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_TbOnly",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined2,assay = "RNA")) 
Idents(mammal.combined2) <- mammal.combined2$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined2$Cells)), Batch =  factor((mammal.combined2$ID3))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_TbOnly",".pdf",sep="") ,width=10,height=10)



#Joint scaling
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$Cells%in%c("Am_d14","EmDisc_d14","Am/EmDisc_d14","Am/EmDisc"))
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_EmDiscOnlyJointScale",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
#subs <- which(mammal.combined$Cells%in%c("Am_d14","EmDisc_d14","Am/EmDisc_d14"))
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$ID3[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_EmDiscOnlyJointScale",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_EmDiscOnlyJointScale",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$ID3[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_EmDiscOnlyJointScale",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_EmDiscOnlyJointScale",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_EmDiscOnlyJointScale",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_EmDiscOnlyJointScale",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_EmDiscOnlyJointScale",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_EmDiscOnlyJointScale",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_EmDiscOnlyJointScale",".pdf",sep="") ,width=10,height=10)




#Joint scaling
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$Cells%in%c("CTB_d14","STB_d14","EVT_d14","putSTB"))
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_TbOnlyJointScale",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_TbOnlyJointScale",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$ID3[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_TbOnlyJointScale",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_TbOnlyJointScale",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_TbOnlyJointScale",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_TbOnlyJointScale",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_TbOnlyJointScale",".pdf",sep="") ,width=10,height=10)




X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_TbOnlyJointScale",".pdf",sep="") ,width=10,height=10)



X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_TbOnlyJointScale",".pdf",sep="") ,width=10,height=10)



X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_TbOnlyJointScale",".pdf",sep="") ,width=10,height=10)

#Now by batch?
#Joint scaling
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C1"))
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch1",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch1",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch1",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch1",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch1",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch1",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch1",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch1",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch1",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch1",".pdf",sep="") ,width=10,height=10)




#Batch 2
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C2"))
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch2",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch2",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch2",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch2",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch2",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch2",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch2",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch2",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch2",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch2",".pdf",sep="") ,width=10,height=10)



#Batch 3
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C3","C4"))
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch3",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch3",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch3",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch3",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch3",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch3",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch3",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch3",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch3",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch3",".pdf",sep="") ,width=10,height=10)

#Batch 6
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C6"))
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch6",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch6",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch6",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch6",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch6",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch6",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch6",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch6",".pdf",sep="") ,width=10,height=35)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch6",".pdf",sep="") ,width=10,height=30)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch6",".pdf",sep="") ,width=10,height=30)




#Batch 7
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C7"))
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch7",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch7",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch7",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch7",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch7",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch7",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch7",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch7",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch7",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch7",".pdf",sep="") ,width=10,height=10)


#Batch 8
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C8"))
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch8",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch8",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch8",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch8",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch8",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch8",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch8",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch8",".pdf",sep="") ,width=10,height=35)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch8",".pdf",sep="") ,width=10,height=30)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch8",".pdf",sep="") ,width=10,height=30)










#Batch and subset by type
#Now by batch?
#Joint scaling
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C1") & mammal.combined$Cells%in%c("CTB_d14","STB_d14","EVT_d14","putSTB") )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch1Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch1Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch1Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch1Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch1Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch1Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch1Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch1Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch1Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch1Tb",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(TbList3,TbList2, TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch1Tb",".pdf",sep="") ,width=10,height=10)



#Batch 2
X <- (GetAssayData(mammal.combined,assay = "RNA")) 

Idents(mammal.combined) <- mammal.combined$Cells
subs <- which(mammal.combined$ID3%in%c("C2") & mammal.combined$Cells%in%c("CTB_d14","STB_d14","EVT_d14","putSTB") )

Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch2Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch2Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch2Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch2Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch2Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch2Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch2Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch2Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch2Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch2Tb",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(TbList3,TbList2,TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch2Tb",".pdf",sep="") ,width=50,height=10)


#Batch 3
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C3","C4") & mammal.combined$Cells%in%c("CTB_d14","STB_d14","EVT_d14","putSTB") )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch3Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch3Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch3Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch3Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch3Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch3Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch3Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch3Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch3Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch3Tb",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(TbList3,TbList2,TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch3Tb",".pdf",sep="") ,width=10,height=10)


#Batch 6
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C6") & mammal.combined$Cells%in%c("CTB_d14","STB_d14","EVT_d14","putSTB") )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch6Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch6Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch6Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch6Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch6Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch6Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch6Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch6Tb",".pdf",sep="") ,width=10,height=35)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch6Tb",".pdf",sep="") ,width=10,height=30)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch6Tb",".pdf",sep="") ,width=10,height=30)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(c(TbList3,TbList2,TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch6Tb",".pdf",sep="") ,width=10,height=30)


#Batch 7
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C7") & mammal.combined$Cells%in%c("CTB_d14","STB_d14","EVT_d14","putSTB") )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch7Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch7Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch7Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch7Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch7Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch7Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch7Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch7Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch7Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch7Tb",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(TbList3,TbList2,TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch7Tb",".pdf",sep="") ,width=10,height=10)

#Batch 8
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C8") & mammal.combined$Cells%in%c("CTB_d14","STB_d14","EVT_d14","putSTB") )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch8Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch8Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch8Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch8Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch8Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch8Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch8Tb",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch8Tb",".pdf",sep="") ,width=10,height=35)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch8Tb",".pdf",sep="") ,width=10,height=30)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch8Tb",".pdf",sep="") ,width=10,height=30)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(c(TbList3,TbList2,TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch8Tb",".pdf",sep="") ,width=10,height=30)









#Batch and subset by type
#Now by batch?
#Joint scaling
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C1") & mammal.combined$Cells%in%c("Am_d14","EmDisc_d14","Am/EmDisc_d14","ExMes_d14","Am/EmDisc","putExMes","Hyp/Am") )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch1EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch1EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch1EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch1EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch1EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch1EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch1EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch1EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch1EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch1EmD",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(TbList3,TbList2, TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch1EmD",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(PrimedNaiveFormList,PrimedNaiveList,hsEAm1,hsEAm2,hsEAm3), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllEmList_Batch1EmD",".pdf",sep="") ,width=10,height=10)




#Batch 2
X <- (GetAssayData(mammal.combined,assay = "RNA")) 

Idents(mammal.combined) <- mammal.combined$Cells
subs <- which(mammal.combined$ID3%in%c("C2") & mammal.combined$Cells%in%c("Am_d14","EmDisc_d14","Am/EmDisc_d14","ExMes_d14","Am/ExMes","putExMes","Hyp/Am")  )

Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch2EmD",".pdf",sep="") ,width=20,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch2EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch2EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch2EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch2EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch2EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch2EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch2EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch2EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch2EmD",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(TbList3,TbList2,TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch2EmD",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(PrimedNaiveFormList,PrimedNaiveList,hsEAm1,hsEAm2,hsEAm3), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllEmList_Batch2EmD",".pdf",sep="") ,width=25,height=10)




#Batch 3
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C3","C4") & mammal.combined$Cells%in%c("Am_d14","EmDisc_d14","Am/EmDisc_d14","ExMes_d14","Am/EmDisc","putExMes","Hyp/Am")  )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch3EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch3EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch3EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch3EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch3EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch3EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch3EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch3EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch3EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch3EmD",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(TbList3,TbList2,TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch3EmD",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(PrimedNaiveFormList,PrimedNaiveList,hsEAm1,hsEAm2,hsEAm3), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllEmList_Batch3EmD",".pdf",sep="") ,width=5,height=10)



#Batch 6
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C6") & mammal.combined$Cells%in%c("Am_d14","EmDisc_d14","Am/EmDisc_d14","ExMes_d14","Am/EmDisc","putExMes","Hyp/Am")  )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch6EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch6EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch6EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch6EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch6EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch6EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch6EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch6EmD",".pdf",sep="") ,width=10,height=35)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch6EmD",".pdf",sep="") ,width=10,height=30)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch6EmD",".pdf",sep="") ,width=10,height=30)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(c(TbList3,TbList2,TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch6EmD",".pdf",sep="") ,width=10,height=30)



X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(PrimedNaiveFormList,PrimedNaiveList,hsEAm1,hsEAm2,hsEAm3), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllEmList_Batch6EmD",".pdf",sep="") ,width=10,height=10)




#Batch 7
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C7") & mammal.combined$Cells%in%c("Am_d14","EmDisc_d14","Am/EmDisc_d14","ExMes_d14","Am/EmDisc","putExMes","Hyp/Am")  )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch7EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch7EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch7EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch7EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch7EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch7EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch7EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch7EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch7EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch7EmD",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(TbList3,TbList2,TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch7EmD",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(PrimedNaiveFormList,PrimedNaiveList,hsEAm1,hsEAm2,hsEAm3), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllEmList_Batch7EmD",".pdf",sep="") ,width=10,height=10)



#Batch 8
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C8") & mammal.combined$Cells%in%c("Am_d14","EmDisc_d14","Am/EmDisc_d14","ExMes_d14","Am/EmDisc","putExMes","Hyp/Am")  )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch8EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch8EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch8EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch8EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch8EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch8EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch8EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch8EmD",".pdf",sep="") ,width=10,height=35)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch8EmD",".pdf",sep="") ,width=10,height=30)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch8EmD",".pdf",sep="") ,width=10,height=30)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(c(TbList3,TbList2,TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch8EmD",".pdf",sep="") ,width=10,height=30)



X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(PrimedNaiveFormList,PrimedNaiveList,hsEAm1,hsEAm2,hsEAm3), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllEmList_Batch8EmD",".pdf",sep="") ,width=10,height=10)






#
subs <- which(mammal.combined$Cells%in%c("CTB_d14","STB_d14","EVT_d14","putSTB")  )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(c(TbList3,TbList2,TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Tb",".pdf",sep="") ,width=30,height=50)


subs <- which(mammal.combined$Cells%in%c("Am_d14","EmDisc_d14","Am/EmDisc_d14","ExMes_d14","Am/EmDisc","Hyp/Am")  )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(PrimedNaiveFormList,PrimedNaiveList,hsEAm1,hsEAm2,hsEAm3), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllEmList_EmD",".pdf",sep="") ,width=30,height=30)











#Batch and subset by type
#Now by batch?
#Joint scaling
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C1") & mammal.combined$Cells%in%c("Am_d14","EmDisc_d14","Am/EmDisc_d14","Am/EmDisc","Hyp/Am") )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch1EmD",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch1EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch1EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch1EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch1EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch1EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch1EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch1EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch1EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch1EmDnoExMes",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(TbList3,TbList2, TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch1EmDnoExMes",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(PrimedNaiveFormList,PrimedNaiveList,hsEAm1,hsEAm2,hsEAm3), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllEmList_Batch1EmDnoExMes",".pdf",sep="") ,width=10,height=10)




#Batch 2
X <- (GetAssayData(mammal.combined,assay = "RNA")) 

Idents(mammal.combined) <- mammal.combined$Cells
subs <- which(mammal.combined$ID3%in%c("C2") & mammal.combined$Cells%in%c("Am_d14","EmDisc_d14","Am/EmDisc_d14","Am/EmDisc","Hyp/Am")  )

Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch2EmDnoExMes",".pdf",sep="") ,width=30,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch2EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch2EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch2EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch2EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch2EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch2EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch2EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch2EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch2EmDnoExMes",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(TbList3,TbList2,TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch2EmDnoExMes",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(PrimedNaiveFormList,PrimedNaiveList,hsEAm1,hsEAm2,hsEAm3), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllEmList_Batch2EmDnoExMes",".pdf",sep="") ,width=30,height=10)




#Batch 3
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C3","C4") & mammal.combined$Cells%in%c("Am_d14","EmDisc_d14","Am/EmDisc_d14","Am/EmDisc","Hyp/Am")  )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch3EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch3EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch3EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch3EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch3EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch3EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch3EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch3EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch3EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch3EmDnoExMes",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(TbList3,TbList2,TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch3EmDnoExMes",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(PrimedNaiveFormList,PrimedNaiveList,hsEAm1,hsEAm2,hsEAm3), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllEmList_Batch3EmDnoExMes",".pdf",sep="") ,width=10,height=10)



#Batch 6
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C6") & mammal.combined$Cells%in%c("Am_d14","EmDisc_d14","Am/EmDisc_d14","Am/EmDisc","Hyp/Am")  )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch6EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch6EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch6EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch6EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch6EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch6EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch6EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch6EmDnoExMes",".pdf",sep="") ,width=10,height=35)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch6EmDnoExMes",".pdf",sep="") ,width=10,height=30)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch6EmDnoExMes",".pdf",sep="") ,width=10,height=30)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(c(TbList3,TbList2,TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch6EmDnoExMes",".pdf",sep="") ,width=10,height=30)



X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(PrimedNaiveFormList,PrimedNaiveList,hsEAm1,hsEAm2,hsEAm3), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllEmList_Batch6EmDnoExMes",".pdf",sep="") ,width=10,height=10)






#Batch 7
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C7") & mammal.combined$Cells%in%c("Am_d14","EmDisc_d14","Am/EmDisc_d14","Am/EmDisc","Hyp/Am")  )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch7EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch7EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch7EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch7EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch7EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch7EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch7EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch7EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch7EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch7EmDnoExMes",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(TbList3,TbList2,TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch7EmDnoExMes",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(PrimedNaiveFormList,PrimedNaiveList,hsEAm1,hsEAm2,hsEAm3), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllEmList_Batch7EmDnoExMes",".pdf",sep="") ,width=10,height=10)



#Batch 8
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C8") & mammal.combined$Cells%in%c("Am_d14","EmDisc_d14","Am/EmDisc_d14","Am/EmDisc","Hyp/Am")  )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveList_Batch8EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(PrimedNaiveFormList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_PrimedNaiveFormList_Batch8EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsTEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])) , Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsTEListList_Batch8EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm1List_Batch8EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm2List_Batch8EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsEAm3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsEAm3List_Batch8EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(hsLAm_TEList, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_hsLAmTEListList_Batch8EmDnoExMes",".pdf",sep="") ,width=10,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList1, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList1_Batch8EmDnoExMes",".pdf",sep="") ,width=10,height=35)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList2, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$ID3[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList2_Batch8EmDnoExMes",".pdf",sep="") ,width=10,height=30)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(TbList3, rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_TbList3_Batch8EmDnoExMes",".pdf",sep="") ,width=10,height=30)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect(c(TbList3,TbList2,TbList1), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllTbList_Batch8EmDnoExMes",".pdf",sep="") ,width=10,height=30)



X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[intersect( c(PrimedNaiveFormList,PrimedNaiveList,hsEAm1,hsEAm2,hsEAm3), rownames(DELIST1)[which(DELIST1$p_val_adj<0.05)]) ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_AllEmList_Batch8EmDnoExMes",".pdf",sep="") ,width=10,height=10)






#Default Assay
DefaultAssay(MatteoRef1) <- "RNA"
DefaultAssay(MatteoRef2) <- "RNA"

HypList1 <- FindMarkers(MatteoRef1,ident.1 = c("Hyp_d12","Hyp_d9","Hyp_d11"), ident.2 = c("EmDisc_d9","EmDisc_d11","EmDisc_d12"), test.use = "MAST", logfc.threshold = log(0))
HypList2 <- FindMarkers(MatteoRef2,ident.1 = c("Hyp_d9","Hyp_d11"), ident.2 = c("EmDisc_d9","EmDisc_d11"), test.use = "MAST", logfc.threshold = log(0))

Idents(D) <- D$Cells
OurList2 <- FindMarkers(D, ident.1 = c("Hyp_d14"),ident.2 = c("Am_d14","EmDisc_d14","Am/EmDisc_d14"),test.use = "MAST", logfc.threshold = log(0))

HypList <- intersect(unique( c( rownames(HypList1)[which(HypList1$p_val_adj<0.05)] ,rownames(HypList2)[which(HypList2$p_val_adj<0.05)] ) ) , rownames(D)) 


#DELIST2


#Batch and subset by type
#Now by batch?
#Joint scaling
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C1") & mammal.combined$Cells%in%c("Hyp_d14","Am_d14","EmDisc_d14","Am/EmDisc_d14","Am/EmDisc","Hyp/Am") )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[rownames(OurList2)[which(OurList2$p_val_adj<0.1)] ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_HypList_Batch1EmDHyp",".pdf",sep="") ,width=10,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C2") & mammal.combined$Cells%in%c("Hyp_d14","Am_d14","EmDisc_d14","Am/EmDisc_d14","Am/EmDisc","Hyp/Am") )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[rownames(OurList2)[which(OurList2$p_val_adj<0.1)] ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_HypList_Batch2EmDHyp",".pdf",sep="") ,width=40,height=10)

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C3","C4") & mammal.combined$Cells%in%c("Hyp_d14","Am_d14","EmDisc_d14","Am/EmDisc_d14","Am/EmDisc","Hyp/Am") )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[rownames(OurList2)[which(OurList2$p_val_adj<0.1)] ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_HypList_Batch3EmDHyp",".pdf",sep="") ,width=5,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C6") & mammal.combined$Cells%in%c("Hyp_d14","Am_d14","EmDisc_d14","Am/EmDisc_d14","Am/EmDisc","Hyp/Am") )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[rownames(OurList2)[which(OurList2$p_val_adj<0.1)] ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_HypList_Batch6EmDHyp",".pdf",sep="") ,width=5,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C7") & mammal.combined$Cells%in%c("Hyp_d14","Am_d14","EmDisc_d14","Am/EmDisc_d14","Am/EmDisc","Hyp/Am") )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[rownames(OurList2)[which(OurList2$p_val_adj<0.1)] ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_HypList_Batch7EmDHyp",".pdf",sep="") ,width=5,height=10)


X <- (GetAssayData(mammal.combined,assay = "RNA")) 
subs <- which(mammal.combined$ID3%in%c("C8") & mammal.combined$Cells%in%c("Hyp_d14","Am_d14","EmDisc_d14","Am/EmDisc_d14","Am/EmDisc","Hyp/Am") )
Idents(mammal.combined) <- mammal.combined$Cells
Xh <- t(scale(t(X)))
Xp <- Xh[rownames(OurList2)[which(OurList2$p_val_adj<0.1)] ,]
Xp <- na.omit(Xp)[,subs]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs])), Batch =  factor((mammal.combined$Genotype2[subs]))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Fullheamtap_All_HypList_Batch8EmDHyp",".pdf",sep="") ,width=10,height=10)



#
TbClustrs <- read_excel("/Users/christopherpenfold/Desktop/EmbryonicDisc.xlsx",sheet = 3) 
#"TP63","KRT7","ELF5"

Idents(mammal.combined) <- mammal.combined$ID3
#cl1 <- c("ACTTTCAAGTACCCTA-1_2","ATCCTATGTAGGCAAC-1_2","ATTATCCCACCCTAGG-1_2","ATTCCCGTCCGTCAAA-1_2","CAACAGTTCGACGACC-1_2","CAACCTCGTCACTTCC-1_2","CAGTTAGCATGAGAAT-1_2","CTAGACATCTCAAAGC-1_2","CTGCTCACACTCACTC-1_2","GTATTGGAGGTTGGAC-1_2","TAATCTCGTTCTCCTG-1_2","TACCCGTCAAACCATC-1_2","TCGCACTAGCCAAGGT-1_2","TTACCATTCTCCAATT-1_2","TTATTGCAGCACAAAT-1_2","TTCACCGGTTCCTAGA-1_2","TTTCAGTGTCAACATC-1_2","TCACACCCACTCCGGA-1_2")
#cl2 <- c("TCTCCGAGTAGGCTCC-1_2","CACTTCGCAGCTGCCA-1_2","AATCGACTCCAGCAAT-1_2","GCCCAGATCGTTGCCT-1_2","GTGAGCCAGCACGTCC-1_2")
#cl3 <- c("AACAACCAGCCTTTCC-1_2","AACACACAGCCTCACG-1_2","ACATTTCAGCTTTCTT-1_2","ACGCACGTCCTACGAA-1_2","ACTCCCACACTCTGCT-1_2","ACTCCCATCAACACGT-1_2","ACTTAGGTCCATTTAC-1_2","AGATAGAGTTCTCTCG-1_2","AGCGATTCACCATTCC-1_2","AGGCATTCAAGACAAT-1_2","AGGGCTCGTGTCACAT-1_2","AGGGTCCGTACAATAG-1_2","ATAGACCTCTGCTTAT-1_2","ATATCCTCAATAGTAG-1_2","ATCACGAGTCCACAGC-1_2","ATCATTCAGGGATGTC-1_2","ATGCGATTCTATTGTC-1_2","ATTCCCGTCCAAGCTA-1_2","ATTTCACTCCTAAACG-1_2","CAACCTCAGTCCTGTA-1_2","CAAGAGGTCCACCCTA-1_2","CACACAACAGGATTCT-1_2","CACGTGGTCAGCTGTA-1_2","CACGTGGTCCCTCATG-1_2","CAGCAGCTCCTCCACA-1_2","CATCCCAAGTCGCGAA-1_2","CCACACTTCTCTATGT-1_2","CCACCATTCCTACCGT-1_2","CCACTTGAGTAAACTG-1_2","CCACTTGCATCGGTTA-1_2","CCATCACCAACAAGTA-1_2","CCCAACTTCTACTTCA-1_2","CCCTCTCAGCCATCCG-1_2","CCGTTCAGTTACGCCG-1_2","CCTACGTAGTCTAGCT-1_2","CCTCTCCAGTCCCTAA-1_2","CGACAGCCATCTTCGC-1_2","CGCCATTTCGTCTACC-1_2","CTACCTGTCCCTTTGG-1_2","CTCAGAATCCCTTCCC-1_2","CTCATGCTCCCAAGTA-1_2","CTCCTTTTCTGCTTAT-1_2","CTGAGGCCACGGGCTT-1_2","CTGCCATTCTTACGTT-1_2","CTTAGGACATATGCGT-1_2","GAGACTTTCCATTCAT-1_2","GAGCTGCCACTCTGCT-1_2","GAGGCAACAATATCCG-1_2","GATTCTTTCCCATAAG-1_2","GCACTAATCTCCGAAA-1_2","GCATCTCTCTTCCCAG-1_2","GCCTGTTAGCCGTTGC-1_2","GCCTGTTTCTATTTCG-1_2","GGAGGATTCGCCACTT-1_2","GGATGTTTCCATCTAT-1_2","GGCTTTCAGCATAGGC-1_2","GGGAAGTCACGTAGTT-1_2","GGTGATTCACCATATG-1_2","GTCACTCTCGGTGTAT-1_2","GTGCTTCTCGGTAAGG-1_2","GTTACCCAGTCTAGCT-1_2","GTTGCTCAGGTACTGG-1_2","GTTTACTTCCTAGCCT-1_2","TACAACGCAACCAATC-1_2","TACGGTATCCTCAGGG-1_2","TATTCCAGTTACACTG-1_2","TCAATCTCACACCGCA-1_2","TCAATTCAGGGTTAAT-1_2","TCACTATTCGCTCTAC-1_2","TCAGCCTTCTTACACT-1_2","TCAGGTATCTGTACAG-1_2","TCAGTCCGTGACTCTA-1_2","TCATCCGGTAAGCGGT-1_2","TCATGCCAGACATCCT-1_2","TCATTTGTCATCTACT-1_2","TCCACCATCCTAAACG-1_2","TCCATGCAGCCATTCA-1_2","TCCATGCGTGCAATGG-1_2","TCCCACACAGTTGTCA-1_2","TCCGTGTGTAAGAACT-1_2","TCCTCCCGTGAGAACC-1_2","TCCTTTCAGACAAGCC-1_2","TCCTTTCTCTGTCCGT-1_2","TCGCTCATCTTTGCTA-1_2","TCTAACTTCTTCCACG-1_2","TCTTAGTTCCAGCCTT-1_2","TCTTCCTCAAGTCATC-1_2","TGATCAGTCTTTGCGC-1_2","TGTTGGATCACTCCGT-1_2","TTACCGCTCTGTAAGC-1_2","TTACGTTTCCATCACC-1_2","TTCAATCGTAGATGTA-1_2","TTCACGCCATTCTCTA-1_2","TTCCTCTTCGCACGAC-1_2","TTCTCTCCAGCAATTC-1_2","TTCTTCCAGTCATCGT-1_2","TTTACCATCATCGCCT-1_2","TTTCCTCTCATCCTAT-1_2")
#cl1 <- str_replace(cl1,"-","-")
#cl2 <- str_replace(cl2,"-","-")
#cl3 <- str_replace(cl3,"-","-")
Ds1 <- subset(mammal.combined,idents="C1")
Idents(Ds1) <- Ds1$CorrectLabel
Idents(Ds1,cells=WhichCells(Ds1,idents=c("putExMes"))) <- "Am/EmDisc_d14"
Idents(Ds1,cells=WhichCells(Ds1,idents=c("Hyp/Am"))) <- "Am/EmDisc_d14"
Idents(Ds1,cells=WhichCells(Ds1,idents=c("ExMes_d14"))) <- "Am/EmDisc_d14"
Idents(Ds1,cells=WhichCells(Ds1,idents=c("Am_d14"))) <- "Am/EmDisc_d14"
Idents(Ds1,cells=WhichCells(Ds1,idents=c("putSTB"))) <- "STB_d14"
Idents(Ds1,cells=WhichCells(Ds1,idents=c("EmDisc_d14"))) <- "Am/EmDisc_d14"

#TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 1" & TbClustrs$Cl=="Cl1")]

Idents(Ds1,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 1" & TbClustrs$Cl=="Cl1")]) <- "STB_Cl1_d14"
Idents(Ds1,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 1" & TbClustrs$Cl=="Cl2")]) <- "STB_Cl2_d14"

#Idents(Ds1,cells=cl1) <- "Cl1"
#Idents(Ds1,cells=cl2) <- "Cl2"
#Idents(Ds1,cells=cl3) <- "Cl3"
#Cl2 Tb
#Hyp/Am EmDisc?

Ds1 <- FindVariableFeatures(Ds1, selection.method = "vst", nfeatures = 3000)
Ds1 <- ScaleData(Ds1, verbose = FALSE)
Ds1 <- RunPCA(Ds1, npcs = 20, verbose = FALSE)
#Ds1 <- RunUMAP(Ds1, reduction = "pca", dims = 1:20)
#Ds1 <- FindNeighbors(Ds1, reduction = "pca", dims = 1:20)#
#colind <- integer( length( levels(Idents(Ds1)) )  )
#for (i in 1:length( levels(Idents(Ds1)) ) ) {
#  colind[i] <- which(cType==levels(Idents(Ds1))[i])
#}
#coluse <- BaseCol[colind]
p<-DimPlot(Ds1,  pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "POU5F1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_OCT4.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "WNT6", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_WNT6.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "DIO2", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_DIO2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "HLA-G", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_HLA-G.pdf",sep=""),width = 10, height = 8,p)



MiniList <- c("DNMT3A","DPPA3","OTX2","POU5F1","PDGFA","PODXL","SFRP1","SOX15","GABRP","VTCN1","WNT6","BAMBI","NOTUM","JAM3",
              "TACSTD2","ENPEP","GATA2","GATA3","TFAP2C","TFAP2A","TEAD3","TP63","OVOL1","CCKBR","CGA","CGB5","CGB8","SDC1","HLA-G",
              "SOX17","GATA4","CER1","LEFTY1","LEFTY2","APOA1","APOB","LINC00261","HHEX","FOXF1","PDGFRA","HAND2","BST2","GATA6","HGF","PRDM1","ASCL2","SNAI1")
X <- GetAssayData(Ds1)
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(Ds1) )), Batch =  factor(( Idents(Ds1)))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch1_markerplot",".pdf",sep="") ,width=10,height=7)

#MiniList <- c("SOX17","TFAP2A","TFAP2C","PRDM1","POU5F1","NANOS3","KIT","DND1","SOX15","DPPA3","CD38","ITGB3","ALPL","KLF4")

MiniList <- c("SOX17","TFAP2A","TFAP2C","PRDM1","POU5F1","NANOS3","KIT","DND1","SOX15","DPPA3","CD38","ITGB3","ALPL","KLF4")
X <- GetAssayData(Ds1)
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(Ds1) )), Batch =  factor(( Idents(Ds1)))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch1_markerplot2",".pdf",sep="") ,width=10,height=3)




p<-FeaturePlot(Ds1, features = "CGA", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_CGA.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "JAM3", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_JAM3.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "VIM", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_VIM.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "HAND2", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_HAND2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "SOX15", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_SOX15.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "PDGFA", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_PDGFA.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "GATA4", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_GATA4.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds1, features = "TTR", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_TTR.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds1, features = "APOA1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_APOA1.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "APOB", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_APOB.pdf",sep=""),width = 10, height = 8,p)

#CTB
p<-FeaturePlot(Ds1, features = "TEAD3", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_TEAD3.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "TP63", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_TP63.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "OVOL1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_OVOL1.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "CCKBR", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_CCKBR.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "SIGLEC6", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_SIGLEC6.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "IFI6", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_IFI6.pdf",sep=""),width = 10, height = 8,p)

#STB
p<-FeaturePlot(Ds1, features = "CGB", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_CGB.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "SDC1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_SDC1.pdf",sep=""),width = 10, height = 8,p)

#
p<-FeaturePlot(Ds1, features = "SOX17", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_SOX17.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "HHEX", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_HHEX.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds1, features = "FOXF1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_HHEX.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds1, features = "ENPEP", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_ENPEP.pdf",sep=""),width = 10, height = 8,p)
#CTb-specific markers (TEAD3, TP63, OVOL1)14, STb-specific markers (CGa, CGb, SDC1)
p<-FeaturePlot(Ds1, features = "TACSTD2", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_TACSTD2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "LINC00261", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_LINC00261.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "FOXF1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_FOXF1.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "GABRP", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_GABRP.pdf",sep=""),width = 10, height = 8,p)

#LINC00261

#cl1 <- c("AACTTCTGTAGCCCTG-1_3","ACTATTCGTCTGTGAT-1_3","ACTTAGGCAGTAACGG-1_3","AGATAGATCTCCCTAG-1_3","AGGATCTAGCGTATAA-1_3","AGGTCATCAGAAGTTA-1_3","ATACTTCTCATGCCAA-1_3","ATCCACCTCATTATCC-1_3","ATTCATCGTGTCTAAC-1_3","CAAGGGACATCTGGGC-1_3","CCACCATAGTGAGTTA-1_3","CGAGAAGCATTGACAC-1_3","CTGGTCTGTATGAGCG-1_3","GATTGGTGTCCAGGTC-1_3","GCCAGGTAGCTTGTGT-1_3","GCCTGTTTCATGGCCG-1_3","GCTACAACAGTAGAAT-1_3","GGAAGTGCATGACAGG-1_3","GGAATCTGTCCTGTCT-1_3","TACTGCCTCGACTCCT-1_3","TGGAACTCAGGTGTTT-1_3","TTGGGATCATGGATCT-1_3","TTTGACTAGGGTCTTT-1_3")
#cl2 <- c("AAACGAAAGCAACTCT-1_3","AACACACCAATGAACA-1_3","AACCAACTCATTGAGC-1_3","AACGTCACAGGTCAAG-1_3","AAGATAGTCTGGTTGA-1_3","AAGCGTTAGCTTGTTG-1_3","AAGTACCAGATTGATG-1_3","AATTCCTTCTAGTTCT-1_3","ACGATCACACCGCTGA-1_3","ACGGTCGGTAGCGTAG-1_3","ACGTACAAGCGCATCC-1_3","ACTGCAAAGACCAACG-1_3","ACTTCGCTCGCGAAGA-1_3","AGCCAGCAGTGACACG-1_3","AGCGTATCAAGGCGTA-1_3","AGGACGAGTTAACAGA-1_3","AGGCATTTCCCGGTAG-1_3","AGGGCTCAGACGTCGA-1_3","AGGTCATCACTCAGAT-1_3","AGGTGTTTCCTTCTTC-1_3","AGTAGCTCAGGACATG-1_3","AGTCAACTCACCATAG-1_3","AGTCTCCGTAGCTTGT-1_3","AGTGCCGAGCTGACAG-1_3","AGTTCGAAGTGCAACG-1_3","ATATCCTTCTGAGATC-1_3","ATCCACCGTGAATGTA-1_3","ATCGGCGGTAGGCAAC-1_3","ATCGGCGTCTGCACCT-1_3","ATGACCACAACCCGCA-1_3","ATGCATGTCTCTTAAC-1_3","ATTATCCGTCCATACA-1_3","ATTCCCGGTTCCACAA-1_3","CAAGACTTCATCCTAT-1_3","CAAGACTTCTTCCGTG-1_3","CAAGAGGAGTAGGCCA-1_3","CAATCGATCCCTTGGT-1_3","CACACAAAGGAGGGTG-1_3","CACATGACAGCAGTTT-1_3","CACATGAGTTCTCCTG-1_3","CACCGTTTCGGTCACG-1_3","CACGTGGAGAAGCGCT-1_3","CAGCAATCATTGTCGA-1_3","CAGCAATGTGTAAACA-1_3","CAGTGCGAGTAGCCAG-1_3","CAGTTAGAGTTTCGGT-1_3","CATACTTAGGGCTAAC-1_3","CATAGACGTACGGATG-1_3","CATGCCTGTCACAGAG-1_3","CATTCTACATCGGATT-1_3","CATTGCCTCCGATAAC-1_3","CCAAGCGAGGCTCACC-1_3","CCACGAGCAACTTGGT-1_3","CCTCAGTAGGGCGAGA-1_3","CCTCAGTGTTCAGCTA-1_3","CCTTCAGCACATCCCT-1_3","CGAATTGCAAGACGGT-1_3","CGAGTGCAGCGCATCC-1_3","CGAGTGCTCGATGCTA-1_3","CGAGTTAGTGTAGTGG-1_3","CGGAACCCAACGACAG-1_3","CGGAGAACACTCAAGT-1_3","CGTTCTGTCCATTTGT-1_3","CTACCCACAGCGTTTA-1_3","CTCAGGGTCAACGTGT-1_3","CTCCACAAGAGGCTGT-1_3","CTCCACAAGTGGATAT-1_3","CTCCCTCCAGACCATT-1_3","CTCGAGGCAGTTTCAG-1_3","CTCTCAGAGTCCGCCA-1_3","CTGCCATCAAAGCGTG-1_3","CTGTGAATCAGACATC-1_3","CTTTCAATCATAGCAC-1_3","GACCCAGAGGGCCAAT-1_3","GACTATGAGCTAGATA-1_3","GAGACCCGTCGTACTA-1_3","GAGTTTGAGGAGTACC-1_3","GATCACATCTAGCATG-1_3","GATCCCTAGACATACA-1_3","GCACGGTAGGGACACT-1_3","GCACGTGAGGTGGCTA-1_3","GCCATTCCATGATGCT-1_3","GCCGTGAGTCAAGCCC-1_3","GCCTGTTAGTCCCGGT-1_3","GCGGAAAGTTGTTGAC-1_3","GCTGCAGTCCGTGCGA-1_3","GCTTCACCATGTCAGT-1_3","GGAACCCTCCAATCTT-1_3","GGAATCTTCCTGGGTG-1_3","GGACGTCAGGGCATGT-1_3","GGATGTTTCATACGGT-1_3","GGCACGTTCATACGGT-1_3","GGGCCATCATGAGGGT-1_3","GGGTAGACATCATTTC-1_3","GGGTATTTCCACGAAT-1_3","GGGTGAAAGTATGAGT-1_3","GGTAATCTCGGCGATC-1_3","GGTGTTAGTTTCGATG-1_3","GGTTAACCAGATCATC-1_3","GTAGTACTCATCGCCT-1_3","GTATTGGTCTCCTGAC-1_3","GTATTTCTCGTAGAGG-1_3","GTCAGCGGTCAAGTTC-1_3","GTCATTTGTATCGCGC-1_3","GTCTCACAGCGGTAGT-1_3","GTGAGTTAGCCTTGAT-1_3","GTGAGTTGTACATACC-1_3","GTTCATTGTTCCGGTG-1_3","TAACTTCGTACGACAG-1_3","TAAGTCGGTCCATAGT-1_3","TACAGGTTCACGGACC-1_3","TACGGGCAGGTAATCA-1_3","TACGGGCCATGGGTCC-1_3","TACTTACAGGTAGTAT-1_3","TAGTGCAAGATACCAA-1_3","TATCGCCTCTCTGGTC-1_3","TATTGCTCACAGAGAC-1_3","TCACAAGAGGACGCAT-1_3","TCATCCGGTGTGGACA-1_3","TCATGTTTCGTGCACG-1_3","TCATTACCAAGTATAG-1_3","TCCGATCCACGCGCAT-1_3","TCCGGGAAGGAATGTT-1_3","TCCTGCAGTGGCTTGC-1_3","TCTACCGTCGATACAC-1_3","TCTTCCTAGGGAGGCA-1_3","TGAATGCGTTCGGACC-1_3","TGAGCATTCATGAGTC-1_3","TGCTCCACACGCTGCA-1_3","TGGGAAGAGAATTGCA-1_3","TGGGAAGAGCTGAGCA-1_3","TGGTGATTCACGGACC-1_3","TGTCAGACATACTGAC-1_3","TTACCATTCATGTCTT-1_3","TTCAGGATCTGTGTGA-1_3","TTCCTCTAGCAATAAC-1_3","TTGAGTGGTACGAAAT-1_3","TTGTTTGCAGACTCTA-1_3","TTTGTTGTCTACAGGT-1_3")
#cl3 <- c("AACACACAGGGCAATC-1_3","AAGTTCGTCCACGTGG-1_3","ACACAGTTCGATACTG-1_3","ACACAGTTCTGCTTAT-1_3","ACACTGATCCATACTT-1_3","ACCAACATCATCACTT-1_3","ACCCTTGAGGTTTACC-1_3","ACGATGTCACAGCCAC-1_3","ACGTACACAGTTGCGC-1_3","ACGTTCCCAATCTGCA-1_3","ACTATCTTCGTTGTTT-1_3","ACTATTCAGTCTACCA-1_3","ACTTATCAGGTAAGGA-1_3","ACTTCGCGTCTAGATC-1_3","AGACCCGAGAGCAACC-1_3","AGCTACATCTCGTGAA-1_3","AGCTTCCTCGCAGTTA-1_3","AGTACCACAACGGCTC-1_3","ATAGGCTCACACGGAA-1_3","ATATCCTAGACTCAAA-1_3","ATCCTATAGACAGTCG-1_3","ATCGGATAGAGGGTAA-1_3","ATCGTGACAGACCAGA-1_3","ATCTCTAGTACCTATG-1_3","ATCTTCACACTTGAGT-1_3","ATGCATGCAGGGTCTC-1_3","ATTACTCCAATACCTG-1_3","ATTCTTGTCAGACCCG-1_3","ATTGTTCCAGCGTGCT-1_3","ATTTCACTCACACCGG-1_3","CACAACAAGTAAACAC-1_3","CATACTTAGGTCGTGA-1_3","CATGCCTCAGGGATAC-1_3","CATGCGGGTGTCTTCC-1_3","CATTGCCTCATCAGTG-1_3","CCCTGATAGCTCTGTA-1_3","CCGGTGAAGCAATTAG-1_3","CCGTAGGGTAAGTTAG-1_3","CCTACGTAGGACTTCT-1_3","CCTATCGAGCACTTTG-1_3","CCTATCGAGTTGCCCG-1_3","CCTCAACTCTTGCAAG-1_3","CCTCTCCTCGTTGTGA-1_3","CCTTTGGAGGTTGTTC-1_3","CGAAGGATCACCACAA-1_3","CGAAGTTAGTGCTCAT-1_3","CGAGGCTAGGACAACC-1_3","CGATGCGCATCATTGG-1_3","CGCAGGTCAGAAACCG-1_3","CGCAGGTTCGGCCTTT-1_3","CGGGTGTAGTTGTAGA-1_3","CGTAGTACAATACCTG-1_3","CTAACCCAGATAGTGT-1_3","CTAACCCAGGCACCAA-1_3","CTCAGGGCACGATTCA-1_3","CTCATCGCATCAACCA-1_3","CTCATGCCACACACGC-1_3","CTCCAACGTTCTCCAC-1_3","CTCCCTCCAGCACAAG-1_3","CTTCAATTCCTCTCTT-1_3","GAAATGAGTCTTGCGG-1_3","GACCCTTTCTAACACG-1_3","GACGCTGTCGAGCCAC-1_3","GAGACCCCACTGCGAC-1_3","GAGACTTCATCAGCGC-1_3","GAGTTACAGCAGCCCT-1_3","GATCACAGTCACTTCC-1_3","GATCATGAGTATCTGC-1_3","GCACATAGTCATCCCT-1_3","GCACGGTCAAAGCGTG-1_3","GCCCGAAGTTGTTTGG-1_3","GGACGTCCAGCCCACA-1_3","GGAGGTATCTAGACCA-1_3","GGCGTCACAATGTTGC-1_3","GTAACCACATTCTGTT-1_3","GTAGCTACAACAGAGC-1_3","GTAGGTTGTCCTCATC-1_3","GTCTACCCAACACGTT-1_3","GTGATGTCACCCTCTA-1_3","GTGGTTATCCGGACTG-1_3","GTGTAACCAAGCAATA-1_3","GTGTGATGTCATATGC-1_3","TAATTCCCACACCAGC-1_3","TACCTGCCATGCCGCA-1_3","TACCTGCGTGTTTCTT-1_3","TACGGTATCATGCAGT-1_3","TACTTACGTACCGGCT-1_3","TACTTCACATGAATCC-1_3","TAGACCAAGGAGACCT-1_3","TCAAGACCAATCTCTT-1_3","TCACAAGTCACGAGGA-1_3","TCACACCCAGGTTACT-1_3","TCACGCTTCAGGAAAT-1_3","TCATACTAGTATAACG-1_3","TCCTCGATCTGCGAGC-1_3","TCGAACACAAGAGCTG-1_3","TCGATTTCAATTTCCT-1_3","TCGCACTCAAATAGCA-1_3","TCGGGTGAGATTCGCT-1_3","TCGGTCTTCCACCTCA-1_3","TCGTGGGCAACGGCCT-1_3","TCTGGCTGTACTCAAC-1_3","TGAACGTAGTTTGAGA-1_3","TGACTCCTCAGTGCGC-1_3","TGAGTCATCATCTACT-1_3","TGATCTTCACACCTGG-1_3","TGATCTTCACACCTGG-1_3","TGCAGGCAGGAGGCAG-1_3","TGCCGAGCATCATTTC-1_3","TGGAGAGAGTATAACG-1_3","TGGTGATCAATCTCGA-1_3","TGTACAGGTTAAGTCC-1_3","TGTCCCAAGTGAATAC-1_3","TTCATGTCACCTTCGT-1_3","TTCTAGTGTAGGAAAG-1_3","TTCTTGACACGAAGAC-1_3","TTGAACGTCTGGCCAG-1_3","TTGACCCAGCTGTCCG-1_3","TTGACCCTCCTTGGAA-1_3","TTTACCATCATCGGGC-1_3","TTTATGCTCCGGACGT-1_3")
#cl4 <- c("AAACCCAAGTAGCAAT-1_3","AACAGGGGTCCAGCCA-1_3","AACCAACCACGTCATA-1_3","AACTTCTAGGACATCG-1_3","AAGGAATTCACTTGTT-1_3","AATCGACAGACGAAGA-1_3","AATCGTGAGCTCGACC-1_3","AATGACCGTTCGGACC-1_3","AATGCCACACCAGTAT-1_3","AATTCCTGTTACGCCG-1_3","ACAAAGACAGGAGGTT-1_3","ACACTGAGTCTGGTTA-1_3","ACAGAAAAGACTCGAG-1_3","ACCTGTCAGTGAGTGC-1_3","ACGGTTAAGGACTATA-1_3","ACGTAGTGTGACGTCC-1_3","ACTATTCAGCACCGTC-1_3","ACTCTCGTCCACTAGA-1_3","ACTGATGAGTCATTGC-1_3","AGACCCGCAGAGGAAA-1_3","AGACCCGCATCGCTAA-1_3","ATACCGAAGGCCCAAA-1_3","ATCCCTGGTACTCGTA-1_3","ATCGCCTTCCGGTTCT-1_3","ATCGGCGGTCCTGAAT-1_3","ATTCACTGTAATTGGA-1_3","ATTCGTTCAGGAACCA-1_3","ATTCGTTTCCTTATGT-1_3","ATTTCACTCTGCACCT-1_3","CAAAGAAGTCATCAGT-1_3","CAAGACTAGGACGCTA-1_3","CAAGAGGAGTTGGAGC-1_3","CACAACACAGGACATG-1_3","CACCGTTCACAATGTC-1_3","CACGGGTCATCCTGTC-1_3","CAGATACCATGACAGG-1_3","CAGCAATAGTTCTACG-1_3","CATCCGTTCTTGGCTC-1_3","CATGCCTCACGTACTA-1_3","CCACCATGTTCAGGTT-1_3","CCCATTGCATCTGCGG-1_3","CCCATTGGTTGCATAC-1_3","CCCTAACCACCGTGAC-1_3","CCGGACAGTATGAGGC-1_3","CCGGGTACAATCGCAT-1_3","CCTAAGACAGGTCCGT-1_3","CCTCAACTCGTCACCT-1_3","CCTCTAGCACCGCTAG-1_3","CCTGCATCAAGTATCC-1_3","CCTTGTGAGCCTGGAA-1_3","CGCCATTCACATTCTT-1_3","CGGGCATCAACCAGAG-1_3","CGGGTCATCAGTCACA-1_3","CGGTCAGGTTGTAAAG-1_3","CGTAATGGTAAGATAC-1_3","CTAACCCAGTAGAGTT-1_3","CTACGGGAGCCGTCGT-1_3","CTATAGGGTGCCTACG-1_3","CTATCCGTCCTGTTGC-1_3","CTCAGAAAGTGATAGT-1_3","CTCCAACTCACCCTCA-1_3","CTCTCAGTCGCCAATA-1_3","CTGAATGTCGCTTGAA-1_3","CTGATCCGTTGCGTAT-1_3","CTGCTCATCAGTGATC-1_3","CTGTATTGTTAAGCAA-1_3","CTGTGGGTCAGCATTG-1_3","CTTCAATAGGATACGC-1_3","CTTCTAACAAGATCCT-1_3","CTTTCAACAACAGTGG-1_3","GAACTGTGTATGGAGC-1_3","GAAGAATAGGAGAGGC-1_3","GACCAATGTGTCATGT-1_3","GACCCTTAGGACTTCT-1_3","GACCTTCCACTTGGCG-1_3","GAGGGTAAGGACTATA-1_3","GAGTCATCACGGGCTT-1_3","GAGTTACTCTCCGCAT-1_3","GAGTTTGAGACTCCGC-1_3","GATAGCTTCAAGCCTA-1_3","GATCCCTCACTACGGC-1_3","GATGAGGAGTTCTACG-1_3","GATGATCAGAAGAGCA-1_3","GCCAGCATCACTACGA-1_3","GCGGAAAAGTTGCGCC-1_3","GCGTTTCAGTATAACG-1_3","GCTGAATCAATGTCTG-1_3","GCTGAATTCTCATAGG-1_3","GGATCTACAAGTATAG-1_3","GGATGTTAGATCGCTT-1_3","GGGACAACAATCGCGC-1_3","GGGTAGAAGTAGGGTC-1_3","GGGTTATCAAATAAGC-1_3","GGTCACGTCTGAGTCA-1_3","GGTGAAGAGGACTTCT-1_3","GGTGAAGAGTGGTGGT-1_3","GGTGATTGTGACTCGC-1_3","GGTGTCGGTACCCGAC-1_3","GGTTAACTCCGGCAGT-1_3","GTAACCAGTCAAGGCA-1_3","GTCAAACTCTCAGTCC-1_3","GTCGAATCAGCTACAT-1_3","GTCGAATCATTGTCGA-1_3","GTCTACCAGAGCCTGA-1_3","GTCTCACAGCGGTATG-1_3","GTCTCACGTATACAGA-1_3","GTGAGCCGTCCTGGTG-1_3","GTGCACGCACAGCATT-1_3","GTTGCGGTCATGCGGC-1_3","GTTTACTAGAGGCGTT-1_3","TAACACGTCACGTAGT-1_3","TAAGTCGGTCACCTTC-1_3","TACTTACCACCCTGAG-1_3","TATACCTCACCAATTG-1_3","TATACCTTCCGTGGGT-1_3","TATATCCGTGGGTCAA-1_3","TCAAGACTCTCCAATT-1_3","TCAAGTGGTCACGTGC-1_3","TCACATTAGGAAGTAG-1_3","TCAGCCTCATCCGGTG-1_3","TCAGTTTTCTCCAATT-1_3","TCATCATGTACCCGCA-1_3","TCATCCGGTGCTATTG-1_3","TCATGGAAGTCTGGAG-1_3","TCATGTTTCCTTCGAC-1_3","TCATGTTTCGAGTCCG-1_3","TCCTCCCCAAACGAGC-1_3","TCGACCTTCTCTCGAC-1_3","TCTACCGGTATTCCGA-1_3","TCTCCGACATTGCAAC-1_3","TCTTAGTGTGGACTAG-1_3","TCTTCCTAGTCTAACC-1_3","TGAGACTAGACCACGA-1_3","TGAGCATGTGGCTTGC-1_3","TGATGGTAGTATAACG-1_3","TGCACGGAGCGAAACC-1_3","TGCATCCGTAGTGGCA-1_3","TGCGACGCACATTGTG-1_3","TGCTTCGCATCGAACT-1_3","TGGATGTAGGGTCAAC-1_3","TGTTACTAGCCTCTCT-1_3","TGTTACTCATAGTCAC-1_3","TGTTCTACAGCGTGAA-1_3","TTAGGGTGTTTCCATT-1_3","TTCACCGCACACCTAA-1_3","TTCATGTGTTGAATCC-1_3","TTCCAATCACATACGT-1_3","TTCCAATGTGCGGTAA-1_3","TTCCTCTGTTCCGTTC-1_3","TTCTCTCAGACTACGG-1_3","TTGCATTAGGGCATGT-1_3","TTGGGCGTCATCTACT-1_3","TTTCATGTCCACACCT-1_3")
#newcl1 <- c("AGATAGATCTCCCTAG-1_3","ACTATTCGTCTGTGAT-1_3","GGAATCTGTCCTGTCT-1_3","AACTTCTGTAGCCCTG-1_3","TGGAACTCAGGTGTTT-1_3","CCACCATAGTGAGTTA-1_3","ATCCACCTCATTATCC-1_3","AGGTCATCAGAAGTTA-1_3","TACTGCCTCGACTCCT-1_3","TTGGGATCATGGATCT-1_3","GATTGGTGTCCAGGTC-1_3","TTTGACTAGGGTCTTT-1_3","ACTTAGGCAGTAACGG-1_3","GCCTGTTTCATGGCCG-1_3","CTGGTCTGTATGAGCG-1_3","GCTACAACAGTAGAAT-1_3","CAAGGGACATCTGGGC-1_3","GGAAGTGCATGACAGG-1_3","AGGATCTAGCGTATAA-1_3","CATTGCCGTAGCTGCC-1_3","GCCAGGTAGCTTGTGT-1_3","ATTCATCGTGTCTAAC-1_3")
#newcl1 <- str_replace(newcl1,"-","-")
#newcl2 <- c("TGTACAGGTTAAGTCC-1_3","GTGTGATGTCATATGC-1_3","TTCTAGTGTAGGAAAG-1_3","CGTAGTACAATACCTG-1_3","TGTCCCAAGTGAATAC-1_3","ATATCCTAGACTCAAA-1_3","ACTTATCAGGTAAGGA-1_3","ATTGTTCCAGCGTGCT-1_3","CCTATCGAGTTGCCCG-1_3","ATGCATGCAGGGTCTC-1_3","TGCAGGCAGGAGGCAG-1_3","ACACAGTTCTGCTTAT-1_3","CCCTGATAGCTCTGTA-1_3","ACGATGTCACAGCCAC-1_3","ACGTACACAGTTGCGC-1_3","GTGTAACCAAGCAATA-1_3","GTAACCACATTCTGTT-1_3","AGTACCACAACGGCTC-1_3","CATTGCCTCATCAGTG-1_3","ATTTCACTCACACCGG-1_3","ATCGGATAGAGGGTAA-1_3","GTGATGTCACCCTCTA-1_3","GAGACTTCATCAGCGC-1_3","CTTCAATTCCTCTCTT-1_3","ACCAACATCATCACTT-1_3","GACCCTTTCTAACACG-1_3","CCTCAACTCTTGCAAG-1_3","GCACGGTCAAAGCGTG-1_3","ATTACTCCAATACCTG-1_3","TCTGGCTGTACTCAAC-1_3","TCAAGACCAATCTCTT-1_3","CGAAGTTAGTGCTCAT-1_3","TGGAGAGAGTATAACG-1_3","GATCATGAGTATCTGC-1_3","TTGACCCTCCTTGGAA-1_3","TCGGTCTTCCACCTCA-1_3","GCCCGAAGTTGTTTGG-1_3","ATCCTATAGACAGTCG-1_3","GACGCTGTCGAGCCAC-1_3","GATCACAGTCACTTCC-1_3","CTCATGCCACACACGC-1_3","TACCTGCCATGCCGCA-1_3","CATACTTAGGTCGTGA-1_3","AGACCCGAGAGCAACC-1_3","TGAACGTAGTTTGAGA-1_3","ATCGTGACAGACCAGA-1_3","CATGCGGGTGTCTTCC-1_3","TCGCACTCAAATAGCA-1_3","GGAGGTATCTAGACCA-1_3","TCGATTTCAATTTCCT-1_3","CTCCAACGTTCTCCAC-1_3","CGCAGGTCAGAAACCG-1_3","TTGAACGTCTGGCCAG-1_3","CCTACGTAGGACTTCT-1_3","TCACACCCAGGTTACT-1_3","AGCTACATCTCGTGAA-1_3","GTAGCTACAACAGAGC-1_3","CGGGTGTAGTTGTAGA-1_3","ACTTCGCGTCTAGATC-1_3","TGGTGATCAATCTCGA-1_3","CGAGGCTAGGACAACC-1_3","GGCGTCACAATGTTGC-1_3","CCGGTGAAGCAATTAG-1_3","CCGTAGGGTAAGTTAG-1_3","AAGTTCGTCCACGTGG-1_3","GTCTACCCAACACGTT-1_3","ATAGGCTCACACGGAA-1_3","GCACATAGTCATCCCT-1_3","AGCTTCCTCGCAGTTA-1_3","CTCCCTCCAGCACAAG-1_3","ACGTTCCCAATCTGCA-1_3","ACCCTTGAGGTTTACC-1_3","TCGGGTGAGATTCGCT-1_3","TACTTACGTACCGGCT-1_3","CTCATCGCATCAACCA-1_3","TCCTCGATCTGCGAGC-1_3","TACCTGCGTGTTTCTT-1_3","AACACACAGGGCAATC-1_3","TCGAACACAAGAGCTG-1_3","CTAACCCAGATAGTGT-1_3","GAAATGAGTCTTGCGG-1_3","TACGGTATCATGCAGT-1_3","GTGGTTATCCGGACTG-1_3","TAGACCAAGGAGACCT-1_3")
#newcl2 <- str_replace(newcl2,"-","-")
#newcl_1 <- c("AGATAGATCTCCCTAG-1_3","ACTATTCGTCTGTGAT-1_3","GGAATCTGTCCTGTCT-1_3","AACTTCTGTAGCCCTG-1_3","TGGAACTCAGGTGTTT-1_3","CCACCATAGTGAGTTA-1_3","ATCCACCTCATTATCC-1_3","AGGTCATCAGAAGTTA-1_3","TACTGCCTCGACTCCT-1_3","TTGGGATCATGGATCT-1_3","GATTGGTGTCCAGGTC-1_3","TTTGACTAGGGTCTTT-1_3","ACTTAGGCAGTAACGG-1_3","GCCTGTTTCATGGCCG-1_3","CTGGTCTGTATGAGCG-1_3","GCTACAACAGTAGAAT-1_3","CAAGGGACATCTGGGC-1_3","GGAAGTGCATGACAGG-1_3","AGGATCTAGCGTATAA-1_3","CATTGCCGTAGCTGCC-1_3","GCCAGGTAGCTTGTGT-1_3","ATTCATCGTGTCTAAC-1_3")
#newcl_2 <- c("GTGTGATGTCATATGC-1_3","TTCTAGTGTAGGAAAG-1_3","CGTAGTACAATACCTG-1_3","TGTCCCAAGTGAATAC-1_3","ATATCCTAGACTCAAA-1_3","ACTTATCAGGTAAGGA-1_3","ATTGTTCCAGCGTGCT-1_3","CCTATCGAGTTGCCCG-1_3","ATGCATGCAGGGTCTC-1_3","TGCAGGCAGGAGGCAG-1_3","ACACAGTTCTGCTTAT-1_3","CCCTGATAGCTCTGTA-1_3","ACGATGTCACAGCCAC-1_3","ACGTACACAGTTGCGC-1_3","GTGTAACCAAGCAATA-1_3","GTAACCACATTCTGTT-1_3","AGTACCACAACGGCTC-1_3","CATTGCCTCATCAGTG-1_3","ATTTCACTCACACCGG-1_3","ATCGGATAGAGGGTAA-1_3","GTGATGTCACCCTCTA-1_3","GAGACTTCATCAGCGC-1_3","CTTCAATTCCTCTCTT-1_3","ACCAACATCATCACTT-1_3","GACCCTTTCTAACACG-1_3","CCTCAACTCTTGCAAG-1_3","GCACGGTCAAAGCGTG-1_3","ATTACTCCAATACCTG-1_3","TCTGGCTGTACTCAAC-1_3","TCAAGACCAATCTCTT-1_3","CGAAGTTAGTGCTCAT-1_3","TGGAGAGAGTATAACG-1_3","GATCATGAGTATCTGC-1_3","TTGACCCTCCTTGGAA-1_3","TCGGTCTTCCACCTCA-1_3","GCCCGAAGTTGTTTGG-1_3","ATCCTATAGACAGTCG-1_3","GACGCTGTCGAGCCAC-1_3","GATCACAGTCACTTCC-1_3","CTCATGCCACACACGC-1_3","TACCTGCCATGCCGCA-1_3","CATACTTAGGTCGTGA-1_3","AGACCCGAGAGCAACC-1_3","TGAACGTAGTTTGAGA-1_3","ATCGTGACAGACCAGA-1_3","CATGCGGGTGTCTTCC-1_3","TCGCACTCAAATAGCA-1_3","GGAGGTATCTAGACCA-1_3","TCGATTTCAATTTCCT-1_3","CTCCAACGTTCTCCAC-1_3","CGCAGGTCAGAAACCG-1_3","TTGAACGTCTGGCCAG-1_3","CCTACGTAGGACTTCT-1_3","TCACACCCAGGTTACT-1_3","AGCTACATCTCGTGAA-1_3","GTAGCTACAACAGAGC-1_3","CGGGTGTAGTTGTAGA-1_3","ACTTCGCGTCTAGATC-1_3","TGGTGATCAATCTCGA-1_3","CGAGGCTAGGACAACC-1_3","GGCGTCACAATGTTGC-1_3","CCGGTGAAGCAATTAG-1_3","CCGTAGGGTAAGTTAG-1_3","AAGTTCGTCCACGTGG-1_3","GTCTACCCAACACGTT-1_3","ATAGGCTCACACGGAA-1_3","GCACATAGTCATCCCT-1_3","AGCTTCCTCGCAGTTA-1_3","CTCCCTCCAGCACAAG-1_3","ACGTTCCCAATCTGCA-1_3","ACCCTTGAGGTTTACC-1_3","TCGGGTGAGATTCGCT-1_3","TACTTACGTACCGGCT-1_3","CTCATCGCATCAACCA-1_3","TCCTCGATCTGCGAGC-1_3","TACCTGCGTGTTTCTT-1_3","AACACACAGGGCAATC-1_3","TCGAACACAAGAGCTG-1_3","CTAACCCAGATAGTGT-1_3","GAAATGAGTCTTGCGG-1_3","TACGGTATCATGCAGT-1_3","TACGGTATCATGCAGT-1_3","GTGGTTATCCGGACTG-1_3","TGTACAGGTTAAGTCC-1_3")
#newcl_1 <- str_replace(newcl_1,"-","-")
#newcl_2 <- str_replace(newcl_2,"-","-")

#Cl1 and Cl3 likely Tb
#Cl2 unknown
Idents(mammal.combined) <- mammal.combined$ID3
Ds2 <- subset(mammal.combined,idents="C2")
Idents(Ds2) <- Ds2$CorrectLabel
#Idents(Ds2,cells=cl1) <- "Cl1"
#Idents(Ds2,cells=cl2) <- "Cl2"
#Idents(Ds2,cells=cl3) <- "Cl3"
#Idents(Ds2,cells=cl4) <- "Cl4" #3 and 4 overlaps?
#Idents(Ds2,cells=newcl_1) <- "newCl1"
#Idents(Ds2,cells=newcl_2) <- "newCl2"
Idents(Ds2,cells=WhichCells(Ds2,idents=c("putExMes"))) <- "Am/EmDisc_d14"
Idents(Ds2,cells=WhichCells(Ds2,idents=c("Am_d14"))) <- "Am/EmDisc_d14"
Idents(Ds2,cells=WhichCells(Ds2,idents=c("putSTB"))) <- "STB_d14"
Idents(Ds2,cells=WhichCells(Ds2,idents=c("Am/EmDisc"))) <- "Am/EmDisc_d14"
Idents(Ds2,cells=WhichCells(Ds2,idents=c("EmDisc_d14"))) <- "Am/EmDisc_d14"


Idents(Ds2,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 2" & TbClustrs$Cl=="Cl1")]) <- "STB_Cl1_d14"
Idents(Ds2,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 2" & TbClustrs$Cl=="Cl2")]) <- "STB_Cl2_d14"
Idents(Ds2,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 2" & TbClustrs$Cl=="Cl3")]) <- "STB_Cl3_d14"
Idents(Ds2,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 2" & TbClustrs$Cl=="Cl4")]) <- "STB_Cl4_d14"
Idents(Ds2,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 2" & TbClustrs$Cl=="Cl5")]) <- "STB_Cl5_d14"
Idents(Ds2,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 2" & TbClustrs$Cl=="Cl6")]) <- "STB_Cl6_d14"
Idents(Ds2,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 2" & TbClustrs$Cl=="Cl7")]) <- "STB_Cl7_d14"
Idents(Ds2,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 2" & TbClustrs$Cl=="Cl8")]) <- "STB_Cl8_d14"
Idents(Ds2,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 2" & TbClustrs$Cl=="Cl9")]) <- "STB_Cl9_d14"
Idents(Ds2,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 2" & TbClustrs$Cl=="Cl10")]) <- "STB_Cl10_d14"
Idents(Ds2,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 2" & TbClustrs$Cl=="Cl11")]) <- "STB_Cl11_d14"
Idents(Ds2,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 2" & TbClustrs$Cl=="Cl12")]) <- "STB_Cl12_d14"
Idents(Ds2,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 2" & TbClustrs$Cl=="Cl13")]) <- "STB_Cl13_d14"
Idents(Ds2,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 2" & TbClustrs$Cl=="Cl14")]) <- "STB_Cl14_d14"
Ds2 <- FindVariableFeatures(Ds2, selection.method = "vst", nfeatures = 3000)
Ds2 <- ScaleData(Ds2, verbose = FALSE)
Ds2 <- RunPCA(Ds2, npcs = 20, verbose = FALSE)
#Ds1 <- RunUMAP(Ds1, reduction = "pca", dims = 1:20)
Ds2 <- FindNeighbors(Ds2, reduction = "pca", dims = 1:20)#
#colind <- integer( length( levels(Idents(Ds2)) )  )
#for (i in 1:length( levels(Idents(Ds2)) ) ) {
#  colind[i] <- which(cType==levels(Idents(Ds2))[i])
#}
#coluse <- BaseCol[colind]
p<-DimPlot(Ds2, pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "POU5F1", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_OCT4.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "WNT6", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_WNT6.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "DIO2", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_DIO2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "HLA-G", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_HLA-G.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds2, features = "CGA", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_CGA.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "JAM3", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_JAM3.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "VIM", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_VIM.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "HAND2", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_HAND2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "SOX15", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_SOX15.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "PDGFA", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_PDGFA.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "GATA4", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_GATA4.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "CER1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_CER1.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds2, features = "BST2", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_BST2.pdf",sep=""),width = 10, height = 8,p)


p<-FeaturePlot(Ds2, features = "NODAL", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_NODAL.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "LEFTY1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_LEFTY1.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "LEFTY2", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_LEFTY2.pdf",sep=""),width = 10, height = 8,p)



p<-FeaturePlot(Ds2, features = "TTR", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_TTR.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds2, features = "APOA1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_APOA1.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "APOB", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_APOB.pdf",sep=""),width = 10, height = 8,p)

#CTB
p<-FeaturePlot(Ds2, features = "TEAD3", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_TEAD3.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "TP63", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_TP63.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "OVOL1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_OVOL1.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "CCKBR", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_CCKBR.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "SIGLEC6", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_SIGLEC6.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "IFI6", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_IFI6.pdf",sep=""),width = 10, height = 8,p)

#STB
p<-FeaturePlot(Ds2, features = "CGB", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_CGB.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "SDC1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_SDC1.pdf",sep=""),width = 10, height = 8,p)

#
p<-FeaturePlot(Ds2, features = "SOX17", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_SOX17.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "HHEX", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_HHEX.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds2, features = "FOXF1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_HHEX.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds2, features = "ENPEP", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_ENPEP.pdf",sep=""),width = 10, height = 8,p)
#CTb-specific markers (TEAD3, TP63, OVOL1)14, STb-specific markers (CGa, CGb, SDC1)
p<-FeaturePlot(Ds2, features = "TACSTD2", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_TACSTD2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "LINC00261", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_LINC00261.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "FOXF1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_FOXF1.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds2, features = "GABRP", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_GABRP.pdf",sep=""),width = 10, height = 8,p)
p<-DimPlot(Ds2, pt.size = 4, reduction = "pca", split.by = "Genotype2",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_GT.pdf",sep=""),width = 20, height = 8,p)


MiniList <- c("DNMT3A","DPPA3","OTX2","POU5F1","PDGFA","PODXL","SFRP1","SOX15","GABRP","VTCN1","WNT6","BAMBI","NOTUM","JAM3","TACSTD2","ENPEP","GATA2","GATA3","TFAP2C","TFAP2A","TEAD3","TP63","OVOL1","CCKBR","CGA","CGB5","CGB8","SDC1","HLA-G","SOX17","GATA4","CER1","LEFTY1","LEFTY2","APOA1","APOB","LINC00261","HHEX","FOXF1","PDGFRA","HAND2","BST2","GATA6","HGF","PRDM1","ASCL2","SNAI1")
X <- GetAssayData(Ds2)
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(Ds2) )), Batch =  factor(( Idents(Ds2)))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch2_markerplot",".pdf",sep="") ,width=10,height=7)

MiniList <- c("SOX17","TFAP2A","TFAP2C","PRDM1","POU5F1","NANOS3")
X <- GetAssayData(Ds2)
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(Ds2) )), Batch =  factor(( Idents(Ds2)))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch2_markerplot2",".pdf",sep="") ,width=50,height=3)



#p<-DimPlot(Ds2, cols = coluse, pt.size = 4, reduction = "pca", split.by = "Genotype2", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic.pdf",sep=""),width = 20, height = 8,p)




Ds3 <- subset(mammal.combined,idents=c("C3","C4") )
Idents(Ds3) <- Ds3$CorrectLabel
Idents(Ds3,cells=WhichCells(Ds3,idents=c("EmDisc_d14"))) <- "Am/EmDisc_d14"
Idents(Ds3,cells=WhichCells(Ds3,idents=c("Am/EmDisc"))) <- "Am/EmDisc_d14"
Idents(Ds3,cells=WhichCells(Ds3,idents=c("ExMes_d14"))) <- "Am/EmDisc_d14"
#Idents(Ds3,cells=WhichCells(Ds3,idents=c("EmDisc_d14"))) <- "Am/EmDisc_d14"
Idents(Ds3,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 3" & TbClustrs$Cl=="Cl1")]) <- "STB_Cl1_d14"
Idents(Ds3,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 3" & TbClustrs$Cl=="Cl2")]) <- "STB_Cl2_d14"
Idents(Ds3,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 3" & TbClustrs$Cl=="Cl3")]) <- "STB_Cl3_d14"
Idents(Ds3,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 3" & TbClustrs$Cl=="Cl4")]) <- "STB_Cl4_d14"
Ds3 <- FindVariableFeatures(Ds3, selection.method = "vst", nfeatures = 3000)
Ds3 <- ScaleData(Ds3, verbose = FALSE)
Ds3 <- RunPCA(Ds3, npcs = 20, verbose = FALSE)
#Ds1 <- RunUMAP(Ds1, reduction = "pca", dims = 1:20)
Ds3 <- FindNeighbors(Ds3, reduction = "pca", dims = 1:20)#
#colind <- integer( length( levels(Idents(Ds3)) )  )
#for (i in 1:length( levels(Idents(Ds3)) ) ) {
#  colind[i] <- which(cType==levels(Idents(Ds3))[i])
#}
#coluse <- BaseCol[colind]
p<-DimPlot(Ds3,  pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "POU5F1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_OCT4.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "WNT6", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_WNT6.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "DIO2", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_DIO2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "HLA-G", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_HLA-G.pdf",sep=""),width = 10, height = 8,p)


p<-FeaturePlot(Ds3, features = "CGA", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_CGA.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "JAM3", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_JAM3.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "VIM", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_VIM.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "HAND2", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_HAND2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "SOX15", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_SOX15.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "PDGFA", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_PDGFA.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "GATA4", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_GATA4.pdf",sep=""),width = 10, height = 8,p)


#CTB
p<-FeaturePlot(Ds3, features = "TEAD3", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_TEAD3.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "TP63", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_TP63.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "OVOL1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_OVOL1.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "CCKBR", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_CCKBR.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "SIGLEC6", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_SIGLEC6.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "IFI6", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_IFI6.pdf",sep=""),width = 10, height = 8,p)

#STB
p<-FeaturePlot(Ds3, features = "CGB", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_CGB.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "SDC1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_SDC1.pdf",sep=""),width = 10, height = 8,p)

#
p<-FeaturePlot(Ds3, features = "SOX17", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_SOX17.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "HHEX", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_HHEX.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds3, features = "FOXF1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_HHEX.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds3, features = "ENPEP", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_ENPEP.pdf",sep=""),width = 10, height = 8,p)
#CTb-specific markers (TEAD3, TP63, OVOL1)14, STb-specific markers (CGa, CGb, SDC1)
p<-FeaturePlot(Ds3, features = "TACSTD2", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_TACSTD2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "LINC00261", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_LINC00261.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "FOXF1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C2Embryonic_FOXF1.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds3, features = "GABRP", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_GABRP.pdf",sep=""),width = 10, height = 8,p)
p<-DimPlot(Ds3, pt.size = 4, reduction = "pca", split.by = "Genotype2",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C3Embryonic_GT.pdf",sep=""),width = 20, height = 8,p)


MiniList <- c("DNMT3A","DPPA3","OTX2","POU5F1","PDGFA","PODXL","SFRP1","SOX15","GABRP","VTCN1","WNT6","BAMBI","NOTUM","JAM3","TACSTD2","ENPEP","GATA2","GATA3","TFAP2C","TFAP2A","TEAD3","TP63","OVOL1","CCKBR","CGA","CGB5","CGB8","SDC1","HLA-G","SOX17","GATA4","CER1","LEFTY1","LEFTY2","APOA1","APOB","LINC00261","HHEX","FOXF1","PDGFRA","HAND2","BST2","GATA6","HGF","PRDM1","ASCL2","SNAI1")
X <- GetAssayData(Ds3)
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(Ds3) )), Batch =  factor(( Idents(Ds3)))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch3_markerplot",".pdf",sep="") ,width=10,height=7)

MiniList <- c("SOX17","TFAP2A","TFAP2C","PRDM1","POU5F1","NANOS3")
X <- GetAssayData(Ds3)
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(Ds3) )), Batch =  factor(( Idents(Ds3)))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch3_markerplot2",".pdf",sep="") ,width=20,height=3)


#cl1 <- c("AACAAGATCTCCAAGA-1_6","ACCTACCAGTTAGAAC-1_6","AGACAGGCAGTGGTGA-1_6","AGGGAGTTCAACCCGG-1_6","ATGCCTCAGTCATCGT-1_6","CGGAATTGTGGGTCAA-1_6","CTGAGCGCAGCTTTCC-1_6","CTGTATTCAAGGCGTA-1_6","GACGCTGGTTCATCGA-1_6","GATTGGTTCGCTAGCG-1_6","GCGGAAATCATCCTGC-1_6","GGAGAACGTGTACGCC-1_6","GGTGTCGCAGCCGTTG-1_6","GTGGTTATCGCCAATA-1_6","TAGTGCACAACTTGCA-1_6","TCGCTCATCCACGGAC-1_6","TCGTAGACACTGATTG-1_6","TCTGCCACATGCCGAC-1_6","TGGATGTGTTAACAGA-1_6","TGGGCTGGTCATCACA-1_6","TTCATTGCACGAAAGC-1_6","TTTCATGTCCGCAACG-1_6","TTTGACTAGGACATCG-1_6")
#cl2 <- c("AAAGAACAGTAAACGT-1_6","AAAGGGCGTCACGACC-1_6","AAGCCATCACACCAGC-1_6","AATCACGAGCGGTAAC-1_6","ACAAAGAAGTTGCATC-1_6","ACAAGCTTCTCATTAC-1_6","ACACGCGAGGAGGGTG-1_6","ACCCTTGGTAGGTTTC-1_6","ACTTCCGTCTGTACAG-1_6","AGAAGTATCCTAGCTC-1_6","AGACACTCATTCTGTT-1_6","AGATAGAAGTCGAAGC-1_6","AGGACGAGTGGCGCTT-1_6","AGTAGCTTCTGGTCAA-1_6","ATACCGAAGGTGCTAG-1_6","ATACCTTCAGCTATTG-1_6","ATCAGGTTCCCTAGGG-1_6","ATCCGTCGTACAAGTA-1_6","ATTCATCAGACGATAT-1_6","CAACAACGTAATGCGG-1_6","CAACGGCTCAAGAAAC-1_6","CAAGACTCAAGATCCT-1_6","CACACAACAGTCTGGC-1_6","CACATGAAGACCACGA-1_6","CACGAATCAATCTGCA-1_6","CAGAGCCCAAGAGTTA-1_6","CAGCCAGCAGATTTCG-1_6","CAGGCCAAGTTTGCTG-1_6","CATAAGCTCGTGCACG-1_6","CATACTTGTCCGGTGT-1_6","CATCGTCTCGTTGTAG-1_6","CCGATCTCACCAAAGG-1_6","CCGGACATCCAAGCAT-1_6","CCGTAGGGTCAAAGAT-1_6","CCTCACACATGCCGCA-1_6","CGCCATTAGTACGAGC-1_6","CGTCCATTCGCACGGT-1_6","CGTGCTTAGTTCCTGA-1_6","CTCGAGGAGCACCTGC-1_6","CTCTCGAGTACTAGCT-1_6","CTGCTCACATGGAATA-1_6","CTGGACGCATGCCGAC-1_6","CTGGACGTCTACCCAC-1_6","CTGGTCTAGTATGTAG-1_6","CTTGATTTCTGGCCGA-1_6","CTTTCGGTCACTAGCA-1_6","GAAGCGACAGATCACT-1_6","GACACGCGTCTGTAGT-1_6","GAGACTTAGCAAATGT-1_6","GAGGCCTAGGATTCCT-1_6","GATCCCTCATTGTAGC-1_6","GATTCGACAAGTGGAC-1_6","GCACGGTTCTTACTGT-1_6","GCATCTCTCGGAATGG-1_6","GCCAACGCATGACGAG-1_6","GCGAGAATCCCGAACG-1_6","GCTGGGTAGGAGACCT-1_6","GGAAGTGAGGCCTTGC-1_6","GGATGTTCACCAGCTG-1_6","GGGACCTCACGGTGTC-1_6","GGGACTCAGCCGGATA-1_6","GGGCTACTCGCAGTTA-1_6","GGGTTATCAAGGACAC-1_6","GGTCTGGTCCACCCTA-1_6","GTACAACGTGAGCCAA-1_6","GTAGCTACAGCTTTCC-1_6","GTAGGTTTCGTACACA-1_6","GTCAAACGTTAGGAGC-1_6","GTCAGCGGTGTTAGCT-1_6","GTCTCACCACTGGATT-1_6","GTCTTTAGTCTTGCTC-1_6","GTGCTGGCACCGGTCA-1_6","GTTGTGAAGAGCAGCT-1_6","TAACTTCTCGCGATCG-1_6","TAAGCACAGAATCGAT-1_6","TAGCACACACCTCGTT-1_6","TAGGAGGTCTATTCGT-1_6","TAGGTTGCAACTCCAA-1_6","TATCTGTCAAGAGAGA-1_6","TATTGCTAGGACGCTA-1_6","TCATGTTGTCGCGGTT-1_6","TCATTCACAACACGTT-1_6","TCCCATGAGTAGCATA-1_6","TCGACCTCATTCAGCA-1_6","TCGCTTGAGCCGATCC-1_6","TCTACATTCGAGAAAT-1_6","TCTGCCAAGTATTAGG-1_6","TCTTTGAAGCAATAGT-1_6","TGAATCGAGTCGTCTA-1_6","TGAGGTTGTTTCCATT-1_6","TGAGGTTTCTACGCAA-1_6","TGATGGTCACTTGACA-1_6","TGCAGGCCACGAAGAC-1_6","TGCAGTACACACACGC-1_6","TGGGCTGCAACCAGAG-1_6","TGGTACATCGCACTCT-1_6","TGTAAGCAGGCTAAAT-1_6","TGTCCACGTAGCGATG-1_6","TGTGCGGCACTAGGTT-1_6","TGTTCCGGTCGAGTGA-1_6","TTACGTTAGATGGTAT-1_6","TTACGTTTCAAGTTGC-1_6","TTAGTCTCAAATGAGT-1_6","TTGCGTCTCTCGCCTA-1_6","TTGGATGAGTTTGTCG-1_6","TTGTTCAGTCGTTATG-1_6","TTTACTGTCACTGTTT-1_6","TATCTGTCAAGAGAGA-1_6")
#cl3 <- c("AAAGTCCAGGGTGAAA-1_6","AACACACTCTCTCTTC-1_6","ACCTGAATCACCGGTG-1_6","AGATAGAAGTCGAAGC-1_6","AGATCGTGTGGTAATA-1_6","AGCATCACAGCACACC-1_6","ATCCGTCGTACAAGTA-1_6","CAACAACGTAATGCGG-1_6","CAACGGCTCAAGAAAC-1_6","CACACAACAGTCTGGC-1_6","CATACTTGTCCGGTGT-1_6","CATTGAGCATAGTCGT-1_6","CCCATTGCACTTGAAC-1_6","CGCGTGACAGCGATTT-1_6","CGTGAATGTTGTAAAG-1_6","CTACGGGTCCTCTAAT-1_6","GACCAATCATTGACAC-1_6","GACTCTCTCCCAGGAC-1_6","GATGACTGTTGCGAAG-1_6","GGTTGTAGTGCCGAAA-1_6","GTCAAACGTTAGGAGC-1_6","GTCTTTAGTCTTGCTC-1_6","GTGCTGGCACCGGTCA-1_6","TAAGCACCAATGGCAG-1_6","TCCCATGTCCCTGGTT-1_6","TCTACATTCGAGAAAT-1_6","TCTACCGCAGTATACC-1_6","TCTTTGAAGCAATAGT-1_6","TGCAGTACACACACGC-1_6","TGCGATACAAAGGCAC-1_6","TGGAGGAGTGTCTTGA-1_6","TTACGTTGTTCCCACT-1_6","TTGGATGAGTTTGTCG-1_6","TTTATGCTCTATGCCC-1_6")


Ds6 <- subset(mammal.combined,idents=c("C6") )
Idents(Ds6) <- Ds6$CorrectLabel

Idents(Ds6,cells=WhichCells(Ds6,idents=c("EmDisc_d14"))) <- "Am/EmDisc_d14"
Idents(Ds6,cells=WhichCells(Ds6,idents=c("Am/EmDisc"))) <- "Am/EmDisc_d14"
Idents(Ds6,cells=WhichCells(Ds6,idents=c("ExMes_d14"))) <- "Am/EmDisc_d14"
#Idents(Ds6,cells=WhichCells(Ds6,idents=c("EmDisc_d14"))) <- "Am/EmDisc_d14"
Idents(Ds6,cells=WhichCells(Ds6,idents=c("Am_d14"))) <- "Am/EmDisc_d14"
Idents(Ds6,cells=WhichCells(Ds6,idents=c("putExMes"))) <- "Am/EmDisc_d14"
Idents(Ds6,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 6" & TbClustrs$Cl=="Cl1")]) <- "STB_Cl1_d14"
Idents(Ds6,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 6" & TbClustrs$Cl=="Cl2")]) <- "STB_Cl2_d14"
Idents(Ds6,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 6" & TbClustrs$Cl=="Cl3")]) <- "STB_Cl3_d14"



#Idents(Ds6,cells=cl1) <- "Cl1"
#Idents(Ds6,cells=cl2) <- "Cl2"
#Idents(Ds6,cells=cl3) <- "Cl3"
#Cl1 3 embryonic
#Cl2 unknown
Ds6 <- FindVariableFeatures(Ds6, selection.method = "vst", nfeatures = 3000)
Ds6 <- ScaleData(Ds6, verbose = FALSE)
Ds6 <- RunPCA(Ds6, npcs = 20, verbose = FALSE)
#Ds1 <- RunUMAP(Ds1, reduction = "pca", dims = 1:20)
Ds6 <- FindNeighbors(Ds6, reduction = "pca", dims = 1:20)#
#colind <- integer( length( levels(Idents(Ds6)) )  )
#for (i in 1:length( levels(Idents(Ds6)) ) ) {
#  colind[i] <- which(cType==levels(Idents(Ds6))[i])
#}
#coluse <- BaseCol[colind]
p<-DimPlot(Ds6,   pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "POU5F1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_OCT4.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "WNT6", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_WNT6.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "DIO2", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_DIO2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "HLA-G", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_HLA-G.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds6, features = "CGA", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_CGA.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "JAM3", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_JAM3.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "VIM", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_VIM.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "HAND2", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_HAND2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "SOX15", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_SOX15.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "PDGFA", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_PDGFA.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "GATA4", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_GATA4.pdf",sep=""),width = 10, height = 8,p)


p<-FeaturePlot(Ds6, features = "BST2", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_BST2.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds6, features = "HAND2", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_HAND2.pdf",sep=""),width = 10, height = 8,p)


#CTB
p<-FeaturePlot(Ds6, features = "TEAD3", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_TEAD3.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "TP63", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_TP63.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "OVOL1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_OVOL1.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "CCKBR", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_CCKBR.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "SIGLEC6", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_SIGLEC6.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "IFI6", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_IFI6.pdf",sep=""),width = 10, height = 8,p)

#STB
p<-FeaturePlot(Ds6, features = "CGB", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_CGB.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "SDC1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_SDC1.pdf",sep=""),width = 10, height = 8,p)

#
p<-FeaturePlot(Ds6, features = "SOX17", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_SOX17.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "HHEX", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_HHEX.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds6, features = "FOXF1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_HHEX.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds6, features = "ENPEP", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_ENPEP.pdf",sep=""),width = 10, height = 8,p)
#CTb-specific markers (TEAD3, TP63, OVOL1)14, STb-specific markers (CGa, CGb, SDC1)
p<-FeaturePlot(Ds6, features = "TACSTD2", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_TACSTD2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "LINC00261", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_LINC00261.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "FOXF1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_FOXF1.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds6, features = "GABRP", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_GABRP.pdf",sep=""),width = 10, height = 8,p)
p<-DimPlot(Ds6, pt.size = 4, reduction = "pca", split.by = "Genotype2",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C6Embryonic_GT.pdf",sep=""),width = 20, height = 8,p)


MiniList <- c("DNMT3A","DPPA3","OTX2","POU5F1","PDGFA","PODXL","SFRP1","SOX15","GABRP","VTCN1","WNT6","BAMBI","NOTUM","JAM3","TACSTD2","ENPEP","GATA2","GATA3","TFAP2C","TFAP2A","TEAD3","TP63","OVOL1","CCKBR","CGA","CGB5","CGB8","SDC1","HLA-G","SOX17","GATA4","CER1","LEFTY1","LEFTY2","APOA1","APOB","LINC00261","HHEX","FOXF1","PDGFRA","HAND2","BST2","GATA6","HGF","PRDM1","ASCL2","SNAI1")
X <- GetAssayData(Ds6)
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(Ds6) )), Batch =  factor(( Idents(Ds6)))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch6_markerplot",".pdf",sep="") ,width=10,height=7)

MiniList <- c("SOX17","TFAP2A","TFAP2C","PRDM1","POU5F1","NANOS3")
X <- GetAssayData(Ds6)
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(Ds6) )), Batch =  factor(( Idents(Ds6)))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch6_markerplot2",".pdf",sep="") ,width=20,height=3)



#cl1 <- c("ACATCGACAGAGGGTT-1_7","ATATCCTAGAGAGCAA-1_7","ATCCCTGAGTGCACCC-1_7","CCTAAGATCCTTCGAC-1_7","GAGGGTACACTCATAG-1_7","GTAGTACCAGCGTTTA-1_7","TGCAGGCGTAGCGAGT-1_7","TTAGGGTGTACCATAC-1_7","AGGGTGATCCGAGATT-1_7")
#cl2 <- c("AACCTGATCCATCCGT-1_7","AGTACTGCATAAGATG-1_7","ATCGTCCGTTGTTTGG-1_7","ATGGAGGCACTTGGCG-1_7","CCTAAGACACTGCGAC-1_7","CCTCACAGTACCTATG-1_7","GACGTTATCCACCTCA-1_7","GAGTTTGAGCATGCAG-1_7","GATTCTTGTATCACGT-1_7","GCAGTTACATTAAGCC-1_7","GGGTCACAGCTACGTT-1_7","GGTTGTAAGGATACGC-1_7","GTAGTACAGCTCGTGC-1_7","GTCACGGGTGTTCCAA-1_7","GTCGTAAGTACGAAAT-1_7","GTTGAACTCCTCAGAA-1_7","TAGGTACCATTCAGCA-1_7","TGTCCTGCACTCTAGA-1_7","TTGTTGTGTAAGCGGT-1_7","TTTACGTCATGATGCT-1_7","TTTAGTCTCATCCCGT-1_7","CCGCAAGTCTGCTAGA-1_7")
#Cl1 tb?
Ds7 <- subset(mammal.combined,idents=c("C7") )
Idents(Ds7) <- Ds7$CorrectLabel

#Idents(Ds7,cells=WhichCells(Ds7,idents=c("EmDisc_d14"))) <- "Am/EmDisc_d14"
Idents(Ds7,cells=WhichCells(Ds7,idents=c("ExMes_d14"))) <- "Am/EmDisc_d14"
Idents(Ds7,cells=WhichCells(Ds7,idents=c("EmDisc_d14"))) <- "Am/EmDisc_d14"
Idents(Ds7,cells=WhichCells(Ds7,idents=c("Am_d14"))) <- "Am/EmDisc_d14"
#Idents(Ds7,cells=WhichCells(Ds7,idents=c("putExMes"))) <- "Am/EmDisc_d14"
Idents(Ds7,cells=WhichCells(Ds7,idents=c("putSTB"))) <- "STB_d14"

Idents(Ds7,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 7" & TbClustrs$Cl=="Cl1")]) <- "STB_Cl1_d14"
Idents(Ds7,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 7" & TbClustrs$Cl=="Cl2")]) <- "STB_Cl2_d14"
Idents(Ds7,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 7" & TbClustrs$Cl=="Cl3")]) <- "STB_Cl3_d14"


#Idents(Ds7,cells=cl1) <- "Cl1"
#Idents(Ds7,cells=cl2) <- "Cl2"
Ds7 <- FindVariableFeatures(Ds7, selection.method = "vst", nfeatures = 3000)
Ds7 <- ScaleData(Ds7, verbose = FALSE)
Ds7 <- RunPCA(Ds7, npcs = 20, verbose = FALSE)
#Ds1 <- RunUMAP(Ds1, reduction = "pca", dims = 1:20)
Ds7 <- FindNeighbors(Ds7, reduction = "pca", dims = 1:20)#
#colind <- integer( length( levels(Idents(Ds7)) )  )
#for (i in 1:length( levels(Idents(Ds7)) ) ) {
#  colind[i] <- which(cType==levels(Idents(Ds7))[i])
#}
#coluse <- BaseCol[colind]
p<-DimPlot(Ds7,   pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "POU5F1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_OCT4.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "WNT6", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_WNT6.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "DIO2", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_DIO2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "HLA-G", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_HLA-G.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds7, features = "CGA", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_CGA.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "JAM3", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_JAM3.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "VIM", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_VIM.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "HAND2", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_HAND2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "SOX15", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_SOX15.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "PDGFA", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_PDGFA.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "GATA4", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_GATA4.pdf",sep=""),width = 10, height = 8,p)



#CTB
p<-FeaturePlot(Ds7, features = "TEAD3", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_TEAD3.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "TP63", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_TP63.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "OVOL1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_OVOL1.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "CCKBR", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_CCKBR.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "SIGLEC6", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_SIGLEC6.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "IFI6", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_IFI6.pdf",sep=""),width = 10, height = 8,p)

#STB
p<-FeaturePlot(Ds7, features = "CGB", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_CGB.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "SDC1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_SDC1.pdf",sep=""),width = 10, height = 8,p)

#
p<-FeaturePlot(Ds7, features = "SOX17", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_SOX17.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "HHEX", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_HHEX.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds7, features = "FOXF1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_HHEX.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds7, features = "ENPEP", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_ENPEP.pdf",sep=""),width = 10, height = 8,p)
#CTb-specific markers (TEAD3, TP63, OVOL1)14, STb-specific markers (CGa, CGb, SDC1)
p<-FeaturePlot(Ds7, features = "TACSTD2", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_TACSTD2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "LINC00261", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_LINC00261.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "FOXF1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_FOXF1.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds7, features = "GABRP", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_GABRP.pdf",sep=""),width = 10, height = 8,p)
p<-DimPlot(Ds7, pt.size = 4, reduction = "pca", split.by = "Genotype2",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C7Embryonic_GT.pdf",sep=""),width = 20, height = 8,p)


MiniList <- c("DNMT3A","DPPA3","OTX2","POU5F1","PDGFA","PODXL","SFRP1","SOX15","GABRP","VTCN1","WNT6","BAMBI","NOTUM","JAM3","TACSTD2","ENPEP","GATA2","GATA3","TFAP2C","TFAP2A","TEAD3","TP63","OVOL1","CCKBR","CGA","CGB5","CGB8","SDC1","HLA-G","SOX17","GATA4","CER1","LEFTY1","LEFTY2","APOA1","APOB","LINC00261","HHEX","FOXF1","PDGFRA","HAND2","BST2","GATA6","HGF","PRDM1","ASCL2","SNAI1")
X <- GetAssayData(Ds7)
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(Ds7) )), Batch =  factor(( Idents(Ds7)))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch7_markerplot",".pdf",sep="") ,width=10,height=7)

MiniList <- c("SOX17","TFAP2A","TFAP2C","PRDM1","POU5F1","NANOS3")
X <- GetAssayData(Ds7)
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(Ds7) )), Batch =  factor(( Idents(Ds7)))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch7_markerplot2",".pdf",sep="") ,width=10,height=3)




Ds8 <- subset(mammal.combined,idents=c("C8") )
Idents(Ds8) <- Ds8$CorrectLabel
Idents(Ds8,cells=WhichCells(Ds8,idents=c("ExMes_d14"))) <- "Am/EmDisc_d14"
Idents(Ds8,cells=WhichCells(Ds8,idents=c("EmDisc_d14"))) <- "Am/EmDisc_d14"

Idents(Ds8,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 8" & TbClustrs$Cl=="Cl1")]) <- "STB_Cl1_d14"
Idents(Ds8,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 8" & TbClustrs$Cl=="Cl2")]) <- "STB_Cl2_d14"
Idents(Ds8,cells= TbClustrs$Barcode[which(TbClustrs$Batch=="Batch 8" & TbClustrs$Cl=="Cl3")]) <- "STB_Cl3_d14"


Ds8 <- FindVariableFeatures(Ds8, selection.method = "vst", nfeatures = 3000)
Ds8 <- ScaleData(Ds8, verbose = FALSE)
Ds8 <- RunPCA(Ds8, npcs = 20, verbose = FALSE)
#Ds1 <- RunUMAP(Ds1, reduction = "pca", dims = 1:20)
Ds8 <- FindNeighbors(Ds8, reduction = "pca", dims = 1:20)#
#colind <- integer( length( levels(Idents(Ds8)) )  )
#for (i in 1:length( levels(Idents(Ds8)) ) ) {
#  colind[i] <- which(cType==levels(Idents(Ds8))[i])
#}
#coluse <- BaseCol[colind]
p<-DimPlot(Ds8,  pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "POU5F1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_OCT4.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "WNT6", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_WNT6.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "DIO2", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_DIO2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "HLA-G", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_HLA-G.pdf",sep=""),width = 10, height = 8,p)


p<-FeaturePlot(Ds8, features = "CGA", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_CGA.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "JAM3", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_JAM3.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "VIM", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_VIM.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "HAND2", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_HAND2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "SOX15", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_SOX15.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "PDGFA", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_PDGFA.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "GATA4", pt.size = 4, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_GATA4.pdf",sep=""),width = 10, height = 8,p)




#CTB
p<-FeaturePlot(Ds8, features = "TEAD3", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_TEAD3.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "TP63", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_TP63.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "OVOL1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_OVOL1.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "CCKBR", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_CCKBR.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "SIGLEC6", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_SIGLEC6.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "IFI6", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_IFI6.pdf",sep=""),width = 10, height = 8,p)

#STB
p<-FeaturePlot(Ds8, features = "CGB", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_CGB.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "SDC1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_SDC1.pdf",sep=""),width = 10, height = 8,p)

#
p<-FeaturePlot(Ds8, features = "SOX17", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_SOX17.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "HHEX", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_HHEX.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds8, features = "FOXF1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_HHEX.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds8, features = "ENPEP", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_ENPEP.pdf",sep=""),width = 10, height = 8,p)
#CTb-specific markers (TEAD3, TP63, OVOL1)14, STb-specific markers (CGa, CGb, SDC1)
p<-FeaturePlot(Ds8, features = "TACSTD2", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_TACSTD2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "LINC00261", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_LINC00261.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "FOXF1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_FOXF1.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds8, features = "GABRP", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_GABRP.pdf",sep=""),width = 10, height = 8,p)
p<-DimPlot(Ds8, pt.size = 4, reduction = "pca", split.by = "Genotype2",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C8Embryonic_GT.pdf",sep=""),width = 20, height = 8,p)


MiniList <- c("DNMT3A","DPPA3","OTX2","POU5F1","PDGFA","PODXL","SFRP1","SOX15","GABRP","VTCN1","WNT6","BAMBI","NOTUM","JAM3","TACSTD2","ENPEP","GATA2","GATA3","TFAP2C","TFAP2A","TEAD3","TP63","OVOL1","CCKBR","CGA","CGB5","CGB8","SDC1","HLA-G","SOX17","GATA4","CER1","LEFTY1","LEFTY2","APOA1","APOB","LINC00261","HHEX","FOXF1","PDGFRA","HAND2","BST2","GATA6","HGF","PRDM1","ASCL2","SNAI1")
X <- GetAssayData(Ds8)
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(Ds8) )), Batch =  factor(( Idents(Ds8)))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch78_markerplot",".pdf",sep="") ,width=10,height=7)

MiniList <- c("SOX17","TFAP2A","TFAP2C","PRDM1")
X <- GetAssayData(Ds8)
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList ,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(Ds8) )), Batch =  factor(( Idents(Ds8)))  )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch8_markerplot2",".pdf",sep="") ,width=10,height=3)


#Panel 1

#list1 <- c("JAM3","GATA2","GATA3","TFAP2A","TFAP2C","CGB8","CGB5","CGA","ERVW-1","PRDM6","TBX3","DIO2","NOTUM","HLA-G","ASCL2","SNAI1",
#           "HGF","SNAI2","HAND1","HAND2","PDGFRA","BST2","GATA6","GATA4","CER1","NODAL","LEFTY1","LEFTY2","APOA1","WNT6","VTCN1","BAMBI","AKAP12","PODXL","SFRP1","SOX15","PDGFA","POU5F1")

#TbList1 <- intersect(intersect(rownames(D),TbMarkers$Gene[which(as.numeric(TbMarkers$'HumanSS2_CTB_vs_STB_adjpval')<0.0005 & as.numeric(TbMarkers$'HumanSS2_CTB_vs_STB_logFC')>1.609438)]),DE1)
#TbList2 <- intersect(intersect(rownames(D),TbMarkers$Gene[which(as.numeric(TbMarkers$'HumanSS2_CTB_vs_EVT_adjpval')<0.0005 & as.numeric(TbMarkers$'HumanSS2_CTB_vs_EVT_logFC')>1.609438)]),DE1)
#TbList3 <- intersect(intersect(rownames(D),TbMarkers$Gene[which(as.numeric(TbMarkers$'HumanSS2_STB_vs_EVT_adjpval')<0.0005 & as.numeric(TbMarkers$'HumanSS2_STB_vs_EVT_logFC')>1.609438)]),DE1)

#TbList4 <- intersect(intersect(rownames(D),TbMarkers$Gene[which(as.numeric(TbMarkers$'HumanSS2_CTB_vs_STB_adjpval')<0.0005 & as.numeric(TbMarkers$'HumanSS2_CTB_vs_STB_logFC')< -1.609438)]),DE1)
#TbList5 <- intersect(intersect(rownames(D),TbMarkers$Gene[which(as.numeric(TbMarkers$'HumanSS2_CTB_vs_EVT_adjpval')<0.0005 & as.numeric(TbMarkers$'HumanSS2_CTB_vs_EVT_logFC')< -1.609438)]),DE1)
#TbList6 <- intersect(intersect(rownames(D),TbMarkers$Gene[which(as.numeric(TbMarkers$'HumanSS2_STB_vs_EVT_adjpval')<0.0005 & as.numeric(TbMarkers$'HumanSS2_STB_vs_EVT_logFC')< -1.609438)]),DE1)


#Here are the heatmaps we want
MiniList1 <- c("JARID2","CASC15","LINC01194","MT1G","MT1H","NLGN4X","RIMS2","CD9","DPPA4","CLU","NAV3","CDH1","CLDN4","PRLR","DMD","CLDN6","EZR","LGR5", "LHX1","OPHN1","SLC3A2","SLC2A3","TANC2","TCIM","COL5A2","COL6A3","TNC","SPARC","COL1A2","COL1A1","COL3A1","DCN","LUM","COL6A1","COL6A2","DPP4","CST3","DNAJC15","DNMT3A","DPPA3","OTX2","POU5F1","PDGFA","PODXL","SFRP1","SOX15","GABRP","VTCN1","WNT6","BAMBI","TFAP2C","TFAP2A","GATA3","CGA","SOX17","GATA4","CER1","LEFTY1","LEFTY2","APOA1","APOB","LINC00261","HHEX","FOXF1","PDGFRA","HAND2","BST2","GATA6","HGF","PRDM1","SNAI1","VIM","FN1")
MiniList2 <- c("DNMT3A","DPPA3","OTX2","NOTUM","JAM3","TACSTD2","ENPEP","GATA2","GATA3","TFAP2C","TFAP2A","TEAD3","TP63","OVOL1","CCKBR","CGA","CGB5","CGB8","SDC1","HLA-C","HLA-G","ASCL2","SNAI1","VIM","FN1","ELF5","MKI67","CDH1","FSTL1","FSTL3","KRT19","KRT7","SPARC","LAMC1","CTNNAL1")

putPGC <- c("CAAGAGGTCCACCCTA-1_2","TTATTGCAGCACAAAT-1_2","TCGCACTAGCCAAGGT-1_2","GTAGATCAGCTGGCCT-1_2","ATCGATGGTAGCGTAG-1_2","GTGTTCCAGGGCAGGA-1_2","CAAGACTGTATTGCCA-1_2","TTCTTCCAGTCATCGT-1_2","CCTGCATCAAGTATCC-1_3","ACACCAACAACCGCTG-1_3","CCTCAACTCTTGCAAG-1_3","TAACTTCGTACAAAGT-1_3","AGCGCCACACGGCACT-1_3","TCATGTTGTCGACGCT-1_3","CAGTGCGAGTAGCCAG-1_3","ATGAGTCAGTCAATCC-1_3","TTACAGGGTCTTCTAT-1_3","CATCGTCTCGTTGTAG-1_6","AAAGAACAGTAAACGT-1_6","ACACGCGAGGAGGGTG-1_6","GCATCTCTCGGAATGG-1_6","GGAAGTGAGGCCTTGC-1_6")
Idents(mammal.combined) <- mammal.combined$Cells
Idents(mammal.combined,cells=putPGC) <- "putPGC"
mammal.combined$Cells <- Idents(mammal.combined)

subs1 <- which(mammal.combined$ID3%in%c("C1") & mammal.combined$Cells%in%c("CTB_d14","STB_d14","EVT_d14","putSTB") )
subs2 <- which(mammal.combined$ID3%in%c("C1") & mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
#cutree_cols = 20, 

EVT
CCGGGTAGTACCTTCC−1_2

CTB 
TCTGTCGCAGTCTACA−1_2
CTCAATTTCCCTCATG−1_2
TACCCACAGCCACTCG−1_2
TCCCACAAGATCACCT−1_2
CGGACACCACACAGCC−1_2
ATGACCACACTACTTT−1_2
TCTCCGAGTAGGCTCC−1_2
TCGGGTGGTTTGAAAG−1_2

CTB to STB
GGGTATTTCATTTGTC−1_2
ACTCTCGCAGTCGAGA−1_2
CTTTCAATCAGCCCAG−1_2
ATCGATGCAATTGCAC−1_2
GTAACACGTAACAGTA−1_2
GATCAGTCAAATCCCA−1_2
TCCACCAAGCAATTCC−1_2
AATGCCAGTAACAGTA−1_2
CATTGTTAGGTGCAGT−1_2
GTGTCCTTCGAACACT−1_2
TCGGGTGGTTTGAAAG−1_2
AATCGACTCCAGCAAT−1_2
ACATCGATCATCTCTA−1_2
CCCTCAAGTCCGAAAG−1_2
TTTACTGCAATCTAGC−1_2
GTCTACCGTTGCTAGT−1_2
ACCAACATCCTACGAA−1_2
GTGAGCCAGCACGTCC−1_2
ACCTACCAGGCAGCTA−1_2
CCCTAACGTAGACAGC−1_2
GCTGGGTCACATGTTG−1_2
TAGCACACATCGATAC−1_2
CTTCAATAGTAGGGTC−1_2
CACTTCGCAGCTGCCA−1_2
GTTACAGTCTTTCCGG−1_2
CCACACTAGTTGCCCG−1_2
GCCCAGATCGTTGCCT−1_2
TCCACCACACCAGTTA−1_2
TCGCTCAAGTCTGCGC−1_2
CCTTGTGCAACTCCAA−1_2
TCATTACGTCTTTATC−1_2

STB2
TGGCGTGTCTGTCGCT−1_2
GAAGAATGTAACGGTG−1_2
TAACACGGTTACCCTC−1_2
ACCTGTCTCGGAATTC−1_2
TGCCGAGGTTACCTGA−1_2
CTCAAGACAAGCGCTC−1_2
TGAATCGAGCCTCTCT−1_2
TGCACGGGTACCGTGC−1_2
AACCTTTTCGTGGACC−1_2
ATCGTGATCCTACCGT−1_2
CATGCAACAGGCATGA−1_2
GTGCGTGGTAACACCT−1_2
TAAGTCGAGTTGCTGT−1_2
TGAGGTTTCTGACCCT−1_2
CCGGTGATCCGTCAAA−1_2
CCTAAGAGTACGTACT−1_2
AGCCACGTCCCTCTAG−1_2
TTGCCTGTCTTTCAGT−1_2
GAACACTTCGCGATCG−1_2
TACGTCCGTCTGTGCG−1_2
CGATGCGGTTTGGCTA−1_2
ATAGAGAGTAGACAAT−1_2
CGTAAGTAGCTATCCA−1_2
ATTCTTGTCCAGCAAT−1_2
CCACAAACACCCGTAG−1_2
TAAGCCAGTTGCGGAA−1_2
GATAGCTGTAGCGAGT−1_2
AAGTCGTTCTTACCGC−1_2
ATTCATCGTCTTACAG−1_2
ATACTTCGTCGTATGT−1_2
GACAGCCTCTGCGGCA−1_2
CAAGGGATCATCACAG−1_2
AATGACCAGTATCTGC−1_2
ACTACGATCTTCTTCC−1_2
ACTCTCGGTTCTCTCG−1_2
CTGATCCTCATACGGT−1_2
GTATTGGTCGAGTCCG−1_2
CCGGTAGTCGGTTGTA−1_2
TAACACGAGAGCGACT−1_2
TTGCATTGTATGGAAT−1_2
GCTTCACGTTTGAAAG−1_2
TTGGGATAGGTTGTTC−1_2
ACTTATCAGTGACACG−1_2
GATTGGTAGTCAACAA−1_2
TTGAACGGTCGAAGCA−1_2

X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList2,]
Xp <- na.omit(Xp)[,subs1]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs1])), Batch =  factor((mammal.combined$Genotype2[subs1])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)



pm1 <- pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 20, filename = paste(saveext,"/DimRed/Tb_Batch1_MinList",".pdf",sep="") ,width=10,height=4)
hc1 <- cutree(pm1$tree_col,20)
saveRDS(hc1,file=paste(saveext,'TbCl_batch1.rds',sep="") )

#lbl <- cutree(hc, 20) # you'll need to change '5' to the number of gene-groups you're interested in


#Clusters1 <- data.frame(cluster = cutree(pm1$tree_row, k = 20))
#x[map$tree_row$order,]


subs1 <- which(mammal.combined$ID3%in%c("C2") & mammal.combined$Cells%in%c("CTB_d14","STB_d14","EVT_d14","putSTB") )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList2,]
Xp <- na.omit(Xp)[,subs1]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs1])), Batch =  factor((mammal.combined$Genotype2[subs1])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm2<-pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Tb_Batch2_MinList",".pdf",sep="") ,width=40,height=4)
#hc2 <- cutree(pm2$tree_col,20)
#saveRDS(hc2,file=paste(saveext,'TbCl_batch2.rds',sep="") )


CTBlist
TCGTGCTAGCAACAGC−1_3
CCGGTGAAGCAATTAG−1_3
CGTTCTGAGACCCTTA−1_3
GGAGCAACATTCTTCA−1_3
CTCATCGCATCAACCA−1_3
GTTAGTGAGCCTAACT−1_3
ATCGATGAGCAAACAT−1_3
TAGGTTGCACGGGTAA−1_3
AGCTACATCTCGTGAA−1_3
CGCAGGTCAGAAACCG−1_3
GCCAGTGTCTATTCGT−1_3
GCGGAAAAGACGGTCA−1_3
AATTCCTAGCCAAGGT−1_3
CCGTAGGGTAAGTTAG−1_3
TACTTCATCGACCAAT−1_3
ATCGTGACAGACCAGA−1_3
GAGGCAAGTTCGGCCA−1_3
TAACTTCTCTACTCAT−1_3
CTGGTCTCAATAAGGT−1_3
TCGAACACAAGAGCTG−1_3
GGAGAACGTCGCGGTT−1_3
CATACTTAGGTCGTGA−1_3
GTCTACCCAACACGTT−1_3
TAGACCAAGGAGACCT−1_3
CTCCAACGTTCTCCAC−1_3
TACTTACGTACCGGCT−1_3
TCGATTTCAATTTCCT−1_3
GAGCCTGGTCGTACAT−1_3
TGTGATGCAACCCTAA−1_3
AGACCCGAGAGCAACC−1_3
TGGATCAGTCCAGCGT−1_3
ACGCACGTCGTTACCC−1_3
AGTGCCGGTTCACCGG−1_3
GTAGGAGAGTCCTACA−1_3
TACGGTATCATGCAGT−1_3
AACAACCAGATTAGCA−1_3
TCGGGTGAGATTCGCT−1_3
CCTACGTAGGACTTCT−1_3
CATTGCCCAACCCGCA−1_3
CTCCCTCCAGCACAAG−1_3
AAGTTCGTCCACGTGG−1_3
ATGACCAAGTTCGGTT−1_3
TCCTCGATCTGCGAGC−1_3
TCTGGCTAGCTGTGCC−1_3
AACCACACAAAGCTAA−1_3
CCACTTGGTTCAGTAC−1_3
AACTTCTAGGTTGTTC−1_3
CGTTCTGTCCAAGAGG−1_3
GAGACTTCATCAGCGC−1_3
CAGGTATTCTGGTCAA−1_3
CTTCAATTCCTCTCTT−1_3
ATCTTCACACTAACCA−1_3
GTCATGACACATATCG−1_3
TGCTCCAGTCATCCCT−1_3
CCGGTAGGTGGCATCC−1_3
GAGGGATTCGTAGTCA−1_3
TGGTGATCAATCTCGA−1_3
AACACACAGGGCAATC−1_3
CTAACCCAGATAGTGT−1_3
GTCTCACAGGCCTGCT−1_3
CATGCGGGTGTCTTCC−1_3
GTCCTCACACAGTCCG−1_3
GTTACGATCTGCTCTG−1_3
TTGAACGTCTGGCCAG−1_3
TGGCGTGAGCGTGAGT−1_3
ACCCTTGAGGTTTACC−1_3
TACCGGGAGCTGTTAC−1_3
CTGTGGGGTCCAGCGT−1_3
GAGTCTATCCCGATCT−1_3
CGAGGCTAGGACAACC−1_3
TGATTTCAGCTCATAC−1_3
ACGGTTAAGAAACTGT−1_3
CCGCAAGGTACCGGCT−1_3
AACAAAGTCCGCCTAT−1_3
GCACATAGTCATCCCT−1_3
TCATTTGCAAGTGGTG−1_3
GTGGTTATCCGGACTG−1_3
TCGCACTCAAATAGCA−1_3
GGTAACTGTGTTCGTA−1_3
CAAGGGAGTACTCGTA−1_3
TGAACGTAGTTTGAGA−1_3
TAATCTCTCTGCGGGT−1_3
CACGTTCCAAGTGGAC−1_3
GTAGCTACAACAGAGC−1_3
CGGGTGTAGTTGTAGA−1_3
TGATTCTGTGTAAATG−1_3
CGGACACAGGTCCCGT−1_3
CTGGCAGCAGGACTTT−1_3
CTCATTATCCCAAGCG−1_3
TACGGTATCCCACAGG−1_3
GACACGCTCCACATAG−1_3
TACCTGCGTGTTTCTT−1_3
GCAGCCAGTATTAAGG−1_3
AGTGACTGTGTAGCAG−1_3
AGGGAGTAGTATCTGC−1_3
GGCGTCACAATGTTGC−1_3
CATACTTGTGAATGAT−1_3
AAACGAATCATTCATC−1_3
ATGGGAGCACGCTGAC−1_3
GCACATACAGCACAAG−1_3
CGTCCATGTACAGTCT−1_3
GAACTGTCAATAGGAT−1_3
TCATGCCAGCCTAGGA−1_3
GCCATGGGTAGCGCCT−1_3
CATCGGGGTAAGGTCG−1_3
GAGCTGCGTTCCATTT−1_3
CACAACAAGATGTTGA−1_3
GATGGAGTCCTAAGTG−1_3
ACGTTCCCAATCTGCA−1_3
AGCTTCCTCGCAGTTA−1_3
ATAGGCTCACACGGAA−1_3
GAAGTAAGTCGTTGGC−1_3
TTCTTGAGTATCGAAA−1_3
GGAGGTATCTAGACCA−1_3
GTGCGTGAGAATCCCT−1_3
TAGTGCATCTCTGGTC−1_3

STB2
GGCTTTCAGCACTTTG−1_3
TGCAGGCAGGAGGCAG−1_3
ACACAGTTCTGCTTAT−1_3
ATGCATGCAGGGTCTC−1_3
CTGTAGAAGAGCACTG−1_3
GCTTGGGTCGCTCTCA−1_3
TACCTGCCATCGAAGG−1_3
ACGTACACAGTTGCGC−1_3
TGCGACGAGCCTATCA−1_3
GATTTCTTCTTGTGCC−1_3
TCCCAGTGTCAAACGG−1_3
GCTTCACTCTTCGATT−1_3
CTCCCAACATCAGTCA−1_3
TACAACGTCCGGTTCT−1_3
TTCGGTCAGGGAACAA−1_3
CGAGTGCCAATACGAA−1_3
GCCATTCTCGACGACC−1_3
ACCAAACCAGAATGTA−1_3
ATAGGCTGTAGGATAT−1_3
AAGACAAGTACGTGAG−1_3
CTTTCAAAGCTGGTGA−1_3
TCCTCTTAGGAACGCT−1_3
ATTCAGGAGCTCGTGC−1_3
TCACTCGTCCACACCT−1_3
CCACACTCACTCAGAT−1_3
CCCTGATAGCTCTGTA−1_3
AGGGTTTTCATTATCC−1_3
GGCGTCATCAGGAAGC−1_3
TCGTGCTGTGGATACG−1_3
CTATCTAAGTGCCCGT−1_3
ATTTCTGCAAGACTGG−1_3
TCCGGGATCCATATGG−1_3
TTTATGCGTGCCGTAC−1_3
GACCCTTTCGTGGTAT−1_3
GCGAGAATCAGTCAGT−1_3
CCAAGCGCACAAGCTT−1_3
TTCCTAACAACAGTGG−1_3
CACTAAGCAGGCTATT−1_3
TACCTGCGTTGCATGT−1_3
TAGGAGGAGTGATCGG−1_3
AAGAACAGTAGTACGG−1_3
AGGCATTGTTGCCTAA−1_3
TTCACCGCAGCTCTGG−1_3
AGGATAATCGGTAAGG−1_3
CTACCTGGTCATAACC−1_3
CTCCAACGTCCTGTCT−1_3
ACTTATCAGGTAAGGA−1_3
ATGGGTTTCGAGCACC−1_3
GATGGAGTCTGGCCGA−1_3
GTGCTTCTCCAACCGG−1_3
CTCAGTCAGGCGTCCT−1_3
CTGCATCAGAGAACCC−1_3
TTTGGAGCAGAGTTGG−1_3
GGTTCTCGTGCAACAG−1_3
TCCTGCAGTTGGCCTG−1_3
GGAGAACCATAGACTC−1_3
GGGACTCAGTCTTCGA−1_3
TCATTCACACTTACAG−1_3
TCGAACAGTCGCACGT−1_3
ATTTACCCACCATTCC−1_3
TTTCATGGTCACGACC−1_3
ATCGCCTCACTAAACC−1_3
GGGCGTTCAACTCCCT−1_3
GTTACCCGTGTTAACC−1_3
AGAAATGCAAGTCCAT−1_3
AAGCGTTAGTAAATGC−1_3
CGCCATTGTGGCTACC−1_3
ACCTGTCGTACTCGTA−1_3
TGTCCACTCGCACGGT−1_3
CTCAAGATCTATCCAT−1_3
GTGATGTCACCCTCTA−1_3
TATCTGTCAGAGGGTT−1_3
TACTGCCCAAAGGGCT−1_3
TTGATGGAGAAATTGC−1_3
TTGCTGCGTTGAGAGC−1_3
ACCTGAAGTTGTGTAC−1_3
TTACCATTCGTGAGAG−1_3
AGACCATAGCATGAAT−1_3
GTAGGTTTCCTTGGAA−1_3
ATGCCTCGTGGGTCAA−1_3
ACTGATGCACTTCAGA−1_3
CACGTGGCACTTGTGA−1_3
GTCTAGAAGCAACAAT−1_3
GCTCAAACATTGCCGG−1_3
AGTCAACTCCATCTGC−1_3
TCAGGGCCATACCATG−1_3
ACCTGTCTCGCTTGCT−1_3
GTATTTCTCAGCAGAG−1_3
TTGTTTGAGATTACCC−1_3
CAACGATTCAACTGGT−1_3
TGTTCCGGTGATATAG−1_3
TATTGCTGTCTGCAAT−1_3
CCGTGAGCAAGTCCCG−1_3
GTAGAAAAGGTACTGG−1_3
GAGGGTAAGAAACCAT−1_3
ATTCATCCAAGTGACG−1_3
AGAACCTTCGGACTTA−1_3
TTACGTTCACGACAGA−1_3
TGCATCCGTGGGTTGA−1_3
AGCCACGCAAGCCTGC−1_3
CGCCAGAAGGCGTTAG−1_3
TCACTATTCACAATGC−1_3
TCTATACAGATACATG−1_3
TCCCACAGTAATTAGG−1_3
TGTGGCGCAACCTAAC−1_3
ATTCCATCACACAGCC−1_3
CTCAATTGTGTCCAAT−1_3
GGGTCACCACCGTACG−1_3
GCTTCACAGGCATCGA−1_3
GGTCTGGAGGCCTGAA−1_3
ATATCCTAGACTCAAA−1_3
TAGATCGGTGTTAAAG−1_3
TCATGAGGTAGTCTGT−1_3
TGTGCGGGTCCTATAG−1_3
ACCTACCTCCCATTTA−1_3
TCTGGCTTCATAGGCT−1_3
CCTCACATCTAACACG−1_3
TTCTAGTGTAGGAAAG−1_3
CGAAGGAAGAATCGCG−1_3
GATAGCTCACTCCGGA−1_3
TAATCTCTCGCAATGT−1_3
CGGAACCCATGGGATG−1_3
ATCCACCAGGCATCAG−1_3
ATGCGATTCGCATTAG−1_3
AGGTCATAGGCCTTGC−1_3
TAGGTTGTCATTACTC−1_3
CAACCAATCAAGAGGC−1_3
TCAGTGACACTCTCGT−1_3
GTCAAACCATCGATGT−1_3
ACACAGTAGCAAGCCA−1_3
GGCTTTCCAACATACC−1_3
TGTCCCAAGTGAATAC−1_3
TTCTTGAGTCCTATAG−1_3
GTCACTCTCGGAGCAA−1_3
GTCTAGACAGTGGTGA−1_3
CTTACCGAGGGTACAC−1_3
CGTAGTACAATACCTG−1_3
CTACCCACATGGGATG−1_3
AGTCATGCAAGAGCTG−1_3
GGTGATTGTTCGTACA−1_3
AAGCCATTCCCGTTGT−1_3
TGGTACAGTGGTCCGT−1_3
ACTGTCCCATGGAACG−1_3
ACTATCTCAAGCGCTC−1_3
CAACAGTGTTTGGCTA−1_3
TACGGTATCACGATAC−1_3
GAGTCTAAGTAGTCCT−1_3
GTAACCACATTCTGTT−1_3
CTGCCTACATCTGGGC−1_3
GTGTGATGTCATATGC−1_3
AAAGTCCCACTGAGGA−1_3
GACCCTTTCTAACACG−1_3
GTGAGCCTCTTACGTT−1_3
GAAATGAGTCTTGCGG−1_3
CGGAATTTCAAGTAAG−1_3
GACTGATAGCGCCATC−1_3
GAGCCTGGTTGTAAAG−1_3
AGAGAGCTCTATTCGT−1_3
GCAGCCATCATCTATC−1_3
GTGGCGTTCGAAGAAT−1_3
ATTGTTCCAGCGTGCT−1_3
GAGACTTAGGTGCCTC−1_3
AGCTTCCGTGATACTC−1_3
ATTTACCTCATGGCCG−1_3
GTAGGAGAGTGCACTT−1_3
TAGACCATCACAATGC−1_3
TGTACAGGTTAAGTCC−1_3
GTGTAACCAAGCAATA−1_3
TCAGCAAGTTGAATCC−1_3
GAAGGACAGGAACATT−1_3
CGTCAAAGTGCTAGCC−1_3
AAAGGATGTTACCCAA−1_3
ACGATGTCACAGCCAC−1_3
ACCACAACATAACAGA−1_3
GGGAAGTCAAAGAGTT−1_3
CTAACTTGTTGTCAGT−1_3
ACGTCCTCACCGTGAC−1_3
CAGATACAGGGCAACT−1_3
CATGCGGTCACAACCA−1_3
CCGGTGACATAACCCA−1_3
ATTATCCCAAATAAGC−1_3
GGGTGAAAGCCGTTAT−1_3
CATTGCCTCATCAGTG−1_3
TTTGATCAGGTAGCCA−1_3
CCGGACAAGTCGAATA−1_3
GCTTTCGAGTTTCGGT−1_3
AGATCCATCATACAGC−1_3
GGCTGTGTCCTTTGAT−1_3
AGACAAACAGGGATAC−1_3
CTCCACAAGCAATAGT−1_3
ACAACCAGTTTATGCG−1_3
TCCTCGACAGACACCC−1_3
GTTGCTCAGCTGCCAC−1_3
TCCTGCATCAGACCTA−1_3
GAAGGACTCGTTCCTG−1_3
TGGAACTCAGGTGTTT−1_3
ATACCTTGTTATGACC−1_3
GTGTGGCCAGTCGCTG−1_3
AGTGCCGTCTTTCTAG−1_3
GTTCATTGTATTTCGG−1_3
CCACTTGTCGCACGAC−1_3
CTCACTGTCACCGGTG−1_3
ACTATTCAGTATTAGG−1_3
GGAATCTGTCCTGTCT−1_3
GTCGCGAGTGGTCCCA−1_3
TGATGCAAGGGTCAAC−1_3
TGCGACGAGATCACCT−1_3
ACCCTCATCATCAGTG−1_3
CATCAAGAGATACTGA−1_3
ATAGGCTGTGTAGTGG−1_3
ACGTTCCGTCAGTTTG−1_3
ATCCACCTCATTATCC−1_3
CTTTCAAAGAGCCATG−1_3
TCTATCAAGAAATCCA−1_3
AGCGCTGCAACCGTGC−1_3
GAGAGGTTCGACACTA−1_3
CCTATCGAGTTGCCCG−1_3
GGGCCATTCTTCCACG−1_3
TGTTTGTGTTAAACAG−1_3
ATCAGGTGTCCTTTGC−1_3
TGACCCTTCGAACCTA−1_3
ATTCGTTTCCTTGACC−1_3
ACCTGAAGTAAGATTG−1_3
GTTATGGCAAAGACGC−1_3
GAGACCCAGCCGCACT−1_3
TTTGACTAGCACTCCG−1_3
TCCTAATTCATTCATC−1_3
ATTATCCAGTCTCTGA−1_3
ATACTTCTCTTGCAGA−1_3
CCTTTGGTCTATGCCC−1_3
TTTGGAGGTGCTCTTC−1_3
CAATGACCACTACCCT−1_3
TAGGGTTTCCTCATAT−1_3
TGCCGAGTCGTTGTTT−1_3
TGCAGTATCGGAATTC−1_3
CATGGTACACACACTA−1_3
TCGGTCTCAGCAGGAT−1_3
GATTGGTTCGCCCAGA−1_3
GCCGATGTCTGAGGCC−1_3
CCGGACATCGAGAGCA−1_3
CAACAACTCCTACGAA−1_3
ACTATTCGTCTGTGAT−1_3
TGCAGATTCGATGGAG−1_3
ATCCCTGGTGGCAGAT−1_3
AGACCATTCTAAACGC−1_3
CATTGTTTCCGATTAG−1_3
TTGGGCGAGTGGTGAC−1_3
AAACGAACACAATTCG−1_3
TCCACCATCATAGACC−1_3
CAGATCACACGCTTAA−1_3
AACCCAAAGCTACTAC−1_3
CACTGGGTCCTATTTG−1_3
GGGACCTTCCAAGCCG−1_3
AGGTTACTCCTTCTGG−1_3
TCAATTCTCACAATGC−1_3
GCATCTCAGACGGTCA−1_3
AACCTTTTCAACTCTT−1_3
TTGGGATCATGGATCT−1_3
CATTGCCGTAGCTGCC−1_3
CATACCCAGCCTAACT−1_3
AGAACAAAGGCTTTCA−1_3
GCAGCTGCACTTTAGG−1_3
ACTATCTGTACACGCC−1_3
TACTGCCTCGACTCCT−1_3
GAGGCAACATGCGGTC−1_3
TTACGCCTCATTGTTC−1_3
ATCCGTCTCAGGTAAA−1_3
CCACCATAGTGAGTTA−1_3
TTCTTCCTCTATCCAT−1_3
GCCCAGATCACGATAC−1_3
TTTACCAAGTGCACTT−1_3
AAGAACAGTATTGGCT−1_3
CCCAACTGTGTTCATG−1_3
TTCTTCCGTTGCGTAT−1_3
ATACTTCTCATGCCAA−1_3
TATACCTTCCAACTAG−1_3
AAGATAGTCTCTGACC−1_3
CACAACAGTAGCGTAG−1_3
TGACAGTCATTGAGCT−1_3
ACAGAAACAGCAGAAC−1_3
TGCTCGTTCGCAATTG−1_3
AGGTCATCAGAAGTTA−1_3
ACTTATCAGCGTGCCT−1_3
ATCACAGAGGCGTTGA−1_3
GTGTGATGTCTTGAGT−1_3
TAATTCCGTTGTTGAC−1_3
CATTCATCATGACTGT−1_3
GTTCGCTGTCATCTAG−1_3
CCTTTGGTCCACAAGT−1_3
GTAGAAATCCGTGGGT−1_3
AAGCATCGTGGTCAAG−1_3
CCTCCTCCAGAAGCGT−1_3
AGCGCCATCGAACGCC−1_3
CATACAGGTCTAACGT−1_3
TCTGGCTTCTCCTGCA−1_3
ACTTCGCGTCTAGATC−1_3
TCACACCCAGGTTACT−1_3
ACCAACATCATCACTT−1_3
CCCTGATGTTGTGGAG−1_3
AGCGCTGAGACCAAAT−1_3
CACCGTTGTATAGCTC−1_3
GCACGGTCAAAGCGTG−1_3
GTGTTAGGTGCCTATA−1_3
CAACGGCTCAGTCCGG−1_3
TCATATCAGTAGTCTC−1_3
TCCATGCGTTGCAACT−1_3
TTTGGAGTCGTAATGC−1_3
ACTTTCACAACATACC−1_3
CACAACAGTCAGTCGC−1_3
CACGGGTTCATTTACC−1_3
AGGTAGGGTTAGCTAC−1_3
TGAGGTTCAGTGTGGA−1_3
TGGTACAAGGAAAGTG−1_3
TCCACGTTCACAAGGG−1_3
TTTCAGTAGCGTCTCG−1_3
TAGGTACTCCCAGCGA−1_3
GTTCATTTCCTAAACG−1_3
CATCGGGGTCAGTTTG−1_3
GTTGCGGTCTCGACGG−1_3
CACGGGTAGTGGTGGT−1_3
TCCCATGAGTCGTCTA−1_3
CTACTATCATCGATGT−1_3
ATTCACTCACGTCTCT−1_3
AACTTCTGTAGCCCTG−1_3
AGTGCCGTCACTGGGC−1_3
AGATAGATCTCCCTAG−1_3
GCTTCACCACCGCTGA−1_3
CTAGGTATCCAGGACC−1_3
AAGACTCTCCAAATGC−1_3
TTTGACTCACACACGC−1_3
ATCTTCACACAAATAG−1_3
GACCAATTCGTACACA−1_3
TCCTTCTAGGCCACCT−1_3
TCAATCTGTGAAAGTT−1_3
TCTATACGTTAACCTG−1_3
TGAGACTCACTCCGGA−1_3
TCATACTAGACGCCCT−1_3
GAGTTTGGTCGAGCAA−1_3
GCCAGGTAGCTTGTGT−1_3
GCTGAATCAAGGTTGG−1_3
TTCCACGGTCTTAGTG−1_3
TCACGGGAGTAAGACT−1_3
GAAGCGACATAACTCG−1_3
GAATCACGTTTCGATG−1_3
TTGACCCTCCTTGGAA−1_3
ATGTCTTCAGAGACTG−1_3
TCGGTCTTCCACCTCA−1_3
ATCAGGTAGAGTAACT−1_3
CATAGACAGGATATAC−1_3
CACTGGGCAGATACTC−1_3
TTTCCTCTCTCCGATC−1_3
TATCTTGCACATAACC−1_3
TTTGATCGTAGCTTGT−1_3
GTAGATCTCGGCTATA−1_3
ACCAAACAGAATCGTA−1_3
ACTGCAATCGCTCTCA−1_3
CAGATCAGTGGAATGC−1_3
CTAGACAAGAAAGTCT−1_3
GGAATCTTCGTTCTAT−1_3
GTTAGACTCGAGTTGT−1_3
TGCGGGTAGACTACGG−1_3
CTTCGGTAGTCTTGGT−1_3
TTCGCTGTCTCTCAAT−1_3
ATTCATCGTGTCTAAC−1_3
GGAGCAAAGACTGAGC−1_3
ACCCTCACAGTTCTAG−1_3
TAGGTTGGTTACCCAA−1_3
CTGATCCGTTCAGCTA−1_3
TATCTGTTCCGCACTT−1_3

EVT
TGGAACTCAGGTGTTT−1_3
CTCCACAAGCAATAGT−1_3
GTAGGAGAGTGCACTT−1_3
GACTGATAGCGCCATC−1_3
ACTATCTCAAGCGCTC−1_3
TAGGAGGAGTGATCGG−1_3
CTGGCAGCAGGACTTT−1_3
CTCATTATCCCAAGCG−1_3
GTAGCTACAACAGAGC−1_3
GTTACGATCTGCTCTG−1_3
TTGAACGTCTGGCCAG−1_3
TGGCGTGAGCGTGAGT−1_3
ACCCTTGAGGTTTACC−1_3
TACCGGGAGCTGTTAC−1_3
CTGTGGGGTCCAGCGT−1_3
GAGTCTATCCCGATCT−1_3
CGAGGCTAGGACAACC−1_3
TGATTTCAGCTCATAC−1_3
ACGGTTAAGAAACTGT−1_3
CCGCAAGGTACCGGCT−1_3
AACAAAGTCCGCCTAT−1_3
GCACATAGTCATCCCT−1_3
TCATTTGCAAGTGGTG−1_3
GTGGTTATCCGGACTG−1_3
TCGCACTCAAATAGCA−1_3
AAGTTCGTCCACGTGG−1_3
ATGACCAAGTTCGGTT−1_3
TCCTCGATCTGCGAGC−1_3
TCTGGCTAGCTGTGCC−1_3
AACCACACAAAGCTAA−1_3
CCACTTGGTTCAGTAC−1_3
AACTTCTAGGTTGTTC−1_3
CGTTCTGTCCAAGAGG−1_3
GAGACTTCATCAGCGC−1_3
AGTGCCGGTTCACCGG−1_3
GTAGGAGAGTCCTACA−1_3
TACGGTATCATGCAGT−1_3
AACAACCAGATTAGCA−1_3
TCGGGTGAGATTCGCT−1_3
CCTACGTAGGACTTCT−1_3
CATTGCCCAACCCGCA−1_3
CTCCCTCCAGCACAAG−1_3
TGGATCAGTCCAGCGT−1_3
TCGTGCTAGCAACAGC−1_3
CCGGTGAAGCAATTAG−1_3
CGTTCTGAGACCCTTA−1_3
GGAGCAACATTCTTCA−1_3
CTCATCGCATCAACCA−1_3
GTTAGTGAGCCTAACT−1_3
ATCGATGAGCAAACAT−1_3
TAGGTTGCACGGGTAA−1_3
AGCTACATCTCGTGAA−1_3
CGCAGGTCAGAAACCG−1_3
GCCAGTGTCTATTCGT−1_3
GCGGAAAAGACGGTCA−1_3
AATTCCTAGCCAAGGT−1_3
CCGTAGGGTAAGTTAG−1_3
TACTTCATCGACCAAT−1_3
ATCGTGACAGACCAGA−1_3
GAGGCAAGTTCGGCCA−1_3
AGCGCCATCGAACGCC−1_3
TGACAGTCATTGAGCT−1_3
CATTGCCGTAGCTGCC−1_3
CATACCCAGCCTAACT−1_3
AGAACAAAGGCTTTCA−1_3
GCAGCTGCACTTTAGG−1_3
ACTATCTGTACACGCC−1_3
CATCGGGGTAAGGTCG−1_3
GAGCTGCGTTCCATTT−1_3
CACAACAAGATGTTGA−1_3
GATGGAGTCCTAAGTG−1_3
ACGTTCCCAATCTGCA−1_3
AGCTTCCTCGCAGTTA−1_3
ATAGGCTCACACGGAA−1_3
GAAGTAAGTCGTTGGC−1_3
TTCTTGAGTATCGAAA−1_3
GGAGGTATCTAGACCA−1_3
GTGCGTGAGAATCCCT−1_3
TAGTGCATCTCTGGTC−1_3
GCACATACAGCACAAG−1_3
CGTCCATGTACAGTCT−1_3
GAACTGTCAATAGGAT−1_3
TCATGCCAGCCTAGGA−1_3
AGCGCTGAGACCAAAT−1_3
ATGTCTTCAGAGACTG−1_3
TCGGTCTTCCACCTCA−1_3
ATCAGGTAGAGTAACT−1_3

strongEVT
CATAGACAGGATATAC−1_3
CACTGGGCAGATACTC−1_3
TTTCCTCTCTCCGATC−1_3
TATCTTGCACATAACC−1_3
TTTGATCGTAGCTTGT−1_3
GTAGATCTCGGCTATA−1_3
ACCAAACAGAATCGTA−1_3
ACTGCAATCGCTCTCA−1_3
CAGATCAGTGGAATGC−1_3
CTAGACAAGAAAGTCT−1_3
GGAATCTTCGTTCTAT−1_3
GTTAGACTCGAGTTGT−1_3
TGCGGGTAGACTACGG−1_3
CTTCGGTAGTCTTGGT−1_3
TTCGCTGTCTCTCAAT−1_3
TCCTTCTAGGCCACCT−1_3
TCAATCTGTGAAAGTT−1_3
TCTATACGTTAACCTG−1_3
TGAGACTCACTCCGGA−1_3

subs1 <- which(mammal.combined$Genotype2%in%c("NC2_1") & mammal.combined$Cells%in%c("CTB_d14","STB_d14","EVT_d14","putSTB") )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList2,]
Xp <- na.omit(Xp)[,subs1]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs1])), Batch =  factor((mammal.combined$Genotype2[subs1])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm2a<-pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 20,  filename = paste(saveext,"/DimRed/Tb_Batch2A_MinList",".pdf",sep="") ,width=25,height=4)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Tb_Batch2_MinList",".pdf",sep="") ,width=40,height=4)
hc2a <- cutree(pm2a$tree_col,20)
saveRDS(hc2a,file=paste(saveext,'TbCl_batch2a.rds',sep="") )



CTB
TTCAGGACATAATCCG−1_3
GAGCTGCTCTTTGCAT−1_3
TGGAGAGAGTATAACG−1_3
ACTTAGGCAGTAACGG−1_3
ACACTGATCCATACTT−1_3
GATCATGAGTATCTGC−1_3
AAAGGATCACACACGC−1_3
GTCAGCGAGCAATAAC−1_3
AGGGTGATCGGCTCTT−1_3
CCTAACCTCTTGATTC−1_3
GAATCGTGTGTCGATT−1_3
ATGACCATCCGTGTCT−1_3
TCAGCCTCACTTCAAG−1_3
ACGTCCTAGTTGTAGA−1_3
TATATCCTCCATTTGT−1_3
GGGACAAAGCGGTAAC−1_3
AAGTCGTTCGAGATAA−1_3
AGCTTCCTCGGACTGC−1_3
AACCCAACATGTTACG−1_3
TCCCAGTTCTGTCGCT−1_3
CGATCGGCAATACAGA−1_3
ATGACCAGTGCGACAA−1_3
GACTATGTCGATGCTA−1_3
CACGTGGGTCGTTGCG−1_3
CCTAAGAGTAACATAG−1_3
TCCACGTCACCACTGG−1_3
GCCTGTTGTAATACCC−1_3
TGTTTGTAGCAGCCTC−1_3
CGTAAGTCAAGCGAGT−1_3
TTCCTCTAGTCATCCA−1_3
ATCTCTAGTACCTATG−1_3
CATGCCTCAGGGATAC−1_3
GAGACCCCACTGCGAC−1_3
GGACGTCCAGCCCACA−1_3
GTGCACGGTGGACTGA−1_3
GTAATGCGTACTGACT−1_3
TCAGTGAAGGCTAGCA−1_3
ATTCAGGGTCTCGACG−1_3
GGGATGAGTACGGCAA−1_3
GTGTGATAGCAGCAGT−1_3
TCTCCGAGTCGATTAC−1_3
GCCTGTTGTACCTAGT−1_3
TCGTGGGCAACGGCCT−1_3
TCATACTAGTATAACG−1_3
TGATCTTCACACCTGG−1_3
TTTACCATCATCGGGC−1_3
ACTATCTTCGTTGTTT−1_3
TCACGCTTCAGGAAAT−1_3
CCTTTGGAGGTTGTTC−1_3
TTCATGTCACCTTCGT−1_3
TTTATGCTCCGGACGT−1_3
ACACAGTTCGATACTG−1_3
TCACAAGTCACGAGGA−1_3
CACAACAAGTAAACAC−1_3
CGAAGGATCACCACAA−1_3
GTAGGTTGTCCTCATC−1_3
AATGAAGGTAGACGTG−1_3
GAGTTACAGCAGCCCT−1_3
CTCAGGGCACGATTCA−1_3
TCATTACAGGACGGAG−1_3
CGCAGGTTCGGCCTTT−1_3
TGAGTCATCATCTACT−1_3
TCATCATCACGCTGAC−1_3
CCTCTCCTCGTTGTGA−1_3
TGACTCCTCAGTGCGC−1_3
ATCTTCACACTTGAGT−1_3
TTGACCCAGCTGTCCG−1_3
ATTCTTGTCAGACCCG−1_3
TTCTTGACACGAAGAC−1_3
GTGCTGGAGACCTGGA−1_3
ACTTCCGAGTTCAACC−1_3
TCTATACCATCCGGTG−1_3
CGAGGAATCGATACGT−1_3
CTGCCATTCTGTCTCG−1_3
TGGTTAGCACAGCGCT−1_3
ATCCACCTCTCAGTCC−1_3
AAACCCAGTGGTTTAC−1_3
AAGCATCCACACGGTC−1_3
GCTGAATAGAAGCTCG−1_3
AATGCCACAATCTCGA−1_3
AGACACTGTTGTTGAC−1_3
ATTGTTCTCCGGCAGT−1_3
CTGTGAAAGTGGTTCT−1_3
TGTAGACCACTAGTAC−1_3
CATGCGGTCTGGGAGA−1_3
CCCTCTCGTAGCGTTT−1_3
TTGGATGAGCTCAGAG−1_3
TCGACGGGTACCGGCT−1_3
AATGGCTGTACTCAAC−1_3
TTGACCCCATCCTTGC−1_3
CTCAACCAGTTCATGC−1_3
GCCAACGTCACACCCT−1_3
TACCCGTTCCAGGACC−1_3
TTTGACTAGCGTGTCC−1_3
TGGTAGTAGTCTCTGA−1_3
CACAACAAGTCTTCGA−1_3
TTTGATCAGGTTATAG−1_3
CGCATGGAGTCTAGCT−1_3
CAGATACGTGCCCGTA−1_3
CGAGAAGCATTGACAC−1_3
GGTTGTAAGACGGATC−1_3
AATGGCTAGGTTGGAC−1_3
CACCAAAAGCGTGTTT−1_3
CTATAGGAGTGATAGT−1_3
CAACGGCTCACATACG−1_3
GACGCTGTCGAGCCAC−1_3
TGCCGAGCATCATTTC−1_3
TGGATGTAGGAGGCAG−1_3
GCCCGAAGTTGTTTGG−1_3
GTTAGACCAGGTTACT−1_3
AAGCCATGTGGCTACC−1_3
ACCATTTCACGGGCTT−1_3
TATTCCAGTGATACAA−1_3
TACCTGCCATGCCGCA−1_3
CTAACCCAGGCACCAA−1_3
TACTTCACATGAATCC−1_3
TTCTAACTCGCCTCTA−1_3
ATTACTCAGCTTTCTT−1_3
CCTCTAGCATGCCGAC−1_3

STB
TCGAAGTAGTCGAAAT−1_3
AGACCATGTACTAACC−1_3
GTAAGTCAGAGTAACT−1_3
CGCATAAGTATCGAGG−1_3
GGGAGATCAACCCGCA−1_3
GTTCCGTCACGACAGA−1_3
GCCCAGATCTCGTGAA−1_3
CACGTTCTCAAAGACA−1_3
CATGCCTTCCACTAGA−1_3
CTCATGCCACACACGC−1_3
CGATGCGCATCATTGG−1_3
TAATTCCCACACCAGC−1_3
GTGCAGCAGGTTAGTA−1_3
CAGGCCAAGGCATCAG−1_3
GTTTGGATCGCACGAC−1_3
AGCGCCACAAATCGTC−1_3
CCGTAGGCAGTCAACT−1_3
TACCCACTCTGTGCAA−1_3
AGGGCCTCATCCCGTT−1_3
CCTCTCCGTGTGATGG−1_3
ATCGCCTAGCACAAAT−1_3
TCACAAGCAGGCGAAT−1_3
CATAGACAGCAACAGC−1_3
CCTATCGAGCACTTTG−1_3
AGTAGTCAGTGTTGTC−1_3
GCTACAACAGAGGTAC−1_3
CCCATTGCAAGTGCAG−1_3
TTCACCGCAGTGTACT−1_3
GCTGAATGTCTTGCTC−1_3
GGTAATCCAACTCCCT−1_3
ACAACCATCTGTACAG−1_3
ATCGTCCCATCGGATT−1_3
CTGCATCTCTCACCCA−1_5
AGTACCACACACCTAA−1_3
GATCACAGTCACTTCC−1_3
CGTAGTATCTCATGCC−1_3
GTGAGGATCTTGGCTC−1_3
CAGTGCGTCGATCCAA−1_3
TTAATCCCAGCGACAA−1_3
GGATGTTCATGCACTA−1_3
GTCGTAAGTGGTCTGC−1_3
TATCGCCCATGGAACG−1_3
AACCATGGTTGTTTGG−1_3
AACAAAGCACAGAGCA−1_3
ATAGGCTGTCCAATCA−1_3
TTACGTTGTCAGATTC−1_3
ATACCGATCCCTAGGG−1_3
TGAGCATTCGCGCTGA−1_3
AAACGCTGTTGGGATG−1_3
TGAGGTTGTGGGATTG−1_3
CAAGCTAGTCCCTCAT−1_3
GAGGGATTCAACGAGG−1_3
CATCAAGCACGCGGTT−1_3
GTCCTCAGTATAGGAT−1_3
TGCATGATCTGAGCAT−1_3
TTACGTTAGCTTTGTG−1_3
CATGCAAAGCTGACTT−1_3
GGTGTTAGTCCCAAAT−1_3
GGGCTACTCAGTCAGT−1_3
TGAGCATGTATGGGAC−1_3
CGTGCTTTCAAGAGGC−1_3
CACCAAATCATCGACA−1_3
CCGCAAGCAAGGTCAG−1_3
TTTGACTAGGGTCTTT−1_3
AGAAATGTCGACGCTG−1_3
TACTTACGTTTGTTCT−1_3
TGATCTTCACTCCTGT−1_3
ATTTCACCATTAGGCT−1_3
TCTCCGACACAGCGCT−1_3
CCCAACTAGGATACCG−1_3
ATTTACCTCACCACAA−1_3
CATTCCGGTCCGGCAT−1_3
CTCCAACAGCGCTTCG−1_3
AGTAGTCGTGCTGATT−1_3
GGGCCATGTAGAATGT−1_3
CCGATCTGTCTTACAG−1_3
GTAGAAACAGACGATG−1_3
TCGGGACGTAGCGTAG−1_3
ATCGTGACAGCAGGAT−1_3
GTCCACTCAGGAGACT−1_3
GGTAATCGTATGAAGT−1_3
ATCGATGCATGTTACG−1_3
GACTTCCGTACGGCAA−1_3
TCGATTTCAATCTCGA−1_3
AGTGTTGGTTTGTTCT−1_3
TTCCGTGTCAGCATTG−1_3
CCACGTTAGTGGTTGG−1_3
AATTTCCCACGCCAGT−1_3
GTAGAGGTCCAGGACC−1_3
GACTCTCCACTGTGTA−1_3
GTAGGAGCATGGGTCC−1_3
TACGGTAAGAAGCGCT−1_3
CCTATCGCATCACCAA−1_3
TGTACAGTCCTGGGAC−1_3
CGACAGCCAAACCGGA−1_3
GTTACAGGTGTGTGTT−1_3
GATTGGTGTCCAGGTC−1_3
GTGCAGCCAGATCATC−1_3
TATATCCTCTAGATCG−1_3
ACTTAGGGTTTGGCTA−1_3
GTAGAAATCGGAATTC−1_3
TGCGGCACACCTGCAG−1_3
ATCTCTAAGTAGGATT−1_3
CATCGGGGTTCTGAGT−1_3
CTGGTCTGTATGAGCG−1_3
TTAGTCTCATGACCCG−1_3
TTGTTTGTCGTCGCTT−1_3
GAGTCTATCAGCTAGT−1_3
CATGGTATCGTAGGGA−1_3
GATGTTGTCTCGACCT−1_3
AGAACCTTCGCCACTT−1_3
CTACCTGTCTCCTACG−1_3
GATGGAGAGTAGGTTA−1_3
GGGCCATGTCAGTCCG−1_3
CATGCAATCCTGTACC−1_3
TCATGAGGTACGATCT−1_3
TAATCTCTCTGCATAG−1_3
CCACCATTCCTTCAGC−1_3
TCCTTTCAGCCAAGCA−1_3
GACTTCCTCGCTCATC−1_3
GTCCACTTCTACTATC−1_3
AGCTCAAAGGTATAGT−1_3
TGTTCATAGCAGGCAT−1_3
TCCTGCATCTGAGTCA−1_3
ACATGCAAGTCATAGA−1_3
CTCATGCTCGCTTTAT−1_3
AACACACAGAAGCCTG−1_3
GACGCTGCATCGCTCT−1_3
AAAGGTATCGCAATGT−1_3
GGAAGTGCATGACAGG−1_3
TGCGGCAGTGACACGA−1_3
CTACTATAGGAGTCTG−1_3
GCCTGTTTCATGGCCG−1_3
GCGTGCAGTTGCATGT−1_3
TCGCAGGGTCTGTTAG−1_3
ACCGTTCCATGATAGA−1_3
GGAGAACAGCGCTGAA−1_3
CAAGGGACATCTGGGC−1_3
CATCGCTTCCGCTAGG−1_3
GCTACAACAGTAGAAT−1_3
ACAGAAAGTTGGACTT−1_3
TCATTGTCAGTTTCAG−1_3
TCTCAGCTCGACTCCT−1_3
TTACAGGGTTGCTCCT−1_3
TCATCCGGTATTGACC−1_3
ATTCTACAGACCATTC−1_3
CCGGTAGTCTCTCCGA−1_3
CACTAAGCACGCTTAA−1_3
GTTGTCCGTCCGGTCA−1_3
TTGTGGAGTTCGAACT−1_3
TCAGGTAGTTTCTATC−1_3
TGTTACTTCAGCCTTC−1_3
GACTCAATCATACGGT−1_3
GTTGTGAAGTCAGGGT−1_3
AGGATCTAGCGTATAA−1_3
GACATCACAATGGCAG−1_3
TGGAGGAAGTTGAAAC−1_3
ATGATCGTCCTTATAC−1_3
CTGGTCTGTTCAAGTC−1_3
AATGGAAGTGACGTCC−1_3
AGGACGATCTTGAACG−1_3
GCCAACGCAACCAGAG−1_3
TGGGCTGTCGCAGTGC−1_3
CTAGGTACAGTCTTCC−1_3
GGGATCCAGGCAATGC−1_3
TGATTCTAGGTTCTAC−1_3
CTGCCATCACTCACTC−1_3
TTCTGTATCATGTCTT−1_3
ATCGGATAGAGGGTAA−1_3
ATCTCTAAGCTCATAC−1_3
CGAAGTTAGTGCTCAT−1_3
TTCGGTCCAAATCAGA−1_3
TCTGGCTGTACTCAAC−1_3
CTTTCAAAGCTTAGTC−1_3
GCAGCTGGTATCTCTT−1_3
AAGGTAATCGGACTGC−1_3
CATCCACGTGACTATC−1_3
TAATTCCAGGACAGTC−1_3
AAGCCATTCGGCCAAC−1_3
AGTACCACAACGGCTC−1_3
CCTCTAGCAGTCTTCC−1_3
GTGCAGCTCTAGCCTC−1_3
CAGTGCGAGGTTCTTG−1_3
CCACCATAGGCAGGGA−1_3
ACCCTCACAATAACCC−1_3
ACTATTCAGTCTACCA−1_3
ACTGTCCGTAACATCC−1_3
CATGCTCCACATAGCT−1_3
TCAAGACCAATCTCTT−1_3
TCCACCAAGTTCCAGT−1_3
CAACCAAGTCGCCACA−1_3
AGTCAACGTCGAACGA−1_3
TACTTACTCAAACCCA−1_3
CATGCGGAGCCTCTGG−1_3
TCACGCTTCAAGCTTG−1_3
GACACGCGTCCAGCCA−1_3
TACGTCCGTATCGTGT−1_3
GTTTACTAGGTACAAT−1_3
CATTCATGTAAGATTG−1_3
TTGACCCTCCATCAGA−1_3
CCACGAGAGCCATCCG−1_3
CTTAGGATCTCATGCC−1_3
TTGCGTCCAACTTGCA−1_3
CTTACCGAGAGGCGTT−1_3
GGTAGAGGTGTTTGCA−1_3
GTAGGAGGTCCTCCAT−1_3
GGTGTCGTCAGACATC−1_3
GACCCAGGTTATTCCT−1_3
AGTCACAAGCTACTGT−1_3
ATTTCACTCACACCGG−1_3
EVTlist
TACCTGCCATGCCGCA−1_3
GGTGTCGTCAGACATC−1_3
GACCCAGGTTATTCCT−1_3
AGTCACAAGCTACTGT−1_3
ATTTCACTCACACCGG−1_3
CCACGAGAGCCATCCG−1_3
CTTAGGATCTCATGCC−1_3
CATTCATGTAAGATTG−1_3
CCTCTAGCAGTCTTCC−1_3
AAGGTAATCGGACTGC−1_3
CATCCACGTGACTATC−1_3
TTCGGTCCAAATCAGA−1_3
TCTGGCTGTACTCAAC−1_3
ATCGGATAGAGGGTAA−1_3
GCCAACGCAACCAGAG−1_3
TGGGCTGTCGCAGTGC−1_3
TTTGACTAGCGTGTCC−1_3
GCCAACGTCACACCCT−1_3
ATCCACCTCTCAGTCC−1_3
CATCGCTTCCGCTAGG−1_3
AAAGGTATCGCAATGT−1_3
GTAGAGGTCCAGGACC−1_3
TTGACCCAGCTGTCCG−1_3
AATGAAGGTAGACGTG−1_3
CACAACAAGTAAACAC−1_3
ATTCAGGGTCTCGACG−1_3
GGGATGAGTACGGCAA−1_3
ATGACCATCCGTGTCT−1_3
ACGTCCTAGTTGTAGA−1_3
GAGCTGCTCTTTGCAT−1_3
TGGAGAGAGTATAACG−1_3
subs1 <- which(mammal.combined$Genotype2%in%c("NC2_2") & mammal.combined$Cells%in%c("CTB_d14","STB_d14","EVT_d14","putSTB") )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList2,]
Xp <- na.omit(Xp)[,subs1]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs1])), Batch =  factor((mammal.combined$Genotype2[subs1])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm2b<-pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 20,  filename = paste(saveext,"/DimRed/Tb_Batch2B_MinList",".pdf",sep="") ,width=20,height=4)
hc2b <- cutree(pm2b$tree_col,20)
saveRDS(hc2b,file=paste(saveext,'TbCl_batch2b.rds',sep="") )



subs1 <- which(mammal.combined$ID3%in%c("C3","C4") & mammal.combined$Cells%in%c("CTB_d14","STB_d14","EVT_d14","putSTB") )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList2,]
Xp <- na.omit(Xp)[,subs1]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs1])), Batch =  factor((mammal.combined$Genotype2[subs1])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm3<-pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  cutree_cols = 20, filename = paste(saveext,"/DimRed/Tb_Batch3_MinList",".pdf",sep="") ,width=10,height=4)
hc3 <- cutree(pm3$tree_col,20)
saveRDS(hc3,file=paste(saveext,'TbCl_batch3.rds',sep="") )

MiniList2b <- c("TACSTD2","HLA-C","KRT19","KRT7",'DNMT3A',"SNAI1","JAM3","CTNNAL1","MKI67","ENPEP","OTX2","ELF5","FN1","FSTL1","SPARC","VIM","CDH1","LAMC1","TP63","CGA","CGB5","CGB8","DPPA3",'CCKBR',"GATA3","OVOL1","TFAP2C","GATA2","SDC1","TFAP2A","TEAD3","NOTUM","ASCL2","HLA-G","FSTL3")


EVTlist <- c(
"ATTCCTACAACAAGAT-1_5",
"TTGGGCGCAGTCGAGA-1_5",
"GACGTTATCAGGTAAA-1_5",
"GTTCGCTCAGGCATTT-1_5",
"ACTTTCATCACACGAT-1_5",
"AGACAAACATGTCGTA-1_5",
"ACGTAGTAGTCTAGCT-1_5",
"GTAAGTCAGCCTCCAG-1_5",
"GCAGCTGCAGTCTCTC-1_5",
"TTCATGTAGCACTTTG-1_5",
"AAAGGGCGTTGCACGC-1_4",
"GCGGAAACATGACGAG-1_5",
"ATTCCTAAGGAAAGGT-1_5",
"AAGACTCCACTTCAAG-1_5",
"ACACAGTTCCGCTTAC-1_5",
"GGTAATCTCCGTGTCT-1_5",
"GTGAGGAGTGGTAATA-1_5",
"TGTGATGCATTCCTCG-1_5")

CTBlist <- c("CTCCTCCGTCAGCTTA-1_5",
"GAGGCCTCACAACGCC-1_5",
"ACCACAACAATTGCAC-1_5",
"CTCCACACAAGCCTGC-1_5",
"TCGACCTCAATGTCTG-1_5",
"CACTGGGTCTCCCTAG-1_5",
"CTGCCTACAACAGTGG-1_5")

STBlist <- c(
"TATCGCCAGGGCGAGA-1_5",
"TCATTGTCATTGCCGG-1_5",
"GATGGAGCATCTTAGG-1_5",
"TAATTCCAGTTGCGCC-1_5",
"ATCGTAGCAGCTTTCC-1_5",
"CCTTCAGGTTAAACCC-1_5",
"GCAGCCAAGCGGTAAC-1_5",
"AGCTACACACTCTGCT-1_5",
"CGTGATACAGATTCGT-1_5",
"GACGTTATCCGTGGGT-1_5",
"CGAGTTAAGCAGGGAG-1_5",
"ATCACAGAGGGTTAAT-1_5",
"TACTGCCGTGCTATTG-1_5",
"CGCATGGTCCCTCTCC-1_4",
"GCTTTCGGTGGACTGA-1_4",
"GTTCTATCATGCCATA-1_4",
"ACATTTCCAATATCCG-1_5",
"CATGCAAAGAAGCTGC-1_5",
"GGTAGAGGTATGGAGC-1_5",
"GTCTCACAGCACGATG-1_5",
"TGCATCCCAAACTAAG-1_5",
"CTTCTCTCACTCAGAT-1_5",
"ATGGGTTCAGACCGCT-1_5",
"TGGTTAGTCAAAGACA-1_5",
"CGGAGAACAGATGCGA-1_5",
"CGTTGGGAGGGAGGGT-1_5",
"TTCATTGAGGTTCATC-1_5",
"CAATTTCTCGAGCCTG-1_5",
"GCTGAATAGAAGTCAT-1_5",
"CTGCATCTCTCACCCA-1_5",
"AAGTGAAAGTACAGAT-1_5",
"ATCTCTAGTCTACAGT-1_5",
"TCGCACTTCATCACCC-1_5",
"AGTAACCAGGTAGACC-1_5",
"GGTTAACAGGTATAGT-1_5")

STB2 <- c(
"CCTGCATAGTTACGAA-1_5",
"CCCATTGGTAACGGTG-1_5",
"GGTGTCGAGTGTTCCA-1_5",
"ATAGGCTAGCCTCATA-1_5",
"ATCCCTGAGGATATAC-1_5",
"AATCGACCAGTAACCT-1_5",
"AGTTCCCTCATTTCCA-1_5",
"AGAAGCGGTCCTACAA-1_5",
"TAACTTCCATCGAAGG-1_5",
"TTCAGGAAGGTATAGT-1_4",
"GACCTTCTCATATGGC-1_5",
"TCATGAGTCGAACTCA-1_5",
"TCTTAGTTCGGAACTT-1_5",
"ATCCTATTCAGCATTG-1_5",
"CTCAGTCAGTAGCAAT-1_5",
"CTGTCGTGTGAGACGT-1_5",
"TCTGTCGAGGTTCTTG-1_5",
"TTCCGTGAGTGTACAA-1_5",
"GTAACACCAAGCTACT-1_5",
"TTCTTGACAGTATACC-1_5")

D3sub <- mammal.combined
Idents(D3sub) <- D3sub$ID3
D3sub <- subset(D3sub,idents =c("C3","C4"))
uL <- as.character(D3sub$CorrectLabel)
uL[1:length(uL)] <- "other"
Idents(D3sub) <- uL

Idents(D3sub,cells=EVTlist) <- "EVT"
Idents(D3sub,cells=CTBlist) <- "CTB1"
Idents(D3sub,cells=STB2) <- "STB"
Idents(D3sub,cells=STBlist) <- "STB"

D <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/Placenta/Placentaharm.rds")
Idents(D) <- D$ID2
#subs1 <- which(D$ID2%in%c("VCT","SCT","EVT"))
AvExp2 <- AverageExpression(D)
Expr2 <- AvExp2$RNA[MiniList2b,c("VCT","SCT","EVT")]
#X <- (GetAssayData(D,assay = "RNA")) 
#Xh <- t(scale(t(X)))
#Xp <- Xh[MiniList2,]
#Xp <- na.omit(Xp)[,subs1]
Expr2 <- t(scale(t(Expr2)))
#Expr2 <- na.omit(Expr2)


AvExp <- AverageExpression(D3sub)
Expr <- AvExp$RNA[MiniList2b,]
Expr <- t(scale(t(Expr)))
#Expr <- na.omit(Expr)

Expr <- Expr[,c("CTB1","STB","EVT")]

Exps <- cbind(Expr,Expr2)

redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm6<-pheatmap(Exps,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,  cluster_rows=TRUE, cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Tb_Batch3_MinListPB",".pdf",sep="") ,width=4,height=4)

EVTlist
ATCGTGAGTTATCTTC−1_6
CAGCCAGTCAGTGTTG−1_6
TTGGGATAGACCTGGA−1_6
TACCCGTTCTACCCAC−1_6
ACGCACGGTAGCGCCT−1_6
CTCCCAATCGTCAAAC−1_6
CAGCGTGAGTGGCGAT−1_6
GTTGCTCAGTAGGGTC−1_6
GCTCAAATCTCACGAA−1_6
GGCGTCAAGCTCGAAG−1_6
ACTGATGGTGCAATGG−1_6
TTGGGTAAGCTCGCAC−1_6
CGGTCAGTCGCTACAA−1_6
CTCCCAAAGGACATCG−1_6
CTCCCAAGTGGATCGA−1_6
CAACGATAGAGGATGA−1_6
CCAAGCGTCATTCACT−1_6
ACCTACCCAATCGCCG−1_6
CCTCATGTCTCCTGTG−1_6
CAACAGTCATAAGCGG−1_6
CTACTATGTAAGCGGT−1_6

STB
TCCTCCCCAAATGGTA−1_6
CTACCCACACAAAGTA−1_6
ATACCGACATCTAACG−1_6
CCGCAAGGTGCACAAG−1_6
AGGCTGCTCGGACCAC−1_6
CATTCTATCATTTCGT−1_6
CCAAGCGTCATTCACT−1_6
CATACTTAGACCAAGC−1_6
TCCGGGATCGCCTCTA−1_6
TTCACGCTCGACCAAT−1_6
AGCATCATCCACGTAA−1_6
GTTGTCCCACTTGAAC−1_6
TAGCACATCAGGCGAA−1_6
GCAGGCTTCTCTCAAT−1_6
TAATCTCGTGGGTATG−1_6
GGGTTATAGGAGATAG−1_6
ACTTTGTAGCCGTTGC−1_6
AGGACGAGTCAGACTT−1_6
AAGTGAATCCGCTAGG−1_6
ACTTTCAGTTAACAGA−1_6
ACCTACCCAATCGCCG−1_6
CCTCATGTCTCCTGTG−1_6
CAACAGTCATAAGCGG−1_6
CTACTATGTAAGCGGT−1_6
GTTGAACCACTGTGAT−1_6
GGCTTGGTCAGACTGT−1_6
GTGGAGAAGGCTCTCG−1_6
GTTAGACTCTTTCGAT−1_6
CGAGTTAGTCGGCCTA−1_6
GAACACTGTTCCGCAG−1_6
GGTAACTGTTTGGGTT−1_6
CGTTCTGGTTCTTGCC−1_6
GACTGATCAGCAGGAT−1_6

CTB
TTCCTCTTCAGTGTTG−1_6
CAACGATAGAGGATGA−1_6
TGGAACTCATCTATCT−1_6
TACACCCGTGACCGAA−1_6
TCTTAGTTCGTGGACC−1_6
TCGGATATCACTCCGT−1_6
ATCTTCAGTTTAGAGA−1_6
CGTAGTAAGAATCGAT−1_6
TTAGTCTAGGAACGCT−1_6
TTCCTCTTCAGTGTTG−1_6
GGGTGAATCTCGGGAC−1_6
TGAGACTGTGCCCTTT−1_6
TTGAGTGAGTATGCAA−1_

CTB to STB2
CTCAGGGCAGGCACTC−1_6
CCGATGGCAATTGTGC−1_6
GATCATGGTGTCTCCT−1_6
AGGGAGTTCCGATTAG−1_6
CTCAGGGTCTGGCCAG−1_6
TATCTGTGTTGCGAAG−1_6
CTGCATCCATTCCTCG−1_6
GGTTAACTCGATACTG−1_6
AGGGTTTGTTCAAGGG−1_6
TTTACGTCAGCTATTG−1_6
CGGGTCAAGAACGCGT−1_6
TGGGAGACAGGATGAC−1_6
TTGCTGCAGCTAAATG−1_6
GGAGGATAGTCGTCTA−1_6
AAGGTAATCTGTCAGA−1_6
GGTAATCTCGGTCTAA−1_6
CAGGCCATCAGGAACG−1_6
CTGAGCGAGACGACGT−1_6

subs1 <- which(mammal.combined$ID3%in%c("C6") & mammal.combined$Cells%in%c("CTB_d14","STB_d14","EVT_d14","putSTB") )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList2,]
Xp <- na.omit(Xp)[,subs1]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs1])), Batch =  factor((mammal.combined$Genotype2[subs1])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm6<-pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 20,  filename = paste(saveext,"/DimRed/Tb_Batch6_MinList",".pdf",sep="") ,width=10,height=4)
hc6 <- cutree(pm6$tree_col,20)

saveRDS(hc6,file=paste(saveext,'TbCl_batch6.rds',sep="") )

CTB
ACATCGACAGAGGGTT−1_7
CCTAAGATCCTTCGAC−1_7
TGCAGGCGTAGCGAGT−1_7
TTAGGGTGTACCATAC−1_7
ATATCCTAGAGAGCAA−1_7
CATGGTAAGGCGAAGG−1_7
GAGGGTACACTCATAG−1_7

STB2
ATCCCTGAGTGCACCC−1_7
TGCTTCGCATAGACTC−1_7
CGAATTGCAGTGCCTG−1_7
CGCCAGATCCTACTGC−1_7
AGGAATAAGGGATCTG−1_7
GTTCCGTGTGACACAG−1_7
TCGCACTCATGTTCGA−1_7
TACCGAAGTTCGGCGT−1_7
TGACGCGAGGCACCAA−1_7
GGTCACGCAGTCGTTA−1_7
GTAGAGGAGCTCCGAC−1_7
TGTAACGGTCTAGGTT−1_7
TTAATCCCACCCATAA−1_7

EVTlist
AGGGTGATCCGAGATT−1_7
GTAGTACCAGCGTTTA−1_7

subs1 <- which(mammal.combined$ID3%in%c("C7") & mammal.combined$Cells%in%c("CTB_d14","STB_d14","EVT_d14","putSTB") )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList2,]
Xp <- na.omit(Xp)[,subs1]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs1])), Batch =  factor((mammal.combined$Genotype2[subs1])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm7<-pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 5,  filename = paste(saveext,"/DimRed/Tb_Batch7_MinList",".pdf",sep="") ,width=6,height=4)
hc7 <- cutree(pm7$tree_col,5)
saveRDS(hc7,file=paste(saveext,'TbCl_batch7.rds',sep="") )

STB
ACAAGCTTCCACACCT−1_8
GAGGCCTAGATGACAT−1_8
CACGGGTGTGTAGCAG−1_8
TTGATGGCAGGCGATA−1_8
ACGTAGTCAGGTGTTT−1_8
AAATGGAAGTCCCAGC−1_8
GACCAATTCGTGCAGC−1_8
GCGGAAAAGTATTAGG−1_8
ACGTAACTCCATTTAC−1_8
AGTCATGCACCTGCAG−1_8
AAATGGATCAACCTCC−1_8
ACGATCATCCACCTCA−1_8
ATTTCACAGAATGTTG−1_8
TAACACGAGGTGAGAA−1_8
CTBlist
CACTAAGGTGAAGCTG−1_8
CCTCCTCTCAAGTGTC−1_8
CGAGGCTCAACCGCCA−1_8
CAGATACGTACTTCCC−1_8
TCAGCCTGTCATACCA−1_8
GATTTCTGTTCCTAGA−1_8
TCACTCGAGTTTGAGA−1_8
AACAAAGAGGCATTTC−1_8
CATCCGTGTTCCTACC−1_8
ATCACTTAGTAGCCAG−1_8
GGGTTTAAGACGAGCT−1_8
CGAGTTAGTATGCTTG−1_8
GGCAGTCAGGTCATCT−1_8
ATCACTTTCAACCTCC−1_8
ATCGATGAGGCCTGCT−1_8
CTGGCAGAGGGTATAT−1_8
TTTACTGGTCGCCTAG−1_8
CCACCATAGCCTGCCA−1_8
TTGCATTCACTAACCA−1_8
CGGACACAGCCTTTGA−1_8
TTTACCAAGACTTCAC−1_8
CAACGGCTCTTGGCTC−1_8
TTTGACTAGCGTCTGC−1_8
CCTCCAAGTAGCTTAC−1_8
CAGTTCCTCTACCAGA−1_8
ACAAGCTTCCCTTGGT−1_8
GGGCGTTCACCCTATC−1_8
GCTACCTAGGCCACCT−1_8
AGATAGATCTATCGGA−1_8
CATTCATGTCGATTCA−1_8
ATCCGTCTCCCTCGTA−1_8
GTAATCGTCTCCTGCA−1_8
GATGCTACAGCAAGAC−1_8
ATTCCTACATCTTCGC−1_8
CATTCCGTCCATCGTC−1_8
CTTCCTTAGACCAAAT−1_8

EVTlist
ACGTAGTCAGGTGTTT−1_8


write.table( as.data.frame(mammal.combined$Genotype), file=paste(saveext,"Genotype.csv",sep=""),sep=",",quote = FALSE)
write.table( as.data.frame(mammal.combined$Genotype2), file=paste(saveext,"Genotype2.csv",sep=""),sep=",",quote = FALSE)


subs1 <- which(mammal.combined$ID3%in%c("C8") & mammal.combined$Cells%in%c("CTB_d14","STB_d14","EVT_d14","putSTB") )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList2,]
Xp <- na.omit(Xp)[,subs1]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs1])), Batch =  factor((mammal.combined$Genotype2[subs1])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm8<-pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 10, filename = paste(saveext,"/DimRed/Tb_Batch8_MinList",".pdf",sep="") ,width=10,height=4)
hc8 <- cutree(pm8$tree_col,10)

saveRDS(hc8,file=paste(saveext,'TbCl_batch8.rds',sep="") )






subs2 <- which(mammal.combined$ID3%in%c("C1") & mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)

pm_1 <- pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 25, filename = paste(saveext,"/DimRed/EmD_Batch1_MinList",".pdf",sep="") ,width=10,height=6)
hc_1 <- cutree(pm_1$tree_col,25)
saveRDS(hc_1,file=paste(saveext,'EmCl_batch1.rds',sep="") )

subs2 <- which(mammal.combined$ID3%in%c("C2") & mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm_2 <- pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 25, filename = paste(saveext,"/DimRed/EmD_Batch2_MinList",".pdf",sep="") ,width=25,height=6)
#hc_2 <- cutree(pm_2$tree_col,25)

#saveRDS(hc_1,file=paste(saveext,'EmCl_batch1.rds',sep="") )

subs2 <- which(mammal.combined$Genotype2%in%c("NC2_1") & mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm_2a <- pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 25, filename = paste(saveext,"/DimRed/EmD_Batch2A_MinList",".pdf",sep="") ,width=25,height=6)
hc_2a <- cutree(pm_2a$tree_col,25)
saveRDS(hc_2a,file=paste(saveext,'EmCl_batch2a.rds',sep="") )


subs2 <- which(mammal.combined$Genotype2%in%c("NC2_2") & mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm_2b <- pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 6, filename = paste(saveext,"/DimRed/EmD_Batch2B_MinList",".pdf",sep="") ,width=10,height=6)
hc_2b <- cutree(pm_2b$tree_col,6)

saveRDS(hc_2b,file=paste(saveext,'EmCl_batch2b.rds',sep="") )

subs2 <- which(mammal.combined$ID3%in%c("C3","C4") & mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm_3 <- pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 4, filename = paste(saveext,"/DimRed/EmD_Batch3_MinList",".pdf",sep="") ,width=6,height=6)
hc_3 <- cutree(pm_3$tree_col,4)

saveRDS(hc_3,file=paste(saveext,'EmCl_batch3.rds',sep="") )



subs2 <- which(mammal.combined$ID3%in%c("C6") & mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm_6 <- pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 15, filename = paste(saveext,"/DimRed/EmD_Batch6_MinList",".pdf",sep="") ,width=14,height=6)
hc_6 <- cutree(pm_6$tree_col,15)

saveRDS(hc_6,file=paste(saveext,'EmCl_batch6.rds',sep="") )


subs2 <- which(mammal.combined$ID3%in%c("C7") & mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm_7 <- pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 5, filename = paste(saveext,"/DimRed/EmD_Batch7_MinList",".pdf",sep="") ,width=5,height=6)
hc_7 <- cutree(pm_7$tree_col,5)

saveRDS(hc_7,file=paste(saveext,'EmCl_batch7.rds',sep="") )


subs2 <- which(mammal.combined$ID3%in%c("C8") & mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pm_8 <- pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 4, filename = paste(saveext,"/DimRed/EmD_Batch8_MinList",".pdf",sep="") ,width=6,height=6)
hc_8 <- cutree(pm_8$tree_col,4)

saveRDS(hc_8,file=paste(saveext,'EmCl_batch8.rds',sep="") )


saveRDS(hc1,file=paste(saveext,"/Batch1_pheatmapclusterTb.rds",sep=""))
saveRDS(hc2a,file=paste(saveext,"/Batch2A_pheatmapclusterTb.rds",sep=""))
saveRDS(hc2b,file=paste(saveext,"/Batch2B_pheatmapclusterTb.rds",sep=""))
saveRDS(hc3,file=paste(saveext,"/Batch3_pheatmapclusterTb.rds",sep=""))
saveRDS(hc6,file=paste(saveext,"/Batch6_pheatmapclusterTb.rds",sep=""))
saveRDS(hc7,file=paste(saveext,"/Batch7_pheatmapclusterTb.rds",sep=""))
saveRDS(hc8,file=paste(saveext,"/Batch8_pheatmapclusterTb.rds",sep=""))

saveRDS(hc_1,file=paste(saveext,"/Batch1_pheatmapclusterembryonic.rds",sep=""))
saveRDS(hc_2a,file=paste(saveext,"/Batch2A_pheatmapclusterembryonic.rds",sep=""))
saveRDS(hc_2b,file=paste(saveext,"/Batch2B_pheatmapclusterembryonic.rds",sep=""))
saveRDS(hc_3,file=paste(saveext,"/Batch3_pheatmapclusterembryonic.rds",sep=""))
saveRDS(hc_6,file=paste(saveext,"/Batch6_pheatmapclusterembryonic.rds",sep=""))
saveRDS(hc_7,file=paste(saveext,"/Batch7_pheatmapclusterembryonic.rds",sep=""))
saveRDS(hc_8,file=paste(saveext,"/Batch8_pheatmapclusterembryonic.rds",sep=""))


Idents(Ds1) <- Ds1$CorrectLabel

Ds1[["percent.mt"]] <- PercentageFeatureSet(Ds1, pattern = "^MT-")
Ds1[["ribo.mt"]] <- PercentageFeatureSet(Ds1, pattern = "^MRPL")


#p1 <- -Ds1$PCl2 + pmax(Ds1$PCl1, Ds1$PCl3, Ds1$PCl4, na.rm = TRUE)
#Ds1$p1 <- p1
p <- VlnPlot(Ds1,feature="percent.mt")
ggsave(filename=paste(saveext,"/DimRed/Vln_C1_mito.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds1,feature="ribo.mt")
ggsave(filename=paste(saveext,"/DimRed/Vln_C1_ribo.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds1,feature="nFeature_RNA")
ggsave(filename=paste(saveext,"/DimRed/Vln_C1_nFeature.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds1,feature="nCount_RNA")
ggsave(filename=paste(saveext,"/DimRed/Vln_C1_nCount.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds1,feature="PCl1")
ggsave(filename=paste(saveext,"/DimRed/Vln_C1_PCl1.pdf",sep=""),width = 10, height = 8,p)


Idents(Ds2) <- Ds2$CorrectLabel

Ds2[["percent.mt"]] <- PercentageFeatureSet(Ds2, pattern = "^MT-")
Ds2[["ribo.mt"]] <- PercentageFeatureSet(Ds2, pattern = "^MRPL")
p <- VlnPlot(Ds2,feature="percent.mt")
ggsave(filename=paste(saveext,"/DimRed/Vln_C2_mito.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds2,feature="ribo.mt")
ggsave(filename=paste(saveext,"/DimRed/Vln_C2_ribo.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds2,feature="nFeature_RNA")
ggsave(filename=paste(saveext,"/DimRed/Vln_C2_nFeature.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds2,feature="nCount_RNA")
ggsave(filename=paste(saveext,"/DimRed/Vln_C2_nCount.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds2,feature="PCl1")
ggsave(filename=paste(saveext,"/DimRed/Vln_C2_PCl1.pdf",sep=""),width = 10, height = 8,p)

Idents(Ds3) <- Ds3$CorrectLabel
Ds3[["percent.mt"]] <- PercentageFeatureSet(Ds3, pattern = "^MT-")
Ds3[["ribo.mt"]] <- PercentageFeatureSet(Ds3, pattern = "^MRPL")
p <- VlnPlot(Ds3,feature="percent.mt")
ggsave(filename=paste(saveext,"/DimRed/Vln_C3_mito.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds3,feature="ribo.mt")
ggsave(filename=paste(saveext,"/DimRed/Vln_C3_ribo.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds3,feature="nFeature_RNA")
ggsave(filename=paste(saveext,"/DimRed/Vln_C3_nFeature.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds3,feature="nCount_RNA")
ggsave(filename=paste(saveext,"/DimRed/Vln_C3_nCount.pdf",sep=""),width = 10, height = 8,p)


Idents(Ds6) <- Ds6$CorrectLabel
Ds6[["percent.mt"]] <- PercentageFeatureSet(Ds6, pattern = "^MT-")
Ds6[["ribo.mt"]] <- PercentageFeatureSet(Ds6, pattern = "^MRPL")
p <- VlnPlot(Ds6,feature="percent.mt")
ggsave(filename=paste(saveext,"/DimRed/Vln_C6_mito.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds6,feature="ribo.mt")
ggsave(filename=paste(saveext,"/DimRed/Vln_C6_ribo.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds6,feature="nFeature_RNA")
ggsave(filename=paste(saveext,"/DimRed/Vln_C6_nFeature.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds6,feature="nCount_RNA")
ggsave(filename=paste(saveext,"/DimRed/Vln_C6_nCount.pdf",sep=""),width = 10, height = 8,p)



Idents(Ds7) <- Ds7$CorrectLabel
Ds7[["percent.mt"]] <- PercentageFeatureSet(Ds7, pattern = "^MT-")
Ds7[["ribo.mt"]] <- PercentageFeatureSet(Ds7, pattern = "^MRPL")
p <- VlnPlot(Ds7,feature="percent.mt")
ggsave(filename=paste(saveext,"/DimRed/Vln_C7_mito.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds7,feature="ribo.mt")
ggsave(filename=paste(saveext,"/DimRed/Vln_C7_ribo.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds7,feature="nFeature_RNA")
ggsave(filename=paste(saveext,"/DimRed/Vln_C7_nFeature.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds7,feature="nCount_RNA")
ggsave(filename=paste(saveext,"/DimRed/Vln_C7_nCount.pdf",sep=""),width = 10, height = 8,p)



Idents(Ds8) <- Ds8$CorrectLabel
Ds8[["percent.mt"]] <- PercentageFeatureSet(Ds8, pattern = "^MT-")
Ds8[["ribo.mt"]] <- PercentageFeatureSet(Ds8, pattern = "^MRPL")
p <- VlnPlot(Ds8,feature="percent.mt")
ggsave(filename=paste(saveext,"/DimRed/Vln_C8_mito.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds8,feature="ribo.mt")
ggsave(filename=paste(saveext,"/DimRed/Vln_C8_ribo.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds8,feature="nFeature_RNA")
ggsave(filename=paste(saveext,"/DimRed/Vln_C8_nFeature.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds8,feature="nCount_RNA")
ggsave(filename=paste(saveext,"/DimRed/Vln_C8_nCount.pdf",sep=""),width = 10, height = 8,p)



p <- VlnPlot(Ds1,feature="PCl2")
ggsave(filename=paste(saveext,"/DimRed/Vln_C1_PCl2.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds1,feature="PCl3")
ggsave(filename=paste(saveext,"/DimRed/Vln_C1_PCl3.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds1,feature="PCl4")
ggsave(filename=paste(saveext,"/DimRed/Vln_C1_PCl4.pdf",sep=""),width = 10, height = 8,p)


p <- VlnPlot(Ds2,feature="PCl1")
ggsave(filename=paste(saveext,"/DimRed/Vln_C2_PCl1.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds2,feature="PCl2")
ggsave(filename=paste(saveext,"/DimRed/Vln_C2_PCl2.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds2,feature="PCl3")
ggsave(filename=paste(saveext,"/DimRed/Vln_C2_PCl3.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds2,feature="PCl4")
ggsave(filename=paste(saveext,"/DimRed/Vln_C2_PCl4.pdf",sep=""),width = 10, height = 8,p)




p <- VlnPlot(Ds6,feature="PCl1")
ggsave(filename=paste(saveext,"/DimRed/Vln_C6_PCl1.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds6,feature="PCl2")
ggsave(filename=paste(saveext,"/DimRed/Vln_C6_PCl2.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds6,feature="PCl3")
ggsave(filename=paste(saveext,"/DimRed/Vln_C6_PCl3.pdf",sep=""),width = 10, height = 8,p)
p <- VlnPlot(Ds6,feature="PCl4")
ggsave(filename=paste(saveext,"/DimRed/Vln_C6_PCl4.pdf",sep=""),width = 10, height = 8,p)



#
GT <- D$Genotype[which(D$ID3=="C1")]
ID <- D$CorrectLabel[which(D$ID3=="C1")]
L1 <- D$PCl1[which(D$ID3=="C1")] #Epi?
L2 <- D$PCl2[which(D$ID3=="C1")] #Embryonic
L3 <- D$PCl3[which(D$ID3=="C1")] #Stromal
L4 <- D$PCl4[which(D$ID3=="C1")] #Epi 3

mean(L1[which(GT%in%c("Epithlial_G4") )])
mean(L4[which(GT%in%c("Epithelial_G3") )])
mean(L3[which(GT%in%c("Embryonic_G1") )])
mean(L2[which(GT%in%c("Stromal_G2") )])


Ds1$Epith1 <- Ds1$PCl1 + Ds1$PCl2
Ds1$Epith2 <- Ds1$PCl3 + Ds1$PCl2
Ds1$Str <- Ds1$PCl4 + Ds1$PCl2
p<-FeaturePlot(Ds1, features = "Epith1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_Epith1.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "Epith2", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_Epith2.pdf",sep=""),width = 10, height = 8,p)
p<-FeaturePlot(Ds1, features = "Str", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1Embryonic_Str.pdf",sep=""),width = 10, height = 8,p)



Idents(D) <- D$ID3
Ds1Full <- subset(D,idents="C1")
Idents(Ds1Full) <- Ds1Full$CorrectLabel
Ds1Full$Delta1 <- Ds1Full$PCl1 - Ds1Full$PCl2
Ds1Full$Delta2 <- Ds1Full$PCl1 - Ds1Full$PCl3
Ds1Full$Delta3 <- Ds1Full$PCl1 - Ds1Full$PCl4

Ds1Full <- FindVariableFeatures(Ds1Full, selection.method = "vst", nfeatures = 4000)
Ds1Full <- ScaleData(Ds1Full)
Ds1Full <- RunPCA(Ds1Full, npcs = 20, verbose = FALSE)
Ds1Full <- RunUMAP(Ds1Full, reduction = "pca", dims = 1:20)
Ds1Full <- FindNeighbors(Ds1Full, reduction = "pca", dims = 1:2)
p<-DimPlot(Ds1Full,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1EmbryoandMat.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p<-DimPlot(Ds1Full,  pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_C1EmbryoandMat.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p<-FeaturePlot(Ds1Full, features = "Delta1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, cols = c("red",  "blue"), min.cutoff = -500, max.cutoff = 500) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1EmbryoandMat_Delta1.pdf",sep=""),width = 10, height = 8,p)

p<-FeaturePlot(Ds1Full, features = "Delta2", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, cols = c("red",  "blue"), min.cutoff = -500, max.cutoff = 500) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1EmbryoandMat_Delta2.pdf",sep=""),width = 10, height = 8,p)


p<-FeaturePlot(Ds1Full, features = "Delta3", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE, cols = c("red", "blue"), blend.threshold=0, min.cutoff = -500, max.cutoff = 500) 
ggsave(filename=paste(saveext,"/DimRed/PCA_C1EmbryoandMat_Delta3.pdf",sep=""),width = 10, height = 8,p)



A <- rbind(Ds1Full$PCl1,Ds1Full$PCl2,Ds1Full$PCl3,Ds1Full$PCl4)
annotation_col = data.frame(Stage = factor(droplevels(Ds1Full$Cells)), Batch =  factor((Ds1Full$Genotype)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
#my_breaks <- seq(-2000, -1000, length.out = 50)
pheatmap(A,color =  redblue1(50), fontsize = 4,  border_color = NA,annotation_col = annotation_col, scale = "column", cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch1_ClusterGenomeLikelihoods",".pdf",sep="") ,width=10,height=3)



Idents(D) <- D$ID3
Ds3Full <- subset(D,idents=c("C4","C3"))
Idents(Ds3Full) <- Ds3Full$CorrectLabel
Ds3Full$Delta1 <- Ds3Full$PCl1 - Ds3Full$PCl2
Ds3Full$Delta2 <- Ds3Full$PCl1 - Ds3Full$PCl3
Ds3Full$Delta3 <- Ds3Full$PCl1 - Ds3Full$PCl4

Ds3Full <- FindVariableFeatures(Ds3Full, selection.method = "vst", nfeatures = 4000)
Ds3Full <- ScaleData(Ds3Full)
Ds3Full <- RunPCA(Ds3Full, npcs = 20, verbose = FALSE)
Ds3Full <- RunUMAP(Ds3Full, reduction = "pca", dims = 1:20)
Ds3Full <- FindNeighbors(Ds3Full, reduction = "pca", dims = 1:2)

A <- rbind(Ds3Full$PCl1,Ds3Full$PCl2,Ds3Full$PCl3,Ds3Full$PCl4,Ds3Full$PCl5,Ds3Full$PCl6)
annotation_col = data.frame(Stage = factor(droplevels(Ds3Full$Cells)), Batch =  factor((Ds3Full$Genotype)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
#my_breaks <- seq(-2000, -1000, length.out = 50)
pheatmap(A,color =  redblue1(50), fontsize = 4,  border_color = NA,annotation_col = annotation_col, scale = "column", cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch3_ClusterGenomeLikelihoods",".pdf",sep="") ,width=10,height=3)



Idents(D) <- D$ID3
Ds6Full <- subset(D,idents=c("C6"))
Idents(Ds6Full) <- Ds6Full$CorrectLabel
Ds6Full$Delta1 <- Ds6Full$PCl1 - Ds6Full$PCl2
Ds6Full$Delta2 <- Ds6Full$PCl1 - Ds6Full$PCl3
Ds6Full$Delta3 <- Ds6Full$PCl1 - Ds6Full$PCl4

Ds6Full <- FindVariableFeatures(Ds6Full, selection.method = "vst", nfeatures = 4000)
Ds6Full <- ScaleData(Ds6Full)
Ds6Full <- RunPCA(Ds6Full, npcs = 20, verbose = FALSE)
Ds6Full <- RunUMAP(Ds6Full, reduction = "pca", dims = 1:20)
Ds6Full <- FindNeighbors(Ds6Full, reduction = "pca", dims = 1:2)

A <- rbind(Ds6Full$PCl1,Ds6Full$PCl2,Ds6Full$PCl3,Ds6Full$PCl4)
annotation_col = data.frame(Stage = factor(droplevels(Ds6Full$Cells)), Batch =  factor((Ds6Full$Genotype)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
#my_breaks <- seq(-2000, -1000, length.out = 50)
pheatmap(A,color =  redblue1(50), fontsize = 4,  border_color = NA,annotation_col = annotation_col, scale = "column", cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch6_ClusterGenomeLikelihoods",".pdf",sep="") ,width=10,height=3)



Idents(D) <- D$ID3
Ds7Full <- subset(D,idents=c("C7"))
Idents(Ds7Full) <- Ds7Full$CorrectLabel
Ds7Full$Delta1 <- Ds7Full$PCl1 - Ds7Full$PCl2
Ds7Full$Delta2 <- Ds7Full$PCl1 - Ds7Full$PCl3
Ds7Full$Delta3 <- Ds7Full$PCl1 - Ds7Full$PCl4

Ds7Full <- FindVariableFeatures(Ds7Full, selection.method = "vst", nfeatures = 4000)
Ds7Full <- ScaleData(Ds7Full)
Ds7Full <- RunPCA(Ds7Full, npcs = 20, verbose = FALSE)
Ds7Full <- RunUMAP(Ds7Full, reduction = "pca", dims = 1:20)
Ds7Full <- FindNeighbors(Ds7Full, reduction = "pca", dims = 1:2)


A <- rbind(Ds7Full$PCl1,Ds7Full$PCl2,Ds7Full$PCl3,Ds7Full$PCl4)
annotation_col = data.frame(Stage = factor(droplevels(Ds7Full$Cells)), Batch =  factor((Ds7Full$Genotype)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
#my_breaks <- seq(-2000, -1000, length.out = 50)
pheatmap(A,color =  redblue1(50), fontsize = 4,  border_color = NA,annotation_col = annotation_col, scale = "column", cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch7_ClusterGenomeLikelihoods",".pdf",sep="") ,width=10,height=3)



Idents(D) <- D$ID3
Ds8Full <- subset(D,idents=c("C8"))
Idents(Ds8Full) <- Ds8Full$CorrectLabel
Ds8Full$Delta1 <- Ds8Full$PCl1 - Ds8Full$PCl2
Ds8Full$Delta2 <- Ds8Full$PCl1 - Ds8Full$PCl3
Ds8Full$Delta3 <- Ds8Full$PCl1 - Ds8Full$PCl4


Ds8Full <- FindVariableFeatures(Ds8Full, selection.method = "vst", nfeatures = 4000)
Ds8Full <- ScaleData(Ds8Full)
Ds8Full <- RunPCA(Ds8Full, npcs = 20, verbose = FALSE)
Ds8Full <- RunUMAP(Ds8Full, reduction = "pca", dims = 1:20)
Ds8Full <- FindNeighbors(Ds8Full, reduction = "pca", dims = 1:2)


A <- rbind(Ds8Full$PCl2,Ds8Full$PCl3,Ds8Full$PCl4,Ds8Full$PCl5)
annotation_col = data.frame(Stage = factor(droplevels(Ds8Full$Cells)), Batch =  factor((Ds8Full$Genotype)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
#my_breaks <- seq(-2000, -1000, length.out = 50)
pheatmap(A,color =  redblue1(50), fontsize = 4,  border_color = NA,annotation_col = annotation_col, scale = "column", cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch8_ClusterGenomeLikelihoods",".pdf",sep="") ,width=10,height=3)


# Stro-1,37 CD90,37 CD133,37 CD146,38 and SUSD2
#ectonucleoside triphosphate diphosphohydrolase-2 (NTPDase2),41 leucine-rich repeat containing G-protein-coupled receptor 5 (LGR5),42 platelet-derived growth factor receptor beta (PDGFR-ß)/CD146,39 side population (SP),1,43 and sushi domain-containing-2 (SUSD2)
# (OCT-4, GDF3, DNMT3B, NANOG, and GABR3)

#N-cadherin+/SSEA-1-

#Epithelial like stem cell OCT4, GDF3, DMRT3B, NANOG, GABR3, NCAD, SSEA-1
#Stromal stem like: Stro-1, CD90, CD133, CD146, SUSD2
# SOX9/SSEA-1

#WNT4 is expressed in the basal lamina and glandular/luminal epithelium of porcine endometrial tissues
#WNT5A is only expressed in the luminal epithelium.





#Now do list enrichment


VE <- FindMarkers(marmoset_dataInVivo2, ident.1 = "VE", ident.2 = c("EmDisc","ExMes","Tb","SYS","Stalk","PGC"), only.pos = TRUE, test.use = "MAST")
SYS <- FindMarkers(marmoset_dataInVivo2, ident.1 = "SYS", ident.2 = c("EmDisc","ExMes","Tb","VE","Stalk","PGC"), only.pos = TRUE, test.use = "MAST")
ExMes <- FindMarkers(marmoset_dataInVivo2, ident.1 = "ExMes", ident.2 = c("EmDisc","VE","Tb","SYS","Stalk","PGC"), only.pos = TRUE, test.use = "MAST")
EmDisc <- FindMarkers(marmoset_dataInVivo2, ident.1 = "EmDisc", ident.2 = c("VE","ExMes","Tb","SYS","Stalk","PGC"), only.pos = TRUE, test.use = "MAST")
Tb <- FindMarkers(marmoset_dataInVivo2, ident.1 = "Tb", ident.2 = c("EmDisc","ExMes","VE","SYS","Stalk","PGC"), only.pos = TRUE, test.use = "MAST")
PGC <- FindMarkers(marmoset_dataInVivo2, ident.1 = "PGC", ident.2 = c("EmDisc","ExMes","Tb","SYS","Stalk","VE"), only.pos = TRUE, test.use = "MAST")
Stalk <- FindMarkers(marmoset_dataInVivo2, ident.1 = "Stalk", ident.2 = c("EmDisc","ExMes","Tb","SYS","VE","PGC"), only.pos = TRUE, test.use = "MAST")

Tr <- FindMarkers(marmoset_dataInVivo2, ident.1 = "Tb_CS3", ident.2 = c("Hyp_CS3","Epi_CS3"), only.pos = TRUE, test.use = "MAST")
Hyp <- FindMarkers(marmoset_dataInVivo2, ident.1 = "Hyp_CS3", ident.2 = c("Epi_CS3","Tb_CS3"), only.pos = TRUE, test.use = "MAST")
Epi <- FindMarkers(marmoset_dataInVivo2, ident.1 = "Epi_CS3", ident.2 = c("Hyp_CS3","Tb_CS3"), only.pos = TRUE, test.use = "MAST")

Ds1 <- AddModuleScore(Ds1, features = list(intersect(rownames(VE),rownames(Ds1))), name = "VE")
p1<-FeaturePlot(Ds1, features = "VE_1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_VE.pdf",sep=""),width = 10, height = 8,p1)


Ds1 <- AddModuleScore(Ds1, features = list(intersect(rownames(Tb),rownames(Ds1))), name = "Tb")
p1<-FeaturePlot(Ds1, features = "Tb1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_Tb.pdf",sep=""),width = 10, height = 8,p1)


Ds1 <- AddModuleScore(Ds1, features = list(intersect(rownames(SYS),rownames(Ds1))), name = "Tb")
p1<-FeaturePlot(Ds1, features = "SYS1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_SYS.pdf",sep=""),width = 10, height = 8,p1)



PrimedNaiveForm <- read_excel("/Users/christopherpenfold/Downloads/mmc3.xls")
PrimedNaive <- read_excel("/Users/christopherpenfold/Downloads/mmc2.xls")
AmnionMarkers <- read_excel("/Users/christopherpenfold/Downloads/AmnionMarkers.xlsx")
TbMarkers <- read_excel("/Users/christopherpenfold/Downloads/AnotatedAeDEupdated2.xlsx",sheet = 2) 
UberList <- c( intersect(PrimedNaiveForm$Gene[1:100], RN),
               intersect(PrimedNaive$Gene[1:100],  RN),
               intersect(PrimedNaiveForm$Gene[1:100],  RN),
               intersect(AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsTE")[1:100]], RN),
               intersect(AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E + hsTE")[1:100]], RN),
               intersect(AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E (+ hsTE)")[1:100]], RN),
               intersect(AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E (+ cyAME-L)")[1:100]], RN),
               intersect(AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E + cyAME-L")[1:100]], RN),
               intersect(AmnionMarkers$Gene[which(AmnionMarkers$marker=="cyAME-L")[1:100]], RN) )



Ds1 <- AddModuleScore(Ds1, features = list(intersect(PrimedNaive$Gene[which( PrimedNaive$logFC>log2(10) & PrimedNaive$FDR<0.05)], rownames(Ds1) )), name = "Primed")
p1<-FeaturePlot(Ds1, features = "Primed1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_Primed.pdf",sep=""),width = 10, height = 8,p1)

Ds1 <- AddModuleScore(Ds1, features = list(intersect(PrimedNaive$Gene[which( PrimedNaive$logFC<-log2(10) & PrimedNaive$FDR<0.05)], rownames(Ds1) )), name = "Naive")
p1<-FeaturePlot(Ds1, features = "Naive1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_Naive.pdf",sep=""),width = 10, height = 8,p1)

Ds1 <- AddModuleScore(Ds1, features = list(intersect(PrimedNaive$Gene[which( PrimedNaiveForm$vsNaive>log2(10) & PrimedNaiveForm$P.Value<0.05)], rownames(Ds1) )), name = "FormN")
p1<-FeaturePlot(Ds1, features = "FormN1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_FormvNaive.pdf",sep=""),width = 10, height = 8,p1)

Ds1 <- AddModuleScore(Ds1, features = list(intersect(PrimedNaive$Gene[which( PrimedNaiveForm$vsPrimed>log2(10) & PrimedNaiveForm$P.Value<0.05)], rownames(Ds1) )), name = "FormP")
p1<-FeaturePlot(Ds1, features = "FormP1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_FormvPrime.pdf",sep=""),width = 10, height = 8,p1)



Ds1 <- AddModuleScore(Ds1, features = list(intersect( AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E (+ hsTE)")], rownames(Ds1) )), name = "AmList")
p1<-FeaturePlot(Ds1, features = "AmList1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_Am.pdf",sep=""),width = 10, height = 8,p1)


Ds1 <- AddModuleScore(Ds1, features = list(intersect( AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E + cyAME-L")], rownames(Ds1) )), name = "LAmList")
p1<-FeaturePlot(Ds1, features = "LAmList1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_LAm.pdf",sep=""),width = 10, height = 8,p1)



Ds1 <- AddModuleScore(Ds1, features = list(intersect(rownames(PGC)[which(PGC$avg_log2FC>0.01 & PGC$p_val_adj<0.05)],rownames(Ds1))), name = "PGC")
p1<-FeaturePlot(Ds1, features = "PGC1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_PGC.pdf",sep=""),width = 10, height = 8,p1)




SCMmarker <- read_excel("/Users/christopherpenfold/Downloads/41586_2023_6604_MOESM3_ESM.xlsx")

SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM")]



Ds1 <- AddModuleScore(Ds1, features = list(intersect( SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM")] ,rownames(Ds1))), name = "ExEm")
p1<-FeaturePlot(Ds1, features = "ExEm1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_ExEm.pdf",sep=""),width = 10, height = 8,p1)


Ds1 <- AddModuleScore(Ds1, features = list(intersect( SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Epiblast")] ,rownames(Ds1))), name = "Epiblast")
p1<-FeaturePlot(Ds1, features = "Epiblast1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_Epiblast.pdf",sep=""),width = 10, height = 8,p1)



Ds1 <- AddModuleScore(Ds1, features = list(intersect( SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Yolk-Sac / Hypoblast")] ,rownames(Ds1))), name = "YS")
p1<-FeaturePlot(Ds1, features = "YS1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_YS.pdf",sep=""),width = 10, height = 8,p1)


Ds1 <- AddModuleScore(Ds1, features = list(intersect( SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Amnion")] ,rownames(Ds1))), name = "Amnion")
p1<-FeaturePlot(Ds1, features = "Amnion1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_Amnion.pdf",sep=""),width = 10, height = 8,p1)



Ds1 <- AddModuleScore(Ds1, features = list(intersect( SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Secondary Yolk-Sac")] ,rownames(Ds1))), name = "SeYS")
p1<-FeaturePlot(Ds1, features = "SeYS1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_SYS.pdf",sep=""),width = 10, height = 8,p1)




Ds1 <- AddModuleScore(Ds1, features = list(intersect( SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Trophoblast")] ,rownames(Ds1))), name = "Troph")
p1<-FeaturePlot(Ds1, features = "Troph1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_Troph.pdf",sep=""),width = 10, height = 8,p1)



Ds1 <- AddModuleScore(Ds1, features = list(intersect( SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Posterior Epiblast")] ,rownames(Ds1))), name = "PostEpi")
p1<-FeaturePlot(Ds1, features = "PostEpi1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_PostEpi.pdf",sep=""),width = 10, height = 8,p1)


Ds1 <- AddModuleScore(Ds1, features = list(intersect( SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Epiblast committed")] ,rownames(Ds1))), name = "ComEpi")
p1<-FeaturePlot(Ds1, features = "ComEpi1", pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Batch1_ComEpi.pdf",sep=""),width = 10, height = 8,p1)


#SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$p_val_adj<0.01)]
L1 <- c(
intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)]),
intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Epiblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)]),
intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Yolk-Sac / Hypoblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)]),
intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Posterior Epiblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)]),
intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Epiblast committed" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)]),
intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Secondary Yolk-Sac" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)]),
intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Amnion" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)]),
intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Trophoblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)]))


subs2 <- which(mammal.combined$ID3%in%c("C1")  )
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=FALSE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch1_List1",".pdf",sep="") ,width=20,height=40)


A <- rbind(Ds1$ExEm1,Ds1$SYS1,Ds1$SeYS1,Ds1$Amnion1,Ds1$Epiblast1,Ds1$PostEpi1,Ds1$ComEpi1,Ds1$Troph1)
annotation_col = data.frame(Stage = factor(droplevels(Ds1$Cells)), Batch =  factor((Ds1$Genotype)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
#my_breaks <- seq(-2000, -1000, length.out = 50)
pheatmap(A,color =  redblue1(50), fontsize = 4,  border_color = NA,labels_row = c("ExMes","VE","SYS","Am","Epi","PEpi","CEpi","Tb"), annotation_col = annotation_col, scale = "column", cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch1_AnoCluster",".pdf",sep="") ,width=10,height=3)



#SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$p_val_adj<0.01)]
subs2 <- which(mammal.combined$ID3%in%c("C1") & mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch1_ExMes",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Yolk-Sac / Hypoblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch1_Hypo",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Amnion" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch1_Am",".pdf",sep="") ,width=15,height=20)


L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Epiblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch1_Epi",".pdf",sep="") ,width=15,height=20)


L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Secondary Yolk-Sac" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch1_SYS",".pdf",sep="") ,width=15,height=20)


hsEAm_TEList1 <- intersect(rownames(D),AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E + hsTE")] )
hsEAm_TEList2 <- intersect(rownames(D),AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E (+ hsTE)")] )
hsEAm1 <- intersect(rownames(D),AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E (+ cyAME-L)")] )
hsEAm2 <- intersect(rownames(D),AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E + cyAME-L")] )
hsEAm3 <- intersect(rownames(D),AmnionMarkers$Gene[which(AmnionMarkers$marker=="cyAME-L")] )
hsLAm_TEList <- intersect(rownames(D),AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsTE + cyAME-L")] )


L1 <- hsEAm1
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch1_Am2",".pdf",sep="") ,width=15,height=20)










#SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$p_val_adj<0.01)]
subs2 <- which(mammal.combined$ID3%in%c("C3","C4") & mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch3_ExMes",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Yolk-Sac / Hypoblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch3_Hypo",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Amnion" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch3_Am",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Epiblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch3_Epi",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Secondary Yolk-Sac" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch3_SYS",".pdf",sep="") ,width=15,height=20)









#SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$p_val_adj<0.01)]
subs2 <- which(mammal.combined$ID3%in%c("C6") & mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch6_ExMes",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Yolk-Sac / Hypoblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch6_Hypo",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Amnion" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch6_Am",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Epiblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch6_Epi",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Secondary Yolk-Sac" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch6_SYS",".pdf",sep="") ,width=15,height=20)







#SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$p_val_adj<0.01)]
subs2 <- which(mammal.combined$ID3%in%c("C7") & mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch7_ExMes",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Yolk-Sac / Hypoblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch7_Hypo",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Amnion" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch7_Am",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Epiblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch7_Epi",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Secondary Yolk-Sac" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch7_SYS",".pdf",sep="") ,width=15,height=20)




#SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$p_val_adj<0.01)]
subs2 <- which(mammal.combined$ID3%in%c("C8") & mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch8_ExMes",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Yolk-Sac / Hypoblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch8_Hypo",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Amnion" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch8_Am",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Epiblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch8_Epi",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Secondary Yolk-Sac" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch8_SYS",".pdf",sep="") ,width=15,height=20)











#SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$p_val_adj<0.01)]
subs2 <- which(mammal.combined$ID3%in%c("C7") & mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch7_ExMes",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Yolk-Sac / Hypoblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch7_Hypo",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Amnion" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch7_Am",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Epiblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch7_Epi",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Secondary Yolk-Sac" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch7_SYS",".pdf",sep="") ,width=15,height=20)




#SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$p_val_adj<0.01)]
subs2 <- which(mammal.combined$Genotype2%in%c("NC2_1") & mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch2A_ExMes",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Yolk-Sac / Hypoblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch2A_Hypo",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Amnion" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch2A_Am",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Epiblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch2A_Epi",".pdf",sep="") ,width=15,height=20)

L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Secondary Yolk-Sac" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE,  filename = paste(saveext,"/DimRed/Batch2A_SYS",".pdf",sep="") ,width=15,height=20)


DELIST1 <- FindMarkers(mammal.combined, ident.1 = c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc"), ident.2 = c("CTB_d14","STB_d14","EVT_d14","putSTB"), test.use = "MAST" )

DE1 <- rownames(DELIST1)[which( DELIST1$p_val_adj<0.05 ) ]

L1 <- c(
  intersect(intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)]), DE1),
  intersect(intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Epiblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)]), DE1),
  intersect(intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Yolk-Sac / Hypoblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)]), DE1),
  intersect(intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Posterior Epiblast" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)]), DE1),
  intersect(intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Epiblast committed" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)]), DE1),
  intersect(intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Secondary Yolk-Sac" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)]), DE1),
  intersect(intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="Amnion" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)]),DE1) )

#epi
"JARID2"
"CASC15"
"LINC01194"
"MT1G"
"MT1H"
"NLGN4X"
"RIMS2"
"CD9"
"DPPA4"
"CLU"
"NAV3"

"CDH1"
"CLDN4"
"PRLR"
"DMD"
"CLDN6"
"EZR"

"LGR5", "LHX1",
"OPHN1"
"SLC3A"
"SLC2A3"
"TANC2"
"TCIM"

"COL5A2"
"COL6A3"
"TNC"
"SPARC"
"COL1A2"
"COL1A1"
"COL3A1"
"BST2"
"DCN"
"LUM"
"COL6A1"
"COL6A2"
"VIM"
"DPP4"
"CST3"
"DNAJC15"

subs2 <- which(mammal.combined$Cells%in%c("putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
L1 <- intersect(rownames(Ds1),SCMmarker$gene[which(SCMmarker$`cluster annotation`=="ExEM" & SCMmarker$avg_log2FC>log2(2) & SCMmarker$p_val_adj<0.05)])
X <- (GetAssayData(mammal.combined,assay = "RNA")) 
Xh <- t(scale(t(X)))
Xp <- Xh[L1,]
Xp <- na.omit(Xp)[,subs2]
annotation_col = data.frame(Stage = factor(droplevels(mammal.combined$Cells[subs2])), Batch =  factor((mammal.combined$Genotype2[subs2])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 50, filename = paste(saveext,"/DimRed/BatchAll_AllMarkers",".pdf",sep="") ,width=45,height=20)







#Redo heatmaps based on clusters

#readRDS(hc1,file=paste(saveext,"/Batch1_pheatmapclusterTb.rds",sep=""))
#readRDS(hc2a,file=paste(saveext,"/Batch2A_pheatmapclusterTb.rds",sep=""))
#readRDS(hc2b,file=paste(saveext,"/Batch2B_pheatmapclusterTb.rds",sep=""))
#readRDS(hc3,file=paste(saveext,"/Batch3_pheatmapclusterTb.rds",sep=""))
#readRDS(hc6,file=paste(saveext,"/Batch6_pheatmapclusterTb.rds",sep=""))
#saveRDS(hc7,file=paste(saveext,"/Batch7_pheatmapclusterTb.rds",sep=""))
#saveRDS(hc8,file=paste(saveext,"/Batch8_pheatmapclusterTb.rds",sep=""))

#saveRDS(hc_1,file=paste(saveext,"/Batch1_pheatmapclusterembryonic.rds",sep=""))
#saveRDS(hc_2a,file=paste(saveext,"/Batch2A_pheatmapclusterembryonic.rds",sep=""))
#saveRDS(hc_2b,file=paste(saveext,"/Batch2B_pheatmapclusterembryonic.rds",sep=""))
#saveRDS(hc_3,file=paste(saveext,"/Batch3_pheatmapclusterembryonic.rds",sep=""))
#saveRDS(hc_6,file=paste(saveext,"/Batch6_pheatmapclusterembryonic.rds",sep=""))
#saveRDS(hc_7,file=paste(saveext,"/Batch7_pheatmapclusterembryonic.rds",sep=""))
#saveRDS(hc_8,file=paste(saveext,"/Batch8_pheatmapclusterembryonic.rds",sep=""))




putPGC <- c("CAAGAGGTCCACCCTA-1_2","TTATTGCAGCACAAAT-1_2","TCGCACTAGCCAAGGT-1_2","GTAGATCAGCTGGCCT-1_2","ATCGATGGTAGCGTAG-1_2","GTGTTCCAGGGCAGGA-1_2","CAAGACTGTATTGCCA-1_2","TTCTTCCAGTCATCGT-1_2","CCTGCATCAAGTATCC-1_3","ACACCAACAACCGCTG-1_3","CCTCAACTCTTGCAAG-1_3","TAACTTCGTACAAAGT-1_3","AGCGCCACACGGCACT-1_3","TCATGTTGTCGACGCT-1_3","CAGTGCGAGTAGCCAG-1_3","ATGAGTCAGTCAATCC-1_3","TTACAGGGTCTTCTAT-1_3","CATCGTCTCGTTGTAG-1_6","AAAGAACAGTAAACGT-1_6","ACACGCGAGGAGGGTG-1_6","GCATCTCTCGGAATGG-1_6","GGAAGTGAGGCCTTGC-1_6")
Idents(mammal.combined) <- mammal.combined$Cells
Idents(mammal.combined,cells=putPGC) <- "putPGC"
mammal.combined$Cells <- Idents(mammal.combined)


mammal.combined <- subset(mammal.combined,idents=c("CTB_d14","STB_d14","EVT_d14","putSTB","putPGC","Hyp_d14","putExMes","ExMes_d14","Am/EmDisc_d14","EmDisc_d14","Am_d14","Hyp/Am","Am/EmDisc") )
Idents(mammal.combined) <- mammal.combined$ID3
D1 <- subset(mammal.combined,idents="C1")
D2 <- subset(mammal.combined,idents="C2")
D3 <- subset(mammal.combined,idents=c("C3","C4") )
D6 <- subset(mammal.combined,idents="C6")
D7 <- subset(mammal.combined,idents="C7")
D8 <- subset(mammal.combined,idents="C8")

saveRDS(D1,file=paste(saveext,"MatteoD1.rds",sep=""))
saveRDS(D2,file=paste(saveext,"MatteoD2.rds",sep=""))
saveRDS(D3,file=paste(saveext,"MatteoD3.rds",sep=""))
saveRDS(D6,file=paste(saveext,"MatteoD6.rds",sep=""))
saveRDS(D7,file=paste(saveext,"MatteoD7.rds",sep=""))
saveRDS(D8,file=paste(saveext,"MatteoD8.rds",sep=""))




