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
D$Cells <- D$CorrectLabel

Idents(D,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="STB")]) <- "STB"
Idents(D,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="putSTB")]) <- "STB"
Idents(D,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="Am")]) <- "Am"
Idents(D,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="Am_2")]) <- "Am"
Idents(D,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="CTB")]) <- "CTB"
Idents(D,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="ExMes")]) <- "ExMes"
Idents(D,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="Hypoblast")]) <- "Hypoblast"
Idents(D,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="STB1")]) <- "STB"
Idents(D,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="SYS")]) <- "SYS"
Idents(D,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="EmDisc")]) <- "EmDisc"
Idents(D,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="EmDisc_2")]) <- "EmDisc"
Idents(D,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="Am/EmDisc")]) <- "Am/EmDisc"
Idents(D,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="EVT")]) <- "EVT"
Idents(D,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="ExMes/SYS")]) <- "ExMes/SYS"
Idents(D,cells=Ano$ID[which(Ano$`Manual neighbour ano`=="pPGC")]) <- "PGC"


Idents(D,WhichCells(D,idents=c("putSTB"))) <- "STB"

Idents(D,WhichCells(D,idents=c("Lumenal","SOX9P","Prolif","SOX9LRG5","Glandular"))) <- "Epithelial"

Idents(D,WhichCells(D,idents=c("SYS"))) <- "YS"
Idents(D,WhichCells(D,idents=c("ExMes/SYS"))) <- "ExMes/YS"
Idents(D,WhichCells(D,idents=c("Stromal fibroblasts"))) <- "Stromal"


pPGC <- c("CGACAGCCATCTTCGC-1_2",
          "CCCAACTTCTACTTCA-1_2",
          "ATAGACCTCTGCTTAT-1_2",
          "AACCTGATCCATCCGT-1_7",
          "TTACAGGGTCTTCTAT-1_3",
          "TAACTTCGTACAAAGT-1_3",
          "CCTCAACTCTTGCAAG-1_3",
          "TCATGTTGTCGACGCT-1_3")

#CCTCAACTCTTGCAAG-1_3
#TAACTTCGTACAAAGT-1_3
#TCATGTTGTCGACGCT-1_3
#TTACAGGGTCTTCTAT-1_3
#WhichCells(D,idents=c("EmDisc","PGC","Am/EmDisc","Am"))

pgcsubs1 <- WhichCells(D,expression = SOX17>0.1)
pgcsubs2 <- WhichCells(D,expression = TFAP2C>0.1)
list1 <- intersect(intersect(pgcsubs1,pgcsubs2),WhichCells(D,idents=c("EmDisc","PGC","Am/EmDisc","Am")))

Idents(D,cells=list1) <- "PGC"


#p<-FeatureScatter(D,feature1 = "SOX17",feature2 = "TFAP2A", cells=WhichCells(D,idents=c("EmDisc","PGC","Am/EmDisc","Am")))
#ggsave(filename=paste(saveext,"/DimRed/tetspdf.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


#pgcsubs3 <- WhichCells(D,expression = TFAP2A>0.5)
#pgcsubs4 <- WhichCells(D,expression = PRDM1>0.5)

#pgcsubs5 <- intersect(pgcsubs1,pgcsubs2)
#pgcsubs6 <- intersect(pgcsubs1,pgcsubs3)
#pgcsubs7 <- intersect(pgcsubs1,pgcsubs4)

#cells2 <- setdiff( colnames(mammal.combinedA) , unique(c(pgcsubs5,pgcsubs6, pgcsubs7)) )
#cells2 <- colnames(mammal.combinedA) 
#cells = unique(c(pgcsubs5,pgcsubs6, pgcsubs7))

#mammal.combinedC <- subset(mammal.combinedA, cells = cells2 ) 



list1 <- WhichCells(D,idents=c("CTB","STB","EVT"))

TbAno <- readRDS("/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/TB_anotatedbyCyno.rds")

list2 <- WhichCells(TbAno,idents=c("CTB","STB","EVT"))

Idents(D,cells=list1) <- "Old"

Idents(D,cells=WhichCells(TbAno,idents=c("CTB"))) <- "CTB"
Idents(D,cells=WhichCells(TbAno,idents=c("STB"))) <- "STB"
Idents(D,cells=WhichCells(TbAno,idents=c("EVT"))) <- "EVT"

Idents(D,cells=WhichCells(D,idents=c("Old"))) <- "STB"


Cols =  c("EmDisc"="#2D7EE0",
                 "Am/EmDisc"="#00B9E3",
                 "Am"="#3CDDDD",
                 "PGC"="#E154E8",
                 "Hypoblast"="#FF8C64",
                 "YS"="#2E892E",
                 "ExMes/YS"="#22BF35",
                 "ExMes"="#9FD13F",
                 "CTB"="#E2C120",
                 "STB"="#D68C34",
                 "EVT"="#BF3737",
                 "Stromal"="#AD71BC",
                 "Epithelial"="#C4807F",
                 "Ciliated"="#7F3C3C")



p<-DimPlot(D, cols = Cols, pt.size = 2, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Embryo_nosplit.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p<-DimPlot(D, cols = Cols, pt.size = 2, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Embryo_nosplit.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


p<-DimPlot(Matteo,  pt.size = 2, reduction = "oldUMAP",  label = TRUE, repel = TRUE) 
ggsave(filename=paste("~/Desktop/UMAP.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)



Cont1 <- colnames(D)[which(D$Cl05%in%c(9,15,11,20,14,21,23) & Idents(D)%in%c("Stromal"))]
Cont2 <- colnames(D)[which(D$Cl05%in%c(9,15,11,20,14,21,23) & Idents(D)%in%c("Epithelial"))]

Idents(D,cells=Cont1) <- "Cont"
Idents(D,cells=Cont2) <- "Cont"

D2 <- subset(D,idents = "Cont", invert=TRUE)

DefaultAssay(D2) <- "RNA"
p<-VlnPlot(D2,features=c("IGFBP3","NPNT"))
ggsave(filename=paste(saveext,"/DimRed/MoreViolin_IGFBP2_NPNT.pdf",sep=""),width = 10, height = 10,p,limitsize = FALSE)


p<-DimPlot(D2, cols = Cols, pt.size = 2, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Embryo_nosplit.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


p<-DimPlot(D2, cols = Cols, pt.size = 2, reduction = "pca",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Embryo_nosplit.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)





p<-DimPlot(D, cols = Cols, pt.size = 2, reduction = "umap",  label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Embryo_nosplit_test.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)



MiniList1 <- c("JARID2","CASC15","LINC01194","MT1G","MT1H","NLGN4X","RIMS2","CD9","DPPA4","CLU","NAV3","CDH1","CLDN4","PRLR","DMD","CLDN6","EZR","LGR5", "LHX1","OPHN1","SLC3A2","SLC2A3","TANC2","TCIM","COL5A2","COL6A3","TNC","SPARC","COL1A2","COL1A1","COL3A1","DCN","LUM","COL6A1","COL6A2","DPP4","CST3","DNAJC15","DNMT3A","DPPA3","OTX2","POU5F1","PDGFA","PODXL","SFRP1","SOX15","GABRP","VTCN1","WNT6","BAMBI","TFAP2C","TFAP2A","GATA3","CGA","SOX17","GATA4","CER1","LEFTY1","LEFTY2","APOA1","APOB","LINC00261","HHEX","FOXF1","PDGFRA","HAND2","BST2","GATA6","HGF","PRDM1","SNAI1","VIM","FN1","SOX17","TFAP2A","TFAP2C","PRDM1","POU5F1","KLF4")
MiniList2 <- c("DNMT3A","DPPA3","OTX2","NOTUM","JAM3","TACSTD2","ENPEP","GATA2","GATA3","TFAP2C","TFAP2A","TEAD3","TP63","OVOL1","CCKBR","CGA","CGB5","CGB8","SDC1","HLA-C","HLA-G","ASCL2","SNAI1","VIM","FN1","ELF5","MKI67","CDH1","FSTL1","FSTL3","KRT19","KRT7","SPARC","LAMC1","CTNNAL1")



AvE <- AverageExpression(D)
Expression1 <- AvE$integrated[,c("CTB","STB","EVT","ExMes","Hypoblast","YS","Am","Am/EmDisc","EmDisc","PGC")]
Expression2 <- AvE$RNA[,c("CTB","STB","EVT","ExMes","Hypoblast","YS","Am","Am/EmDisc","EmDisc","PGC")]

MiniList2 <- c("DNMT3A","DPPA3","OTX2","NOTUM","JAM3","TACSTD2","ENPEP","GATA2","GATA3","TFAP2C","TFAP2A","TEAD3","TP63","OVOL1","CCKBR","CGA","CGB5","CGB8","SDC1","HLA-C","HLA-G","ASCL2","SNAI1","VIM","FN1","ELF5","MKI67","CDH1","FSTL1","FSTL3","KRT19","KRT7","SPARC","LAMC1","CTNNAL1")
MiniList1 <- c("JARID2","CASC15","LINC01194","MT1G","MT1H","NLGN4X","RIMS2","CD9","DPPA4","CLU","NAV3","CDH1","CLDN4","PRLR","DMD","CLDN6","EZR","LGR5", "LHX1","OPHN1","SLC3A2","SLC2A3","TANC2","TCIM","COL5A2","COL6A3","TNC","SPARC","COL1A2","COL1A1","COL3A1","DCN","LUM","COL6A1","COL6A2","DPP4","CST3","DNAJC15","DNMT3A","DPPA3","OTX2","POU5F1","PDGFA","PODXL","SFRP1","SOX15","GABRP","VTCN1","WNT6","BAMBI","TFAP2C","TFAP2A","GATA3","CGA","SOX17","GATA4","CER1","LEFTY1","LEFTY2","APOA1","APOB","LINC00261","HHEX","FOXF1","PDGFRA","HAND2","BST2","GATA6","HGF","PRDM1","SNAI1","VIM","FN1")


Markers1 <- c("OPHN1",
              "TEAD3",
              "GATA2",
              "GATA3",
              "CD34",
              "EZR",
              "JAM3",
              "CGB8",
              "CGB5",
              "CGA",
              "SDC1",
              "ERVW-1",
              "PRDM6",
              "TBX3",
              "DIO2",
              "NOTUM",
              "ASCL2",
              "HLA-G",
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
pheatmap(Expression2[Markers1,],color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_selectmarkers_alttbAno2.pdf",sep=""),scale="row",width=10,height=32)





mat_breaks <- seq(-3, 3, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap(Expression2[MiniList1,],color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_selectmarkers_alttbAno2_MinList1.pdf",sep=""),scale="row",width=10,height=32)



mat_breaks <- seq(-3, 3, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFD185","#FF0B07"))
pheatmap(Expression2[MiniList2,],color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/Markers/LineageHeatmap_selectmarkers_alttbAno2_MinList2.pdf",sep=""),scale="row",width=10,height=32)

