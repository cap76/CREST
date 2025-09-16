library(Seurat) 
library(ggplot2) 
library(cowplot) 
library(Matrix) 
library(dplyr) 
library("Matrix")
#library(harmony)

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

#Load in the previous run to get informatio about clusters and anotations
mammal.combined <- readRDS(file = paste(saveext,"Seurat_combined_ano1.rds",sep=""))

mammal.combined$ano <- Idents(mammal.combined)

#Load in the data
D1 <- readRDS("../Mole_experimentRv5_1.rds")
D2 <- readRDS("../Mole_experimentRv5_2.rds")
D3 <- readRDS("../Mole_experimentRv5_3.rds")
D3B <- readRDS("../Mole_experimentRv5_3B.rds")
D6 <- readRDS("../Mole_experimentRv5_6.rds")
D7 <- readRDS("../Mole_experimentRv5_7.rds")
D8 <- readRDS("../Mole_experimentRv5_8.rds")

#Load in the genotypes
BC1 <- read.table("Genotype/D1_4Cl.tsv",sep="\t",header=TRUE)
BC2 <- read.table("Genotype/D2_4Cl.tsv",sep="\t",header=TRUE)
BC3 <- read.table("Genotype/D3_6Cl.tsv",sep="\t",header=TRUE)
BC3B <- read.table("Genotype/D3B_6Cl.tsv",sep="\t",header=TRUE)
BC6 <- read.table("Genotype/D6_4Cl.tsv",sep="\t",header=TRUE)
BC7 <- read.table("Genotype/D7_4Cl.tsv",sep="\t",header=TRUE)
BC8 <- read.table("Genotype/D8_4Cl.tsv",sep="\t",header=TRUE)
d1 <- as.data.frame(colnames(D1))
d2 <- as.data.frame(colnames(D2))
d3 <- as.data.frame(colnames(D3))
d3B <- as.data.frame(colnames(D3B))
d6 <- as.data.frame(colnames(D6))
d7 <- as.data.frame(colnames(D7))
d8 <- as.data.frame(colnames(D8))
colnames(d1) <- "barcode"
colnames(d2) <- "barcode"
colnames(d3) <- "barcode"
colnames(d3B) <- "barcode"
colnames(d6) <- "barcode"
colnames(d7) <- "barcode"
colnames(d8) <- "barcode"
m1 <- merge(x = d1, y = BC1, by = "barcode", all = TRUE)
m2 <- merge(x = d2, y = BC2, by = "barcode", all = TRUE)
m3 <- merge(x = d3, y = BC3, by = "barcode", all = TRUE)
m3B <- merge(x = d3B, y = BC3B, by = "barcode", all = TRUE)
m6 <- merge(x = d6, y = BC6, by = "barcode", all = TRUE)
m7 <- merge(x = d7, y = BC7, by = "barcode", all = TRUE)
m8 <- merge(x = d8, y = BC8, by = "barcode", all = TRUE)


#Now update labels in Seurat objects
#D1 <- readRDS("../Mole_experimentRv5_1.rds")
Idents(D1) <- "NAssigned"
Idents(D1,cells=BC1$barcode[which(BC1$assignment==0)]) <- "G1"
Idents(D1,cells=BC1$barcode[which(BC1$assignment==1)]) <- "G2"
Idents(D1,cells=BC1$barcode[which(BC1$assignment==2)]) <- "G3"
Idents(D1,cells=BC1$barcode[which(BC1$assignment==3)]) <- "G4"
D1$PCl1 <- m1[,"cluster0"]
D1$PCl2 <- m1[,"cluster1"]
D1$PCl3 <- m1[,"cluster2"]
D1$PCl4 <- m1[,"cluster3"]

#D2 <- readRDS("../Mole_experimentRv5_2.rds")
Idents(D2) <- "NAssigned"
Idents(D2,cells=BC2$barcode[which(BC2$assignment==0)]) <- "G1"
Idents(D2,cells=BC2$barcode[which(BC2$assignment==1)]) <- "G2"
Idents(D2,cells=BC2$barcode[which(BC2$assignment==2)]) <- "G3"
Idents(D2,cells=BC2$barcode[which(BC2$assignment==3)]) <- "G4"
D2$PCl1 <- m2[,"cluster0"]
D2$PCl2 <- m2[,"cluster1"]
D2$PCl3 <- m2[,"cluster2"]
D2$PCl4 <- m2[,"cluster3"]

#D3 <- readRDS("../Mole_experimentRv5_3.rds")
Idents(D3) <- "NAssigned"
Idents(D3,cells=BC3$barcode[which(BC3$assignment==0)]) <- "G1"
Idents(D3,cells=BC3$barcode[which(BC3$assignment==1)]) <- "G2"
#Idents(D3,cells=BC3$barcode[which(BC3$assignment==2)]) <- "G3"
Idents(D3,cells=BC3$barcode[which(BC3$assignment==3)]) <- "G4"
Idents(D3,cells=BC3$barcode[which(BC3$assignment==4)]) <- "G5"
Idents(D3,cells=BC3$barcode[which(BC3$assignment==5)]) <- "G6"
D3$PCl1 <- m3[,"cluster0"]
D3$PCl2 <- m3[,"cluster1"]
D3$PCl3 <- m3[,"cluster2"]
D3$PCl4 <- m3[,"cluster3"]
D3$PCl5 <- m3[,"cluster4"]
D3$PCl6 <- m3[,"cluster5"]

#D3B <- readRDS("../Mole_experimentRv5_3B.rds")
Idents(D3B) <- "NAssigned"
Idents(D3B,cells=BC3B$barcode[which(BC3B$assignment==0)]) <- "G1"
Idents(D3B,cells=BC3B$barcode[which(BC3B$assignment==1)]) <- "G2"
Idents(D3B,cells=BC3B$barcode[which(BC3B$assignment==2)]) <- "G3"
Idents(D3B,cells=BC3B$barcode[which(BC3B$assignment==3)]) <- "G4"
Idents(D3B,cells=BC3B$barcode[which(BC3B$assignment==4)]) <- "G5"
Idents(D3B,cells=BC3B$barcode[which(BC3B$assignment==5)]) <- "G6"

D3B$PCl1 <- m3B[,"cluster0"]
D3B$PCl2 <- m3B[,"cluster1"]
D3B$PCl3 <- m3B[,"cluster2"]
D3B$PCl4 <- m3B[,"cluster3"]
D3B$PCl5 <- m3B[,"cluster4"]
D3B$PCl6 <- m3B[,"cluster5"]

#D6 <- readRDS("../Mole_experimentRv5_6.rds")
Idents(D6) <- "NAssigned"
Idents(D6,cells=BC6$barcode[which(BC6$assignment==0)]) <- "G1"
Idents(D6,cells=BC6$barcode[which(BC6$assignment==1)]) <- "G2"
Idents(D6,cells=BC6$barcode[which(BC6$assignment==2)]) <- "G3"
Idents(D6,cells=BC6$barcode[which(BC6$assignment==3)]) <- "G4"

D6$PCl1 <- m6[,"cluster0"]
D6$PCl2 <- m6[,"cluster1"]
D6$PCl3 <- m6[,"cluster2"]
D6$PCl4 <- m6[,"cluster3"]

#D7 <- readRDS("../Mole_experimentRv5_7.rds")
Idents(D7) <- "NAssigned"
Idents(D7,cells=BC7$barcode[which(BC7$assignment==0)]) <- "G1"
Idents(D7,cells=BC7$barcode[which(BC7$assignment==1)]) <- "G2"
Idents(D7,cells=BC7$barcode[which(BC7$assignment==2)]) <- "G3"
Idents(D7,cells=BC7$barcode[which(BC7$assignment==3)]) <- "G4"

D7$PCl1 <- m7[,"cluster0"]
D7$PCl2 <- m7[,"cluster1"]
D7$PCl3 <- m7[,"cluster2"]
D7$PCl4 <- m7[,"cluster3"]

#D8 <- readRDS("../Mole_experimentRv5_8.rds")
Idents(D8) <- "NAssigned"
Idents(D8,cells=BC8$barcode[which(BC8$assignment==0)]) <- "G1"
Idents(D8,cells=BC8$barcode[which(BC8$assignment==1)]) <- "G2"
Idents(D8,cells=BC8$barcode[which(BC8$assignment==2)]) <- "G3"
Idents(D8,cells=BC8$barcode[which(BC8$assignment==3)]) <- "G4"

D8$PCl1 <- m8[,"cluster0"]
D8$PCl2 <- m8[,"cluster1"]
D8$PCl3 <- m8[,"cluster2"]
D8$PCl4 <- m8[,"cluster3"]

D1$Genotype <- Idents(D1)
D2$Genotype <- Idents(D2)
D3$Genotype <- Idents(D3)
D3B$Genotype <- Idents(D3B)
D6$Genotype <- Idents(D6)
D7$Genotype <- Idents(D7)
D8$Genotype <- Idents(D8)

D1$ID3 <- "C1"
D2$ID3 <- "C2"
D3$ID3 <- "C3"
D3B$ID3 <- "C3B"
D6$ID3 <- "C6"
D7$ID3 <- "C7"
D8$ID3 <- "C8"

#Get cell labels for specific subsets of the data: EmD, Hyp, and Epith
List1 <-  which(mammal.combined$Cl15%in%c(31,32,37,24,5,8,0,19) & mammal.combined$ID3=="C1")
List2 <-  which(mammal.combined$Cl15%in%c(31,32,37,24,5,8,0,19) & mammal.combined$ID3=="C2")
List3 <-  which(mammal.combined$Cl15%in%c(31,32,37,24,5,8,0,19) & mammal.combined$ID3=="C3")
List4 <-  which(mammal.combined$Cl15%in%c(31,32,37,24,5,8,0,19) & mammal.combined$ID3=="C4")
List6 <-  which(mammal.combined$Cl15%in%c(31,32,37,24,5,8,0,19) & mammal.combined$ID3=="C6")
List7 <-  which(mammal.combined$Cl15%in%c(31,32,37,24,5,8,0,19) & mammal.combined$ID3=="C7")
List8 <-  which(mammal.combined$Cl15%in%c(31,32,37,24,5,8,0,19) & mammal.combined$ID3=="C8")

Idents(D1,cells=gsub("-1_7","-1",names(which(mammal.combined$ano%in%c("EmDisc_d14") & mammal.combined$ID3=="C1"))) ) <- "EmDisc_d14"
Idents(D2,cells=gsub("-1_8","-1",names(which(mammal.combined$ano%in%c("EmDisc_d14") & mammal.combined$ID3=="C2"))) ) <- "EmDisc_d14"
Idents(D3,cells=gsub("-1_9","-1",names(which(mammal.combined$ano%in%c("EmDisc_d14") & mammal.combined$ID3=="C3"))) ) <- "EmDisc_d14"
Idents(D3B,cells=gsub("-1_10","-1",names(which(mammal.combined$ano%in%c("EmDisc_d14") & mammal.combined$ID3=="C4"))) ) <- "EmDisc_d14"
Idents(D6,cells=gsub("-1_11","-1",names(which(mammal.combined$ano%in%c("EmDisc_d14") & mammal.combined$ID3=="C6"))) ) <- "EmDisc_d14"
Idents(D7,cells=gsub("-1_12","-1",names(which(mammal.combined$ano%in%c("EmDisc_d14") & mammal.combined$ID3=="C7"))) ) <- "EmDisc_d14"
Idents(D8,cells=gsub("-1_13","-1",names(which(mammal.combined$ano%in%c("EmDisc_d14") & mammal.combined$ID3=="C8"))) ) <- "EmDisc_d14"

Idents(D1,cells=gsub("-1_7","-1",names(which(mammal.combined$ano%in%c("Hyp_d14") & mammal.combined$ID3=="C1"))) ) <- "Hyp_d14"
Idents(D2,cells=gsub("-1_8","-1",names(which(mammal.combined$ano%in%c("Hyp_d14") & mammal.combined$ID3=="C2"))) ) <- "Hyp_d14"
Idents(D3,cells=gsub("-1_9","-1",names(which(mammal.combined$ano%in%c("Hyp_d14") & mammal.combined$ID3=="C3"))) ) <- "Hyp_d14"
Idents(D3B,cells=gsub("-1_10","-1",names(which(mammal.combined$ano%in%c("Hyp_d14") & mammal.combined$ID3=="C4"))) ) <- "Hyp_d14"
Idents(D6,cells=gsub("-1_11","-1",names(which(mammal.combined$ano%in%c("Hyp_d14") & mammal.combined$ID3=="C6"))) ) <- "Hyp_d14"
Idents(D7,cells=gsub("-1_12","-1",names(which(mammal.combined$ano%in%c("Hyp_d14") & mammal.combined$ID3=="C7"))) ) <- "Hyp_d14"
Idents(D8,cells=gsub("-1_13","-1",names(which(mammal.combined$ano%in%c("Hyp_d14") & mammal.combined$ID3=="C8"))) ) <- "Hyp_d14"

Idents(D1,cells=gsub("-1_7","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 1") & mammal.combined$ID3=="C1"))) ) <- "Unciliated epithelia 1"
Idents(D2,cells=gsub("-1_8","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 1") & mammal.combined$ID3=="C2"))) ) <- "Unciliated epithelia 1"
Idents(D3,cells=gsub("-1_9","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 1") & mammal.combined$ID3=="C3"))) ) <- "Unciliated epithelia 1"
Idents(D3B,cells=gsub("-1_10","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 1") & mammal.combined$ID3=="C4"))) ) <- "Unciliated epithelia 1"
Idents(D6,cells=gsub("-1_11","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 1") & mammal.combined$ID3=="C6"))) ) <- "Unciliated epithelia 1"
Idents(D7,cells=gsub("-1_12","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 1") & mammal.combined$ID3=="C7"))) ) <- "Unciliated epithelia 1"
Idents(D8,cells=gsub("-1_13","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 1") & mammal.combined$ID3=="C8"))) ) <- "Unciliated epithelia 1"

Idents(D1,cells=gsub("-1_7","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 2") & mammal.combined$ID3=="C1"))) ) <- "Unciliated epithelia 2"
Idents(D2,cells=gsub("-1_8","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 2") & mammal.combined$ID3=="C2"))) ) <- "Unciliated epithelia 2"
Idents(D3,cells=gsub("-1_9","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 2") & mammal.combined$ID3=="C3"))) ) <- "Unciliated epithelia 2"
Idents(D3B,cells=gsub("-1_10","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 2") & mammal.combined$ID3=="C4"))) ) <- "Unciliated epithelia 2"
Idents(D6,cells=gsub("-1_11","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 2") & mammal.combined$ID3=="C6"))) ) <- "Unciliated epithelia 2"
Idents(D7,cells=gsub("-1_12","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 2") & mammal.combined$ID3=="C7"))) ) <- "Unciliated epithelia 2"
Idents(D8,cells=gsub("-1_13","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 2") & mammal.combined$ID3=="C8"))) ) <- "Unciliated epithelia 2"

Idents(D1,cells=gsub("-1_7","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 3") & mammal.combined$ID3=="C1"))) ) <- "Unciliated epithelia 3"
Idents(D2,cells=gsub("-1_8","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 3") & mammal.combined$ID3=="C2"))) ) <- "Unciliated epithelia 3"
Idents(D3,cells=gsub("-1_9","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 3") & mammal.combined$ID3=="C3"))) ) <- "Unciliated epithelia 3"
Idents(D3B,cells=gsub("-1_10","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 3") & mammal.combined$ID3=="C4"))) ) <- "Unciliated epithelia 3"
Idents(D6,cells=gsub("-1_11","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 3") & mammal.combined$ID3=="C6"))) ) <- "Unciliated epithelia 3"
Idents(D7,cells=gsub("-1_12","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 3") & mammal.combined$ID3=="C7"))) ) <- "Unciliated epithelia 3"
Idents(D8,cells=gsub("-1_13","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 3") & mammal.combined$ID3=="C8"))) ) <- "Unciliated epithelia 3"

Idents(D1,cells=gsub("-1_7","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 4") & mammal.combined$ID3=="C1"))) ) <- "Unciliated epithelia 4"
Idents(D2,cells=gsub("-1_8","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 4") & mammal.combined$ID3=="C2"))) ) <- "Unciliated epithelia 4"
Idents(D3,cells=gsub("-1_9","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 4") & mammal.combined$ID3=="C3"))) ) <- "Unciliated epithelia 4"
Idents(D3B,cells=gsub("-1_10","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 4") & mammal.combined$ID3=="C4"))) ) <- "Unciliated epithelia 4"
Idents(D6,cells=gsub("-1_11","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 4") & mammal.combined$ID3=="C6"))) ) <- "Unciliated epithelia 4"
Idents(D7,cells=gsub("-1_12","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 4") & mammal.combined$ID3=="C7"))) ) <- "Unciliated epithelia 4"
Idents(D8,cells=gsub("-1_13","-1",names(which(mammal.combined$ano%in%c("Unciliated epithelia 4") & mammal.combined$ID3=="C8"))) ) <- "Unciliated epithelia 4"

Idents(D1,cells=gsub("-1_7","-1",names(which(mammal.combined$ano%in%c("Stromal fibroblasts") & mammal.combined$ID3=="C1"))) ) <- "Stromal fibroblasts"
Idents(D2,cells=gsub("-1_8","-1",names(which(mammal.combined$ano%in%c("Stromal fibroblasts") & mammal.combined$ID3=="C2"))) ) <- "Stromal fibroblasts"
Idents(D3,cells=gsub("-1_9","-1",names(which(mammal.combined$ano%in%c("Stromal fibroblasts") & mammal.combined$ID3=="C3"))) ) <- "Stromal fibroblasts"
Idents(D3B,cells=gsub("-1_10","-1",names(which(mammal.combined$ano%in%c("Stromal fibroblasts") & mammal.combined$ID3=="C4"))) ) <- "Stromal fibroblasts"
Idents(D6,cells=gsub("-1_11","-1",names(which(mammal.combined$ano%in%c("Stromal fibroblasts") & mammal.combined$ID3=="C6"))) ) <- "Stromal fibroblasts"
Idents(D7,cells=gsub("-1_12","-1",names(which(mammal.combined$ano%in%c("Stromal fibroblasts") & mammal.combined$ID3=="C7"))) ) <- "Stromal fibroblasts"
Idents(D8,cells=gsub("-1_13","-1",names(which(mammal.combined$ano%in%c("Stromal fibroblasts") & mammal.combined$ID3=="C8"))) ) <- "Stromal fibroblasts"

Idents(D1,cells=gsub("-1_7","-1",names(which(mammal.combined$ano%in%c("CTB_d14") & mammal.combined$ID3=="C1"))) ) <- "CTB_d14"
Idents(D2,cells=gsub("-1_8","-1",names(which(mammal.combined$ano%in%c("CTB_d14") & mammal.combined$ID3=="C2"))) ) <- "CTB_d14"
Idents(D3,cells=gsub("-1_9","-1",names(which(mammal.combined$ano%in%c("CTB_d14") & mammal.combined$ID3=="C3"))) ) <- "CTB_d14"
Idents(D3B,cells=gsub("-1_10","-1",names(which(mammal.combined$ano%in%c("CTB_d14") & mammal.combined$ID3=="C4"))) ) <- "CTB_d14"
Idents(D6,cells=gsub("-1_11","-1",names(which(mammal.combined$ano%in%c("CTB_d14") & mammal.combined$ID3=="C6"))) ) <- "CTB_d14"
Idents(D7,cells=gsub("-1_12","-1",names(which(mammal.combined$ano%in%c("CTB_d14") & mammal.combined$ID3=="C7"))) ) <- "CTB_d14"
Idents(D8,cells=gsub("-1_13","-1",names(which(mammal.combined$ano%in%c("CTB_d14") & mammal.combined$ID3=="C8"))) ) <- "CTB_d14"

Idents(D1,cells=gsub("-1_7","-1",names(which(mammal.combined$ano%in%c("STB_d14") & mammal.combined$ID3=="C1"))) ) <- "STB_d14"
Idents(D2,cells=gsub("-1_8","-1",names(which(mammal.combined$ano%in%c("STB_d14") & mammal.combined$ID3=="C2"))) ) <- "STB_d14"
Idents(D3,cells=gsub("-1_9","-1",names(which(mammal.combined$ano%in%c("STB_d14") & mammal.combined$ID3=="C3"))) ) <- "STB_d14"
Idents(D3B,cells=gsub("-1_10","-1",names(which(mammal.combined$ano%in%c("STB_d14") & mammal.combined$ID3=="C4"))) ) <- "STB_d14"
Idents(D6,cells=gsub("-1_11","-1",names(which(mammal.combined$ano%in%c("STB_d14") & mammal.combined$ID3=="C6"))) ) <- "STB_d14"
Idents(D7,cells=gsub("-1_12","-1",names(which(mammal.combined$ano%in%c("STB_d14") & mammal.combined$ID3=="C7"))) ) <- "STB_d14"
Idents(D8,cells=gsub("-1_13","-1",names(which(mammal.combined$ano%in%c("STB_d14") & mammal.combined$ID3=="C8"))) ) <- "STB_d14"

#Idents(D1,cells=gsub("-1_7","-1",names(which(mammal.combined$ano%in%c("EVT") & mammal.combined$ID3=="C1"))) ) <- "EVT"
#Idents(D2,cells=gsub("-1_8","-1",names(which(mammal.combined$ano%in%c("EVT") & mammal.combined$ID3=="C2"))) ) <- "EVT"
#Idents(D3,cells=gsub("-1_9","-1",names(which(mammal.combined$ano%in%c("EVT") & mammal.combined$ID3=="C3"))) ) <- "EVT"
#Idents(D3B,cells=gsub("-1_10","-1",names(which(mammal.combined$ano%in%c("EVT") & mammal.combined$ID3=="C4"))) ) <- "EVT"
#Idents(D6,cells=gsub("-1_11","-1",names(which(mammal.combined$ano%in%c("EVT") & mammal.combined$ID3=="C6"))) ) <- "EVT"
#Idents(D7,cells=gsub("-1_12","-1",names(which(mammal.combined$ano%in%c("EVT") & mammal.combined$ID3=="C7"))) ) <- "EVT"
#Idents(D8,cells=gsub("-1_13","-1",names(which(mammal.combined$ano%in%c("EVT") & mammal.combined$ID3=="C8"))) ) <- "EVT"

#Subset data base on these cells
D1 <- subset(D1,cells=gsub("-1_7","-1",names(List1)))
D2 <- subset(D2,cells=gsub("-1_8","-1",names(List2)))
D3 <- subset(D3,cells=gsub("-1_9","-1",names(List3)))
D3B <- subset(D3B,cells=gsub("-1_10","-1",names(List4)))
D6 <- subset(D6,cells=gsub("-1_11","-1",names(List6)))
D7 <- subset(D7,cells=gsub("-1_12","-1",names(List7)))
D8 <- subset(D8,cells=gsub("-1_13","-1",names(List8)))

#Load in human ref
MergedData1 <- readRDS(file="Merged1.rds")
MergedData2 <- readRDS(file="Merged2.rds")
MergedData3 <- readRDS(file="Merged3.rds")
MergedData4 <- readRDS(file="Merged4.rds")
MergedData5 <- readRDS(file="Merged5.rds")
MergedData6 <- readRDS(file="Merged6.rds")
#Give old dataset a unique ID
MergedData1$species1 <- "Old"
MergedData2$species1 <- "Old"
MergedData3$species1 <- "Old"
MergedData4$species1 <- "Old"
MergedData5$species1 <- "Old"
MergedData6$species1 <- "Old"
#Get subset of data
List1 <-  which(mammal.combined$Cl15%in%c(31,32,37,24,5,8,0,19) & mammal.combined$ID3=="B1")
List2 <-  which(mammal.combined$Cl15%in%c(31,32,37,24,5,8,0,19) & mammal.combined$ID3=="B2")
List3 <-  which(mammal.combined$Cl15%in%c(31,32,37,24,5,8,0,19) & mammal.combined$ID3=="B3")
List4 <-  which(mammal.combined$Cl15%in%c(31,32,37,24,5,8,0,19) & mammal.combined$ID3=="B4")
List5 <-  which(mammal.combined$Cl15%in%c(31,32,37,24,5,8,0,19) & mammal.combined$ID3=="B5")
List6 <-  which(mammal.combined$Cl15%in%c(31,32,37,24,5,8,0,19) & mammal.combined$ID3=="B6")
Dp1 <- subset(MergedData1,cells=gsub("_1","",names(List1)))
Dp2 <- subset(MergedData2,cells=gsub("_2","",names(List2)))
Dp3 <- subset(MergedData3,cells=gsub("_3","",names(List3)))
Dp4 <- subset(MergedData4,cells=gsub("_4","",names(List4)))
Dp5 <- subset(MergedData5,cells=gsub("_5","",names(List5)))
Dp6 <- subset(MergedData6,cells=gsub("_6","",names(List6)))


#Now do itegration
mammal.anchors <- FindIntegrationAnchors(object.list = list(D1,D2,D3B,D6,D8,Dp1,Dp2,Dp3,Dp5), dims = 1:20, anchor.features = 4000, k.filter=20)
mammal.combined1 <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined1) <- "integrated"
mammal.combined1 <- ScaleData(mammal.combined1)
mammal.combined1 <- RunPCA(mammal.combined1, npcs = 20, verbose = FALSE)
mammal.combined1 <- RunUMAP(mammal.combined1, reduction = "pca", dims = 1:20)
mammal.combined1 <- FindNeighbors(mammal.combined1, reduction = "pca", dims = 1:20)

#saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined_rerun.rds",sep=""))

#sdadasd

#mammal.combined$Cells <- Idents(mammal.combined)


DimPlot(mammal.combined1,   reduction = "umap", split.by = "ID3", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat_split_rerun",".pdf",sep=""),width = 70, height = 10, limitsize = FALSE, useDingbats=FALSE)
sadsadsa

DimPlot(mammal.combined,  reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat_rerun",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE)

#saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined_.rds",sep=""))


mammal.combined <- FindClusters(mammal.combined, resolution = 1.0)
DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-Cl10",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-Cl10",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
mammal.combined$Cl10 <- Idents(mammal.combined)

mammal.combined <- FindClusters(mammal.combined, resolution = 1.5)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-Cl15",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
mammal.combined$Cl15 <- Idents(mammal.combined)
DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-Cl15",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)


mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-Cl05",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-Cl05",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
mammal.combined$Cl05 <- Idents(mammal.combined)


#PutaEmb1 <- colnames(mammal.combined)[which(mammal.combined$Cl15=="7" & mammal.combined$species2=="D")]
#PutaEmb2 <- colnames(mammal.combined)[which(mammal.combined$Cl15%in%c(13,19,24) & mammal.combined$species2=="D")]
#saveRDS(PutaEmb1,file="PutativeEmb1.rds")
#saveRDS(PutaEmb2,file="PutativeEmb2.rds")



DefaultAssay(mammal.combined) <- "RNA"
#AvExp <- FindAllMarkers(mammal.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#write.csv(as.data.frame(AvExp$RNA),file="AverageExpComb.csv",quote=FALSE)

FeaturePlot(mammal.combined, features = "SOX17", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_SOX17.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "SOX2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_SOX2.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "POU5F1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_POU5F1.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NANOG", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_NANOG.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "DAPP1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_DAPP1.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "KLF17", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_KLF17.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "TBXT", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_T.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "MIXL1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_MIXL1.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "EOMES", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_EOMES.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "HNF1B", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_HF1B.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "GSG1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GSG1.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "GATM", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GATM.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "NOTO", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_NOTO.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "FXYD3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_FXYD3.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)


FeaturePlot(mammal.combined, features = "PADI1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_PADI1.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "SP100", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_SP100.pdf",sep=""),width = 550, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "MST1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_MST1.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)



FeaturePlot(mammal.combined, features = "RPS6KA5", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "GADD45G", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GADD45G.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "CGA", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_CGA.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)


FeaturePlot(mammal.combined, features = "GATA2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GATA2.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "GATA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GATA3.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "TFAP2C", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_TFAP2C.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "TFAP2A", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_TFAP2A.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)


FeaturePlot(mammal.combined, features = "GATA6", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GATA6.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "SOX17", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_SOX17.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "APOA1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_APOA1.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "PDGFRA", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_PDGFRA.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "OTX2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_OTX2.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "GSC", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GSC.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "LEFTY2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_LEFTY2.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "TTR", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_TTR.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "SNAI2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_SNAI2.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "HGF", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_HGF.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "HAND1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_HAND1.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "ISL1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_ISL1.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "VTCN1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_VTCN1.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "MESP1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_MESP1.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "NANOS3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_NANOS3.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "IGFBP1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_IGFBP1.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "IGFBP2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_IGFBP2.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "IGFBP6", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_IGFBP6.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "PRL", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_PRL.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

#FeaturePlot(mammal.combined, features = "CD73", split.by = "species1", cols =  c("lightgrey", "darkblue"))
#ggsave(filename=paste(saveext,"/Markers/Marker_CD73.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)


FeaturePlot(mammal.combined, features = "EPCAM", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_EPCAM.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "UCA1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_UCA1.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "FOXJ1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_FOXJ1.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)





mammal.anchors <- FindIntegrationAnchors(object.list = list(D1,D2,D3B,D6,D8,Dp1,Dp2,Dp3,Dp4,Dp5,Dp6), dims = 1:20, anchor.features = 4000, k.filter=20)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined_2.rds",sep=""))

#sdadasd

mammal.combined$Cells <- Idents(mammal.combined)

DimPlot(mammal.combined,   reduction = "umap", split.by = "species1", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat_splitII","2.pdf",sep=""),width = 110, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined,  reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_SeuratII","2.pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE)

#saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined_.rds",sep=""))
mammal.combined <- FindClusters(mammal.combined, resolution = 1.0)
DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-Cl10","2.pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-Cl10","2.pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
mammal.combined$Cl10 <- Idents(mammal.combined)
mammal.combined <- FindClusters(mammal.combined, resolution = 1.5)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-Cl15","2.pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
mammal.combined$Cl15 <- Idents(mammal.combined)
DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-Cl15","2.pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-Cl05","2.pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-Cl05","2.pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
mammal.combined$Cl05 <- Idents(mammal.combined)


#Idents(mammal.combined) <- mammal.combined$ID3

#saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined.rds",sep=""))
Idents(mammal.combined) <- mammal.combined$Cells

newID <- as.character(mammal.combined$Cells)
#newID[which(mammal.combined$Cl05%in%c(0,10,13,15) & mammal.combined$species2=="D")] <- "Stromal fibroblasts"
newID[which(mammal.combined$Cl05%in%c(7,8) & mammal.combined$species2=="D")] <- "Ciliated"

newID[which(mammal.combined$Cl05%in%c(0,1,5) & mammal.combined$species2=="D")] <- "Unciliated epithelia 2"
newID[which(mammal.combined$Cl05%in%c(2,3,4) & mammal.combined$species2=="D")] <- "Unciliated epithelia 1"
newID[which(mammal.combined$Cl15%in%c(6) & mammal.combined$species2=="D")] <- "Unciliated epithelia 2"


newID[which(mammal.combined$Cl05%in%c(9) )] <- "Unciliated epithelia 3"
#newID[which(mammal.combined$Cl05%in%c(19) & mammal.combined$species2=="D")] <- "Hyp_d14"
#newID[which(mammal.combined$Cl05%in%c(6,8,16) & mammal.combined$species2=="D")] <- "CTB_d14"
#newID[which(mammal.combined$Cl05%in%c(5,9,12) & mammal.combined$species2=="D")] <- "STB_d14"
#newID[which(mammal.combined$Cl05%in%c(9) & mammal.combined$species2=="D")] <- "EVT_d14"
newID[which(mammal.combined$Cl05%in%c(6,10) & mammal.combined$species2=="D")] <- "EmDisc_d14"

Idents(mammal.combined) <- newID

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- cols[colind]


species3 <- as.character(mammal.combined$species1)

species3[which(species3=="Ref1")] <- "Ref"
species3[which(species3=="Ref2")] <- "Ref"
species3[which(species3=="Ref3")] <- "Ref"
species3[which(species3=="Ref4")] <- "Ref"
species3[which(species3=="Ref5")] <- "Ref"
species3[which(species3=="Ref6")] <- "Ref"

mammal.combined$species3 <- species3

saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined_ano1.rds",sep=""))
	


#DimPlot(mammal.combined,   reduction = "umap", split.by = "species1", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"UMAP_Seurat_splitIIICol",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)
#DimPlot(mammal.combined,  reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"UMAP_SeuratIIICol",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, cols = coluse,  reduction = "umap", split.by = "species3", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_SeuratIIICol",".pdf",sep=""),width = 26, height = 10, limitsize = FALSE, useDingbats=FALSE)

DimPlot(mammal.combined, cols = coluse,  reduction = "pca", split.by = "species3", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_SeuratIIICol",".pdf",sep=""),width = 26, height = 10, limitsize = FALSE, useDingbats=FALSE)


asdsada

SIGNAL<-read.table("../Thorsten/Other/LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]
TF<- read.table("../Thorsten/Other/TF.txt",header = F)
TF <- TF$V1

EF<- read.table("../Thorsten/Other/epifactors.txt",header = F)
EF <- EF$V1

DefaultAssay(mammal.combined) <- "RNA"

mammal.combined$newID <- newID
Idents(mammal.combined) <- paste(mammal.combined$species2,mammal.combined$newID,sep="")


Ae <- AverageExpression(mammal.combined)
Ae <- Ae$RNA
#Ae$gene <- rownames(Ae)

#Ae <- readRDS(file = paste(saveext,"/Ae.rds",sep=""))

M1 <- FindMarkers(mammal.combined, ident.1=c("RefEmDisc_d9","RefEmDisc_d11","RefEmDisc_d12"),ident.2=c("RefUnciliated epithelia 2"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1=c("DEmDisc_d14"),ident.2 = c("DUnciliated epithelia 2"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Cl1 <- M1
Cl2 <- M2

Ae5 <- as.data.frame(Ae)
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_log2FC) > log2(1.2))], "Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_log2FC) > log2(1.2))], "Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_log2FC) > log2(1.2))],rownames(Cl2)[which(abs(Cl2$avg_log2FC) > log2(1.2))]),"Indi"] <- 3

Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_log2FC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_log2FC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'DEmDisc_d14'>0.1 | Ae5$'RefEmDisc_d9'>0.1 | Ae5$'RefUnciliated epithelia 2'>0.1 | Ae5$'DUnciliated epithelia 2'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_log2FC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_log2FC)>log(1.5))]),TF))
genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_log2FC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_log2FC)>log(1.5))]),EF))
genes.to.label3 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_log2FC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_log2FC)>log(1.5))]),SIGNAL1))
genes.to.label4 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_log2FC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_log2FC)>log(1.5))]),SIGNAL2))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_log2FC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_log2FC)>log(1.5))]),SIGNAL3))

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  
scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + 
geom_hline(yintercept=log2(1.2),linetype="dashed") + geom_hline(yintercept=-log2(1.2),linetype="dashed")  + 
geom_vline(xintercept = log2(1.2), linetype="dashed") + geom_vline(xintercept = -log2(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"EmDisc_d14_vs_d9-12_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +
scale_color_manual(values=c('lightgrey','#777777','#555555','black')) +
geom_hline(yintercept=log2(1.2),linetype="dashed") + geom_hline(yintercept=-log2(1.2),linetype="dashed")  +
geom_vline(xintercept = log2(1.2), linetype="dashed") + geom_vline(xintercept = -log2(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"EmDisc_d14_vs_d9-12_EF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +
scale_color_manual(values=c('lightgrey','#777777','#555555','black')) +
geom_hline(yintercept=log2(1.2),linetype="dashed") + geom_hline(yintercept=-log2(1.2),linetype="dashed")  +
geom_vline(xintercept = log2(1.2), linetype="dashed") + geom_vline(xintercept = -log2(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label3, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"EmDisc_d14_vs_d9-12_LIG1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +
scale_color_manual(values=c('lightgrey','#777777','#555555','black')) +
geom_hline(yintercept=log2(1.2),linetype="dashed") + geom_hline(yintercept=-log2(1.2),linetype="dashed")  +
geom_vline(xintercept = log2(1.2), linetype="dashed") + geom_vline(xintercept = -log2(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label4, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"EmDisc_d14_vs_d9-12_LIG3.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +
scale_color_manual(values=c('lightgrey','#777777','#555555','black')) +
geom_hline(yintercept=log2(1.2),linetype="dashed") + geom_hline(yintercept=-log2(1.2),linetype="dashed")  +
geom_vline(xintercept = log2(1.2), linetype="dashed") + geom_vline(xintercept = -log2(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"EmDisc_d14_vs_d9-12_LIG3.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



M1 <- FindMarkers(mammal.combined, ident.1=c("RefCiliated"),ident.2=c("RefUnciliated epithelia 2"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
M2 <- FindMarkers(mammal.combined, ident.1=c("DCiliated"),ident.2 = c("DUnciliated epithelia 2"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

Cl1 <- M1
Cl2 <- M2

Ae5 <- as.data.frame(Ae)
Ae5$Pval1 <- NA
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(Cl1)[which(abs(Cl1$avg_log2FC) > log2(1.2))], "Indi"] <- 1
Ae5[rownames(Cl2)[which(abs(Cl2$avg_log2FC) > log2(1.2))], "Indi"] <- 2
Ae5[intersect(rownames(Cl1)[which(abs(Cl1$avg_log2FC) > log2(1.2))],rownames(Cl2)[which(abs(Cl2$avg_log2FC) > log2(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(Cl1),"X"] <- Cl1$avg_log2FC
Ae5[rownames(Cl2),"Y"] <- Cl2$avg_log2FC
Ae5$Keep <- 0
Ae5[unique(c(rownames(Cl1)[which(Cl1$p_val_adj<0.1)],rownames(Cl2)[which(Cl2$p_val_adj<0.1)])),"Keep"] <- 1
#Filter out lowly expressed genes. In this case if log CP10k > 0.1 in any of the comparison set keep it
Ae5 <- Ae5[which((Ae5$'DCiliated'>0.1 | Ae5$'RefCiliated'>0.1 | Ae5$'RefUnciliated epithelia 2'>0.1 | Ae5$'DUnciliated epithelia 2'>0.1)  ) ,]
#Create a set of genes we want to label e.g., TFs
genes.to.label1 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_log2FC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_log2FC)>log(1.5))]),TF))
genes.to.label2 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_log2FC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_log2FC)>log(1.5))]),EF))
genes.to.label3 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_log2FC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_log2FC)>log(1.5))]),SIGNAL1))
genes.to.label4 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_log2FC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_log2FC)>log(1.5))]),SIGNAL2))
genes.to.label5 = c( intersect(c(rownames(Cl1)[which(abs(Cl1$avg_log2FC)>log(1.5))],rownames(Cl2)[which(abs(Cl2$avg_log2FC)>log(1.5))]),SIGNAL3))

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +
scale_color_manual(values=c('lightgrey','#777777','#555555','black')) +
geom_hline(yintercept=log2(1.2),linetype="dashed") + geom_hline(yintercept=-log2(1.2),linetype="dashed")  +
geom_vline(xintercept = log2(1.2), linetype="dashed") + geom_vline(xintercept = -log2(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cill_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +
scale_color_manual(values=c('lightgrey','#777777','#555555','black')) +
geom_hline(yintercept=log2(1.2),linetype="dashed") + geom_hline(yintercept=-log2(1.2),linetype="dashed")  +
geom_vline(xintercept = log2(1.2), linetype="dashed") + geom_vline(xintercept = -log2(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cil_EF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +
scale_color_manual(values=c('lightgrey','#777777','#555555','black')) +
geom_hline(yintercept=log2(1.2),linetype="dashed") + geom_hline(yintercept=-log2(1.2),linetype="dashed")  +
geom_vline(xintercept = log2(1.2), linetype="dashed") + geom_vline(xintercept = -log2(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label3, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cil_LIG1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +
scale_color_manual(values=c('lightgrey','#777777','#555555','black')) +
geom_hline(yintercept=log2(1.2),linetype="dashed") + geom_hline(yintercept=-log2(1.2),linetype="dashed")  +
geom_vline(xintercept = log2(1.2), linetype="dashed") + geom_vline(xintercept = -log2(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label4, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cil_LIG3.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +
scale_color_manual(values=c('lightgrey','#777777','#555555','black')) +
geom_hline(yintercept=log2(1.2),linetype="dashed") + geom_hline(yintercept=-log2(1.2),linetype="dashed")  +
geom_vline(xintercept = log2(1.2), linetype="dashed") + geom_vline(xintercept = -log2(1.2), linetype="dashed")
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "log FC expression dataset 1", y = "log FC expresssion dataset 2")
ggsave(filename=paste(saveext,"Cil_LIG3.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



sdasdas

DimPlot(mammal.combined,cols = coluse,  reduction = "pca", split.by = "species1", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat_splitIICol",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, cols = coluse, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_SeuratIICol",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE)

saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined_ano1.rds",sep=""))



#mammal.combined <- readRDS(file = paste(saveext,"Seurat_combined.rds",sep=""))

#DefaultAssay(mammal.combined) <- "RNA"


DimPlot(mammal.combined, reduction = "pca", split.by = "species1", label = TRUE, repel = TRUE,dims=c(1,3))
ggsave(filename=paste(saveext,"PCA_13_Seurat_splitII",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)

DimPlot(mammal.combined, reduction = "pca", split.by = "species1", label = TRUE, repel = TRUE,dims=c(2,3))
ggsave(filename=paste(saveext,"PCA_23_Seurat_splitII",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "pca", split.by = "species1", label = TRUE, repel = TRUE,dim=c(1,4))
ggsave(filename=paste(saveext,"PCA_14_Seurat_splitII",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "pca", split.by = "species1", label = TRUE, repel = TRUE,dims=c(2,4))
ggsave(filename=paste(saveext,"PCA_24_Seurat_splitII",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "pca", split.by = "species1", label = TRUE, repel = TRUE,dims=c(3,4))
ggsave(filename=paste(saveext,"PCA_34_Seurat_splitII",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)


DimPlot(mammal.combined, reduction = "umap", split.by = "species1", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat_splitII",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_SeuratII",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE)

DimPlot(mammal.combined, reduction = "pca", split.by = "species1", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat_splitII",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_SeuratII",".pdf",sep=""),width = 20, height = 10, limitsize = FALSE, useDingbats=FALSE)


#mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20, n.components = 3L)

#DimPlot(mammal.combined, reduction = "pca", split.by = "species1", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"PCA_Seurat_splitII",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)

#DimPlot(mammal.combined, reduction = "umap", split.by = "species1", label = TRUE, repel = TRUE, dims=c(1,2))
#ggsave(filename=paste(saveext,"UMAP_12_Seurat_splitII",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)
#DimPlot(mammal.combined, reduction = "umap", split.by = "species1", label = TRUE, repel = TRUE,dims=c(1,3))
#ggsave(filename=paste(saveext,"UMAP_13_Seurat_splitII",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)
#DimPlot(mammal.combined, reduction = "umap", split.by = "species1", label = TRUE, repel = TRUE,dims=c(2,3))
#ggsave(filename=paste(saveext,"UMAP_23_Seurat_splitII",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)



DefaultAssay(mammal.combined) <- "RNA"


#AvExp <- FindAllMarkers(mammal.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#write.csv(as.data.frame(AvExp$RNA),file="AverageExpComb.csv",quote=FALSE)

FeaturePlot(mammal.combined, features = "SOX17", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_SOX17.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "SOX2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_SOX2.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "POU5F1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_POU5F1.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NANOG", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_NANOG.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "DAPP1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_DAPP1.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "KLF17", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_KLF17.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "T", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_T.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "MIXL1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_MIXL1.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "EOMES", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_EOMES.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "HNF1B", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_HF1B.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "GSG1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GSG1.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "GATM", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GATM.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "NOTO", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_NOTO.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "FXYD3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_FXYD3.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)


FeaturePlot(mammal.combined, features = "PADI1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_PADI1.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "SP100", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_SP100.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "MST1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_MST1.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)


FeaturePlot(mammal.combined, features = "RPS6KA5", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "GADD45G", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GADD45G.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "CGA", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_CGA.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)


FeaturePlot(mammal.combined, features = "GATA2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GATA2.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "GATA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GATA3.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "TFAP2C", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_TFAP2C.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "TFAP2A", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_TFAP2A.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)


FeaturePlot(mammal.combined, features = "GATA6", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GATA6.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "SOX17", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_SOX17.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "APOA1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_APOA1.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "PDGFRA", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_PDGFRA.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "OTX2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_OTX2.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "GSC", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GSC.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "LEFTY2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_LEFTY2.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "TTR", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_TTR.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "SNAI2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_SNAI2.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "HGF", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_HGF.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "HAND1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_HAND1.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "ISL1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_ISL1.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "VTCN1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_VTCN1.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "MESP1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_MESP1.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "NANOS3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_NANOS3.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "IGFBP1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_IGFBP1.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "IGFBP2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_IGFBP2.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "IGFBP6", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_IGFBP6.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "PRL", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_PRL.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)

#FeaturePlot(mammal.combined, features = "CD73", split.by = "species1", cols =  c("lightgrey", "darkblue"))
#ggsave(filename=paste(saveext,"/Markers/Marker_CD73.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)


FeaturePlot(mammal.combined, features = "EPCAM", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_EPCAM.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "UCA1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_UCA1.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "FOXJ1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_FOXJ1.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)



sadsadsadsada

DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat",".pdf",sep=""),width = 22, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat",".pdf",sep=""),width = 22, height = 10, limitsize = FALSE, useDingbats=FALSE)
DefaultAssay(mammal.combined) <- "integrated"


mammal.combined <- FindClusters(mammal.combined, resolution = 1.0)

DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_SeuratCl",".pdf",sep=""),width = 22, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_SeuratCl",".pdf",sep=""),width = 22, height = 10, limitsize = FALSE, useDingbats=FALSE)

DefaultAssay(mammal.combined) <- "RNA"
Markers <- FindAllMarkers(mammal.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(as.data.frame(Markers),file="AverageExpComb.csv",quote=FALSE)



#Dp1 <- subset(Dp1,idents=c("dP1","dP2","dS2","dS3","Epi1","Epi2"),invert=TRUE)
#Dp2 <- subset(Dp2,idents=c("dP1","dP2","dS2","dS3","Epi1","Epi2"),invert=TRUE)
#Dp3 <- subset(Dp3,idents=c("dP1","dP2","dS2","dS3","Epi1","Epi2"),invert=TRUE)
#Dp4 <- subset(Dp4,idents=c("dP1","dP2","dS2","dS3","Epi1","Epi2"),invert=TRUE)
#Dp5 <- subset(Dp5,idents=c("dP1","dP2","dS2","dS3","Epi1","Epi2"),invert=TRUE)
#Dp6 <- subset(Dp6,idents=c("dP1","dP2","dS2","dS3","Epi1","Epi2"),invert=TRUE)


#mammal.anchors <- FindIntegrationAnchors(object.list = list(D1,D2,D3,D3B,D6,Dp1,Dp2,Dp3,Dp4,Dp5,Dp6), dims = 1:20, anchor.features = 4000)
#mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
#DefaultAssay(mammal.combined) <- "integrated"
#mammal.combined <- ScaleData(mammal.combined)
#mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
#mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

#saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined-EmbOnly.rds",sep=""))

#DefaultAssay(mammal.combined) <- "RNA"


#D1 <- subset(mammal.combined,idents=c("D"))
#D2 <- subset(mammal.combined,idents=c("Ref"))
#Idents(D1) <- D1$Cells
#Idents(D2) <- D2$Cells
#DefaultAssay(D1) <- "RNA"
#DefaultAssay(D2) <- "RNA"



mammal.combined <- readRDS(file = paste(saveext,"Seurat_combined-EmbOnly.rds",sep=""))


DimPlot(mammal.combined, reduction = "pca", split.by = "species1", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat_split-EmbOnly",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "umap", split.by = "species1", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat_split-EmbOnly",".pdf",sep=""),width = 120, height = 10, limitsize = FALSE, useDingbats=FALSE)

DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-EmbOnly",".pdf",sep=""),width = 22, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-EmbOnly",".pdf",sep=""),width = 22, height = 10, limitsize = FALSE, useDingbats=FALSE)

DefaultAssay(mammal.combined) <- "integrated"

mammal.combined <- FindClusters(mammal.combined, resolution = 1.0)

DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-EmbOnlyCl",".pdf",sep=""),width = 22, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-EmbOnlyCl",".pdf",sep=""),width = 22, height = 10, limitsize = FALSE, useDingbats=FALSE)





#mammal.anchors <- FindIntegrationAnchors(object.list = list(D1,D2,D3,D3B,D6), dims = 1:20, anchor.features = 4000)
#mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
#DefaultAssay(mammal.combined) <- "integrated"
#mammal.combined <- ScaleData(mammal.combined)
#mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
#mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

#saveRDS(mammal.combined, file = paste(saveext,"Seurat_combined-NewDataOnly.rds",sep=""))


mammal.combined <- readRDS(file = paste(saveext,"Seurat_combined-NewDataOnly.rds",sep=""))


DimPlot(mammal.combined, reduction = "pca", split.by = "species1", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat_split-NewOnly",".pdf",sep=""),width = 50, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "umap", split.by = "species1", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat_split-NewOnly",".pdf",sep=""),width = 50, height = 10, limitsize = FALSE, useDingbats=FALSE)


DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-NewOnly",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-NewOnly",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)

DefaultAssay(mammal.combined) <- "integrated"

mammal.combined <- FindClusters(mammal.combined, resolution = 1.0)

DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"PCA_Seurat-NewOnlyCl",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"UMAP_Seurat-NewOnlyCl",".pdf",sep=""),width = 12, height = 10, limitsize = FALSE, useDingbats=FALSE)


DefaultAssay(mammal.combined) <- "RNA"


Markers <- FindAllMarkers(mammal.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(as.data.frame(Markers),file="AverageExp.csv",quote=FALSE)

FeaturePlot(mammal.combined, features = "SOX17", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_SOX17-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "SOX2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_SOX2-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "POU5F1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_POU5F1-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NANOG", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_NANOG-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "DAPP1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_DAPP1-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "KLF17", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_KLF17-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "T", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_T-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "MIXL1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_MIXL1-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "EOMES", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_EOMES-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "HNF1B", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_HF1B-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "GSG1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GSG1-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "GATM", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GATM-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "NOTO", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_NOTO-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "FXYD3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_FXYD3-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)


FeaturePlot(mammal.combined, features = "PADI1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_PADI1-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "SP100", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_SP100-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "MST1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_MST1-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)


FeaturePlot(mammal.combined, features = "RPS6KA5", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "GADD45G", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GADD45G-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "NCOA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_RPS6KA5-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "CGA", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_CGA-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)


FeaturePlot(mammal.combined, features = "GATA2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GATA2-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "GATA3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GATA3-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "TFAP2C", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_TFAP2C-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "TFAP2A", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_TFAP2A-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)


FeaturePlot(mammal.combined, features = "GATA6", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GATA6-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "SOX17", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_SOX17-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "APOA1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_APOA1-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "PDGFRA", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_PDGFRA-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "OTX2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_OTX2-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "GSC", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_GSC-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "LEFTY2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_LEFTY2-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "TTR", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_TTR-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "SNAI2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_SNAI2-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "HGF", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_HGF-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "HAND1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_HAND1-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "ISL1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_ISL1-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "VTCN1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_VTCN1-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "MESP1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_MESP1-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "NANOS3", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_NANOS3-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
#DefaultAssay(mammal.combined) <- "RNA"

FeaturePlot(mammal.combined, features = "IGFBP1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_IGFBP1-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "IGFBP2", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_IGFBP2-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "IGFBP6", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_IGFBP6-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

FeaturePlot(mammal.combined, features = "PRL", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_PRL-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)

#FeaturePlot(mammal.combined, features = "CD73", split.by = "species1", cols =  c("lightgrey", "darkblue"))
#ggsave(filename=paste(saveext,"/Markers/Marker_CD73.pdf",sep=""),width = 120, height = 10,limitsize = FALSE)


FeaturePlot(mammal.combined, features = "EPCAM", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_EPCAM-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "UCA1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_UCA1-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)
FeaturePlot(mammal.combined, features = "FOXJ1", split.by = "species1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Markers/Marker_FOXJ1-NewOnlyCl.pdf",sep=""),width = 50, height = 10,limitsize = FALSE)


#D1 <- subset(mammal.combined,idents=c("D"))
#D2 <- subset(mammal.combined,idents=c("Ref"))
#Idents(D1) <- D1$Cells
#Idents(D2) <- D2$Cells
#DefaultAssay(D1) <- "RNA"
#DefaultAssay(D2) <- "RNA"


