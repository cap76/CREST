
cType <- c("Hyp/Am","Hyp_d14","Am/EmDisc","ExMes","EmDisc/Am1","EmDisc/Am2","EmDisc/Am3","EmDisc/Am4","SYS_CS6/7","EmDisc_CS6/7","PGC_CS6/7","Am","EmDisc","PGC","VE","SYS","ExMes","None","Other","Gland_CS5","Gland_CS6","Gland","Stroma_CS5","Stroma","Remodelled_CS5","Remodelled_CS6","Remodelled","Stalk_CS6","2307","2308","Other","Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS2","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","Stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","Am","Amnoid_bead","BMP_MEF","BMP_noMEF","EmD","EmDisc","ActA_MEF","ActA_noMEF","SB43_MEF","CHIR_MEF","FGF_noMEF","Am_CS5_PGC","Am_CS5_ExMesPGC","EmDisc_CS5_Am","PGC_CS6","Stalk_CS7_PGC","EmDisc_CS7_PGC","EmDisc_CS7_Am","Am_CS6_EmDisc","EmDisc_CS5_Gast","EmDisc_CS6_Am","EmDisc_CS6_Gast","EmDisc_CS7_Gast","Stalk_CS7_PGC","Gland_CS6_","ReStroma_CS5","ReStroma_CS6","ReStroma_CS7","Stroma_CS5","Stroma_CS6","Stroma_CS7","VE_CS6","VE_CS7","Am_CS7_EmDisc","EmDisc_Gast_CS6","EmDisc_Gast_CS7","Other")
BaseCol <- c("black","#F04C04","#0233BF","#967700","#0c9cf5","#877bd6","#0767DA","#5F54C7","#d17600","#0767DA","#E6E600","#877bd6","#0c9cf5","#E6E600","#F04C04","#E68600","#e6c800","#B3B2B2","#f5f2d0","#B3B2B2","#A5A4A3","#969593","#DFDFDF","#D0D0D0","#878684","#797775","#6A6866","#754C24","lightgrey","lightgrey","lightgrey","#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#603813","#E6E600","#BFBF04","#999903","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#7b3294","#c2a5cf","#a6dba0","#008837","#ca0020","#f4a582","#92c5de","#0571b0","#d8b365","#5ab4ac","#4d4d4d","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#0233BF","#5F54C7","#0c9cf5","#0767DA","#0767DA","#0233BF","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0",'#D74404',"black","#5F54C7","#0767DA","#0233BF","#754C24,","#f5f2d0")



#mammal.combined <- readRDS(file=paste(saveext,"/BatchAll_marmosetintegration_wMapre.rds",sep=""))
#saveRDS(mammal.combined,file=paste(saveext,"/BatchAll_marmosetintegration_wpre.rds",sep=""))
mammal.combined <- readRDS(file=paste(saveext,"/BatchAll_marmosetintegration_wMa.rds",sep=""))
#mammal.combined <- readRDS(file=paste(saveext,"/BatchAll_marmosetintegration.rds",sep=""))
#mammal.combined <- readRDS(file=paste(saveext,"/BatchAll_marmosetintegration_wMa_EmDAm.rds",sep="") )
#Hyp_d14       putExMes      ExMes_d14     Am/EmDisc_d14 EmDisc_d14    Am_d14        Hyp/Am        putPGC        Am/EmDisc    
#0c9cf5	#0767DA
EmDAmList1 <- readRDS(file=paste(saveext,"/DimRed/EmDiscList1.rds",sep=""))
EmDAmList2 <- readRDS(file=paste(saveext,"/DimRed/EmDiscList2.rds",sep=""))
EmDAmList3 <- readRDS(file=paste(saveext,"/DimRed/EmDiscList3.rds",sep=""))
EmDAmList4 <- readRDS(file=paste(saveext,"/DimRed/EmDiscList4.rds",sep=""))
ExMes <- readRDS(file=paste(saveext,"/DimRed/ExMes.rds",sep=""))
ExMes2 <- readRDS(file=paste(saveext,"/DimRed/ExMes2.rds",sep=""))

PGC <- WhichCells(mammal.combined,idents=c("putPGC"))

Hyp <- readRDS(file=paste(saveext,"/DimRed/Hyp.rds",sep=""))
SYS <- readRDS(file=paste(saveext,"/DimRed/SYS.rds",sep=""))

Idents(mammal.combined,cells = WhichCells(mammal.combined,idents="putPGC") ) <- "PGC"
Idents(mammal.combined,cells = WhichCells(mammal.combined,idents="ExMes_d14") ) <- "SYS"
Idents(mammal.combined,cells = WhichCells(mammal.combined,idents="putExMes") ) <- "SYS"

Idents(mammal.combined,cells = EmDAmList1) <- "EmDisc/Am1"
Idents(mammal.combined,cells = EmDAmList2) <- "EmDisc/Am2"
Idents(mammal.combined,cells = EmDAmList3) <- "EmDisc/Am3"
Idents(mammal.combined,cells = EmDAmList4) <- "EmDisc/Am4"
Idents(mammal.combined,cells = ExMes) <- "ExMes"
Idents(mammal.combined,cells = ExMes2) <- "ExMes"
Idents(mammal.combined,cells = Hyp) <- "Hyp_d14"
Idents(mammal.combined,cells = SYS) <- "SYS"

mammal.combined$ID <- Idents(mammal.combined)

Idents(mammal.combined) <- mammal.combined$ID # <- Idents(mammal.combined)

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]


p<-DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Nakamura_umap.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)

p<-DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Nakamura_pca.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)

mammal.combined <- FindClusters(mammal.combined,resolution = 0.5)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Nakamura_umap_Cl.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)
mammal.combined$Cl <- Idents(mammal.combined)

EmDAmList1 <- colnames(mammal.combined)[which(mammal.combined$Cl=="3" & mammal.combined$Dataset=="Matteo")]
EmDAmList2 <- colnames(mammal.combined)[which(mammal.combined$Cl=="1" & mammal.combined$Dataset=="Matteo")]
EmDAmList3 <- colnames(mammal.combined)[which(mammal.combined$Cl=="4" & mammal.combined$Dataset=="Matteo")]
EmDAmList4 <- colnames(mammal.combined)[which(mammal.combined$Cl=="5" & mammal.combined$Dataset=="Matteo")]

#saveRDS(EmDAmList1,file=paste(saveext,"/DimRed/EmDiscList1.rds",sep=""))
#saveRDS(EmDAmList2,file=paste(saveext,"/DimRed/EmDiscList2.rds",sep=""))
#saveRDS(EmDAmList3,file=paste(saveext,"/DimRed/EmDiscList3.rds",sep=""))
#saveRDS(EmDAmList4,file=paste(saveext,"/DimRed/EmDiscList4.rds",sep=""))
#saveRDS(ExMes,file=paste(saveext,"/DimRed/ExMes.rds",sep=""))

#saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)=="0")],file=paste(saveext,"/DimRed/ExMes2.rds",sep=""))



saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)=="3" & mammal.combined$ID=="Hyp/Am" )],file=paste(saveext,"/DimRed/Hyp.rds",sep=""))
saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)=="2" & mammal.combined$ID=="Hyp/Am" )],file=paste(saveext,"/DimRed/SYS.rds",sep=""))

saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)%in%c("4","0") )],file=paste(saveext,"/DimRed/ExMes_Cl4_0.rds",sep=""))
saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)%in%c("3") )],file=paste(saveext,"/DimRed/Hyp_Cl3.rds",sep=""))
saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)%in%c("2") )],file=paste(saveext,"/DimRed/SYS_Cl2.rds",sep=""))
saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)%in%c("5","6","7") )],file=paste(saveext,"/DimRed/EmDisc_Cl5_6_7.rds",sep=""))
saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)%in%c("9") )],file=paste(saveext,"/DimRed/EmDisc_Cl9.rds",sep=""))






saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)%in%c("4","0") )],file=paste(saveext,"/DimRed/ExMes_Cl4_0.rds",sep=""))
saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)%in%c("3") )],file=paste(saveext,"/DimRed/Hyp_Cl3.rds",sep=""))
saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)%in%c("2") )],file=paste(saveext,"/DimRed/SYS_Cl2.rds",sep=""))

saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)%in%c("0") )],file=paste(saveext,"/DimRed/ExMes_Cl0.rds",sep=""))
saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)%in%c("4") )],file=paste(saveext,"/DimRed/ExMes_Cl4.rds",sep=""))

saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)%in%c("5") )],file=paste(saveext,"/DimRed/EmDisc_Cl5.rds",sep=""))
saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)%in%c("6") )],file=paste(saveext,"/DimRed/EmDisc_Cl6.rds",sep=""))
saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)%in%c("7") )],file=paste(saveext,"/DimRed/EmDisc_Cl7.rds",sep=""))
saveRDS(colnames(mammal.combined)[which(mammal.combined$Dataset=="Matteo" & Idents(mammal.combined)%in%c("9") )],file=paste(saveext,"/DimRed/EmDisc_Cl9.rds",sep=""))



ExMes <- c("AACACACCAATGAACA-1_3", "AAGCGTTAGCTTGTTG-1_3", "AGCGTATCAAGGCGTA-1_3",
"AGGCATTTCCCGGTAG-1_3", "AGGTCATAGATGCCGA-1_3", "AGTAGCTCAGGACATG-1_3",
"ATGCATGTCTCTTAAC-1_3", "ATTATCCGTCCATACA-1_3", "CAAGACTTCTTCCGTG-1_3",
"CACGTGGAGAAGCGCT-1_3", "CAGTGCGAGTAGCCAG-1_3", "CAGTTAGAGTTTCGGT-1_3",
"CATGCCTGTCACAGAG-1_3", "CATTGCCTCCGATAAC-1_3", "CCAAGCGAGGCTCACC-1_3",
"CCACGAGCAACTTGGT-1_3", "CCTCAGTGTTCAGCTA-1_3", "CGAATTGCAAGACGGT-1_3",
"CTCGAGGCAGTTTCAG-1_3", "CTGAATGAGAGCCCAA-1_3", "GACCCAGAGGGCCAAT-1_3",
"GATCCCTAGACATACA-1_3", "GCCGTGAGTCAAGCCC-1_3", "GCCTGTTAGTCCCGGT-1_3",
"GGAATCTTCCTGGGTG-1_3", "GGGAGATTCATGAAAG-1_3", "GGGTGAAAGTATGAGT-1_3",
"GGTAATCTCGGCGATC-1_3", "GTCATTTGTATCGCGC-1_3", "GTCTCACAGCGGTAGT-1_3",
"TAACTTCGTACGACAG-1_3", "TACTTACAGGTAGTAT-1_3", "TCATCCGGTGTGGACA-1_3",
"TCATTACCAAGTATAG-1_3", "TGAATGCGTTCGGACC-1_3", "TGGGAAGAGCTGAGCA-1_3",
"TTCAGGATCTGTGTGA-1_3", "TGGGCTGGTCATCACA-1_6")

#[1] "ACTTAGGTCCATTTAC-1_2" "CACGTGGTCCCTCATG-1_2" "CCACTTGAGTAAACTG-1_2" "CCATCACCAACAAGTA-1_2" "CTGCCATTCTTACGTT-1_2"
#[6] "TCGCTCATCTTTGCTA-1_2" "GTAACCATCATCTGTT-1_5" "AAAGAACAGTAAACGT-1_6" "AGACACTCATTCTGTT-1_6" "ATACCTTCAGCTATTG-1_6"
#[11] "CCTCACACATGCCGCA-1_6" "CGCGTGACAGCGATTT-1_6" "CGTCCATTCGCACGGT-1_6" "CTTTCGGTCACTAGCA-1_6" "GACGCTGGTTCATCGA-1_6"
#[16] "GACTCTCTCCCAGGAC-1_6" "GAGGCCTAGGATTCCT-1_6" "GTACAACGTGAGCCAA-1_6" "GTAGCTACAGCTTTCC-1_6" "TGAGGTTGTTTCCATT-1_6"
#[21] "TGCGATACAAAGGCAC-1_6" "TGGAGGAGTGTCTTGA-1_6" "TGGGCTGGTCATCACA-1_6" "ATGGAGGCACTTGGCG-1_7" "GACGTTATCCACCTCA-1_7"
#[26] "GAGTTTGAGCATGCAG-1_7" "GTAGTACAGCTCGTGC-1_7" "GTCGTAAGTACGAAAT-1_7" "TTGTTGTGTAAGCGGT-1_7" "TTTACGTCATGATGCT-1_7"
#[31] "AAGATAGAGAACCGCA-1_8" "AATAGAGGTCACTTAG-1_8" "ACCAAACCAAGCTGCC-1_8" "ATCGCCTTCTGAGTCA-1_8" "ATGAGTCTCCTTATAC-1_8"
#[36] "CACCAAACACCCATAA-1_8" "CCACACTCACTACTTT-1_8" "GCCCGAACAAACGGCA-1_8" "GGTAACTGTATGCGTT-1_8" "GTCTCACGTTGGCTAT-1_8"
#[41] "TACACCCAGTTTGTCG-1_8" "TGTACAGGTCTAGGTT-1_8" "TTCCAATCATCTGTTT-1_8"

#Amnion clustr

#[1] "ACTTTCAAGTACCCTA-1_2" "ATCCTATGTAGGCAAC-1_2" "ATCGATGGTAGCGTAG-1_2" "ATTATCCCACCCTAGG-1_2" "ATTCCCGTCCGTCAAA-1_2"
#[6] "CAAGACTGTATTGCCA-1_2" "CAGTTAGCATGAGAAT-1_2" "CCCTTAGCACGACTAT-1_2" "CCGTTCAGTTACGCCG-1_2" "CTGCTCACACTCACTC-1_2"
#[11] "GCATCTCTCTTCCCAG-1_2" "GTAGATCAGCTGGCCT-1_2" "GTATTGGAGGTTGGAC-1_2" "TAATCTCGTTCTCCTG-1_2" "TCACTATTCGCTCTAC-1_2"
#[16] "TCGCACTAGCCAAGGT-1_2" "TTACCATTCTCCAATT-1_2" "TTCACCGGTTCCTAGA-1_2" "TTTCAGTGTCAACATC-1_2" "AAACCCAAGTAGCAAT-1_3"
#[21] "AACCAACCACGTCATA-1_3" "AACGTCACAGGTCAAG-1_3" "AACTTCTAGGACATCG-1_3" "AAGGAATTCACTTGTT-1_3" "AATCGACAGACGAAGA-1_3"
#[26] "AATCGTGAGCTCGACC-1_3" "AATGACCGTTCGGACC-1_3" "AATGCCACACCAGTAT-1_3" "ACACCAACAACCGCTG-1_3" "ACACTGAGTCTGGTTA-1_3"
#[31] "ACCTGTCAGTGAGTGC-1_3" "ACGGTCGGTAGCGTAG-1_3" "ACGGTTAAGGACTATA-1_3" "ACGTAGTGTGACGTCC-1_3" "ACTATTCAGCACCGTC-1_3"
#[36] "ACTCTCGTCCACTAGA-1_3" "ACTGATGAGTCATTGC-1_3" "AGACCCGCAGAGGAAA-1_3" "AGACCCGCATCGCTAA-1_3" "AGCGCCACACGGCACT-1_3"
#[41] "ATATCCTTCTGAGATC-1_3" "ATCGCCTTCCGGTTCT-1_3" "ATGAGTCAGTCAATCC-1_3" "ATTCACTGTAATTGGA-1_3" "ATTCGTTCAGGAACCA-1_3"
#[46] "CAAAGAAGTCATCAGT-1_3" "CAACAACGTATTAAGG-1_3" "CAAGACTAGGACGCTA-1_3" "CAAGAGGAGTTGGAGC-1_3" "CACAACACAGGACATG-1_3"
#[51] "CACGGGTCATCCTGTC-1_3" "CAGCAATAGTTCTACG-1_3" "CATCCGTTCTTGGCTC-1_3" "CATGCCTCACGTACTA-1_3" "CCCTAACCACCGTGAC-1_3"
#[56] "CCGGACAGTATGAGGC-1_3" "CCGGGTACAATCGCAT-1_3" "CCTAAGACAGGTCCGT-1_3" "CCTCAACTCGTCACCT-1_3" "CCTCAACTCTTGCAAG-1_3"
#[61] "CCTGCATCAAGTATCC-1_3" "CCTTGTGAGCCTGGAA-1_3" "CGCCATTCACATTCTT-1_3" "CGGGCATCAACCAGAG-1_3" "CGGGTCATCAGTCACA-1_3"
#[66] "CGGTCAGGTTGTAAAG-1_3" "CGTAATGGTAAGATAC-1_3" "CTAACCCAGTAGAGTT-1_3" "CTACGGGAGCCGTCGT-1_3" "CTATCCGTCCTGTTGC-1_3"
#[71] "CTCAGAAAGTGATAGT-1_3" "CTCCAACTCACCCTCA-1_3" "CTCTCAGTCGCCAATA-1_3" "CTGATCCGTTGCGTAT-1_3" "CTGTATTGTTAAGCAA-1_3"
#[76] "CTTCAATAGGATACGC-1_3" "CTTTCAACAACAGTGG-1_3" "GAAGAATAGGAGAGGC-1_3" "GACCAATGTGTCATGT-1_3" "GACCCTTAGGACTTCT-1_3"
#[81] "GACCTTCCACTTGGCG-1_3" "GAGGGTAAGGACTATA-1_3" "GAGTCATCACGGGCTT-1_3" "GAGTTACTCTCCGCAT-1_3" "GATAGCTCACTTGAAC-1_3"
#[86] "GATCCCTCACTACGGC-1_3" "GATGAGGAGTTCTACG-1_3" "GATGATCAGAAGAGCA-1_3" "GCCAGCATCACTACGA-1_3" "GCGGAAAAGTTGCGCC-1_3"
#[91] "GCGTTTCAGTATAACG-1_3" "GCTGAATCAATGTCTG-1_3" "GGATCTACAAGTATAG-1_3" "GGGACAACAATCGCGC-1_3" "GGTGAAGAGGACTTCT-1_3"
#[96] "GGTGAAGAGTGGTGGT-1_3" "GGTGATTGTGACTCGC-1_3" "GGTTAACTCCGGCAGT-1_3" "GTAACCAGTCAAGGCA-1_3" "GTCAAACTCTCAGTCC-1_3"
#[101] "GTCGAATCAGCTACAT-1_3" "GTCGAATCATTGTCGA-1_3" "GTCTCACGTATACAGA-1_3" "GTGAGCCGTCCTGGTG-1_3" "GTGCACGCACAGCATT-1_3"
#[106] "GTTGCGGTCATGCGGC-1_3" "GTTTACTAGAGGCGTT-1_3" "TAACACGTCACGTAGT-1_3" "TAACTTCGTACAAAGT-1_3" "TAAGTCGGTCACCTTC-1_3"
#[111] "TATACCTCACCAATTG-1_3" "TATACCTTCCGTGGGT-1_3" "TATATCCGTGGGTCAA-1_3" "TCAAGACTCTCCAATT-1_3" "TCAAGTGGTCACGTGC-1_3"
#[116] "TCACATTAGGAAGTAG-1_3" "TCAGCCTCATCCGGTG-1_3" "TCAGTTTTCTCCAATT-1_3" "TCATCCGGTGCTATTG-1_3" "TCATGGAAGTCTGGAG-1_3"
#[121] "TCATGTTGTCGACGCT-1_3" "TCATGTTTCGAGTCCG-1_3" "TCCTCCCCAAACGAGC-1_3" "TCTACCGGTATTCCGA-1_3" "TCTCCGACATTGCAAC-1_3"
#[126] "TCTTAGTGTGGACTAG-1_3" "TGAGACTAGACCACGA-1_3" "TGCGACGCACATTGTG-1_3" "TGGATGTAGGGTCAAC-1_3" "TGTTACTAGCCTCTCT-1_3"
#[131] "TGTTACTCATAGTCAC-1_3" "TGTTCTACAGCGTGAA-1_3" "TTACAGGGTCTTCTAT-1_3" "TTAGGGTGTTTCCATT-1_3" "TTCACCGCACACCTAA-1_3"
#[136] "TTCATGTGTTGAATCC-1_3" "TTCCAATGTGCGGTAA-1_3" "TTCCTCTGTTCCGTTC-1_3" "TTCTCTCAGACTACGG-1_3" "TTGCATTAGGGCATGT-1_3"
#[141] "TTGGGCGTCATCTACT-1_3" "AGCGCTGTCGCGGACT-1_5" "CAGCCAGCAAGCGAAC-1_5" "CCTCTAGTCATTGCTT-1_5" "TATTTCGCACCTGAAT-1_5"
#[146] "CAACAACGTAATGCGG-1_6" "CTGAGCGCAGCTTTCC-1_6" "GATGACTGTTGCGAAG-1_6" "GATTGGTTCGCTAGCG-1_6" "TCTTTGAAGCAATAGT-1_6"
#[151] "TTTCATGTCCGCAACG-1_6" "GATTCTTGTATCACGT-1_7"

#Odd EmDisc

#[1] "AAACGAAAGCAACTCT-1_3" "AACCAACTCATTGAGC-1_3" "AAGATAGTCTGGTTGA-1_3" "AAGTACCAGATTGATG-1_3" "AATTCCTTCTAGTTCT-1_3"
#[6] "ACGATCACACCGCTGA-1_3" "ACGTACAAGCGCATCC-1_3" "ACTGCAAAGACCAACG-1_3" "ACTTCGCTCGCGAAGA-1_3" "AGCCAGCAGTGACACG-1_3"
#[11] "AGGACGAGTTAACAGA-1_3" "AGGGCTCAGACGTCGA-1_3" "AGGTCATCACTCAGAT-1_3" "AGGTGTTTCCTTCTTC-1_3" "AGTCAACTCACCATAG-1_3"
#[16] "AGTCTCCGTAGCTTGT-1_3" "AGTGCCGAGCTGACAG-1_3" "AGTTCGAAGTGCAACG-1_3" "ATCCACCGTGAATGTA-1_3" "ATCGGCGGTAGGCAAC-1_3"
#[21] "ATCGGCGTCTGCACCT-1_3" "ATGACCACAACCCGCA-1_3" "ATTCCCGGTTCCACAA-1_3" "CAAGACTTCATCCTAT-1_3" "CAAGAGGAGTAGGCCA-1_3"
#[26] "CAATCGATCCCTTGGT-1_3" "CACACAAAGGAGGGTG-1_3" "CACATGACAGCAGTTT-1_3" "CACATGAGTTCTCCTG-1_3" "CACCGTTTCGGTCACG-1_3"
#[31] "CAGCAATCATTGTCGA-1_3" "CAGCAATGTGTAAACA-1_3" "CATACTTAGGGCTAAC-1_3" "CATAGACGTACGGATG-1_3" "CATTCTACATCGGATT-1_3"
#[36] "CCTCAGTAGGGCGAGA-1_3" "CCTTCAGCACATCCCT-1_3" "CGAGTGCAGCGCATCC-1_3" "CGAGTGCTCGATGCTA-1_3" "CGAGTTAGTGTAGTGG-1_3"
#[41] "CGGAACCCAACGACAG-1_3" "CGGAGAACACTCAAGT-1_3" "CGTTCTGTCCATTTGT-1_3" "CTACCCACAGCGTTTA-1_3" "CTCAGGGTCAACGTGT-1_3"
#[46] "CTCCACAAGAGGCTGT-1_3" "CTCCACAAGTGGATAT-1_3" "CTCCCTCCAGACCATT-1_3" "CTCTCAGAGTCCGCCA-1_3" "CTGCCATCAAAGCGTG-1_3"
#[51] "CTGTGAATCAGACATC-1_3" "CTTTCAATCATAGCAC-1_3" "GACTATGAGCTAGATA-1_3" "GAGACCCGTCGTACTA-1_3" "GAGTTTGAGGAGTACC-1_3"
#[56] "GATCACATCTAGCATG-1_3" "GCACGGTAGGGACACT-1_3" "GCACGTGAGGTGGCTA-1_3" "GCCATTCCATGATGCT-1_3" "GCGGAAAGTTGTTGAC-1_3"
#[61] "GCTGCAGTCCGTGCGA-1_3" "GCTTCACCATGTCAGT-1_3" "GGAACCCTCCAATCTT-1_3" "GGACGTCAGGGCATGT-1_3" "GGATGTTTCATACGGT-1_3"
#[66] "GGCACGTTCATACGGT-1_3" "GGGCCATCATGAGGGT-1_3" "GGGTAGACATCATTTC-1_3" "GGGTATTTCCACGAAT-1_3" "GGTGTTAGTTTCGATG-1_3"
#[71] "GGTTAACCAGATCATC-1_3" "GTAGTACTCATCGCCT-1_3" "GTATTGGTCTCCTGAC-1_3" "GTATTTCTCGTAGAGG-1_3" "GTCAGCGGTCAAGTTC-1_3"
#[76] "GTGAGTTAGCCTTGAT-1_3" "GTGAGTTGTACATACC-1_3" "GTTCATTGTTCCGGTG-1_3" "TAAGTCGGTCCATAGT-1_3" "TACAGGTTCACGGACC-1_3"
#[81] "TACGGGCAGGTAATCA-1_3" "TACGGGCCATGGGTCC-1_3" "TAGTGCAAGATACCAA-1_3" "TATCGCCTCTCTGGTC-1_3" "TATTGCTCACAGAGAC-1_3"
#[86] "TCACAAGAGGACGCAT-1_3" "TCATGTTTCGTGCACG-1_3" "TCCGATCCACGCGCAT-1_3" "TCCGGGAAGGAATGTT-1_3" "TCCTGCAGTGGCTTGC-1_3"
#[96] "TGGTGATTCACGGACC-1_3" "TGTCAGACATACTGAC-1_3" "TTACCATTCATGTCTT-1_3" "TTCCTCTAGCAATAAC-1_3" "TTGAGTGGTACGAAAT-1_3"
#[91] "TCTACCGTCGATACAC-1_3" "TCTTCCTAGGGAGGCA-1_3" "TGAGCATTCATGAGTC-1_3" "TGCTCCACACGCTGCA-1_3" "TGGGAAGAGAATTGCA-1_3"
#[101] "TTGTTTGCAGACTCTA-1_3" "TTTGTTGTCTACAGGT-1_3" "AGTAGCTCACTACCGG-1_5" "CTGCCTACAACAAGTA-1_5" "AGACCATCAGTGGTGA-1_6"
#[106] "ATGCCTCAGTCATCGT-1_6" "CGGAATTGTGGGTCAA-1_6" "CTGTATTCAAGGCGTA-1_6" "TCGTAGACACTGATTG-1_6" "TGGATGTGTTAACAGA-1_6"
#[111] "TTTGACTAGGACATCG-1_6"












#mammal.combined <- readRDS(file=paste(saveext,"/BatchAll_marmosetintegration_wMapre.rds",sep=""))
#saveRDS(mammal.combined,file=paste(saveext,"/BatchAll_marmosetintegration_wpre.rds",sep=""))
mammal.combined <- readRDS(file=paste(saveext,"/BatchAll_marmosetintegration_wMa.rds",sep=""))


ExMes <- readRDS(file=paste(saveext,"/DimRed/ExMes_Cl4_0.rds",sep=""))
Hyp <- readRDS(file=paste(saveext,"/DimRed/Hyp_Cl3.rds",sep=""))
SYS <- readRDS(file=paste(saveext,"/DimRed/SYS_Cl2.rds",sep=""))
EmDAmList1 <- readRDS(file=paste(saveext,"/DimRed/EmDisc_Cl5.rds",sep=""))
EmDAmList2 <- readRDS(file=paste(saveext,"/DimRed/EmDisc_Cl6.rds",sep=""))
EmDAmList3 <- readRDS(file=paste(saveext,"/DimRed/EmDisc_Cl7.rds",sep=""))

EmDAmList4 <- readRDS(file=paste(saveext,"/DimRed/EmDisc_Cl9.rds",sep=""))


#PGC <- WhichCells(mammal.combined,idents="putPGC") ) 

Idents(mammal.combined,cells = ExMes) <- "ExMes"
Idents(mammal.combined,cells = Hyp) <- "Hyp_d14"
Idents(mammal.combined,cells = SYS) <- "SYS"
Idents(mammal.combined,cells = EmDAmList1) <- "EmDisc/Am1"
Idents(mammal.combined,cells = EmDAmList2) <- "EmDisc/Am2"
Idents(mammal.combined,cells = EmDAmList3) <- "EmDisc/Am3"

Idents(mammal.combined,cells = EmDAmList4) <- "EmDisc/Am4"
Idents(mammal.combined,cells = PGC) <- "PGC"


mammal.combined$ID <- Idents(mammal.combined)

Idents(mammal.combined) <- mammal.combined$ID # <- Idents(mammal.combined)

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]


p<-DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Nakamura_umap.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)

p<-DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Nakamura_pca.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)



p<-DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Nakamura_pca_ID.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)


mammal.combined$AllID <- Idents(mammal.combined)
Idents(mammal.combined) <- mammal.combined$Dataset
mammal.combined <- subset(mammal.combined,idents="Matteo")
Idents(mammal.combined) <- mammal.combined$AllID
DefaultAssay(mammal.combined) <- "RNA"

AllCl <- FindAllMarkers(mammal.combined, test.use = "MAST", only.pos = TRUE)
write.table(AllCl,file=paste(saveext,'AllCl.csv',sep=""),sep=",")


#########


#mammal.combined <- readRDS(file=paste(saveext,"/BatchAll_marmosetintegration_wMapre.rds",sep=""))
#saveRDS(mammal.combined,file=paste(saveext,"/BatchAll_marmosetintegration_wpre.rds",sep=""))
mammal.combined <- readRDS(file=paste(saveext,"/BatchAll_marmosetintegration_wMa_EmDAm.rds",sep=""))


ExMes <- readRDS(file=paste(saveext,"/DimRed/ExMes_Cl4_0.rds",sep=""))
Hyp <- readRDS(file=paste(saveext,"/DimRed/Hyp_Cl3.rds",sep=""))
SYS <- readRDS(file=paste(saveext,"/DimRed/SYS_Cl2.rds",sep=""))
EmDAmList1 <- readRDS(file=paste(saveext,"/DimRed/EmDisc_Cl5.rds",sep=""))
EmDAmList2 <- readRDS(file=paste(saveext,"/DimRed/EmDisc_Cl6.rds",sep=""))
EmDAmList3 <- readRDS(file=paste(saveext,"/DimRed/EmDisc_Cl7.rds",sep=""))

EmDAmList4 <- readRDS(file=paste(saveext,"/DimRed/EmDisc_Cl9.rds",sep=""))


#PGC <- WhichCells(mammal.combined,idents="putPGC") ) 

Idents(mammal.combined,cells = ExMes) <- "ExMes"
Idents(mammal.combined,cells = Hyp) <- "Hyp_d14"
Idents(mammal.combined,cells = SYS) <- "SYS"
Idents(mammal.combined,cells = EmDAmList1) <- "EmDisc/Am1"
Idents(mammal.combined,cells = EmDAmList2) <- "EmDisc/Am2"
Idents(mammal.combined,cells = EmDAmList3) <- "EmDisc/Am3"

Idents(mammal.combined,cells = EmDAmList4) <- "EmDisc/Am4"
Idents(mammal.combined,cells = PGC) <- "PGC"


mammal.combined$ID <- Idents(mammal.combined)

Idents(mammal.combined) <- mammal.combined$ID # <- Idents(mammal.combined)

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]


p<-DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Nakamura_umap.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)

p<-DimPlot(mammal.combined,  cols = coluse, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Nakamura_pca.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)

mammal.combined <- FindClusters(mammal.combined,resolution = 0.5)
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Nakamura_pca_Cl.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)

p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Nakamura_pca_Cl_ID.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)

MiniList2 <- c("DNMT3A","DPPA3","OTX2","NOTUM","JAM3","TACSTD2","ENPEP","GATA2","GATA3","TFAP2C","TFAP2A","TEAD3","TP63","OVOL1","CCKBR","CGA","CGB5","CGB8","SDC1","HLA-C","HLA-G","ASCL2","SNAI1","VIM","FN1","ELF5","MKI67","CDH1","FSTL1","FSTL3","KRT19","KRT7","SPARC","LAMC1","CTNNAL1")
MiniList1 <- c("JARID2","CASC15","LINC01194","MT1G","MT1H","NLGN4X","RIMS2","CD9","DPPA4","CLU","NAV3","CDH1","CLDN4","PRLR","DMD","CLDN6","EZR","LGR5", "LHX1","OPHN1","SLC3A2","SLC2A3","TANC2","TCIM","COL5A2","COL6A3","TNC","SPARC","COL1A2","COL1A1","COL3A1","DCN","LUM","COL6A1","COL6A2","DPP4","CST3","DNAJC15","DNMT3A","DPPA3","OTX2","POU5F1","PDGFA","PODXL","SFRP1","SOX15","GABRP","VTCN1","WNT6","BAMBI","TFAP2C","TFAP2A","GATA3","CGA","SOX17","GATA4","CER1","LEFTY1","LEFTY2","APOA1","APOB","LINC00261","HHEX","FOXF1","PDGFRA","HAND2","BST2","GATA6","HGF","PRDM1","SNAI1","VIM","FN1")

subs1 <- which(mammal.combined$Dataset%in%c("Matteo") )

X <- GetAssayData(mammal.combined, assay = "RNA")
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList1,]
Xp <- na.omit(Xp)[,subs1]
annotation_col = data.frame(Stage = factor(droplevels( Idents(mammal.combined)[subs1])), ID = factor(( mammal.combined$ID3[subs1])) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-2, 2, length.out = 50)
pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 20,cutree_rows = 20, filename = paste(saveext,"/DimRed/LineageLabelledExtendedList",".pdf",sep="") ,width=40,height=10)




mammal.combined$Cl <- Idents(mammal.combined)

Idents(mammal.combined) <- mammal.combined$Dataset
mammal.combined <- subset(mammal.combined,idents="Matteo")
Idents(mammal.combined) <- mammal.combined$Cl
DefaultAssay(mammal.combined) <- "RNA"

AllCl <- FindAllMarkers(mammal.combined, test.use = "MAST", only.pos = TRUE)
write.table(AllCl,file=paste(saveext,'AllCl.csv',sep=""),sep=",")




mammal.combined$ID <- Idents(mammal.combined)
Idents(mammal.combined) <- mammal.combined$Dataset
mammal.combinedA <- subset(mammal.combined,idents="Matteo")
DefaultAssay(mammal.combinedA) <- "RNA"
p<-FeaturePlot(mammal.combinedA, features = "POU5F1", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Nakamura_umap_POU5F1.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p<-FeaturePlot(mammal.combinedA, features = "POU5F1", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Nakamura_pca_POU5F1.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

Idents(mammal.combinedA) <- mammal.combinedA$ID
p<-DimPlot(mammal.combinedA,  pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca__ID.pdf",sep=""),width = 8, height = 8,p,limitsize = FALSE)



#
DefaultAssay(mammal.combinedA) <- "RNA"
PrimedNaiveForm <- read_excel("/Users/christopherpenfold/Downloads/mmc3.xls")
PrimedNaive <- read_excel("/Users/christopherpenfold/Downloads/mmc2.xls")
AmnionMarkers <- read_excel("/Users/christopherpenfold/Downloads/AmnionMarkers.xlsx")
TbMarkers <- read_excel("/Users/christopherpenfold/Downloads/AnotatedAeDEupdated2.xlsx",sheet = 2) 

#TbList1 <- intersect(rownames(Dsubset2),TbMarkers$Gene[which(as.numeric(TbMarkers$'HumanSS2_CTB_vs_STB_adjpval')<0.05)])
#TbList2 <- intersect(rownames(Dsubset2),TbMarkers$Gene[which(as.numeric(TbMarkers$'HumanSS2_CTB_vs_EVT_adjpval')<0.05)])
#TbList3 <- intersect(rownames(Dsubset2),TbMarkers$Gene[which(as.numeric(TbMarkers$'HumanSS2_STB_vs_EVT_adjpval')<0.05)])
PrimedNaiveFormList <- PrimedNaiveForm$Gene[which(PrimedNaiveForm$P.Value<0.0001 & PrimedNaiveForm$vsNaive > 0)] 
PrimedNaiveFormList2 <- PrimedNaiveForm$Gene[which(PrimedNaiveForm$P.Value<0.0001 & PrimedNaiveForm$vsPrimed > 0)]                                                                                                                                       

PrimedNaiveList1 <- PrimedNaive$Gene[which(PrimedNaive$FDR<0.0001 & PrimedNaive$logFC > 0)] 
PrimedNaiveList2 <- PrimedNaive$Gene[which(PrimedNaive$FDR<0.0001 & PrimedNaive$logFC < 0)] 

hsTEList <- AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsTE")] 
hsEAm_TEList1 <- AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E + hsTE")] 
hsEAm_TEList2 <- AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E (+ hsTE)")] 
hsEAm1 <- AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E (+ cyAME-L)")] 
hsEAm2 <- AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsAME-E + cyAME-L")] 
hsEAm3 <- AmnionMarkers$Gene[which(AmnionMarkers$marker=="cyAME-L")] 
hsLAm_TEList <- AmnionMarkers$Gene[which(AmnionMarkers$marker=="hsTE + cyAME-L")] 


mammal.combinedA <- AddModuleScore(mammal.combinedA, features = list(intersect(hsEAm_TEList1,rownames(mammal.combined))), name = "AmnionTE")
p1<-FeaturePlot(mammal.combinedA, features = "AmnionTE1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/AmnionTE.pdf",sep=""),width = 8, height = 8,p1)

mammal.combinedA <- AddModuleScore(mammal.combinedA, features = list(intersect(hsEAm_TEList2,rownames(mammal.combined))), name = "AmnTE")
p1<-FeaturePlot(mammal.combinedA, features = "AmnTE1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/AmnionTE2.pdf",sep=""),width = 8, height = 8,p1)

mammal.combinedA <- AddModuleScore(mammal.combinedA, features = list(intersect(hsLAm_TEList,rownames(mammal.combined))), name = "AmTE")
p1<-FeaturePlot(mammal.combinedA, features = "AmTE1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/AmnionTE3.pdf",sep=""),width = 8, height = 8,p1)



mammal.combinedA <- AddModuleScore(mammal.combinedA, features = list(intersect(PrimedNaiveFormList,rownames(mammal.combined))), name = "FormvNaive")
p1<-FeaturePlot(mammal.combinedA, features = "FormvNaive1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/FormNaive1.pdf",sep=""),width = 8, height = 8,p1)


mammal.combinedA <- AddModuleScore(mammal.combinedA, features = list(intersect(PrimedNaiveFormList2,rownames(mammal.combined))), name = "FormvPrimed")
p1<-FeaturePlot(mammal.combinedA, features = "FormvPrimed1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/FormPrimed1.pdf",sep=""),width = 8, height = 8,p1)


mammal.combinedA <- AddModuleScore(mammal.combinedA, features = list(intersect(PrimedNaiveList1,rownames(mammal.combined))), name = "NaivevPrimed")
p1<-FeaturePlot(mammal.combinedA, features = "NaivevPrimed1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/NaivevPrimed1.pdf",sep=""),width = 8, height = 8,p1)


mammal.combinedA <- AddModuleScore(mammal.combinedA, features = list(intersect(PrimedNaiveList2,rownames(mammal.combined))), name = "PrimedvNaive")
p1<-FeaturePlot(mammal.combinedA, features = "PrimedvNaive1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PrimedvNaive1.pdf",sep=""),width = 8, height = 8,p1)


p<-FeaturePlot(mammal.combinedA, features = "MALAT1", reduction = "pca", split.by = "ID3", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_MALAT1.pdf",sep=""),width = 50, height = 8,p,limitsize = FALSE)


p<-FeaturePlot(mammal.combinedA, features = "VTCN1", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_VTCN1.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p<-FeaturePlot(mammal.combinedA, features = "WNT6", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_WNT6.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


p<-FeaturePlot(mammal.combinedA, features = "WNT6", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_WNT6.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


p<-FeaturePlot(mammal.combinedA, features = "PMAIP1", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_PMAIP1.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


p<-FeaturePlot(mammal.combinedA, features = "TDGF1", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_TDGF1.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p<-FeaturePlot(mammal.combinedA, features = "SCG3", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_SCG3.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p<-FeaturePlot(mammal.combinedA, features = "PMAIP1", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_PMAIP1.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)
p<-FeaturePlot(mammal.combinedA, features = "GPRC5B", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_GPRC5B.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p<-FeaturePlot(mammal.combinedA, features = "MT1H", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_MT1H.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p<-FeaturePlot(mammal.combinedA, features = "PTPRZ1", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_PTPRZ1.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p<-FeaturePlot(mammal.combinedA, features = "MT1G", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_MT1G.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


p<-FeaturePlot(mammal.combinedA, features = "SFRP2", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_SFRP2.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


p<-FeaturePlot(mammal.combinedA, features = "SFRP2", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_SFRP2.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p<-FeaturePlot(mammal.combinedA, features = "SLC7A3", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_SLC7A3.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p<-FeaturePlot(mammal.combinedA, features = "DPPA4", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_DPPA4.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


p<-FeaturePlot(mammal.combinedA, features = "EPCAM", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_EPCAM.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


p<-FeaturePlot(mammal.combinedA, features = "SEPHS1", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_SEPHS1.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p<-FeaturePlot(mammal.combinedA, features = "RARRES2", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_RARRES2.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p<-FeaturePlot(mammal.combinedA, features = "USP44", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_USP44.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


p<-FeaturePlot(mammal.combinedA, features = "POLR3G", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_POLR3G.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


p<-FeaturePlot(mammal.combinedA, features = "CLDN6", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_CLDN6.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


p<-FeaturePlot(mammal.combinedA, features = "CLU", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_CLU.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)



p<-FeaturePlot(mammal.combinedA, features = "PDGFA", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_PDGFA.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


p<-FeaturePlot(mammal.combinedA, features = "SOX15", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_SOX15.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


p<-FeaturePlot(mammal.combinedA, features = "SOX2", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_SOX2.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)


p<-FeaturePlot(mammal.combinedA, features = "NANOG", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_NANOG.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p<-FeaturePlot(mammal.combinedA, features = "TFAP2A", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_TFAP2A.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p<-FeaturePlot(mammal.combinedA, features = "TFAP2C", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_TFAP2C.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p<-FeaturePlot(mammal.combinedA, features = "KRT7", reduction = "pca", split.by = "Dataset", pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_KRT7.pdf",sep=""),width = 10, height = 8,p,limitsize = FALSE)

p<-FeatureScatter(mammal.combinedA,feature1 = "POU5F1",feature2 = "WNT6",pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_OCT4_WNT6.pdf",sep=""),width = 8, height = 8,p,limitsize = FALSE)


p<-FeatureScatter(mammal.combinedA,feature1 = "POU5F1",feature2 = "WNT6",pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_OCT4_WNT6.pdf",sep=""),width = 8, height = 8,p,limitsize = FALSE)


p<-FeatureScatter(mammal.combinedA,feature1 = "POU5F1",feature2 = "TFAP2A",pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_OCT4_TFAP2A.pdf",sep=""),width = 8, height = 8,p,limitsize = FALSE)


p<-FeatureScatter(mammal.combinedA,feature1 = "POU5F1",feature2 = "TFAP2C",pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_OCT4_TFAP2C.pdf",sep=""),width = 8, height = 8,p,limitsize = FALSE)


p<-FeatureScatter(mammal.combinedA,feature1 = "POU5F1",feature2 = "GATA2",pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_OCT4_GATA2.pdf",sep=""),width = 8, height = 8,p,limitsize = FALSE)

p<-FeatureScatter(mammal.combinedA,feature1 = "POU5F1",feature2 = "GATA3",pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_OCT4_GATA3.pdf",sep=""),width = 8, height = 8,p,limitsize = FALSE)


p<-FeatureScatter(mammal.combinedA,feature1 = "DPPA4",feature2 = "WNT6",pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_DPPA4_WNT6.pdf",sep=""),width = 8, height = 8,p,limitsize = FALSE)

p<-FeatureScatter(mammal.combinedA,feature1 = "PMAIP1",feature2 = "WNT6",pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_PMAIP1_WNT6.pdf",sep=""),width = 8, height = 8,p,limitsize = FALSE)


p<-FeatureScatter(mammal.combinedA,feature1 = "PMAIP1",feature2 = "TFAP2C",pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_PMAIP1_TFAP2C.pdf",sep=""),width = 8, height = 8,p,limitsize = FALSE)

p<-FeatureScatter(mammal.combinedA,feature1 = "PMAIP1",feature2 = "TFAP2A",pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_PMAIP1_TFAP2A.pdf",sep=""),width = 8, height = 8,p,limitsize = FALSE)



p<-FeatureScatter(mammal.combinedA,feature1 = "SOX17",feature2 = "TFAP2C",pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_SOX17_TFAP2C.pdf",sep=""),width = 8, height = 8,p,limitsize = FALSE)

p<-FeatureScatter(mammal.combinedA,feature1 = "SOX17",feature2 = "TFAP2A",pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_SOX17_TFAP2A.pdf",sep=""),width = 8, height = 8,p,limitsize = FALSE)


pgcsubs1 <- WhichCells(mammal.combinedA,expression = SOX17>0.5)
pgcsubs2 <- WhichCells(mammal.combinedA,expression = TFAP2C>0.5)
pgcsubs3 <- WhichCells(mammal.combinedA,expression = TFAP2A>0.5)
pgcsubs4 <- WhichCells(mammal.combinedA,expression = PRDM1>0.5)

pgcsubs5 <- intersect(pgcsubs1,pgcsubs2)
pgcsubs6 <- intersect(pgcsubs1,pgcsubs3)
pgcsubs7 <- intersect(pgcsubs1,pgcsubs4)

cells2 <- setdiff( colnames(mammal.combinedA) , unique(c(pgcsubs5,pgcsubs6, pgcsubs7)) )
cells2 <- colnames(mammal.combinedA) 
cells = unique(c(pgcsubs5,pgcsubs6, pgcsubs7))

mammal.combinedC <- subset(mammal.combinedA, cells = cells2 ) 

mammal.combinedD <- mammal.combinedA
Idents(mammal.combinedD,cells=WhichCells(mammal.combinedD,idents=c("EmDisc/Am1","EmDisc/Am2","EmDisc/Am3", "EmDisc/Am4", "EmDisc/Am5", "PGC"))) <- "EmDisc"
Idents(mammal.combinedD,cells= cells) <- "PGC"
Cl1 <- FindAllMarkers(mammal.combinedD,only.pos = TRUE, test.use = "MAST")

write.table(Cl1,file=paste(saveext,'AllPGCmarkers.csv',sep=""),sep=",")

mammal.combinedB <- subset(mammal.combinedA, cells = unique(c(pgcsubs5,pgcsubs6, pgcsubs7)) )
MiniList1 <- c("DPPA3","POU5F1","PDGFA","SFRP1","SOX15","GABRP","VTCN1","WNT6","BAMBI","TFAP2C","TFAP2A","GATA3","SOX17","GATA4","NANOS3","CER1","LEFTY1","LEFTY2","APOA1","APOB","PDGFRA","HAND2","BST2","GATA6","HGF","PRDM1","SNAI1","VIM","FN1")
MiniList1 <- c("POU5F1","PDGFA","VTCN1","WNT6","BAMBI","TFAP2C","TFAP2A","SOX17","GATA4","NANOS3","CER1","LEFTY1","LEFTY2","HAND2","BST2","PRDM1","SNAI1")
MiniList1 <- c("TFAP2C","TFAP2A","SOX17","NANOS3","PRDM1")
MiniList1 <- c("TFAP2C","TFAP2A","SOX17","NANOS3","PRDM1","POU5F1", "ALPL", "KLF4", "LIN28A", "KIT", "FUT4",  "DPPA3", "ZFP42", "CD38", "EPCAM", "ITGA6","ITGB3", "FGFR3", "KIT", "UTF1", "TET1", "TET2", "TET3", "DND1","SOX15", "GATA4", "PRDM14", "SALL4","CDX2", "GATA3","PRMT5")

X <- GetAssayData(mammal.combinedB, assay = "RNA")
Xh <- t(scale(t(X)))
Xp <- X[MiniList1,]
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(mammal.combinedB))), ID = factor(( mammal.combinedB$ID3)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(0, 1.2, length.out = 50)
pm <- pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 3, filename = paste(saveext,"/DimRed/LineageLabelledExtendedListPGC",".pdf",sep="") ,width=10,height=4)

X <- GetAssayData(mammal.combinedB, assay = "RNA")
#Xh <- t(scale(t(X)))
Xp <- exp(X[MiniList1,]-1)*100
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(mammal.combinedB))), ID = factor(( mammal.combinedB$ID3)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(0, 70, length.out = 50)
pm <- pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=TRUE, cluster_cols=TRUE, cutree_cols = 4, filename = paste(saveext,"/DimRed/LineageLabelledExtendedListPGC2",".pdf",sep="") ,width=10,height=4)

MiniList1 <- c("SOX17","TFAP2A",
"PRDM1",
"TFAP2C",
"GATA3",
"TET2",
"TET1",
"TET3",
"POU5F1",
"PRMT5",
"KLF4",
"FUT4",
"SOX15",
"ZFP42",
"DPPA3",
"DND1",
"CDX2",
"PRDM14",
"GATA4",
"UTF1",
"KIT",
"CD38",
"KIT",
"ALPL",
"LIN28A",
"FGFR3",
"NANOS3",
"ITGB3",
"SALL4",
"EPCAM",
"ITGA6")

X <- GetAssayData(mammal.combinedC, assay = "RNA")
Xh <- t(scale(t(X)))
Xp <- Xh[MiniList1,] #exp(X[MiniList1,]-1)*100
Xp <- na.omit(Xp)
annotation_col = data.frame(Stage = factor(droplevels( Idents(mammal.combinedC))), ID = factor(( mammal.combinedC$ID3)) )
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
my_breaks <- seq(-1, 1, length.out = 50)
pm <- pheatmap(Xp,color =  redblue1(50), fontsize = 4, breaks = my_breaks, border_color = NA,annotation_col = annotation_col,  cluster_rows=FALSE, cluster_cols=TRUE, cutree_cols = 12, filename = paste(saveext,"/DimRed/LineageLabelledExtendedListPGC3",".pdf",sep="") ,width=20,height=4)


#PMAIP1 TDGF1 SCG3 FOXD3-AS1 HPAT5 L1TD1 POU5F1 GPRC5B MT1H PTPRZ1 MT1G SFRP2 SLC7A3 DPPA4 EPCAM SEPHS1 USP44 RARRES2 POLR3G CLDN6 CLU

#mammal.combinedA <- FindNeighbors(mammal.combinedA, reduction = "pca", dims = 1:20)#
#mammal.combinedA <- FindClusters(mammal.combinedA,resolution = 0.5)

#p<-DimPlot(mammal.combinedA,  pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
#ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_Cl.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)



mammal.combined <- AddModuleScore(mammal.combinedA, features = list(intersect(PrimedNaiveList2,rownames(mammal.combined))), name = "PrimedvNaive")
p1<-FeaturePlot(mammal.combined, features = "PrimedvNaive", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PrimedvNaive.pdf",sep=""),width = 16, height = 8,p1)



mammal.combined <- AddModuleScore(mammal.combinedA, features = list(intersect(PrimedNaiveList1,rownames(mammal.combined))), name = "NaivevPrimed")
p1<-FeaturePlot(mammal.combined, features = "NaivevPrimed1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/NaivevPrimed.pdf",sep=""),width = 16, height = 8,p1)

mammal.combined <- AddModuleScore(mammal.combinedA, features = list(intersect(PrimedNaiveList2,rownames(mammal.combined))), name = "PrimedvNaive")
p1<-FeaturePlot(mammal.combined, features = "PrimedvNaive1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PrimedvNaive.pdf",sep=""),width = 8, height = 8,p1)

mammal.combined <- AddModuleScore(mammal.combinedA, features = list(intersect(hsEAm1,rownames(mammal.combined))), name = "Amnion")
p1<-FeaturePlot(mammal.combined, features = "Amnion1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Amnion.pdf",sep=""),width = 8, height = 8,p1)

mammal.combined <- AddModuleScore(mammal.combinedA, features = list(intersect(hsEAm2,rownames(mammal.combined))), name = "Amn")
p1<-FeaturePlot(mammal.combined, features = "Amn1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Amnion2.pdf",sep=""),width = 8, height = 8,p1)

mammal.combined <- AddModuleScore(mammal.combined, features = list(intersect(hsEAm3,rownames(mammal.combined))), name = "Am")
p1<-FeaturePlot(mammal.combined, features = "Am1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Amnion3.pdf",sep=""),width = 8, height = 8,p1)


marmoset_dataInVivo2$oldID <- Idents(marmoset_dataInVivo2)
Idents(marmoset_dataInVivo2,cells=WhichCells(marmoset_dataInVivo2,idents=c("SYS_CS5","SYS_CS6","SYS_CS7"))) <- "SYS"
Idents(marmoset_dataInVivo2,cells=WhichCells(marmoset_dataInVivo2,idents=c("VE_CS5","VE_CS6","VE_CS7"))) <- "VE"
Idents(marmoset_dataInVivo2,cells=WhichCells(marmoset_dataInVivo2,idents=c("Am_CS5","Am_CS6","Am_CS7"))) <- "Am"
Idents(marmoset_dataInVivo2,cells=WhichCells(marmoset_dataInVivo2,idents=c("EmDisc_CS5","EmDisc_CS6","EmDisc_CS7"))) <- "EmDisc"
Idents(marmoset_dataInVivo2,cells=WhichCells(marmoset_dataInVivo2,idents=c("PGC_CS5","PGC_CS6","PGC_CS7"))) <- "PGC"
Idents(marmoset_dataInVivo2,cells=WhichCells(marmoset_dataInVivo2,idents=c("ExMes_CS5","ExMes_CS6","ExMes_CS7","Stalk_CS6","Stalk_CS7"))) <- "ExMes"


Mrk2 <- FindAllMarkers(marmoset_dataInVivo2, test.use = "MAST", only.pos = TRUE)
LIST3 <- rownames(Mrk2)[which(Mrk2$cluster=="Am")]  
LIST4 <- rownames(Mrk2)[which(Mrk2$cluster=="VE")]  
LIST5 <- rownames(Mrk2)[which(Mrk2$cluster=="SYS")]  
LIST6 <- rownames(Mrk2)[which(Mrk2$cluster=="EmDisc")]  
LIST7 <- rownames(Mrk2)[which(Mrk2$cluster=="ExMes")]  
LIST8 <- rownames(Mrk2)[which(Mrk2$cluster=="PGC")] 

DefaultAssay(mammal.combined) <- "RNA"

mammal.combined <- AddModuleScore(mammal.combined, features = list(intersect(LIST3,rownames(mammal.combined))), name = "MarmAmnion")
p1<-FeaturePlot(mammal.combined, features = "MarmAmnion1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/DEMarmAmnion.pdf",sep=""),width = 8, height = 8,p1)


mammal.combined <- AddModuleScore(mammal.combined, features = list(intersect(LIST6,rownames(mammal.combined))), name = "MarmEmDisc")
p1<-FeaturePlot(mammal.combined, features = "MarmEmDisc1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/DEMMarmEmDisc.pdf",sep=""),width = 8, height = 8,p1)


mammal.combined <- AddModuleScore(mammal.combined, features = list(intersect(LIST7,rownames(mammal.combined))), name = "MarmExMes")
p1<-FeaturePlot(mammal.combined, features = "MarmExMes1", pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/DEMarmExMes.pdf",sep=""),width = 8, height = 8,p1)


Idents(mammal.combined) <- paste(mammal.combined$ID3,mammal.combined$Genotype,sep="_")
p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Nakamura_pca_GT1.pdf",sep=""),width = 60, height = 8,p,limitsize = FALSE)


mammal.combined <- readRDS(file=paste(saveext,"/BatchAll_marmosetintegration_wMa.rds",sep=""))


Hyp <- c(
"ATGAGTCAGTCAATCC-1_3",
"ACCAAACCAAGCTGCC-1_8",
"CTCACTGGTAATGATG-1_2",
"GTTACCCAGTCTAGCT-1_2",
"ATGCGATTCTATTGTC-1_2",
"GTTTACTTCCTAGCCT-1_2",
"GTTGCTCTCACCTCGT-1_2",
"GTAACACCACACAGAG-1_7",
"ATCACGACAAGAGGCT-1_7",
"GGTTAACAGAGCTGCA-1_7",
"AGGTAGGTCCCTCTTT-1_5",
"CATACTTGTCCGGTGT-1_6",
"ACTTTGTGTGGCTGCT-1_2",
"TATCAGGCACAAGTTC-1_2",
"CCTGCATCATTAAGCC-1_2",
"TCTTCCTTCAACCCGG-1_2",
"TTCCTAACACTATCGA-1_2",
"CTTCTCTCAGGAATAT-1_2",
"CAAGAGGTCCACCCTA-1_2",
"CACGTGGTCAGCTGTA-1_2",
"CCCTCAACACCCTTAC-1_2",
"CTCCATGGTATCTTCT-1_2",
"TCTCCGACAGGGTCTC-1_2",
"GTCACTCTCGGTGTAT-1_2",
"TATCGCCGTAACACCT-1_2",
"ATTCACTCATACTTTC-1_2",
"CACTAAGCAATCTAGC-1_2",
"GTAATCGTCAGTGCGC-1_2",
"CCTCACATCCCGAACG-1_2",
"TCCTCCCGTGAGAACC-1_2",
"AAAGAACAGCCGTAAG-1_2",
"TCCATGCGTGCAATGG-1_2",
"GATTCTTCATGGAAGC-1_2",
"GGATCTATCATCACAG-1_2",
"TCAGGTATCTGTACAG-1_2",
"AACGTCATCTTACCGC-1_2",
"ATGGGAGAGGGACAGG-1_2",
"ATTACCTGTTCCACGG-1_2",
"TCGTGCTCACGTTCGG-1_2",
"CCATCACCAACATACC-1_2",
"TGCGGCATCTCATTGT-1_2",
"GCTTGGGGTCCGAAGA-1_2",
"AAAGTCCGTCCAGCGT-1_2",
"AGGGAGTTCTAATTCC-1_2",
"CCTCAACTCGTGAGAG-1_2",
"AGCGCCAGTGGATCAG-1_2",
"TGGATCAGTACACGCC-1_2")

SYS <- c(
"AGCATCACAGCACACC-1_6",
"GACGTTAAGCAGTAAT-1_2",
"TACGGTATCCTCAGGG-1_2",
"TCTACCGCAGTATACC-1_6",
"GTCGTAAGTACGAAAT-1_7",
"TTCCTCTTCGCACGAC-1_2",
"ACAAGCTTCTCATTAC-1_6",
"CTACCTGTCCCTTTGG-1_2",
"ACCTGAATCACCGGTG-1_6",
"GCAGTTACATTAAGCC-1_7",
"CCACTTGAGTAAACTG-1_2",
"CCTAAGACACTGCGAC-1_7",
"GATTCTTGTATCACGT-1_7",
"CAACAGTTCGACGACC-1_2",
"TCCCACACAGTTGTCA-1_2",
"GTTGAACTCCTCAGAA-1_7",
"GAGCTGCCACTCTGCT-1_2",
"AGAAGTATCCTAGCTC-1_6",
"TAGGTACCATTCAGCA-1_7",
"TTCTCTCCAGCAATTC-1_2",
"CTCCTTTTCTGCTTAT-1_2",
"CCACACTCACTACTTT-1_8",
"AACAACCAGCCTTTCC-1_2",
"TGTCCTGCACTCTAGA-1_7",
"CTGCCATTCTTACGTT-1_2",
"TTGTTGTGTAAGCGGT-1_7",
"TCCACCATCCTAAACG-1_2",
"TTTACCATCATCGCCT-1_2",
"TTACGTTGTTCCCACT-1_6",
"GATGACTGTTGCGAAG-1_6",
"GGGCTACTCGCAGTTA-1_6",
"TCACACCCACTCCGGA-1_2",
"TTCACGCCATTCTCTA-1_2",
"ATTTCACTCCTAAACG-1_2",
"TCCTTTCTCTGTCCGT-1_2",
"CAAGACTGTATTGCCA-1_2",
"GTAGATCAGCTGGCCT-1_2",
"TCAATTCAGGGTTAAT-1_2",
"ACTCCCATCAACACGT-1_2",
"GTGTTCCAGGGCAGGA-1_2",
"ACATTTCAGCTTTCTT-1_2",
"GCATCTCTCTTCCCAG-1_2",
"GCAACCGTCGCTTACC-1_8",
"GGGAAGTCACGTAGTT-1_2",
"TTGTTCAGTCGTTATG-1_6",
"TTTATGCTCTATGCCC-1_6",
"TCTCTGGCAAGCAATA-1_2",
"TTCAATCGTAGATGTA-1_2",
"AGTACTGCATAAGATG-1_7",
"AGATAGAGTTCTCTCG-1_2",
"ACCCTTGGTAGGTTTC-1_6",
"GATTCTTTCCCATAAG-1_2",
"CCGATCTCACCAAAGG-1_6",
"GATCCCTCATTGTAGC-1_6",
"AACACACAGCCTCACG-1_2",
"TCAGCCTTCTTACACT-1_2",
"GGGTTATCAAGGACAC-1_6",
"CCGCAAGTCTGCTAGA-1_7",
"TAGGTTGCAACTCCAA-1_6",
"TTGCGTCTCTCGCCTA-1_6",
"ATCACGAGTGCATACT-1_4",
"CTACGGGTCCTCTAAT-1_6",
"TCGCTCATCTTTGCTA-1_2",
"TTACCGCTCTGTAAGC-1_2",
"ATTCCCGTCCAAGCTA-1_2",
"GCGAGAATCCCGAACG-1_6",
"TAATCTCGTTCTCCTG-1_2",
"CATCCCAAGTCGCGAA-1_2",
"GAGACTTTCCATTCAT-1_2",
"GTCTCACCACTGGATT-1_6",
"CTCATGCTCCCAAGTA-1_2",
"ATCATTCAGGGATGTC-1_2",
"TCCATGCAGCCATTCA-1_2",
"TCAGTCCGTGACTCTA-1_2",
"TCCGTGTGTAAGAACT-1_2",
"GGGTCACAGCTACGTT-1_7",
"GGTTGTAAGGATACGC-1_7",
"TCTTAGTTCCAGCCTT-1_2",
"CCCATTGCACTTGAAC-1_6",
"CCTCTCCAGTCCCTAA-1_2",
"TCAATCTCACACCGCA-1_2",
"TTACGTTTCCATCACC-1_2",
"CCTCACAGTACCTATG-1_7",
"TGATCAGTCTTTGCGC-1_2",
"AAGCCATCACACCAGC-1_6",
"TGTTGGATCACTCCGT-1_2",
"TTCTTCCAGTCATCGT-1_2",
"CAGTTAGTCCACCTCA-1_2",
"AAAGTCCAGGGTGAAA-1_6",
"CACACAACAGGATTCT-1_2",
"ACGCACGTCCTACGAA-1_2",
"AGCGATTCACCATTCC-1_2",
"CTCAGAATCCCTTCCC-1_2",
"GTCACGGGTGTTCCAA-1_7",
"CAGGCCAAGTTTGCTG-1_6",
"GCACGGTTCTTACTGT-1_6",
"TTTCCTCTCATCCTAT-1_2",
"TTTAGTCTCATCCCGT-1_7",
"TTAGTCTCAAATGAGT-1_6",
"GCCTGTTAGCCGTTGC-1_2",
"ATCAGGTTCCCTAGGG-1_6",
"AACACACTCTCTCTTC-1_6",
"TAAGCACCAATGGCAG-1_6",
"AAGATAGAGAACCGCA-1_8",
"CACCAAACACCCATAA-1_8",
"AGTAGCTTCTGGTCAA-1_6",
"CGCCATTAGTACGAGC-1_6",
"GAAGCGACAGATCACT-1_6",
"AAAGGGCGTCACGACC-1_6",
"GTCAGCGGTGTTAGCT-1_6",
"TAACTTCTCGCGATCG-1_6",
"GGGACCTCACGGTGTC-1_6",
"TGTAAGCAGGCTAAAT-1_6",
"CTGGACGTCTACCCAC-1_6",
"TCCCATGAGTAGCATA-1_6",
"TCTGCCAAGTATTAGG-1_6",
"TGTTCCGGTCGAGTGA-1_6",
"CAAGACTCAAGATCCT-1_6",
"TAGGAGGTCTATTCGT-1_6",
"GCCAACGCATGACGAG-1_6",
"GGTCTGGTCCACCCTA-1_6",
"CCGGACATCCAAGCAT-1_6",
"GACACGCGTCTGTAGT-1_6",
"CCGTTCAGTTACGCCG-1_2",
"CTTAGGACATATGCGT-1_2",
"TGCAGGCCACGAAGAC-1_6",
"TGGTACATCGCACTCT-1_6",
"AGGGCTCGTGTCACAT-1_2",
"GGAAGTGAGGCCTTGC-1_6",
"CGTGCTTAGTTCCTGA-1_6",
"CTCGAGGAGCACCTGC-1_6",
"TGTCCACGTAGCGATG-1_6",
"TGAGGTTTCTACGCAA-1_6",
"TGTGCGGCACTAGGTT-1_6",
"CTCTCGAGTACTAGCT-1_6",
"ATACCGAAGGTGCTAG-1_6",
"TATCTGTCAAGAGAGA-1_6",
"GAGGCAACAATATCCG-1_2",
"CAGAGCCCAAGAGTTA-1_6",
"CAATGACAGAGTACCG-1_6",
"CACGTGGTCCCTCATG-1_2",
"GGAGGATTCGCCACTT-1_2",
"GGGACTCAGCCGGATA-1_6",
"GGCTTTCAGCATAGGC-1_2",
"AGACACTCATTCTGTT-1_6",
"CCATCACCAACAAGTA-1_2",
"TCATTCACAACACGTT-1_6",
"TTACGTTTCAAGTTGC-1_6",
"GACTCTCTCCCAGGAC-1_6",
"CGTGAATGTTGTAAAG-1_6",
"TTATTGCAGCACAAAT-1_2",
"TCGACCTCATTCAGCA-1_6",
"GCACTAATCTCCGAAA-1_2",
"TCATTTGTCATCTACT-1_2",
"CCACCATTCCTACCGT-1_2",
"TCTAACTTCTTCCACG-1_2",
"CGTCCATTCGCACGGT-1_6",
"CCGTAGGGTCAAAGAT-1_6",
"CTTTCGGTCACTAGCA-1_6",
"ACTCCCACACTCTGCT-1_2",
"ACTTCCGTCTGTACAG-1_6",
"CACGAATCAATCTGCA-1_6",
"TTTACTGTCACTGTTT-1_6",
"TATTCCAGTTACACTG-1_2",
"TCACTATGTGAGAACC-1_2",
"CTCAACCAGTAGAGTT-1_8",
"GCCTGTTTCTATTTCG-1_2",
"GCCCGAACAAACGGCA-1_8",
"GTAGCTACAGCTTTCC-1_6",
"GAATCGTAGGGTGAGG-1_8",
"CCACTTGCATCGGTTA-1_2",
"TCCTCCCTCGGTCTAA-1_8",
"CGCCATTTCGTCTACC-1_2",
"AATCACGAGCGGTAAC-1_6",
"GAGACTTAGCAAATGT-1_6",
"GTTGCTCAGGTACTGG-1_2",
"CCACACTTCTCTATGT-1_2",
"TTACGTTAGATGGTAT-1_6",
"CTTGATTTCTGGCCGA-1_6",
"TATTGCTAGGACGCTA-1_6",
"CAACCTCAGTCCTGTA-1_2",
"ACAAAGAAGTTGCATC-1_6",
"ATTCATCAGACGATAT-1_6",
"CACATGAAGACCACGA-1_6",
"GATTCGACAAGTGGAC-1_6",
"ATATCCTCAATAGTAG-1_2",
"TAAGCACAGAATCGAT-1_6",
"ACGCACGCACTCCGAG-1_8",
"GGATGTTTCCATCTAT-1_2",
"CTAGACATCTCAAAGC-1_2",
"TATCGCCTCCACTAGA-1_2",
"GGTGATTCACCATATG-1_2",
"ACACGCGAGGAGGGTG-1_6",
"CTGGACGCATGCCGAC-1_6",
"CCTACGTAGTCTAGCT-1_2",
"GTTGTGAAGAGCAGCT-1_6",
"TAGCACACACCTCGTT-1_6",
"TCCCATGTCCCTGGTT-1_6",
"TACAACGCAACCAATC-1_2",
"CTGAGGCCACGGGCTT-1_2",
"CATAAGCTCGTGCACG-1_6",
"ATCGATGGTAGCGTAG-1_2",
"AGGGTCCGTACAATAG-1_2",
"GGATGTTCACCAGCTG-1_6",
"ATCACGAGTCCACAGC-1_2",
"TCATCCGGTAAGCGGT-1_2",
"GTAGGTTTCGTACACA-1_6",
"TCATGTTGTCGCGGTT-1_6",
"AGGCATTCAAGACAAT-1_2",
"CCCTCTCAGCCATCCG-1_2",
"GTCTCACGTTGGCTAT-1_8",
"GGTAACTGTATGCGTT-1_8",
"TACACCCAGTTTGTCG-1_8",
"GTAACCATCATCTGTT-1_5",
"CGCGTGACAGCGATTT-1_6",
"CCTCACACATGCCGCA-1_6",
"CAGCAGCTCCTCCACA-1_2",
"CAGCCAGCAGATTTCG-1_6",
"TGAATCGAGTCGTCTA-1_6",
"AGATGAAGTTCGTAAC-1_5",
"TTACTGTTCGTTCAGA-1_8",
"ATGGAGGCACTTGGCG-1_7",
"GAGGCCTAGGATTCCT-1_6",
"TGTACAGGTCTAGGTT-1_8",
"TCTTCCTCAAGTCATC-1_2",
"AGATCGTGTGGTAATA-1_6",
"ATCGCCTTCTGAGTCA-1_8",
"ATACCTTCAGCTATTG-1_6",
"GAGTTTGAGCATGCAG-1_7",
"TTTACGTCATGATGCT-1_7",
"ACTTAGGTCCATTTAC-1_2",
"AAAGAACAGTAAACGT-1_6",
"TGAGGTTGTTTCCATT-1_6",
"GTACAACGTGAGCCAA-1_6",
"ATCGTCCGTTGTTTGG-1_7",
"AATAGAGGTCACTTAG-1_8",
"ATGAGTCTCCTTATAC-1_8",
"TCGCTTGAGCCGATCC-1_6",
"TGGGCTGCAACCAGAG-1_6",
"CATCGTCTCGTTGTAG-1_6",
"CTGCTCACATGGAATA-1_6",
"TGCGATACAAAGGCAC-1_6",
"TGGAGGAGTGTCTTGA-1_6",
"TCATGCCAGACATCCT-1_2",
"AGGACGAGTGGCGCTT-1_6")

ExMes <- c(
"CATTGAGCATAGTCGT-1_6",
"GACCAATCATTGACAC-1_6",
"GGTTGTAGTGCCGAAA-1_6",
"GACGTTATCCACCTCA-1_7",
"GTAGTACAGCTCGTGC-1_7",
"GGGAGATTCATGAAAG-1_3",
"CAGTTAGAGTTTCGGT-1_3",
"GGGTGAAAGTATGAGT-1_3",
"TGGGAAGAGCTGAGCA-1_3",
"TACTTACAGGTAGTAT-1_3",
"AGTAGCTCAGGACATG-1_3",
"TGAATGCGTTCGGACC-1_3",
"GCCGTGAGTCAAGCCC-1_3",
"ATGCATGTCTCTTAAC-1_3",
"TAACTTCGTACGACAG-1_3",
"TCATCCGGTGTGGACA-1_3",
"CATGCCTGTCACAGAG-1_3",
"TTCAGGATCTGTGTGA-1_3",
"CATTGCCTCCGATAAC-1_3",
"CCAAGCGAGGCTCACC-1_3",
"CACGTGGAGAAGCGCT-1_3",
"CTGAATGAGAGCCCAA-1_3",
"GACCCAGAGGGCCAAT-1_3",
"GGTAATCTCGGCGATC-1_3",
"CCACGAGCAACTTGGT-1_3",
"GTCATTTGTATCGCGC-1_3",
"CGAATTGCAAGACGGT-1_3",
"CAAGACTTCTTCCGTG-1_3",
"ATTATCCGTCCATACA-1_3",
"CCTCAGTGTTCAGCTA-1_3",
"CTCGAGGCAGTTTCAG-1_3",
"GATCCCTAGACATACA-1_3",
"AAGCGTTAGCTTGTTG-1_3",
"GGAATCTTCCTGGGTG-1_3",
"AGGTCATAGATGCCGA-1_3",
"GCCTGTTAGTCCCGGT-1_3",
"GTCTCACAGCGGTAGT-1_3",
"AACACACCAATGAACA-1_3",
"TCATTACCAAGTATAG-1_3")

#SYS or PGC?
pPGC <- c("CGACAGCCATCTTCGC-1_2",
"CCCAACTTCTACTTCA-1_2",
"ATAGACCTCTGCTTAT-1_2",
"AACCTGATCCATCCGT-1_7",
"TTACAGGGTCTTCTAT-1_3",
"TAACTTCGTACAAAGT-1_3",
"CCTCAACTCTTGCAAG-1_3",
"TCATGTTGTCGACGCT-1_3")


AmEmD1 <- c(
"GTCGAATCATTGTCGA-1_3",
"AACTTCTAGGACATCG-1_3",
"ATTCGTTCAGGAACCA-1_3",
"GAGTCATCACGGGCTT-1_3",
"TCTACCGGTATTCCGA-1_3",
"CGGAGAACACTCAAGT-1_3",
"TTACCATTCATGTCTT-1_3",
"CACATGAGTTCTCCTG-1_3",
"CAGCAATCATTGTCGA-1_3",
"CAAGACTTCATCCTAT-1_3",
"TTTGTTGTCTACAGGT-1_3",
"CATAGACGTACGGATG-1_3",
"CAATCGATCCCTTGGT-1_3",
"GCTGCAGTCCGTGCGA-1_3",
"CTCCACAAGTGGATAT-1_3",
"GCTGGGTAGGAGACCT-1_6",
"CAACAACGTATTAAGG-1_3",
"GATCAGTCACTAAACC-1_3")

AmEmD2 <- c(
"CTCTCGAAGCCTCAGC-1_4",
"CAGTGCGAGTAGCCAG-1_3",
"AGGCATTTCCCGGTAG-1_3",
"AGCGTATCAAGGCGTA-1_3",
"CCTGCATCAAGTATCC-1_3",
"TGCTTCGCATCGAACT-1_3",
"AGTGCCGAGCTGACAG-1_3",
"CATGCCTCACGTACTA-1_3",
"GTCAAACTCTCAGTCC-1_3",
"AGACCATCAGTGGTGA-1_6",
"ATGCCTCAGTCATCGT-1_6",
"ACACCAACAACCGCTG-1_3",
"AGCGCCACACGGCACT-1_3",
"TTCCTCTGTTCCGTTC-1_3",
"AACGTCACAGGTCAAG-1_3",
"GATTGGTTCGCTAGCG-1_6",
"TATACCTCACCAATTG-1_3",
"CTGAGCGCAGCTTTCC-1_6",
"TTCATTGCACGAAAGC-1_6",
"ACAGAAAAGACTCGAG-1_3",
"GAACTGTGTATGGAGC-1_3",
"CCACCATGTTCAGGTT-1_3",
"GGTGTCGGTACCCGAC-1_3",
"ATCCCTGGTACTCGTA-1_3",
"CCCATTGCATCTGCGG-1_3",
"GGATGTTAGATCGCTT-1_3",
"ATACCGAAGGCCCAAA-1_3",
"CTATAGGGTGCCTACG-1_3",
"CTGTGGGTCAGCATTG-1_3",
"GGTCACGTCTGAGTCA-1_3",
"CTGCTCATCAGTGATC-1_3",
"GTCTACCAGAGCCTGA-1_3",
"ACCTACCAGTTAGAAC-1_6",
"CTATCTAAGCGGACAT-1_3",
"TCTGCCACATGCCGAC-1_6",
"ATCGGCGGTCCTGAAT-1_3",
"TCATCATGTACCCGCA-1_3",
"TCAGTTTTCTCCAATT-1_3",
"TTACCATTCTCCAATT-1_2",
"CAGTTAGCATGAGAAT-1_2",
"CAACCTCGTCACTTCC-1_2",
"TTTCAGTGTCAACATC-1_2",
"TAGTGCACAACTTGCA-1_6",
"TCGTAGACACTGATTG-1_6",
"CCTCTAGTCATTGCTT-1_5",
"TTTCATGTCCGCAACG-1_6",
"TGGATGTGTTAACAGA-1_6",
"ATTCCCGTCCGTCAAA-1_2",
"TACCCGTCAAACCATC-1_2",
"CGGAATTGTGGGTCAA-1_6",
"AGTAGCTCACTACCGG-1_5",
"TTTGACTAGGACATCG-1_6",
"GGGTTATCAAATAAGC-1_3",
"CAAAGAAGTCATCAGT-1_3",
"AATCGTGAGCTCGACC-1_3",
"CTCAGAAAGTGATAGT-1_3",
"TCGGGCAGTAGCTCGC-1_5",
"GTGGTTATCGCCAATA-1_6",
"CTGCCTACAACAAGTA-1_5",
"GACGCTGGTTCATCGA-1_6",
"ACGTAGTGTGACGTCC-1_3",
"TGTTCTACAGCGTGAA-1_3",
"TGCGACGCACATTGTG-1_3",
"CTAACCCAGTAGAGTT-1_3",
"GACCCTTAGGACTTCT-1_3",
"AGACCCGCATCGCTAA-1_3",
"CAAGACTAGGACGCTA-1_3",
"TCTTAGTGTGGACTAG-1_3",
"GGGACAACAATCGCGC-1_3",
"GTTGCGGTCATGCGGC-1_3",
"TAACACGTCACGTAGT-1_3",
"GGATCTACAAGTATAG-1_3",
"AATCGACAGACGAAGA-1_3",
"TCATGGAAGTCTGGAG-1_3",
"CCGGACAGTATGAGGC-1_3",
"TGTTACTCATAGTCAC-1_3",
"GAGGGTAAGGACTATA-1_3",
"GGTGAAGAGTGGTGGT-1_3",
"GTCTCACAGCGGTATG-1_3",
"AGCGCTGTCGCGGACT-1_5",
"CTACGGGAGCCGTCGT-1_3",
"GTAACCAGTCAAGGCA-1_3",
"GCTGAATTCTCATAGG-1_3",
"GATAGCTTCAAGCCTA-1_3",
"TGCACGGAGCGAAACC-1_3",
"GAGTTTGAGACTCCGC-1_3",
"GTCTGTCAGCGACATG-1_3",
"CCCATTGGTTGCATAC-1_3",
"CTGAATGTCGCTTGAA-1_3",
"TTCCAATCACATACGT-1_3",
"GGGTAGAAGTAGGGTC-1_3",
"CCTCTAGCACCGCTAG-1_3",
"TCATGTTTCCTTCGAC-1_3",
"TTGGGCGTCATCTACT-1_3",
"CTGATCCGTTGCGTAT-1_3",
"TCATGTTTCGAGTCCG-1_3",
"GGTGATTGTGACTCGC-1_3",
"GGTTAACTCCGGCAGT-1_3",
"AATGACCGTTCGGACC-1_3",
"CGGTCAGGTTGTAAAG-1_3",
"GAGTTACTCTCCGCAT-1_3",
"GTGCACGCACAGCATT-1_3",
"TGGATGTAGGGTCAAC-1_3",
"ATATCCTTCTGAGATC-1_3",
"CTTTCAATCATAGCAC-1_3",
"GACCAATGTGTCATGT-1_3",
"GACCTTCCACTTGGCG-1_3",
"TGTTACTAGCCTCTCT-1_3",
"CCTAAGACAGGTCCGT-1_3",
"TCAAGTGGTCACGTGC-1_3",
"CTTTCAACAACAGTGG-1_3",
"ACGGTTAAGGACTATA-1_3",
"TCCTCCCCAAACGAGC-1_3",
"ATTCACTGTAATTGGA-1_3",
"GCGGAAAAGTTGCGCC-1_3",
"GATGAGGAGTTCTACG-1_3",
"GCTGAATCAATGTCTG-1_3",
"TTCATGTGTTGAATCC-1_3",
"GTTTACTAGAGGCGTT-1_3",
"AGACCCGCAGAGGAAA-1_3",
"GAAGAATAGGAGAGGC-1_3",
"CCTTGTGAGCCTGGAA-1_3",
"GTGAGCCGTCCTGGTG-1_3",
"TCTCCGACATTGCAAC-1_3",
"CACAACACAGGACATG-1_3",
"TTCACCGCACACCTAA-1_3",
"GATCCCTCACTACGGC-1_3",
"TCATCCGGTGCTATTG-1_3",
"CGTAATGGTAAGATAC-1_3",
"TATACCTTCCGTGGGT-1_3",
"AAGGAATTCACTTGTT-1_3",
"AATTCCTGTTACGCCG-1_3",
"ATTTCACTCTGCACCT-1_3",
"AACAGGGGTCCAGCCA-1_3",
"CGGGCATCAACCAGAG-1_3",
"TGCATCCGTAGTGGCA-1_3",
"CAAGAGGAGTTGGAGC-1_3",
"GTCGAATCAGCTACAT-1_3",
"TCACATTAGGAAGTAG-1_3",
"CAGCAATAGTTCTACG-1_3",
"CCGGGTACAATCGCAT-1_3",
"ACACTGAGTCTGGTTA-1_3",
"CGGGTCATCAGTCACA-1_3",
"TATATCCGTGGGTCAA-1_3",
"ACCTGTCAGTGAGTGC-1_3",
"GGTGAAGAGGACTTCT-1_3",
"GCCAGCATCACTACGA-1_3",
"ACTATTCAGCACCGTC-1_3",
"ACTGATGAGTCATTGC-1_3",
"CTATCCGTCCTGTTGC-1_3",
"TGAGACTAGACCACGA-1_3",
"CTTCAATAGGATACGC-1_3",
"TCAGCCTCATCCGGTG-1_3",
"AAACCCAAGTAGCAAT-1_3",
"CCTCAACTCGTCACCT-1_3",
"CATCCGTTCTTGGCTC-1_3",
"CTCCAACTCACCCTCA-1_3",
"TCAAGACTCTCCAATT-1_3",
"TAAGTCGGTCACCTTC-1_3",
"GCGTTTCAGTATAACG-1_3",
"TTAGGGTGTTTCCATT-1_3",
"CCCTAACCACCGTGAC-1_3",
"TTCCAATGTGCGGTAA-1_3",
"AACCAACCACGTCATA-1_3",
"TTCTCTCAGACTACGG-1_3")

AmEmD3 <- c(
"TCGCACTAGCCAAGGT-1_2",
"AGATAGAAGTCGAAGC-1_6",
"CACACAACAGTCTGGC-1_6",
"ACTTTCAAGTACCCTA-1_2",
"ATCCTATGTAGGCAAC-1_2",
"AGACAGGCAGTGGTGA-1_6",
"GGAGAACGTGTACGCC-1_6",
"TGGGCTGGTCATCACA-1_6",
"GGTGTCGCAGCCGTTG-1_6",
"TCACTATTCGCTCTAC-1_2",
"TTCACCGGTTCCTAGA-1_2",
"AGGGAGTTCAACCCGG-1_6",
"CTGTATTCAAGGCGTA-1_6",
"ACGGTCGGTAGCGTAG-1_3",
"TTGAGTGGTACGAAAT-1_3")

AmEmD4 <- c(
"CTCCACAAGAGGCTGT-1_3",
"GCACGGTAGGGACACT-1_3",
"GTATTTCTCGTAGAGG-1_3",
"AGTCTCCGTAGCTTGT-1_3",
"ATCCACCGTGAATGTA-1_3",
"GCTTCACCATGTCAGT-1_3",
"GTAGTACTCATCGCCT-1_3",
"GCGGAAAGTTGTTGAC-1_3",
"TAAGTCGGTCCATAGT-1_3",
"CGAGTTAGTGTAGTGG-1_3",
"AACCAACTCATTGAGC-1_3",
"TCTTCCTAGGGAGGCA-1_3",
"AGTTCGAAGTGCAACG-1_3",
"CTGTGAATCAGACATC-1_3",
"ATCGGCGTCTGCACCT-1_3",
"CATACTTAGGGCTAAC-1_3",
"TAGTGCAAGATACCAA-1_3",
"CTCAGGGTCAACGTGT-1_3",
"AGGTGTTTCCTTCTTC-1_3",
"CGTTCTGTCCATTTGT-1_3",
"GATCACATCTAGCATG-1_3",
"TATCGCCTCTCTGGTC-1_3",
"AGCCAGCAGTGACACG-1_3",
"CGGAACCCAACGACAG-1_3",
"GCCATTCCATGATGCT-1_3",
"TTGTTTGCAGACTCTA-1_3",
"CCTTCAGCACATCCCT-1_3",
"TTGCATTAGGGCATGT-1_3",
"CTACCCACAGCGTTTA-1_3",
"TGCTCCACACGCTGCA-1_3",
"GTGAGTTGTACATACC-1_3",
"CACATGACAGCAGTTT-1_3",
"TACGGGCAGGTAATCA-1_3",
"TCCTGCAGTGGCTTGC-1_3",
"ATGACCACAACCCGCA-1_3",
"TGTCAGACATACTGAC-1_3",
"TCGACCTTCTCTCGAC-1_3",
"AGGGCTCAGACGTCGA-1_3",
"GGTGTTAGTTTCGATG-1_3",
"AAGTACCAGATTGATG-1_3",
"TCTACCGTCGATACAC-1_3",
"AATTCCTTCTAGTTCT-1_3",
"CTCCCTCCAGACCATT-1_3",
"CACGGGTCATCCTGTC-1_3",
"GGCACGTTCATACGGT-1_3",
"GAGACCCGTCGTACTA-1_3",
"GTGAGTTAGCCTTGAT-1_3",
"AGGTCATCACTCAGAT-1_3",
"GGGTATTTCCACGAAT-1_3",
"ACGCACGGTGATATAG-1_3",
"ATTATCCCACCCTAGG-1_2",
"CACCGTTCACAATGTC-1_3",
"AATGCCACACCAGTAT-1_3",
"CGCCATTCACATTCTT-1_3",
"TCACAAGAGGACGCAT-1_3",
"TGGTGATTCACGGACC-1_3",
"AAACGAAAGCAACTCT-1_3",
"CTCTCAGAGTCCGCCA-1_3",
"GAGTTTGAGGAGTACC-1_3",
"CATTCTACATCGGATT-1_3",
"ACGTACAAGCGCATCC-1_3",
"CGAGTGCAGCGCATCC-1_3",
"TCCGATCCACGCGCAT-1_3",
"TGAGCATTCATGAGTC-1_3",
"CAAGAGGAGTAGGCCA-1_3",
"CCTCAGTAGGGCGAGA-1_3",
"GGAACCCTCCAATCTT-1_3",
"ACTGCAAAGACCAACG-1_3",
"GGATGTTTCATACGGT-1_3",
"AAGATAGTCTGGTTGA-1_3",
"GCACGTGAGGTGGCTA-1_3",
"TCATGTTTCGTGCACG-1_3",
"ACGATCACACCGCTGA-1_3",
"CACCGTTTCGGTCACG-1_3",
"GGGTAGACATCATTTC-1_3",
"GGGCCATCATGAGGGT-1_3",
"GTTCATTGTTCCGGTG-1_3",
"ATCGGCGGTAGGCAAC-1_3",
"CTGCCATCAAAGCGTG-1_3",
"ATCGCCTTCCGGTTCT-1_3",
"TCCGGGAAGGAATGTT-1_3",
"AAGGTAACAGTTAGGG-1_3",
"ACTCTCGTCCACTAGA-1_3",
"GTCAGCGGTCAAGTTC-1_3",
"TACAGGTTCACGGACC-1_3",
"GTATTGGTCTCCTGAC-1_3",
"TGGGAAGAGAATTGCA-1_3",
"GACTATGAGCTAGATA-1_3",
"ACTTCGCTCGCGAAGA-1_3",
"TATTGCTCACAGAGAC-1_3",
"CAGCAATGTGTAAACA-1_3",
"GGTTAACCAGATCATC-1_3",
"ACAAAGACAGGAGGTT-1_3",
"TTGCGTCCAGTTAAAG-1_3",
"CAGATACCATGACAGG-1_3",
"CTTCGGTTCCGAAATC-1_3",
"ATTCGTTTCCTTATGT-1_3",
"GGGATCCTCTGGACTA-1_3",
"TGATGGTAGTATAACG-1_3",
"CTTCTAACAAGATCCT-1_3",
"CAGGTATGTTCACGAT-1_3",
"TCTTCCTAGTCTAACC-1_3",
"AACAAGATCTCCAAGA-1_6",
"GCGGAAATCATCCTGC-1_6",
"CTAAGTGTCCGTAATG-1_5",
"GTAATGCTCGGTGTAT-1_5",
"TGAGCATGTGGCTTGC-1_3",
"CACACAAAGGAGGGTG-1_3",
"CAGCCAGCAAGCGAAC-1_5",
"TATTTCGCACCTGAAT-1_5",
"CTGTATTGTTAAGCAA-1_3",
"CAACAACGTAATGCGG-1_6",
"CTCATTACAAACACGG-1_5",
"TCTACATTCGAGAAAT-1_6",
"ATCCGTCGTACAAGTA-1_6",
"CTGCTCACACTCACTC-1_2",
"TTGGATGAGTTTGTCG-1_6",
"GTCAAACGTTAGGAGC-1_6",
"GTCTTTAGTCTTGCTC-1_6",
"CAACGGCTCAAGAAAC-1_6",
"GTGCTGGCACCGGTCA-1_6",
"CTCTCAGTCGCCAATA-1_3",
"TCTTTGAAGCAATAGT-1_6",
"TCGCTCATCCACGGAC-1_6",
"TGCAGTACACACACGC-1_6",
"AGTCAACTCACCATAG-1_3",
"TTTCATGTCCACACCT-1_3")

EmDisc1 <- c(
"TCCTTTCAGACAAGCC-1_2",
"TTCCAATCATCTGTTT-1_8",
"GTATTGGAGGTTGGAC-1_2",
"GCATCTCTCGGAATGG-1_6",
"TGATGGTCACTTGACA-1_6",
"CGAGTGCTCGATGCTA-1_3",
"GGACGTCAGGGCATGT-1_3",
"AGGACGAGTTAACAGA-1_3",
"TACGGGCCATGGGTCC-1_3",
"ATTCCCGGTTCCACAA-1_3",
"TTCCTCTAGCAATAAC-1_3",
"TACTTACCACCCTGAG-1_3",
"CCCTTAGCACGACTAT-1_2",
"GATGATCAGAAGAGCA-1_3",
"GTCTCACGTATACAGA-1_3",
"GTGCTTCTCGGTAAGG-1_2",
"CTGGTCTAGTATGTAG-1_6",
"GACCCAGTCTCGAGTA-1_3",
"GGAGGTACATAATCGC-1_3",
"CTTACCGCATGACGTT-1_3",
"GTACAGTGTAACGGTG-1_3",
"ACTATGGTCTCTCCGA-1_3",
"GATAGCTCACTTGAAC-1_3")

EmDisc2 <- c(
"CTTTCAACAACAGTGG-1_3",
"ACGGTTAAGGACTATA-1_3",
"TCCTCCCCAAACGAGC-1_3",
"ATTCACTGTAATTGGA-1_3",
"GCGGAAAAGTTGCGCC-1_3",
"GATGAGGAGTTCTACG-1_3",
"GCTGAATCAATGTCTG-1_3",
"TTCATGTGTTGAATCC-1_3",
"GTTTACTAGAGGCGTT-1_3",
"AGACCCGCAGAGGAAA-1_3",
"GAAGAATAGGAGAGGC-1_3",
"TCATCATGTACCCGCA-1_3",
"TTTGACTAGGACATCG-1_6",
"TGCGACGCACATTGTG-1_3",
"GGATCTACAAGTATAG-1_3",
"TCATGGAAGTCTGGAG-1_3",
"CTATCCGTCCTGTTGC-1_3",
"AAACCCAAGTAGCAAT-1_3",
"CATCCGTTCTTGGCTC-1_3",
"CTCCAACTCACCCTCA-1_3",
"GCGTTTCAGTATAACG-1_3",
"TTAGGGTGTTTCCATT-1_3")

Idents(mammal.combined,cells = Hyp) <- "Hyp_d14"
Idents(mammal.combined,cells = ExMes) <- "ExMes"
Idents(mammal.combined,cells = SYS) <- "SYS"

Idents(mammal.combined,cells = EmDisc1) <- "EmDisc1"
Idents(mammal.combined,cells = EmDisc2) <- "EmDisc2"
Idents(mammal.combined,cells = AmEmD3) <- "AmEmD3"
Idents(mammal.combined,cells = AmEmD2) <- "AmEmD2"

Idents(mammal.combined,cells = AmEmD1) <- "AmEmD1"

Idents(mammal.combined,cells = pPGC) <- "putPGC1"
Idents(mammal.combined,cells = PGC) <- "putPGC2"


p<-DimPlot(mammal.combined,  pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_umap_heatmapano.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)

p<-DimPlot(mammal.combined,   pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca_heatmapano.pdf",sep=""),width = 30, height = 8,p,limitsize = FALSE)


p<-DimPlot(mammal.combined,   pt.size = 4, reduction = "pca", split.by = "ID3", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/Marmoset_BatchAll_Ma_pca2_heatmapano.pdf",sep=""),width = 50, height = 8,p,limitsize = FALSE)


