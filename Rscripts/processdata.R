library(Seurat)

Alt_0 <- Read10X(data.dir = "SRR11869229_aligned/outs/filtered_gene_bc_matrices/GRCh38/")
Alt1 <- CreateSeuratObject(counts = Alt_0, assay = "RNA",min.cells = 0, min.features = 0)
Alt1 <- NormalizeData(Alt1, verbose = FALSE)
Alt1$ID1 <- "d20"
Alt1$Dataset <- "Early-mid-sec"
saveRDS(Alt1,file="Endo10_.rds")


Alt_0 <- Read10X(data.dir = "SRR11869227_aligned/outs/filtered_gene_bc_matrices/GRCh38/")
Alt1 <- CreateSeuratObject(counts = Alt_0, assay = "RNA",min.cells = 0, min.features = 0)
Alt1 <- NormalizeData(Alt1, verbose = FALSE)
Alt1$ID1 <- "d22"
Alt1$Dataset <- "Mid-sec"

saveRDS(Alt1,file="Endo1.rds")
Alt_0 <- Read10X(data.dir = "SRR11869223_aligned/outs/filtered_gene_bc_matrices/GRCh38/")
Alt1 <- CreateSeuratObject(counts = Alt_0, assay = "RNA",min.cells = 0, min.features = 0)
Alt1 <- NormalizeData(Alt1, verbose = FALSE)
Alt1$ID1 <- "d17"
Alt1$Dataset <- "Early"
saveRDS(Alt1,file="Endo2_.rds")

Alt_0 <- Read10X(data.dir = "SRR11869225_aligned/outs/filtered_gene_bc_matrices/GRCh38/")
Alt1 <- CreateSeuratObject(counts = Alt_0, assay = "RNA",min.cells = 0, min.features = 0)
Alt1 <- NormalizeData(Alt1, verbose = FALSE)
Alt1$ID1 <- "d22"
Alt1$Dataset <- "Prolif"
saveRDS(Alt1,file="Endo3_.rds")

Alt_0 <- Read10X(data.dir = "SRR11869238_aligned/outs/filtered_gene_bc_matrices/GRCh38/")
Alt1 <- CreateSeuratObject(counts = Alt_0, assay = "RNA",min.cells = 0, min.features = 0)
Alt1 <- NormalizeData(Alt1, verbose = FALSE)
Alt1$ID1 <- "d26"
Alt1$Dataset <- "Early-sec"
saveRDS(Alt1,file="Endo4_.rds")

Alt_0 <- Read10X(data.dir = "SRR11869236_aligned/outs/filtered_gene_bc_matrices/GRCh38/")
Alt1 <- CreateSeuratObject(counts = Alt_0, assay = "RNA",min.cells = 0, min.features = 0)
Alt1 <- NormalizeData(Alt1, verbose = FALSE)
Alt1$ID1 <- "d20"
Alt1$Dataset <- "Early-sec"
saveRDS(Alt1,file="Endo5_.rds")

Alt_0 <- Read10X(data.dir = "SRR11869235_aligned/outs/filtered_gene_bc_matrices/GRCh38/")
Alt1 <- CreateSeuratObject(counts = Alt_0, assay = "RNA",min.cells = 0, min.features = 0)
Alt1 <- NormalizeData(Alt1, verbose = FALSE)
Alt1$ID1 <- "d26"
Alt1$Dataset <- "Late-sec"
saveRDS(Alt1,file="Endo6_.rds")

Alt_0 <- Read10X(data.dir = "SRR11869233_aligned/outs/filtered_gene_bc_matrices/GRCh38/")
Alt1 <- CreateSeuratObject(counts = Alt_0, assay = "RNA",min.cells = 0, min.features = 0)
Alt1 <- NormalizeData(Alt1, verbose = FALSE)
Alt1$ID1 <- "d19"
Alt1$Dataset <- "Early"
saveRDS(Alt1,file="Endo7_.rds")

Alt_0 <- Read10X(data.dir = "SRR11869231_aligned/outs/filtered_gene_bc_matrices/GRCh38/")
Alt1 <- CreateSeuratObject(counts = Alt_0, assay = "RNA",min.cells = 0, min.features = 0)
Alt1 <- NormalizeData(Alt1, verbose = FALSE)
Alt1$ID1 <- "d23"
Alt1$Dataset <- "Early"
saveRDS(Alt1,file="Endo8_.rds")

Alt_0 <- Read10X(data.dir = "SRR11869240_aligned/outs/filtered_gene_bc_matrices/GRCh38/")
Alt1 <- CreateSeuratObject(counts = Alt_0, assay = "RNA",min.cells = 0, min.features = 0)
Alt1 <- NormalizeData(Alt1, verbose = FALSE)
Alt1$ID1 <- "d16"
Alt1$Dataset <- "Early"
saveRDS(Alt1,file="Endo9_.rds")

Alt_0 <- Read10X(data.dir = "SRR11869229_aligned/outs/filtered_gene_bc_matrices/GRCh38/")
Alt1 <- CreateSeuratObject(counts = Alt_0, assay = "RNA",min.cells = 0, min.features = 0)
Alt1 <- NormalizeData(Alt1, verbose = FALSE)
Alt1$ID1 <- "d20"
Alt1$Dataset <- "Early-mid-sec"
saveRDS(Alt1,file="Endo10_.rds")

