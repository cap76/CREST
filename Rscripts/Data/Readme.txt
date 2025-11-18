A processed Seurat object containing the gene expression and dimensionality reduction, necessary to replot Figure 6A/B can be found at https://drive.google.com/file/d/1MYrjAxGue18yibqu1eg91ePOFfaYeA6U/view?usp=share_link

A processeed series of processed Seurat objects can be downloaded here: https://drive.google.com/drive/folders/1uV__ZgSRPSG03T-9CCFQggu6I-RVPOcS?usp=share_link. This includes an aligned dataset containing all cells that passed QC (CREST_anotated.rds, https://drive.google.com/file/d/1MYrjAxGue18yibqu1eg91ePOFfaYeA6U/view?usp=share_link). This file holds 15673 cells with the following meta-information:
1) Idents = inferred anotation
2) ID3 = processed batch

It has the following dimensitonality reductions:
1) oldUMAP is an aligned dimensionality reduction to a composite reference (human cycling endometrium and human in vitro cultured embryos. Used for Figure 3. 
2) harmony and harmumap are batch harmonised alignment of the CREST data only using harmony
3) pca, umap, are batch alignment of the CREST data only using Seurats integrate functions.

This file can be donwloaded manually and placed in the Data folder or downloaded using python gdown package:

gdown --id '1MYrjAxGue18yibqu1eg91ePOFfaYeA6U'

Processed data for trophoblast cells aligned to trophoblast cells in cells to human, rhesus, and marmoset for replotting Figure 6D/E can be downloaded from (https://drive.google.com/file/d/1PHpCaxXt54X7zorPuzuf3WyBny5VSaU5/view?usp=share_link) or via gdown

gdown --id '1PHpCaxXt54X7zorPuzuf3WyBny5VSaU5'

Metadata in this file includes Idents=inferred anotation, Species=dataset/species.

The file CREST_epithelia.rds (~1.5GB; https://drive.google.com/file/d/1VphkBKRnWA_2MBEpvqioVGqd9Gr1efvN/view?usp=share_link) is a Seurat object of the CREST epithelial cells aligned to a reference human endometrial dataset. With the following metadata (EpithelialAno = cell anotation, Idents = dataset/stage). Metadata (EpithelialAno = anotations, ID6 = datasets/stage)

The file CREST_stroma.rds (https://drive.google.com/file/d/1BiAa8r-tcEeiEEe0FRwf3bqLOWGYlY7l/view?usp=share_link) contains the aligned stromal cells from CREST to an endometrial reference (Idents = anotation, ID3 = batch).  


Finally, the CREST_aligned_to_human_rhesus_marmoset.rds (https://drive.google.com/file/d/1ZopPHsXsp6YcRO2k5b--V9ZJ3J8ezJxi/view?usp=share_link) contains an alignment of CREST embryonic lineages to human, rhesus, and marmoset reference datasets. 
