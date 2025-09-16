# CREST
Cell-engineered receptive endometrial scaffold technology (CREST) by Molè, Elderkin, Zorzan, Penfold, Horsley et al.

Here we provide R scripts and intermediate processed data for the paper Modelling human embryo implantation in vitro (Rscripts folder) and MATLAB scripts for projecting processed data onto CS6 marmoset embryos (MATLAB folder).

A processeed series of processed Seurat objects can be downloaded here: https://drive.google.com/drive/folders/1uV__ZgSRPSG03T-9CCFQggu6I-RVPOcS?usp=share_link

This includes an aligned dataset containing all cells that passed QC (CREST_anotated.rds, https://drive.google.com/file/d/1MYrjAxGue18yibqu1eg91ePOFfaYeA6U/view?usp=share_link). This file holds 15673 cells with the following meta-information:
1) Idents = inferred anotation
2) ID3 = processed batch

It has the following dimensitonality reductions:
1) oldUMAP is an aligned dimensionality reduction to a composite reference (human cycling endometrium and human in vitro cultured embryos. Used for Figure 3. 
2) harmony and harmumap are batch harmonised alignment of the CREST data only using harmony
3) pca, umap, are batch alignment of the CREST data only using Seurats integrate functions.

The file CREST_epithelia.rds (~1.5GB; https://drive.google.com/file/d/1VphkBKRnWA_2MBEpvqioVGqd9Gr1efvN/view?usp=share_link) is a Seurat object of the CREST epithelial cells aligned to a reference human endometrial dataset. With the following metadata (EpithelialAno = cell anotation, Idents = dataset/stage). Metadata (EpithelialAno = anotations, ID6 = datasets/stage)

The file CREST_stroma.rds (https://drive.google.com/file/d/1BiAa8r-tcEeiEEe0FRwf3bqLOWGYlY7l/view?usp=share_link) contains the aligned stromal cells from CREST to an endometrial reference (Idents = anotation, ID3 = batch).  

The file CREST_trophoblast_aligned.rds (https://drive.google.com/file/d/1PHpCaxXt54X7zorPuzuf3WyBny5VSaU5/view?usp=share_link) contains aligned trophoblast cells to human, rhesus, and marmoset references. Metadata includes Idents=inferred anotation, Species=dataset/species).

Finally, the CREST_aligned_to_human_rhesus_marmoset.rds (https://drive.google.com/file/d/1ZopPHsXsp6YcRO2k5b--V9ZJ3J8ezJxi/view?usp=share_link) contains an alignment of CREST embryonic lineages to human, rhesus, and marmoset reference datasets. 

