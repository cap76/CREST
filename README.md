# CREST
Code for analysis of data from cell-engineered receptive endometrial scaffold technology (CREST) by Molè, Elderkin, Zorzan, Penfold, Horsley et al.

Here we provide R scripts and intermediate processed data for the paper Modelling human embryo implantation in vitro (Rscripts folder) and MATLAB scripts for projecting processed data onto CS6 marmoset embryos (MATLAB folder).

A processeed Seurat object (~1GB) can be downloaded from here: https://drive.google.com/file/d/1MYrjAxGue18yibqu1eg91ePOFfaYeA6U/view?usp=share_link

This file holds 15673 cells with the following meta-information:
1) Idents = inferred anotation
2) ID3 = processed batch

Spreadsheets of this information can also be found in the Files folder (Meta.csv).

It has the following dimensitonality reductions:
1) oldUMAP is an aligned dimensionality reduction to a composite reference (human cycling endometrium and human in vitro cultured embryos. Used for Figure 3. 
2) harmony and harmumap are batch harmonised alignment of the CREST data only using harmony
3) pca, umap, are batch alignment of the CREST data only using Seurats integrate functions.


