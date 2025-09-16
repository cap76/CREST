Matlab files for mapping CREST data to CS6 embryo marmoset. Requires SpatialModelling codebase from https://github.com/Boroviak-Lab/SpatialModelling

MapHuman2Marmoset.m
Script loading in processed data and plotting cells projected onto the CS6 mamrmoset embryo.

ProjectData.m
Function for MapHuman2Marmoset.m that does the spatial mapping

Calculated indexes for SNNs for each cells (and a control permited set) are saved in the Files folder (C1_set1_withPreImpl_byCl.csv, C1perm_set1_withPreImpl_byCl.csv, S1_set1_all_withPreImpl_byCl.csv). These can be generated from an aligned Seurat object using the Align3D_withPreImpl.R script (example aligned datasets can be found here: https://drive.google.com/file/d/1ZopPHsXsp6YcRO2k5b--V9ZJ3J8ezJxi/view?usp=share_link). 

Additional processed intermediate datasets are available from: https://drive.google.com/drive/folders/1uV__ZgSRPSG03T-9CCFQggu6I-RVPOcS?usp=share_link
