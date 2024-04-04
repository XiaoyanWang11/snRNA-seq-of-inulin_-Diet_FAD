# -*- coding: utf-8 -*-
#Author: Xiaoyan Wang
#Copyright (c) 2024 __CarlosLab@CCMU/CIBR__. All rights reserved.

### cell type markers for cerebellum
CB_markers <- c("Gabra6",   # Granule
             "Lypd6", "Prkcd",   # molecular layer interneurons (MLIs)
             "Nxph1",   # MLI2
             "Mobp", "Mbp", "Mog",     # Oligo/ODC
             "Gdf10",    # Bergmann
             "Aldh1l1", "Aqp4",     # Astrocyte
             "Cx3cr1",     # Microglia
             "Lgi2",     # Golgi
             "Eomes",    # excitatory unipolar brush cells (UBCs)
             "Ppp1r17", # Purkinje
             "Klhl1",   # Purkinje layer interneurons (PLIs), 
             "Pdgfra",   # OPC
             "Ttr")      # Choroid
## cerebellum
gg_title <- "cerebellum"
Obj_Integrate <- RenameIdents(Obj_Integrate, '0' = "Granule", '1' = "Granule", '2' =  "Granule", 
                              '3' = "MLI1", '4' = "Oligo", '5' = "Bergmann", '6' = "MLI2", 
                              '7' = "Astrocyte", '8' = "Microglia", '9' = "Golgi", '10' = "UBC", 
                              '11' = "Purkinje", '12' = "OPC", '13' = "Choroid")
Obj_Integrate$celltype <- Idents(Obj_Integrate)
Idents(Obj_Integrate) <- "celltype"

Obj_Integrate@active.ident <- factor(x = Obj_Integrate@active.ident, 
                                     levels = c("Granule", "MLI1", "MLI2", "Oligo", "Bergmann", 
                                                "Astrocyte", "Microglia", "Golgi", "UBC", "Purkinje", "OPC", "Choroid"))

DimPlot(Obj_Integrate, reduction = "umap", label = T, label.color = "black", cols = CB_colors, raster=FALSE) + ggtitle(gg_title)
ggsave(paste0(opt$opd, gg_title, "_annotated_DimPlot_",Sys.Date(),".png"), width = 10, height = 7)

### subset each cell type and re-clustering, keep all clusters after reclustering, resolution=0.1
