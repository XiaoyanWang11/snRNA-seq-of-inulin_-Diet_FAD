# -*- coding: utf-8 -*-
#Author: Xiaoyan Wang
#Copyright (c) 2024 __CarlosLab@CCMU/CIBR__. All rights reserved.

IB.BS_colors <- c("Neuron" = "palevioletred2", "Astrocyte" = "olivedrab2", "Microglia" = "mediumpurple1",
            "Oligo" = "skyblue", "OPC" = "chocolate1", "Pericyte" = "royalblue", "Unknown" = "turquoise")

### cell type markers for forebrain, interbrain and brainstem
markers <- c("Grin2b", "Rbfox3", "Grin1",    # neuron
             "Aldh1l1", "Aqp4",  # astrocyte
             "Cx3cr1",                    # microglia
             "Mbp", "Mog",        # Oligo
             "Pdgfra", "Vcan",    # OPC
             "Vtn")                       # pericyte

## brainstem
gg_title <- "brainstem"
Obj_Integrate <- RenameIdents(BS_Integrated, '0' = "Oligo", '1' = "Neuron", 
                              '2' =  "Neuron", '3' = "Astrocyte", '4' = "Microglia", '5' = "Neuron", 
                              '6' = "Neuron", '7' = "OPC", '8' = "Neuron", '9' = "Neuron", 
                              '10' = "Neuron", '11' = "Neuron", '12' = "Neuron", '13' = "Pericyte", '14' = "Unknown")
Obj_Integrate$celltype <- Idents(Obj_Integrate)
Idents(Obj_Integrate) <- "celltype"

Obj_Integrate@active.ident <- factor(x = Obj_Integrate@active.ident, 
                                     levels = c("Neuron", "Astrocyte", "Microglia", "Oligo", "OPC", "Pericyte", "Unknown"))

DimPlot(Obj_Integrate, reduction = "umap", label = T, label.color = "black", cols = IB.BS_colors, raster=FALSE) + ggtitle(gg_title)
ggsave(paste0(opt$opd, gg_title, "_annotated_DimPlot_",Sys.Date(),".png"), width = 10, height = 7)


### subset each cell type and re-clustering
## Neuron
gg_title <- "Neuron"
Obj_Neuron <- find_cluster(Obj_Integrate, idents = "Neuron", res = 0.1, gg_title = gg_title)
Obj_N <- find_cluster2(Obj_Neuron, idents = 7, res = 0.1, gg_title = gg_title)   # remove cluster 7

## Astrocyte
gg_title <- "Astrocyte"
Obj_Astrocyte <- find_cluster(Obj_Integrate, idents = "Astrocyte", res = 0.3, gg_title = gg_title)
Obj_A <- find_cluster2(Obj_Astrocyte, idents = c(3,5,7), res = 0.1, gg_title = gg_title)     # remove cluster 3,5,7

## Microglia
gg_title <- "Microglia"
Obj_Microglia <- find_cluster(Obj_Integrate, idents = "Microglia", res = 0.3, gg_title = gg_title)
Obj_Microglia2 <- find_cluster2(Obj_Microglia, idents = c(3,4,6), res = 0.1, gg_title = gg_title)     # remove cluster 3,4,6

## Oligodendrocyte
gg_title <- "Oligodendrocyte"
Obj_Oligo <- find_cluster(Obj_Integrate, idents = "Oligo", res = 0.3, gg_title = gg_title)
Obj_Ol <- find_cluster1(Obj_Oligo, gg_title = gg_title)

## OPC
gg_title <- "OPC"
Obj_OPC <- find_cluster(Obj_Integrate, idents = "OPC", res = 0.3, gg_title = gg_title)
Obj_OP <- find_cluster2(Obj_OPC, idents = c(4,5,6), res = 0.1, gg_title = gg_title)    # remove cluster 4,5,6
