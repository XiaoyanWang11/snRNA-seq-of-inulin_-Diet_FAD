# -*- coding: utf-8 -*-
#Author: Xiaoyan Wang
#Copyright (c) 2024 __CarlosLab@CCMU/CIBR__. All rights reserved.

FB_colors <- c("Neuron" = "palevioletred2", "Astrocyte" = "olivedrab2", "Microglia" = "mediumpurple1",
            "Oligo" = "skyblue", "OPC" = "chocolate1")

### cell type markers for forebrain, interbrain and brainstem
markers <- c("Grin2b", "Rbfox3", "Grin1",    # neuron
             "Aldh1l1", "Aqp4",  # astrocyte
             "Cx3cr1",                    # microglia
             "Mbp", "Mog",        # Oligo
             "Pdgfra", "Vcan")    # OPC
            
## forebrain
gg_title <- "forebrain"
Obj_Integrate <- RenameIdents(FB_Integrated, '0' = "Oligo", '1' = "Neuron", '2' =  "Microglia",         
                              '3' = "Neuron", '4' = "Neuron", '5' =   "Neuron", '6' = "Neuron", 
                              '7' = "Neuron", '8' = "Neuron", '9' = "Neuron", '10' = "Astrocyte", 
                              '11' = "OPC", '12' = "Neuron", '13' = "Neuron", '14' = "Neuron", 
                              '15' = "Neuron", '16' = "Neuron")
Obj_Integrate$celltype <- Idents(Obj_Integrate)
Idents(Obj_Integrate) <- "celltype"
Obj_Integrate@active.ident <- factor(x = Obj_Integrate@active.ident, levels = c("Neuron", "Astrocyte", "Microglia", "Oligo", "OPC"))

DimPlot(Obj_Integrate, reduction = "umap", label = T, label.color = "black", label.size = 8, cols = FB_colors, raster=FALSE)
ggsave(paste0(opt$opd, gg_title, "_annotated_DimPlot_", Sys.Date(),".png"), width = 10, height = 7)

VlnPlot(Obj_Integrate, features = markers, stack = T, flip=T, pt.size = 0,
        raster=FALSE, fill.by = "ident", split.by = "group_AD", split.plot = T, assay = "RNA")
ggsave(paste0(opt$opd, gg_title, "_annotation.celltype.VlnPlot_", Sys.Date(), ".png"), width = 6, height = 6)

N <- table(Obj_Integrate$orig.ident, Obj_Integrate@active.ident)
write.csv(N, paste0(opt$opd, gg_title, "_nucleiNumbers_", Sys.Date(), ".csv"))
  Quan <- N %>% as.data.frame.matrix()
  Quan$Sum <- rowSums(Quan)
  Quan <- Quan/Quan$Sum
  Quan <- Quan[,-ncol(Quan)]
  write.csv(Quan, paste0(opt$opd, gg_title, "_cellproportion_", Sys.Date(), ".csv"))

### subset each cell type and re-clustering
if(T){
  # subset celltype from integrated dataset
  find_cluster <- function(Seurat_Obj, nPCs = 15, idents, res, gg_title){
    Obj <- subset(Seurat_Obj, idents = idents)
    Obj <- ScaleData(Obj)  # gene 
    Obj <- RunPCA(Obj)
    Obj <- RunUMAP(Obj, reduction = "pca", dims = 1:nPCs)
    Obj <- FindNeighbors(Obj, dims = 1:nPCs)
    Obj <- FindClusters(Obj, resolution = res)
    return(Obj)
  }
  
  # save re-clustering result
  find_cluster1 <- function(Obj, gg_title){
    saveRDS(Obj, file = paste0(opt$opd, gg_title, "_", Sys.Date(), ".rds"))
    DimPlot(Obj, reduction = "umap", group.by = "ident", label = T, repel = T,
            label.color = "black", raster=FALSE) + ggtitle(gg_title)
    ggsave(paste0(opt$opd, gg_title, "_DimPlot_", Sys.Date(),".png"), width = 6, height = 4)
    
    Quan <- table(Obj$orig.ident,Obj@active.ident) %>% as.data.frame.matrix()
    Quan$Sum <- rowSums(Quan)
    Quan <- Quan/Quan$Sum
    Quan <- Quan[,-ncol(Quan)]
    write.csv(Quan, paste0(opt$opd, gg_title, "_cellproportion_", Sys.Date(), ".csv"))
    return(Obj)
  }
  
  # remove mix/unknown cluster after subset
  find_cluster2 <- function(Seurat_Obj, nPCs = 15, idents, res, gg_title){
    Obj <- subset(Seurat_Obj, idents = idents, invert = TRUE)
    Obj <- ScaleData(Obj)
    Obj <- RunPCA(Obj)
    Obj <- RunUMAP(Obj, reduction = "pca", dims = 1:nPCs)
    Obj <- FindNeighbors(Obj, dims = 1:nPCs)
    Obj <- FindClusters(Obj, resolution = res)
    saveRDS(Obj, file = paste0(opt$opd, gg_title, "_", Sys.Date(), ".rds"))
    DimPlot(Obj, reduction = "umap", group.by = "ident", label = T, repel = T,
            label.color = "black", raster=FALSE) + ggtitle(gg_title)
    ggsave(paste0(opt$opd, gg_title, "_DimPlot_", Sys.Date(),".png"), width = 6, height = 4)
    
    Quan <- table(Obj$orig.ident,Obj@active.ident) %>% as.data.frame.matrix()
    Quan$Sum <- rowSums(Quan)
    Quan <- Quan/Quan$Sum
    Quan <- Quan[,-ncol(Quan)]
    write.csv(Quan, paste0(opt$opd, gg_title, "_subtype.proportion_", Sys.Date(), ".csv"))
    return(Obj)
  }
}

## subset neuron
gg_title <- "Neuron"
Obj_Neuron <- find_cluster(Obj_Integrated, idents = "Neuron", res = 0.1, gg_title = gg_title)
Obj_N <- find_cluster2(Obj_Neuron, idents = 18, res = 0.1, gg_title = gg_title)    # remove cluster 18
Obj_N <- RenameIdents(Obj_N, '0' = "Ex", '1' = "Ex", '2' =  "In", 
                              '3' = "In", '4' = "Ex", '5' =   "Ex", '6' = "In", 
                              '7' = "Ex", '8' = "In", '9' = "Ex", '10' = "In", 
                              '11' = "Ex", '12' = "Ex", '13' = "Ex", '14' = "Ex", 
                              '15' = "In", '16' = "Ex", '17' = "In", '18' = "Ex")
Obj_N$subtype <- Idents(Obj_N)
Idents(Obj_N) <- "subtype"

## subset Astrocyte
gg_title <- "Astrocyte"
Obj_Astrocyte <- find_cluster(Obj_Integrate, idents = "Astrocyte", res = 0.3, gg_title = gg_title)
Obj_A <- find_cluster1(Obj_Astrocyte, gg_title = gg_title)

## subset Microglia
gg_title <- "Microglia"
Obj_Microglia <- find_cluster(Obj_Integrate, idents = "Microglia", res = 0.3, gg_title = gg_title)
Obj_Microglia2 <- find_cluster2(Obj_Microglia, idents = 6, res = 0.1, gg_title = gg_title)  # remove cluster 6

## subset Oligodendrocyte
gg_title <- "Oligodendrocyte"
Obj_Oligo <- find_cluster(Obj_Integrate, idents = "Oligo", res = 0.3, gg_title = gg_title)
Obj_Ol <- find_cluster1(Obj_Oligo, gg_title = gg_title)

## subset OPC
gg_title <- "OPC"
Obj_OPC <- find_cluster(Obj_Integrate, idents = "OPC", res = 0.3, gg_title = gg_title)
Obj_OP <- find_cluster2(Obj_OPC, idents = 6, res = 0.1, gg_title = gg_title)  # remove cluster 6

### differential expression genes(DEGs) between diet groups in each celltype  
markers <- FindMarkers(Obj_A, ident.1 = "FBLI", ident.2 = "FB93", assay = "RNA") %>% filter(p_val_adj < 0.05)
deg <- markers
up_DEG <- row.names(deg[deg$avg_log2FC > 0.3 & deg$p_val_adj < 0.05,]) 
down_DEG <- row.names(deg[deg$avg_log2FC < -0.3 & deg$p_val_adj < 0.05,]) 

UP_GO <- enrichGO(gene = up_DEG, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = T)    
DOWN_GO <- enrichGO(gene = down_DEG, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                    pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = T)

# select relevant GO-BP pathways and plot top 10 pathways for each celltype
up_down <- rbind(Neuron_UP_GO, Astrocyte_UP_GO, Astrocyte_DOWN_GO,
                 Microglia_UP_GO, Microglia_DOWN_GO, Oligo_UP_GO, 
                 Oligo_DOWN_GO, OPC_UP_GO, OPC_DOWN_GO)  
up_down <- up_down %>% filter(up_down$p.adjust < 0.05) 
ggplot(up_down, aes(Description, -log10(p.adjust), fill = change)) +
     geom_bar(position = position_dodge(), stat = "identity", width = 0.6) +  
      scale_fill_manual(
        values = c("up" = "#a32a31", "down" = "#3665a6"),
        labels = c("Up-regulated pathways in Inulin_AD group", "Down-regulated pathways in Inulin_AD group")) +
      labs(title = "GO pathways", y = "-log10(p.adjust)")
ggsave(paste0(opt$opd, gg_title, "_DEGs.GO-BP_pathway_", Sys.Date(), ".png"), width = 14, height = 8)

