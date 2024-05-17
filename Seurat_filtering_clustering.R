# -*- coding: utf-8 -*-
#Author: Xiaoyan Wang
#Copyright (c) 2024 __CarlosLab@CCMU/CIBR__. All rights reserved.

# Generate seurat object, filtering doublets, integrate and clustering in FAD fiber project

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
})

opt <- list(input = "/../outs/CellBender_feature_bc_matrix_filtered.h5",
            opd = "/diretory_to_output/")

## forebrain, the same procedure for other 3 regions
samples <- list.dirs(opt$input, recursive = F, full.names = F)

## read h5 file and create Seurat object
FAD_FB.list <- lapply(samples, function(i){
  print(i)
  counts <- Read10X_h5(paste0(opt$input, i, "/outs/CellBender_feature_bc_matrix_filtered.h5"))
  sce <- CreateSeuratObject(counts, project=i, min.cells=10, min.features = 500)
  return(sce)
})
names(FAD_FB.list) =  samples

## Filtering low-quality nuclei and doublets
RunSeurat <- function(ObjList, Filter_nFeature_RNA = c(500, 6000),
                      Filter_nCount_RNA = c(1000, 50000), MitoPercent = 0.2,
                      nPCs = 15, DoubletPercent = 0.08) {
  ObjList <- lapply(X = ObjList, FUN = function(Obj) {
    #1. Low-quality cells or empty droplets will often have very few genes. 
    #2. Cell doublets or multiplets may exhibit an aberrantly high gene count. 
    #3. Low-quality / dying cells often exhibit extensive mitochondrial contamination
    Obj[["percent_mt"]] <- PercentageFeatureSet(Obj, pattern = "^mt-")
    Obj <- subset(Obj, subset =
                    nFeature_RNA > Filter_nFeature_RNA[1] & nFeature_RNA < Filter_nFeature_RNA[2]
                  & nCount_RNA > Filter_nCount_RNA[1] & nCount_RNA < Filter_nCount_RNA[2]
                  & percent_mt < MitoPercent)
    Obj <- NormalizeData(Obj, normalization.method = "LogNormalize", scale.factor = 10000)
    Obj <- FindVariableFeatures(Obj, selection.method = "vst", nfeatures = 2000)
    Obj <- ScaleData(Obj)
    Obj <- RunPCA(Obj)
    Obj <- RunUMAP(Obj, dims = 1:nPCs, reduction = "pca", verbose = T)
  })
  for(i in 1:length(ObjList)){
    # find optimal_pK
    sweep.res <- paramSweep_v3(ObjList[[i]], PCs = 1:nPCs, sct = FALSE)
    sweep <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- find.pK(sweep)
    optimal_pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    
    # find nExp_poi
    nExp_poi <- round(DoubletPercent*nrow(ObjList[[i]]@meta.data))
    ObjList[[i]] <- doubletFinder_v3(ObjList[[i]], PCs = 1:15, pN = 0.25, pK = optimal_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    Idents(ObjList[[i]]) <- ObjList[[i]]@meta.data[,paste0("DF.classifications_0.25_",optimal_pK,"_",nExp_poi)]
    p <- DimPlot(ObjList[[i]], reduction = "umap", label = T, repel = T, label.box =T, label.color = "white")
    plot(p)
    # Obj@meta.data[,ncol(Obj@meta.data)] 即 Obj@meta.data[,9] 
     ToRemove <- rownames( ObjList[[i]]@meta.data[ObjList[[i]]@meta.data[,ncol(ObjList[[i]]@meta.data)] == "Doublet",] )
    ObjList[[i]] <- ObjList[[i]][ ,!colnames(ObjList[[i]]) %in% ToRemove]
    cat(paste0(length(ToRemove)," doublets were removed !!!\n"))
  }
  return(ObjList)
}
FB_DoubletRD <- RunSeurat(FB.list)

# integrate data
ObjList <- lapply(X = FB_DoubletRD, FUN = function(Obj) {
  Obj <- FindVariableFeatures(Obj, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = ObjList, nfeatures = 2000) # nfeatures = Obj@assays$RNA@var.features
anchors <- FindIntegrationAnchors(object.list = ObjList, anchor.features = features)
Obj_Integrated <- IntegrateData(anchorset = anchors)
DefaultAssay(Obj_Integrated) <- "integrated"

# clustering
nPCs = 15
Obj_Integrated <- ScaleData(Obj_Integrated)
Obj_Integrated <- RunPCA(Obj_Integrated, npcs = nPCs)
Obj_Integrated <- RunUMAP(Obj_Integrated, dims = 1:nPCs, reduction = "pca", verbose = T)
Obj_Integrated <- FindNeighbors(Obj_Integrated, reduction = "pca", dims = 1:nPCs)
Obj_Integrated <- FindClusters(Obj_Integrated, resolution = 0.1)
# group by diets
control_FB <- c("Mu93FB_1", "Mu93FB_2", "Mu93FB_3")
Obj_Integrated$group <- sapply(Obj_Integrated$orig.ident, function(x) ifelse (x %in% control_FB, "FB93", "FBLI"))
saveRDS(Obj_Integrated, file = paste0(opt$opd, "FB_Integrated_res0.1_", Sys.Date(), ".rds"))

# markers of each cluster
gg_title <- "forebrain"
AllMarkers <- FindAllMarkers(Obj_Integrated, only.pos = TRUE, min.pct = 0.25,  # upregulation
                               logfc.threshold = 0.25, assay = "RNA") %>% filter(p_val_adj < 0.05)
write.csv(AllMarkers, file = paste0(opt$opd, gg_title, "_Allmarkers_",Sys.Date(),".csv"), quote = F)   #  Idents(Obj) <- "seurat_clusters"

                               TopMarkersList <- AllMarkers %>% group_by(cluster) %>% slice_max(n = NMarkers, order_by = avg_log2FC) %>% as.data.frame()
write.csv(TopMarkersList, file = paste0(opt$opd, gg_title, "_Top5markers_",Sys.Date(),".csv"), quote = F)

                               





    
