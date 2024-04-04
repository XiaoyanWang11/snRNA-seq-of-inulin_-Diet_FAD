
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(DoubletFinder)
library(sctransform)

opt <- list(input_FB = "/../FAD_FB/",
            input_MB = "/../FAD_MB/",
            input_BS = "/../FAD_BS/",
            opd = "/../FAD_FB.MB.BS/")
### import rds file which has been removed doublets by DoubletFinder function
CellBender_FB <- readRDS(paste0(opt$input_FB, "FAD_FB_DoubletRD.rds"))
CellBender_MB <- readRDS(paste0(opt$input_MB, "FAD_MB_DoubletRD.rds"))
CellBender_BS <- readRDS(paste0(opt$input_BS, "FAD_BS_DoubletRD.rds"))

# Integrate 3 regions

FAD_DoubletRD <- c(CellBender_FB, CellBender_MB, CellBender_BS)

ObjList <- lapply(X = FAD_DoubletRD, FUN = function(Obj) {
  Obj <- FindVariableFeatures(Obj, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = ObjList, nfeatures = 2000) # nfeatures = Obj@assays$RNA@var.features
anchors <- FindIntegrationAnchors(object.list = ObjList, anchor.features = features)
Obj_Integrated <- IntegrateData(anchorset = anchors)
DefaultAssay(Obj_Integrated) <- "integrated"

nPCs = 15
Obj_Integrated <- ScaleData(Obj_Integrated)
Obj_Integrated <- RunPCA(Obj_Integrated, npcs = nPCs)
Obj_Integrated <- RunUMAP(Obj_Integrated, dims = 1:nPCs, reduction = "pca", verbose = T)
Obj_Integrated <- FindNeighbors(Obj_Integrated, reduction = "pca", dims = 1:nPCs)
Obj_Integrated <- FindClusters(Obj_Integrated, resolution = 0.3)

# add diet group
group_diet <- data.frame(orig.ident = c('Mu93FB_1', 'Mu93FB_2', 'Mu93FB_3', 'Mu93MB_1', 'Mu93MB_2', 'Mu93MB_3',
                                        'Mu93BS_1', 'Mu93BS_2', 'Mu93BS_3',
                                        'MuLIFB_1' , 'MuLIFB_2' , 'MuLIFB_3', 'MuLIMB_1' , 'MuLIMB_2' , 'MuLIMB_3',
                                        'MuLIBS_1' , 'MuLIBS_2' , 'MuLIBS_3'),
                         group = c('Ctrl_AD','Ctrl_AD','Ctrl_AD', 'Ctrl_AD','Ctrl_AD','Ctrl_AD',
                                   'Ctrl_AD','Ctrl_AD','Ctrl_AD', 
                                   'Inulin_AD','Inulin_AD','Inulin_AD', 'Inulin_AD','Inulin_AD','Inulin_AD',
                                   'Inulin_AD','Inulin_AD','Inulin_AD'))

orig <- data.frame(orig.ident = Obj_Integrated@meta.data$orig.ident)
new1 <- merge(orig, group_diet, by = "orig.ident")

Obj_Integrated$diet <- new1$group
Obj_Integrated[['diet']] <- new1$group[match(Obj_Integrated$orig.ident, new1$orig.ident)]

# add region group
group_region <- data.frame(orig.ident = c('Mu93FB_1', 'Mu93FB_2', 'Mu93FB_3', 'MuLIFB_1' , 'MuLIFB_2' , 'MuLIFB_3',
                                          'Mu93MB_1', 'Mu93MB_2', 'Mu93MB_3', 'MuLIMB_1' , 'MuLIMB_2' , 'MuLIMB_3',
                                          'Mu93BS_1', 'Mu93BS_2', 'Mu93BS_3', 'MuLIBS_1' , 'MuLIBS_2' , 'MuLIBS_3'),
                           group = c('forebrain','forebrain','forebrain', 'forebrain','forebrain','forebrain',
                                     'interbrain','interbrain','interbrain', 'interbrain','interbrain','interbrain',
                                     'brainstem','brainstem','brainstem', 'brainstem','brainstem','brainstem'))

orig <- data.frame(orig.ident = Obj_Integrated@meta.data$orig.ident)
new1 <- merge(orig, group_region, by = "orig.ident")
Obj_Integrated$region <- new1$group
Obj_Integrated[['region']] <- new1$group[match(Obj_Integrated$orig.ident, new1$orig.ident)]
saveRDS(Obj_Integrated, file = paste0(opt$opd, "FAD_FB.MB.BS_res0.3_", Sys.Date(), ".rds"))

combine_colors <- c("Neuron" = "palevioletred2", "Astrocyte" = "olivedrab2", "Microglia" = "mediumpurple1",
            "Oligo" = "skyblue", "OPC" = "chocolate1", "Pericyte" = "royalblue",
            "Endothelium" = "goldenrod", "Unknown" = "turquoise")

### annotate each cluster
gg_title <- "FAD_region_combine"
markers <- c("Grin2b", "Rbfox3", "Grin1",          # Neuron
             "Mbp", "Mog",               # Oligo
             "Aldh1l1", "Aqp4",         # Astrocyte
             "Cx3cr1",                    # Microglia
             "Pdgfra", "Vcan",        # OPC
             "Vtn",                               # Pericyte
             "Flt1", "Pecam1")                    # EC
             
VlnPlot(Obj_Integrated, features = markers, stack = T, flip=T, pt.size = 0, 
        raster=FALSE, fill.by = "ident", split.by = "diet", split.plot = TRUE, assay = "RNA")

Obj_Integrated <- RenameIdents(Obj_Integrated, '0' = "Neuron", '1' = "Oligo", '2' = "Neuron", '3' = "Astrocyte", 
                    '4' = "Oligo", '5' = "Microglia", '6' = "Neuron", '7' = "Microglia", 
                    '8' = "Neuron", '9' = "OPC", '10' = "Neuron", '11' = "Neuron",
                    '12' = "Neuron", '13' = "Neuron", '14' = "Neuron", '15' = "Neuron",
                    '16' = "Neuron", '17' = "Pericyte", '18' = "Neuron", '19' = "Neuron", 
                    '20' = "Unknown", '21' = "Unknown", '22' = "Endothelium", 
                    '23' = "Unknown", '24' = "Unknown", '25' = "Unknown")
Obj_Integrated$celltype <- Idents(Obj_Integrated)
Idents(Obj_Integrated) <- "celltype"

DimPlot(Obj_Integrated, reduction = "umap", label = T, repel = T, label.color = "black", 
        label.size = 6, group.by = "ident", split.by = "region",  raster=FALSE)
ggsave(paste0(opt$opd, gg_title, "_DimPlot_",Sys.Date(),".png"))

## subset each cell type and re-clustering
# subset Astrocyte
gg_title <- "Astrocyte"
Obj_Astrocyte <- find_cluster(Obj_Integrated, idents = "Astrocyte", res = 0.3, gg_title = gg_title)
Obj_A <- find_cluster2(Obj_Astrocyte, idents = c(6, 7), res = 0.3, gg_title = gg_title)

# subset Microglia
gg_title <- "Microglia"
Obj_Microglia <- find_cluster(Obj_Integrated, idents = "Microglia", res = 0.1, gg_title = gg_title)
Obj_Microglia2 <- find_cluster1(Obj_Microglia, gg_title = gg_title)


## GO-BP pathway analysis

# find markers for each cluster
AllMarkers <- FindAllMarkers(Obj, only.pos = TRUE, min.pct = 0.25,  # upregulation
                             logfc.threshold = 0.3, assay = "RNA") %>% filter(p_val_adj < 0.05)
# GO-BP pathway for each cluster
DEG_GO <- function(Seurat_DEG_file, nGO =10, gg_title){
  DEG_list <- split(Seurat_DEG_file$gene, f = Seurat_DEG_file$cluster)
  
  # convert gene symbol to Entrez id
  GeneID_list_ID <- DEG_list %>% map(~{
    gene.df <- AnnotationDbi::select(org.Mm.eg.db, keys = .x, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
    gene <- gene.df$ENTREZID
    gene <- gene[which(!is.na(gene))]
    gene <- unique(gene)
    return(gene)
  })
  
  # GO analysis
  compGO <- compareCluster(geneCluster = DEG_list, fun= "enrichGO", pvalueCutoff = 0.05, keyType = "SYMBOL",
                           pAdjustMethod = "BH", OrgDb = org.Mm.eg.db, readable = T, ont = "BP")
  write.csv(compGO@compareClusterResult, file = paste0(opt$opd, gg_title,"_GO-BP_", Sys.Date(), ".csv"))
  TopGOList <- compGO@compareClusterResult %>% group_by(Cluster) %>% slice_min(n = 10, order_by = p.adjust)
  write.csv(TopGOList, file = paste0(opt$opd, gg_title,"_TopGO-BP_", Sys.Date(), ".csv"))
  
  TopGOList$Description <- factor(TopGOList$Description, 
                                   levels = unique(TopGOList$Description[order(TopGOList$Cluster)]))
  ggplot(data = TopGOList, 
         aes(x = Cluster, y = Description, 
             color = Cluster, size = p.adjust)) + 
    geom_point(alpha=0.8) +
    guides(colour = guide_legend(order = 1), 
           size = guide_legend(order = 2)) +  
    theme(axis.title.y.left = element_blank())
  ggsave(paste0(opt$opd, gg_title,"_TopGO-BP_", Sys.Date(), ".png"), width = 8, height = 10)
}     
gg_title <- "Astrocyte"; Obj <- Obj_A; DEG <- DEG_GO(Astrocyte_AllMarkers, gg_title = gg_title)












