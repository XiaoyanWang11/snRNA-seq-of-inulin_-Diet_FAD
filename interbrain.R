IB.BS_colors <- c("Neuron" = "palevioletred2", "Astrocyte" = "olivedrab2", "Microglia" = "mediumpurple1",
            "Oligo" = "skyblue", "OPC" = "chocolate1", "Pericyte" = "royalblue", "Unknown" = "turquoise")

### cell type markers for forebrain, interbrain and brainstem
markers <- c("Grin2b", "Rbfox3", "Grin1",    # neuron
             "Aldh1l1", "Aqp4",  # astrocyte
             "Cx3cr1",                    # microglia
             "Mbp", "Mog",        # Oligo
             "Pdgfra", "Vcan",    # OPC
             "Vtn")                       # pericyte

gg_title <- "interbrain"
Obj_Integrate <- RenameIdents(IB_Integrated, '0' = "Oligo", '1' = "Neuron", 
                                  '2' =  "Neuron", '3' = "Astrocyte", '4' = "Microglia", '5' = "Oligo", '6' = "OPC", 
                                  '7' = "Neuron", '8' = "Neuron", '9' = "Neuron", '10' = "Pericyte", 
                                  '11' = "Unknown", '12' = "Unknown")
Obj_Integrate$celltype <- Idents(Obj_Integrate)
Idents(Obj_Integrate) <- "celltype"

Obj_Integrate@active.ident <- factor(x = Obj_Integrate@active.ident, levels = c("Neuron", "Astrocyte", "Microglia", "Oligo", "OPC", "Pericyte", "Unknown"))
DimPlot(Obj_Integrate, reduction = "umap", label = T, label.color = "black", label.size = 8, cols = IB.BS_colors, raster=FALSE)
ggsave(paste0(opt$opd, gg_title, "_annotated_DimPlot_", Sys.Date(),".png"), width = 10, height = 7)


## Neuron
gg_title <- "Neuron"
Obj_Neuron <- find_cluster(Obj_Integrate, idents = "Neuron", res = 0.1, gg_title = gg_title)
Obj_N <- find_cluster1(Obj_Neuron, gg_title = gg_title)

## Astrocyte
gg_title <- "Astrocyte"
Obj_Astrocyte <- find_cluster(Obj_Integrate, idents = "Astrocyte", res = 0.3, gg_title = gg_title)
Obj_A <- find_cluster2(Obj_Astrocyte, idents = c(4,6), res = 0.3, gg_title = gg_title)       # remove cluster4,6

## Microglia
gg_title <- "Microglia"
Obj_Microglia <- find_cluster(Obj_Integrate, idents = "Microglia", res = 0.3, gg_title = gg_title)
Obj_Microglia2 <- find_cluster2(Obj_Microglia, idents = 4, res = 0.1, gg_title = gg_title)  # remove cluster4

## Oligodendrocyte
gg_title <- "Oligodendrocyte"
Obj_Oligo <- find_cluster(Obj_Integrate, idents = "Oligo", res = 0.3, gg_title = gg_title)
Obj_Ol <- find_cluster1(Obj_Oligo, gg_title = gg_title)

## OPC
gg_title <- "OPC"
Obj_OPC <- find_cluster(Obj_Integrate, idents = "OPC", res = 0.3, gg_title = gg_title)
Obj_OP <- find_cluster2(Obj_OPC, idents = c(4,5), res = 0.1, gg_title = gg_title)    # remove cluster 4,5


### DEGs of cluster0 vs cluster2
markers <- FindMarkers(Obj, ident.1 = "0", ident.2 = "2", assay = "RNA") %>% filter(p_val_adj < 0.05)
deg <- markers
up_DEG <- row.names(deg[deg$avg_log2FC > 0.3 & deg$p_val_adj < 0.05,]) 
down_DEG <- row.names(deg[deg$avg_log2FC < -0.3 & deg$p_val_adj < 0.05,])
UP_GO <- enrichGO(gene = up_DEG, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                      pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = T)
DOWN_GO <- enrichGO(gene = down_DEG, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                      pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = T)
# select relevant pathways and plot top 10
up_down <- rbind(head(up_go, 10L), head(down_go, 10L))
up_down <- up_down %>% 
   filter(up_down$p.adjust < 0.05) 
ggplot(up_down, aes(-log10(p.adjust), Description, fill = change)) +
      geom_bar(position = position_dodge(), stat = "identity", width = 0.6) +  
      scale_fill_manual(
        values = c("up" = "#a32a31", "down" = "#3665a6"),
        labels = c("Up-regulated pathways", "Down-regulated pathways")) +
      ggtitle("GO pathways")
ggsave(paste0(opt$opd, gg_title, "_DEGs.GO-BP_", Sys.Date(), ".png"), width = 11, height = 8)

