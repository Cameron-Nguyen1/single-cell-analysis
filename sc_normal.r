#!/bin/env Rscript
#Normal workflow
library(pacman,lib="")
dep = c("plyr","data.table","devtools","Seurat","tidyverse","miQC","SeuratWrappers","flexmix","SingleCellExperiment","SummarizedExperiment","RColorBrewer","ggsankey",
"ggplot2","cowplot","SingleR","scran","celldex","ComplexHeatmap","pheatmap","circlize","GGally","forcats","dplyr","patchwork","pals","harmony",
"ggpubr","paletteer","ggridges","fgsea","UCell","FactoMineR","factoextra","zeallot","gridExtra","grid","optparse")
#lapply(dep,install.packages,character.only=TRUE,repos='http://cran.us.r-project.org',INSTALL_opts='--no-lock',lib="")
p_load(dep,character.only=TRUE,lib="",install=FALSE,update=FALSE)
source("/sc_functions.r")

S_Combined = readRDS(paste0(getwd(),"/S_Combined_postQC_stage1.rds"))
DefaultAssay(S_Combined) = "RNA"
S_Combined = Seurat::NormalizeData(S_Combined, normalization.method = "LogNormalize", scale.factor = 10000)
S_Combined = Seurat::FindVariableFeatures(S_Combined, selection.method = "vst", nfeatures = 4000)
S_Combined = Seurat::ScaleData(S_Combined)
S_Combined=RunPCA(S_Combined,npcs=50) #Test 50 dimensions
bestdims_c = bdim(S_Combined)
S_Combined = RunUMAP(S_Combined, reduction = "pca", dims = 1:bestdims_c, verbose = F) ##UMAP Cell Cycle Scores
S_Combined = FindNeighbors(S_Combined, reduction = "pca", dims = 1:bestdims_c)
S_Combined = FindClusters(S_Combined, resolution = 0.2)

S_Combined@meta.data$CC.Diff = S_Combined@meta.data$S.Score - S_Combined@meta.data$G2M.Score ##PCA Cell Cycle Scores
S_Combined@meta.data$Full = "All CMOs"
pc.df = cbind(S_Combined@reductions$pca@cell.embeddings, S_Combined@meta.data)
p1 = ggplot(data = pc.df, aes(x = PC_1, y = PC_2, color = S.Score)) + geom_point(size = 0.1) +
  scale_color_viridis_c("S.Score", option = "inferno", direction = -1) + facet_grid(.~Full)
p2 = ggplot(data = pc.df, aes(x = PC_1, y = PC_2, color = S.Score)) + geom_point(size = 0.1) +
  scale_color_viridis_c("S.Score", option = "inferno", direction = -1) + facet_grid(.~orig.ident)
p3 = ggplot(data = pc.df, aes(x = PC_1, y = PC_2, color = G2M.Score)) + geom_point(size = 0.1) +
  scale_color_viridis_c("G2M.Score", option = "inferno", direction = -1) + facet_grid(.~Full)
p4 = ggplot(data = pc.df, aes(x = PC_1, y = PC_2, color = G2M.Score)) + geom_point(size = 0.1) +
  scale_color_viridis_c("G2m.Score", option = "inferno", direction = -1) + facet_grid(.~orig.ident)
pdf("Cell_Cycle_Scores_PCA.pdf",width=21,height=12)
p1 + p2 + patchwork::plot_layout(ncol = 2, widths = c(1, 2), guides="collect")
p3 + p4 + patchwork::plot_layout(ncol = 2, widths = c(1, 2), guides="collect")
dev.off()

umap.df = cbind(S_Combined@reductions$umap@cell.embeddings, S_Combined@meta.data)
p1 = ggplot(data = umap.df, aes(x = UMAP_1, y = UMAP_2, color = S.Score)) + geom_point(size = 0.1) +
  scale_color_viridis_c("S.Score", option = "inferno", direction = -1) + facet_grid(.~Full)
p2 = ggplot(data = umap.df, aes(x = UMAP_1, y = UMAP_2, color = S.Score)) + geom_point(size = 0.1) +
  scale_color_viridis_c("S.Score", option = "inferno", direction = -1) + facet_grid(.~orig.ident)
p3 = ggplot(data = umap.df, aes(x = UMAP_1, y = UMAP_2, color = G2M.Score)) + geom_point(size = 0.1) +
  scale_color_viridis_c("G2M.Score", option = "inferno", direction = -1) + facet_grid(.~Full)
p4 = ggplot(data = umap.df, aes(x = UMAP_1, y = UMAP_2, color = G2M.Score)) + geom_point(size = 0.1) +
  scale_color_viridis_c("G2m.Score", option = "inferno", direction = -1) + facet_grid(.~orig.ident)
pdf("Cell_Cycle_Scores_UMAP.pdf",width=21,height=12)
p1 + p2 + patchwork::plot_layout(ncol = 2, widths = c(1, 2), guides="collect")
p3 + p4 + patchwork::plot_layout(ncol = 2, widths = c(1, 2), guides="collect")
dev.off()

pdf("Cells_By_Phase_Barplot.pdf",width=8)
cell_cycle.cols = c(`G1`="#003049", `G2M`="#d62828", `S`="#f77f00")
ggplot(data = S_Combined@meta.data, aes(x = orig.ident, fill = Phase)) +
  geom_bar(position = "fill", color = "black") + scale_fill_manual("", values = cell_cycle.cols) + xlab("") + ylab("Proportion of cells")
dev.off()

p1 = DimPlot(S_Combined, group.by = "Phase", cols = cell_cycle.cols, reduction = "pca", split.by = "Full") + ggtitle("")
p2 = DimPlot(S_Combined, group.by = "Phase", cols = cell_cycle.cols, reduction = "pca", split.by = "orig.ident") + ggtitle("")
a1 = p1 + p2 + patchwork::plot_layout(ncol = 2, widths = c(1, 2), guides = "collect")
p3 = DimPlot(S_Combined, group.by = "Phase", cols = cell_cycle.cols, split.by = "Full") + ggtitle("")
p4 = DimPlot(S_Combined, group.by = "Phase", cols = cell_cycle.cols, split.by = "orig.ident") + ggtitle("")
a2 = p3 + p4 + patchwork::plot_layout(ncol = 2, widths = c(1, 2), guides = "collect")
ggsave(filename="Cells_By_Phase_PCA.pdf",plot=a1,width=16)
ggsave(filename="Cells_By_Phase_UMAP.pdf",plot=a2,width=16)
saveRDS(S_Combined,"S_Combined_stage2.rds")