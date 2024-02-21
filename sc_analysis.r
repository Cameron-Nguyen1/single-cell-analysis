#!/bin/env Rscript
library(pacman,lib="/your/R/Library")
library(optparse,lib="/your/R/Library")
dep = c("plyr","data.table","devtools","Seurat","tidyverse","miQC","SeuratWrappers","flexmix","SingleCellExperiment","SummarizedExperiment","RColorBrewer","ggsankey",
"ggplot2","cowplot","SingleR","scran","celldex","ComplexHeatmap","pheatmap","circlize","GGally","forcats","dplyr","patchwork","pals","harmony",
"ggpubr","paletteer","ggridges","fgsea","UCell","FactoMineR","factoextra","zeallot","gridExtra","grid","SeuratDisk")
#lapply(dep,install.packages,character.only=TRUE,repos='http://cran.us.r-project.org',INSTALL_opts='--no-lock',lib="/your/R/Library")
p_load(dep,character.only=TRUE,lib="/your/R/Library",install=FALSE,update=FALSE)
library(Azimuth,lib="/your/R/Library")

option_list = list(
  make_option(c("--wd"), type="character", default=NULL,help="Input folder, where all the cellranger files are being stored.", metavar="character"),
  make_option(c("--samples"),type="character", default=NULL,help="Location of per sample out folders from 10x AGGR/Count pipeline. Contains count matrices.", metavar="character")
  )
p = parse_args(OptionParser(option_list=option_list))

source("/your/R/Library/sc_functions.r") # Source the sc_functions.r file, where it is located on your system.
base=paste0(as.character(p[2]),"/outs/per_sample_outs/") 
parent=as.character(p[1])
S_Combined.INT = readRDS(paste0(p[1],"/S_Combined.INT.stage2.rds"))
S_Combined.harm = readRDS(paste0(p[1],"/S_Combined.harm.stage2.rds"))
S_Combined.SCT = readRDS(paste0(p[1],"/S_Combined.SCT.stage2.rds"))
DefaultAssay(S_Combined.INT) = "integrated"
DefaultAssay(S_Combined.harm) = "RNA"
DefaultAssay(S_Combined.SCT) = "SCT"

S_Combined.SCT = aziPredict(dataset=S_Combined.SCT,
        reference="/your/reference/location",
        assay="SCT",hs_mm="mm",outfile="SCT_AP_SCORES.csv")
S_Combined.harm = aziPredict(dataset=S_Combined.harm,
        reference="/your/reference/location",
        assay="RNA",hs_mm="mm",outfile="Harmony_AP_SCORES.csv")

saveRDS(S_Combined.SCT,"S_Combined.SCTP.rds")
saveRDS(S_Combined.harm,"S_Combined.harmP.rds")

#COLORS#
sample.cols = paletteer_d("ggsci::category20_d3")[1:length(list.files(base))]
heatmap.cols <- circlize::colorRamp2(c(0, 1), c("white", "black"))
heatmap_ann.cols <- circlize::colorRamp2(c(0, 0.1, 1), c("#EDD9A3", "#EA4F88", "#4B2991"))
cluster.cols <- pals::polychrome(n = 30) #Supports up to 36 cluster colors
names(cluster.cols) = as.character(0:29)
c(cluster.cols["1"],cluster.cols["4"],cluster.cols["8"],cluster.cols["9"]) %<-% c("#CDC1C5","#43CD80","#97FFFF","#9AFF9A")
cell_cycle.cols = c(`G1`="#003049", `G2M`="#d62828", `S`="#f77f00")

sig.cols <- c(`FDR`="#FF00CC", `Nominal`="#ffb3fd", `NS`="#C8C8CD")
heatmap_data.cols = circlize::colorRamp2(breaks=c(0:5),hcl_palette="Inferno")
heatmap_scale.cols = circlize::colorRamp2(breaks=c(-3.5,0,3.5),hcl_palette = "Lisbon")
heatmap_marker.cols = circlize::colorRamp2(breaks=c(0:5), hcl_palette="Inferno")
#END COLORS#

#UMAP by CMO
p1 = DimPlot(S_Combined.INT, reduction = "umap", group.by = "orig.ident", cols = alpha(sample.cols, 0.7)) + ggtitle("Integrated SCT")
p2 = DimPlot(S_Combined.harm, reduction = "umap", group.by = "orig.ident", cols = alpha(sample.cols, 0.7)) + ggtitle("logNorm + PCA + Harmony")
p3 = DimPlot(S_Combined.SCT, reduction = "umap", group.by = "orig.ident", cols = alpha(sample.cols, 0.7)) + ggtitle("SCT + CCA")
g1 = p2 + p3 + p1 + patchwork::plot_layout(ncol = 3, guides="collect")
ggsave(plot=g1,file="UMAP_CMOs.pdf",width=15,height=11)

#PCA by CMO
p1.1 = DimPlot(S_Combined.INT, reduction = "pca", group.by = "orig.ident", cols = alpha(sample.cols, 0.7)) + ggtitle("Integrated SCT")
p2.1 = DimPlot(S_Combined.harm, reduction = "pca", group.by = "orig.ident", cols = alpha(sample.cols, 0.7)) + ggtitle("logNorm + PCA + Harmony")
p2.2 = DimPlot(S_Combined.SCT, reduction = "pca", group.by = "orig.ident", cols = alpha(sample.cols, 0.7)) + ggtitle("SCT + CCA")
g1 = p2.1 +  p2.2 + p1.1 + patchwork::plot_layout(ncol = 3, guides="collect")
ggsave(filename="PCA_CMOs.pdf",plot=g1,width=15,height=11)

#CMO Barplot
tab <- as.data.frame(table(S_Combined.SCT@meta.data$SCT_snn_res.0.2, S_Combined.INT@meta.data$orig.ident))
p1.2 <- ggplot(tab, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual("", values = sample.cols) + xlab("Clusters") + ylab("Proportion of cells") + geom_hline(yintercept = 0) + ggtitle("Merged + log-Norm")
tab2 <- as.data.frame(table(S_Combined.harm@meta.data$RNA_snn_res.0.2, S_Combined.harm@meta.data$orig.ident))
p2.2 <- ggplot(tab2, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual("", values = sample.cols) + xlab("Clusters") + ylab("Proportion of cells") + geom_hline(yintercept = 0) + ggtitle("log-Norm + Harmony")
tab3 <- as.data.frame(table(S_Combined.INT@meta.data$integrated_snn_res.0.2, S_Combined.INT@meta.data$orig.ident))
p3.2 <- ggplot(tab3, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual("", values = sample.cols) + xlab("Clusters") + ylab("Proportion of cells") + geom_hline(yintercept = 0) + ggtitle("Integrated SCT")
g1 = p1.2 + p2.2 + p3.2 + patchwork::plot_layout(ncol = 3, guides="collect")
ggsave("CMO_Barplot.pdf",width=15,height=11,plot=g1)

#Heatmap Proportions by CMO
hm1 = hm_prop(S_Combined.INT)
hm2 = hm_prop(S_Combined.harm)
hm3 = hm_prop(S_Combined.SCT)
z = grid.arrange(hm1,hm2,hm3,ncol=1)
ggsave(plot=z,file="CMO_Cluster_Proportions.pdf",height=7,width=7)

#By PCA Dimensions
objs = c(S_Combined.INT,S_Combined.SCT,S_Combined.harm)
pdf("By_PC_Axis.pdf",width=12,height=12)
lapply(objs,function(x){
    qr = names(x@commands)
    if(length(qr[grepl("harmony",qr)])>=1){title_wow="Harmony"}
    if(length(qr[grepl("integrated",qr)])>=1){title_wow="Integrated"}
    if(length(qr[grepl("SCT",qr)])>=1){title_wow="SCT"}
    pc <- as.data.frame(x[["pca"]]@cell.embeddings)
    pc$orig.ident <- x@meta.data$orig.ident
    q = GGally::ggpairs(pc, columns = 1:4, upper = "blank",
                            mapping = ggplot2::aes(color = orig.ident, alpha = 0.5),
                            progress = FALSE) +
        scale_fill_manual("", values = sample.cols) + scale_color_manual("", values = sample.cols) + ggtitle(title_wow)
    q
})
dev.off()

#Louvain Clustering UMAP
a = cowplot::plot_grid(ncol = 2,
    DimPlot(S_Combined.harm, reduction = "umap", label = T, group.by = "RNA_snn_res.0.2") + ggtitle("Louvain Harm - res 0.2"),
    DimPlot(S_Combined.harm, reduction = "umap", label = T, group.by = "RNA_snn_res.0.5") + ggtitle("Louvain Harm - res 0.5")
)
b =cowplot::plot_grid(ncol = 2,
    DimPlot(S_Combined.SCT, reduction = "umap", label = T, group.by = "SCT_snn_res.0.2") + ggtitle("Louvain SCT - res 0.2"),
    DimPlot(S_Combined.SCT, reduction = "umap", label = T, group.by = "SCT_snn_res.0.5") + ggtitle("Louvain SCT - res 0.5")
)
c =cowplot::plot_grid(ncol = 2,
    DimPlot(S_Combined.INT, reduction = "umap", label = T, group.by = "integrated_snn_res.0.2") + ggtitle("Louvain Integrated - res 0.2"),
    DimPlot(S_Combined.INT, reduction = "umap", label = T, group.by = "integrated_snn_res.0.5") + ggtitle("Louvain Integrated - res 0.5")
)
z = a+b+c+plot_layout(ncol=1)
ggsave(plot=z,file="Louvain_UMAP_Clustering.pdf",width=25,height=30)

#Louvain Heatmap
Louvain_Annotation_HM(S_Combined.INT,"integrated_snn_res.0.2","Annot_Cluster_Prop_INT.pdf","mm")
Louvain_Annotation_HM(S_Combined.harm,"RNA_snn_res.0.2","Annot_Cluster_Prop_Harm.pdf","mm")
Louvain_Annotation_HM(S_Combined.SCT,"SCT_snn_res.0.2","Annot_Cluster_Prop_SCT.pdf","mm")

#UMAP Annotations
make_UMAP_Annotation(S_Combined.harm,"RNA_snn_res.0.2","UMAP_ANNOT_HARM.pdf","mm")
make_UMAP_Annotation(S_Combined.INT,"integrated_snn_res.0.2","UMAP_ANNOT_INT.pdf","mm")
make_UMAP_Annotation(S_Combined.SCT,"SCT_snn_res.0.2","UMAP_ANNOT_SCT.pdf","mm")

#Azimuth annotations
make_UMAP_Azimuth(S_Combined.harm,"AziAnnot_harmony.pdf")
make_UMAP_Azimuth(S_Combined.SCT,"AziAnnot_SCT.pdf")

find_plot_markers(S_Combined.harm,"RNA","HARM_LFC","HARM_ANNOT_HM","Dotplot_Features_Harm","mm")
S_Combined.SCT = PrepSCTFindMarkers(S_Combined.SCT)
find_plot_markers(S_Combined.SCT,"SCT","SCT_LFC","SCT_ANNOT_HM","Dotplot_Features_SCT","mm")
#find_plot_markers(S_Combined.INT,"integrated","INT_LFC","INT_ANNOT_HM","Dotplot_Features_INT") not available for "integrated" assay data

feature_plots(S_Combined.harm,"Annotation_Cluster_Barplot_Harm","Ridgeplot_Features_Harm","RNA","mm")
feature_plots(S_Combined.SCT,"Annotation_Cluster_Barplot_SCT","Ridgeplot_Features_SCT","SCT","mm")
feature_plots(S_Combined.INT,"Annotation_Cluster_Barplot_INT","Ridgeplot_Features_INT","integrated","mm")

tx = list("Mock"=c("CMO_1","CMO_2","CMO_3","CMO_4"),"Infected"=c("CMO_5","CMO_6","CMO_7","CMO_8"))
subcluster(S_Combined.harm,tx=tx,integ="SubCluster",parent,"mm")
subcluster(S_Combined.SCT,tx=tx,integ="SubCluster",parent,"mm")