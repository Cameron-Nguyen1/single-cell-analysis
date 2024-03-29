#!/bin/env Rscript
#DECLARE FUNCTIONS / LOAD LIBRARIES#
library(pacman,lib="/your/R/Library")
dep = c("plyr","data.table","devtools","Seurat","tidyverse","miQC","SeuratWrappers","flexmix","SingleCellExperiment","SummarizedExperiment","RColorBrewer","ggsankey",
"ggplot2","cowplot","SingleR","scran","celldex","ComplexHeatmap","pheatmap","circlize","GGally","forcats","dplyr","patchwork","pals","harmony",
"ggpubr","paletteer","ggridges","fgsea","UCell","FactoMineR","factoextra","zeallot","gridExtra","grid","optparse")
#lapply(dep,install.packages,character.only=TRUE,repos='http://cran.us.r-project.org',INSTALL_opts='--no-lock',lib="/your/R/Library")
p_load(dep,character.only=TRUE,lib="/your/R/Library",install=FALSE,update=FALSE)

option_list = list(
  make_option(c("--hs_or_mm"),type="character", default=NULL,help="Either HS or MM for homo sapiens and mus musculus respectively.", metavar="character")
)
p = parse_args(OptionParser(option_list=option_list))

source("/your/R/Library/sc_functions.r") # Source the sc_functions.r file, where it is located on your system.

#BEGIN SCT#
S_Combined = readRDS(paste0(getwd(),"/S_Combined_postQC_stage1.rds"))

myl = SplitObject(S_Combined,split.by="orig.ident")
myl = lapply(myl,SCTransform,method='glmGamPoi',vst.flavor='v2',vars.to.regress=c('percent.mt','S.Score','G2M.Score'),return.only.var.genes=FALSE,verbose=FALSE)
features = SelectIntegrationFeatures(object.list = myl, nfeatures = 4000)
myl = PrepSCTIntegration(object.list = myl, anchor.features = features)
anchors = FindIntegrationAnchors(object.list = myl, normalization.method = "SCT", anchor.features = features)
S_Combined.SCT = IntegrateData(anchorset = anchors, normalization.method = "SCT")  ###This can bug out if cells per sample is < 100. Apply k.weight = min(sample)
S_Combined.INT = S_Combined.SCT
DefaultAssay(S_Combined.SCT) = "SCT"
DefaultAssay(S_Combined.INT) = "integrated"

VariableFeatures(S_Combined.SCT[["SCT"]]) = rownames(S_Combined.SCT[["SCT"]]@scale.data)
S_Combined.SCT= RunPCA(S_Combined.SCT,npcs=50)
bestdims_sct = bdim(S_Combined.SCT)
S_Combined.SCT = RunUMAP(S_Combined.SCT, dims = 1:bestdims_sct, verbose = FALSE)
S_Combined.SCT = Seurat::FindNeighbors(S_Combined.SCT,reduction="pca",dims=1:bestdims_sct)
S_Combined.SCT = Seurat::FindClusters(S_Combined.SCT, algorithm = 1, resolution = c(0.2, 0.5))

S_Combined.INT= RunPCA(S_Combined.INT,npcs=50)
bestdims_int = bdim(S_Combined.INT)
S_Combined.INT = RunUMAP(S_Combined.INT, dims = 1:bestdims_int, verbose = FALSE)
S_Combined.INT = Seurat::FindNeighbors(S_Combined.INT,reduction="pca",dims=1:bestdims_int)
S_Combined.INT = Seurat::FindClusters(S_Combined.INT, algorithm = 1, resolution = c(0.2, 0.5))

if (p[2] == "hs"){
  c(S_Combined.INT,immgen.ref,mrsd.ref) %<-% annotate_human(S_Combined.INT,"integrated")[c(1,2,3)]
  c(S_Combined.SCT,immgen.ref,mrsd.ref) %<-% annotate_human(S_Combined.SCT,"SCT")[c(1,2,3)]
}else{
  c(S_Combined.INT,immgen.ref,mrsd.ref) %<-% annotate_mouse(S_Combined.INT,"integrated")[c(1,2,3)]
  c(S_Combined.SCT,immgen.ref,mrsd.ref) %<-% annotate_mouse(S_Combined.SCT,"SCT")[c(1,2,3)]
}

S_Combined.SCT = clean(S_Combined.SCT)
S_Combined.INT = clean(S_Combined.INT)

saveRDS(S_Combined.SCT,"S_Combined.SCT.stage2.rds")
saveRDS(S_Combined.INT,"S_Combined.INT.stage2.rds")