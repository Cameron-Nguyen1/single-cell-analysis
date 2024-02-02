#!/bin/env Rscript
#Harmony workflow
library(pacman,lib="")
dep = c("plyr","data.table","devtools","Seurat","tidyverse","miQC","SeuratWrappers","flexmix","SingleCellExperiment","SummarizedExperiment","RColorBrewer","ggsankey",
"ggplot2","cowplot","SingleR","scran","celldex","ComplexHeatmap","pheatmap","circlize","GGally","forcats","dplyr","patchwork","pals","harmony",
"ggpubr","paletteer","ggridges","fgsea","UCell","FactoMineR","factoextra","zeallot","gridExtra","grid","optparse")
#lapply(dep,install.packages,character.only=TRUE,repos='http://cran.us.r-project.org',INSTALL_opts='--no-lock',lib="")
p_load(dep,character.only=TRUE,lib="",install=FALSE,update=FALSE)

option_list = list(
  make_option(c("--wd"), type="character", default=NULL,help="Input folder, where all the cellranger files are being stored.", metavar="character")
 # make_option(c("--output"), type="character", default=NULL,help="Decides what folder is the root for output.", metavar="character")
  )
p = parse_args(OptionParser(option_list=option_list))

source("/sc_functions.r")
S_Combined = readRDS(paste0(getwd(),"/S_Combined_postQC_stage1.rds"))
S_Combined = NormalizeData(S_Combined, normalization.method = "LogNormalize", scale.factor = 10000)
S_Combined = FindVariableFeatures(S_Combined, selection.method = "vst", nfeatures = 4000)
S_Combined = ScaleData(S_Combined)
S_Combined=RunPCA(S_Combined,npcs=50) #Test 50 dimensions
S_Combined.harm = RunHarmony(S_Combined, group.by.vars = "orig.ident", plot_convergence = FALSE) #this is probably what makes "Rplots.pdf"
bestdims_c = bdim(S_Combined)
S_Combined.harm = RunUMAP(S_Combined.harm, reduction = "harmony", dims = 1:bdim(S_Combined.harm))
S_Combined.harm = FindNeighbors(S_Combined.harm, reduction = "harmony", dims = 1:bdim(S_Combined.harm))
S_Combined.harm = FindClusters(S_Combined.harm,algorithm = 1, resolution = c(0.2, 0.5))

c(S_Combined.harm,immgen.ref,mrsd.ref) %<-% annotate_human(S_Combined.harm,"RNA")[c(1,2,3)]
S_Combined.harm =  clean(S_Combined.harm)

#c(preds.harm.immgen,preds.harm.mrsd) %<-% c(marker_predict(S_Combined.harm,immgen.ref,"RNA")[1],marker_predict(S_Combined.harm,mrsd.ref,"RNA")[1])
#saveRDS(preds.harm.immgen,"preds.harm.immgen.rds")
#saveRDS(preds.harm.mrsd,"preds.harm.mrsd.rds")

saveRDS(S_Combined.harm,"S_Combined.harm.stage2.rds")