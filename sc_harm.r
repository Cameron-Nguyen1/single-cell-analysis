#!/bin/env Rscript
#Harmony workflow
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
p[1] = tolower(as.character(p[1]))

source("/your/R/Library/sc_functions.r") # Source the sc_functions.r file, where it is located on your system.
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

if (p[2] == "hs"){
  c(S_Combined.harm,immgen.ref,mrsd.ref) %<-% annotate_human(S_Combined.harm,"RNA")[c(1,2,3)]
}else{
  c(S_Combined.harm,immgen.ref,mrsd.ref) %<-% annotate_mouse(S_Combined.harm,"RNA")[c(1,2,3)]
}
S_Combined.harm =  clean(S_Combined.harm)

saveRDS(S_Combined.harm,"S_Combined.harm.stage2.rds")