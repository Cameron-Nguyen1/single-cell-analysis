#!/bin/env Rscript
library(pacman,lib="/your/R/Library")
dep = c("plyr","data.table","devtools","Seurat","tidyverse","miQC","SeuratWrappers","flexmix","SingleCellExperiment","SummarizedExperiment","RColorBrewer","ggsankey",
"ggplot2","cowplot","SingleR","scran","celldex","ComplexHeatmap","pheatmap","circlize","GGally","forcats","dplyr","patchwork","pals","harmony",
"ggpubr","paletteer","ggridges","fgsea","UCell","FactoMineR","factoextra","zeallot","gridExtra","grid","optparse")
#lapply(dep,install.packages,character.only=TRUE,repos='http://cran.us.r-project.org',INSTALL_opts='--no-lock',lib="/your/R/Library")
p_load(dep,character.only=TRUE,lib="/your/R/Library",install=FALSE,update=FALSE)
option_list = list(
  make_option(c("--hs_or_mm"),type="character", default=NULL,help="Either HS or MM for homo sapiens and mus musculus respectively.", metavar="character"),
  make_option(c("--samples"),type="character", default=NULL,help="Location of per sample out folders from 10x AGGR/Mutli/Count pipeline. Contains count matrices.", metavar="character")
  )
p = parse_args(OptionParser(option_list=option_list))

source("/your/R/Library/sc_functions.r") # Source the sc_functions.r file, where it is located on your system.
p3t = tolower(as.character(p[1]))
base = as.character(p[2])




#base=paste0(as.character(p[2]),"/")
sample.cols = paletteer_d("ggsci::category20_d3")[1:length(list.files(base))]
for (sample in list.files(base)){ #Generate Seurat objects and assign percent.mt
  print(paste0("S",sample,"=ReadMtx(cells='",base,as.character(sample),"/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz',features='",base,as.character(sample),"/count/sample_filtered_feature_bc_matrix/features.tsv.gz',mtx='",base,as.character(sample),"/count/sample_filtered_feature_bc_matrix/matrix.mtx.gz'"))
  eval(parse(text=paste0("S",sample,"=ReadMtx(cells='",base,as.character(sample),"/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz',features='",base,as.character(sample),"/count/sample_filtered_feature_bc_matrix/features.tsv.gz',mtx='",base,as.character(sample),"/count/sample_filtered_feature_bc_matrix/matrix.mtx.gz')")))
  eval(parse(text=paste0("S",sample,"=rem_cmo(S",sample,")")))
  eval(parse(text=paste0("S",sample,"_Seurat=CreateSeuratObject(S",sample,",project='CMO_",sample,"',assay='RNA')")))
  eval(parse(text=paste0("S",sample,"_Seurat=PercentageFeatureSet(S",sample,"_Seurat,pattern = '^mt-', col.name='percent.mt')")))
}
### Add a function here to find variable features, good to know before we normalize everything ###

store = objects()[objects() %like% "Seurat$"]  ###Works for now, fix later, cell.ids (little tags at beginning of each sequence) aren't being set. Orig.ID is handled.
iter = 0
for (obj in store){
  iter = iter+1
  if (iter > (length(store))){
    break
  }
  if (iter >= 3){
    S_Combined = merge(S_Combined,y=get(store[iter]),project="REPC")
  }
  if (iter == 1){
    S_Combined = merge(get(store[1]),y=get(store[2]),project="REPC")
  }
}
DefaultAssay(S_Combined) = "RNA"
S_Combined = PercentageFeatureSet(S_Combined, pattern = "^mt-|^MT-", col.name = "percent.mt")  # Mitochondria % for merged set
pdf("Merged_preQC_Violin.pdf",width=15)
VlnPlot(object = S_Combined, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), cols=sample.cols, pt.size=0, ncol = 3) & geom_jitter(alpha = 0.25, size = 0.1) #preQC stats
dev.off()

nUMI.low = 500                   #QC SECTION
nUMI.high = 40000
percent_mito.low = -Inf
percent_mito.high = 10
nGene.low = 250
nGene.high = 5000

S_Combined@meta.data$Filter = ifelse(S_Combined@meta.data$nFeature_RNA > nGene.high | S_Combined@meta.data$nFeature_RNA < nGene.low, "Remove", "Keep")
S_Combined@meta.data$Filter = ifelse(S_Combined@meta.data$percent.mt > percent_mito.high, "Remove", S_Combined@meta.data$Filter)
S_Combined@meta.data$Filter = ifelse(S_Combined@meta.data$nCount_RNA > nUMI.high | S_Combined@meta.data$nCount_RNA < nUMI.low, "Remove", S_Combined@meta.data$Filter)

filter.cols = c(`Keep`="grey25", `Remove`="firebrick")
S_Combined@meta.data$orig.ident = factor(S_Combined@meta.data$orig.ident)
S_Combined@meta.data$orig.ident = factor(S_Combined@meta.data$orig.ident)
p1 = ggplot(data = S_Combined@meta.data,
       aes(x = orig.ident, y = nCount_RNA, fill = orig.ident)) +
  geom_jitter(aes(color = Filter), size = 0.5, position = position_jitter(0.1)) + geom_violin(aes(alpha = 0.1)) +
  scale_fill_manual("", values = sample.cols) + scale_color_manual("", values = filter.cols) +
  xlab("Identity") + theme(legend.position = "none")
p2 = ggplot(data = S_Combined@meta.data,
       aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident)) +
  geom_jitter(aes(color = Filter), size = 0.5, position = position_jitter(0.1)) + geom_violin(aes(alpha = 0.1)) +
  scale_fill_manual("", values = sample.cols) + scale_color_manual("", values = filter.cols) +
  xlab("Identity") + theme(legend.position = "none")
p3 = ggplot(data = S_Combined@meta.data, aes(x = orig.ident, y = percent.mt, fill = orig.ident)) +
  geom_jitter(aes(color = Filter), size = 0.5, position = position_jitter(0.1)) + geom_violin(aes(alpha = 0.1)) +
  scale_fill_manual("", values = sample.cols) + scale_color_manual("", values = filter.cols) +
  xlab("Identity") + theme(legend.position = "none")
g1 = cowplot::plot_grid(p1, p2, p3, ncol = 3)
p1 = ggplot(data = S_Combined@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = Filter)) + geom_point() +
  geom_hline(yintercept = c(nGene.high, nGene.low)) + geom_vline(xintercept = c(nUMI.high, nUMI.low)) +
  facet_grid(.~orig.ident) + scale_color_manual(values = filter.cols) + theme(legend.position = "none")
p2 = ggplot(data = S_Combined@meta.data, aes(x = percent.mt, y = log10(nFeature_RNA), color = Filter)) + geom_point() +
  geom_vline(xintercept = c(percent_mito.high)) + geom_hline(yintercept = c(log10(nGene.high), log10(nGene.low))) +
  facet_grid(.~orig.ident) + scale_color_manual(values = filter.cols) + theme(legend.position ="none")
p3 = ggplot(data = S_Combined@meta.data, aes(x = percent.mt, y = log10(nCount_RNA), color = Filter)) + geom_point() +
  geom_vline(xintercept = c(percent_mito.high)) + geom_hline(yintercept = c(log10(nUMI.high), log10(nUMI.low))) +
  facet_grid(.~orig.ident) + scale_color_manual(values = filter.cols) + theme(legend.position = "none")
g2 = cowplot::plot_grid(p1, p2, p3, ncol = 1)
ggsave(filename="Violin_Features.pdf",plot=g1,height=10,width=18)
ggsave(filename="Scatter_Features.pdf",plot=g2,height=10,width=18)

S_Combined = subset(x = S_Combined, subset = nFeature_RNA > nGene.low & nFeature_RNA < nGene.high & nCount_RNA > nUMI.low & nCount_RNA < nUMI.high & percent.mt < percent_mito.high)
pdf("Merged_postQC_Violin.pdf",width=15) #POSTQC VIOLIN
VlnPlot(object = S_Combined, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), cols=sample.cols, pt.size=0, ncol = 3) & geom_jitter(alpha = 0.25, size = 0.1) #preQC stats
dev.off()

if (p3t == "mm"){
  h2m.df = read.table(paste0(getwd(),"/HMD_HumanPhenotype.rpt.txt"), sep = "\t") #Cell Cycle Scoring
  h2m.df = h2m.df[, c("V1","V3")]
  colnames(h2m.df) = c("Gene.name", "ortholog_name")
  h2m.df = h2m.df[!duplicated(h2m.df$Gene.name), ]
  h2m.df = h2m.df[!duplicated(h2m.df$ortholog_name), ]
  mmus_s = h2m.df[which(h2m.df$Gene.name%in%cc.genes$s.genes), ]$ortholog_name
  mmus_g2m = h2m.df[which(h2m.df$Gene.name%in%cc.genes$g2m.genes), ]$ortholog_name
  S_Combined = CellCycleScoring(S_Combined, s.features = mmus_s, g2m.features = mmus_g2m, set.ident = FALSE) #Species specific
}else{
  h2m.df = read.table(paste0(getwd(),"/HMD_HumanPhenotype.rpt.txt"), sep = "\t") #Cell Cycle Scoring
  h2m.df = h2m.df[,"V1"]
  h2m.df = h2m.df[!duplicated(h2m.df)]
  #h2m.df = h2m.df[!duplicated(h2m.df$ortholog_name), ]
  hs_s = h2m.df[which(h2m.df%in%cc.genes$s.genes)]
  hs_g2m = h2m.df[which(h2m.df%in%cc.genes$g2m.genes)]
  S_Combined = CellCycleScoring(S_Combined, s.features = hs_s, g2m.features = hs_g2m, set.ident = FALSE,nbin=23) #Species specific
}

saveRDS(S_Combined,"S_Combined_postQC_stage1.rds")