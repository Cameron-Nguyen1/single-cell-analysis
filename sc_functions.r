rem_cmo = function(x){
    x = x[!grepl("CMO",rownames(x)),]
    return(x)
}
generate_annot_cols = function(x){
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vec = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  vec = colnames(x@meta.data)[colnames(x@meta.data) %like% "pruned" | colnames(x@meta.data) %like% "^predicted.*[0-9]$"]
  myl = list()
  for (var in vec){
    ctypes = unique(x[[var]][,1])
    ctypes.cols = col_vec[1:length(ctypes)]
    ctype_vector = unlist(as.list(setNames(ctypes.cols, ctypes)))
    myl[[var]] = ctype_vector
  }
  return(myl)
}
aziPredict = function(dataset,assay,reference,hs_mm,outfile){
    options(Azimuth.map.ndims = 50)
    if (hs_mm == "hs"){
      levels = c("ann_level_3","ann_level_4")
    }else{
      levels =  c("celltype_level1","celltype_level2","celltype_level3")
    }
    dataset =  RunAzimuth(query=dataset,reference=reference,assay=assay,annotation.levels=levels)
    mzl = list("Map_Score.50"=table(dataset$mapping.score > .50),"Map_Score.75"=table(dataset$mapping.score > .75))
    ids = paste0("predicted.",levels,".score")
    newframe = dataset[[ids]]
    for (i in 1:length(names(newframe))){
        mz_50 = list(name50=table(newframe[[i]] >= .50))
        mz_75 = list(name75=table(newframe[[i]] >= .75))
        mzl = c(mzl, mz_50, mz_75)
    }
    lnum = 1
    for (i in 3:length(names(mzl))){
        if (i %% 2 != 0 && i != 3){
            lnum = lnum + 1
        }
        if (i %% 2 == 0){
            names(mzl)[i] = paste0("L",lnum,"_Pred.Score_75")
        }else{
            names(mzl)[i] = paste0("L",lnum,"_Pred.Score_50")
        }
    }
    df_pred = data.frame(mzl)
    df_pred = data.frame(t(df_pred[,!grepl("Var1",colnames(df_pred))]))
    df_pred = cbind(df_pred,list(pct.True=df_pred[,2]/(df_pred[,2]+df_pred[,1])))
    colnames(df_pred) = c("Below.Scr","Above.Scr","Pct.Above")
    write.csv(x=df_pred,file=outfile)
    return(dataset)
}
annotate_mouse = function(x,assay){
  immgen.ref = celldex::ImmGenData(ensembl=FALSE)
  mrsd.ref = celldex::MouseRNAseqData(ensembl=FALSE)

  # perform predictions
  immgen.preds = SingleR::SingleR(test = x@assays[[assay]]@data, ref = immgen.ref, labels = immgen.ref$label.main, fine.tune = TRUE,
                  de.method = "classic", de.n = 10, clusters = NULL, assay.type.test = "logcounts", assay.type.ref = "logcounts")
  mrsd.preds = SingleR::SingleR(test = x@assays[[assay]]@data, ref = mrsd.ref, labels = mrsd.ref$label.main, fine.tune = TRUE,
                  de.method = "classic", de.n = 10, clusters = NULL, assay.type.test = "logcounts", assay.type.ref = "logcounts")
  myl2 = list()
  myl3 = list()
  for (obj in objects()[grepl("preds",objects())]){
    namae = str_split_1(obj,"\\.")
    myl2[[obj]] = SingleR::plotScoreHeatmap(get(obj),show.pruned=TRUE,main=paste(namae[1],namae[2]))[[4]]
    myl3[[obj]] = SingleR::plotDeltaDistribution(get(obj)) + ggtitle(paste(namae[1],namae[2]))
  }
  matr = matrix(c(1,3,5,2,4,6),nrow=3,ncol=2)
  g=grid.arrange(arrangeGrob(grobs= myl2,layout_matrix=matr))
  g2=grid.arrange(arrangeGrob(grobs= myl3,ncol=2,layout_matrix=matr))
  ggsave(paste0(assay,"_HM_CONF.pdf"),plot=g,width=16,height=15,device='pdf')
  ggsave(paste0(assay,"_Delta_CONF.pdf"),plot=g2,width=16,height=15,device='pdf')
  x@meta.data$immgen.pred <- immgen.preds$labels
  x@meta.data$immgen.pred_pruned <- immgen.preds$pruned.labels # contains NA
  x@meta.data$mrsd.pred <- mrsd.preds$labels
  x@meta.data$mrsd.pred_pruned <- mrsd.preds$pruned.labels # contains NA
  return(c(x,immgen.ref,mrsd.ref))
}

annotate_human = function(x,assay){
  HPCA.ref = celldex::HumanPrimaryCellAtlasData(ensembl=FALSE)
  DICE.ref = celldex::DatabaseImmuneCellExpressionData(ensembl=FALSE)

  # perform predictions
  DICE.preds = SingleR::SingleR(test = x@assays[[assay]]@data, ref = DICE.ref, labels = DICE.ref$label.main, fine.tune = TRUE,
                  de.method = "classic", de.n = 10, clusters = NULL, assay.type.test = "logcounts", assay.type.ref = "logcounts")
  HPCA.preds = SingleR::SingleR(test = x@assays[[assay]]@data, ref = HPCA.ref, labels = HPCA.ref$label.main, fine.tune = TRUE,
                  de.method = "classic", de.n = 10, clusters = NULL, assay.type.test = "logcounts", assay.type.ref = "logcounts")
  myl2 = list()
  myl3 = list()
  for (obj in objects()[grepl("preds",objects())]){
    namae = str_split_1(obj,"\\.")
    myl2[[obj]] = SingleR::plotScoreHeatmap(get(obj),show.pruned=TRUE,main=paste(namae[1],namae[2]))[[4]]
    myl3[[obj]] = SingleR::plotDeltaDistribution(get(obj)) + ggtitle(paste(namae[1],namae[2]))
  }
  matr = matrix(c(1,3,5,2,4,6),nrow=3,ncol=2)
  g=grid.arrange(arrangeGrob(grobs= myl2,layout_matrix=matr))
  g2=grid.arrange(arrangeGrob(grobs= myl3,ncol=2,layout_matrix=matr))
  ggsave(paste0(assay,"_HM_CONF.pdf"),plot=g,width=16,height=15,device='pdf')
  ggsave(paste0(assay,"_Delta_CONF.pdf"),plot=g2,width=16,height=15,device='pdf')
  x@meta.data$DICE.pred <- DICE.preds$labels
  x@meta.data$DICE.pred_pruned <- DICE.preds$pruned.labels # contains NA
  x@meta.data$HPCA.pred <- HPCA.preds$labels
  x@meta.data$HPCA.pred_pruned <- HPCA.preds$pruned.labels # contains NA
  return(c(x,DICE.ref,HPCA.ref))
}

clean = function(x){
    for (i in colnames(x[[]])[colnames(x[[]]) %like% "pred_pruned"]){
        ref = eval(parse(text=paste0("x@meta.data[i]$",colnames(x@meta.data[i])))) #Finds prediction column in seurat object
        ref[is.na(ref)] = "Others"
       # x@meta.data[i][,1] = forcats::fct_collapse(ref,`B cells` = "B cells, pro")     #Uses assay label to clean individual assignments
        freq <- as.data.frame(table(x@meta.data[i]))
        freq <- freq %>% mutate("prop" = prop.table(Freq))
        colnames(freq) <- c("cells", "n", "proportion")
        freq$percent <- freq$proportion * 100
        freq <- freq[order(freq$n, decreasing = TRUE), ]
        others <- c(as.vector(freq$cells[which(freq$n < 50)]),NA)
        x@meta.data[i][,1] <- forcats::fct_collapse(ref, Others = others,`B cells` = "B cells, pro")
    }
    return(x)
}
marker_predict = function(x,reference,assay){ ##Probably redundant function
#Marker Gene Detection
#Classic approach / Cell Level / Bulk RNAseq Reference
classic.pred <- SingleR::SingleR(test = x@assays[[assay]]@data,ref = reference, labels = reference$label.main,fine.tune = TRUE,
                  de.method = "classic", de.n = 10, clusters = NULL, assay.type.test = "logcounts", assay.type.ref = "logcounts")

#Pairwise Wilcoxon ranked sum test / Cell Level / Single-cell Reference
wilcox.pred <- SingleR::SingleR(test = x@assays[[assay]]@data, ref = reference,
                labels = reference$label.main, fine.tune = TRUE, de.method = "wilcox",
                de.n = 10, clusters = NULL, assay.type.test = "logcounts", assay.type.ref = "logcounts")
#Cluster Level / Not recommended
#cluster.pred <- SingleR::SingleR(test = x@assays[["RNA"]]@data,
#                ref = get(reference), labels = get(reference)$label.main,
#                fine.tune = TRUE, de.method = "classic",  de.n = 10, clusters = x@meta.data$integrated_snn_res.0.2,
#                assay.type.test = "logcounts", assay.type.ref = "logcounts")
return(list(classic.pred,wilcox.pred))
}
make_UMAP_Annotation = function(x,full_assay,pdf_name,hs_or_mm){
  if (hs_or_mm == "mm"){gb=list("mrsd.pred_pruned","immgen.pred_pruned")
  }else{gb=list("HPCA.pred_pruned","DICE.pred_pruned")}
  annot_colors = generate_annot_cols(x)
  p1 <- DimPlot(x, reduction = "umap", label=T, group.by = full_assay, cols = cluster.cols) + ggtitle(paste0("Louvain clusters - resolution ",substr(full_assay,nchar(full_assay)-3,nchar(full_assay))))
  p3 <- DimPlot(x, reduction = "umap", label=T, group.by = gb[[1]], cols = annot_colors[unlist(gb[1])][[1]]) + ggtitle(gb[1])
  p4 <- DimPlot(x, reduction = "umap", label=T, group.by = gb[[2]], cols = annot_colors[unlist(gb[2])][[1]]) + ggtitle(gb[2])
  g1 = p1 + p3 + p4 + patchwork::plot_layout(ncol = 2, nrow = 2)
  ggsave(pdf_name,width=20,height=24,plot=g1)
}

make_UMAP_Azimuth = function(x,pdf_name){
  annot_colors = generate_annot_cols(x)
  ref = colnames(x@meta.data)
  plst = list()
  for (i in ref[grepl("^predicted",ref) & !grepl("score",ref)]){
      plst[[i]] = DimPlot(x,reduction="umap", label=T, repel=T, label.size=2.5, group.by=i, cols = annot_colors[[i]]) + ggtitle(paste0("Azimuth "),i)
  }
  g=marrangeGrob(grobs=plst,ncol=2,nrow=2)
  ggsave(paste0(pdf_name,".pdf"),plot=g,width=24,height=16,device="pdf", limitsize = FALSE)
}
hm_prop = function(x){
    strang = deparse(substitute(x))
    if (strang %like% "harm"){
        c_title = "log-Norm + Harmony"
        tab <- prop.table(table(x@meta.data$orig.ident, x@meta.data$RNA_snn_res.0.2), margin = 2)
        nc <- length(levels(as.factor(x@meta.data$RNA_snn_res.0.2)))}
    else if(strang %like% "SCT"){
        c_title = "SCT + CCA"
        tab <- prop.table(table(x@meta.data$orig.ident, x@meta.data$SCT_snn_res.0.2), margin = 2)
        nc <- length(levels(as.factor(x@meta.data$SCT_snn_res.0.2)))}
    else{
        c_title = "Integrated SCT"
        tab <- prop.table(table(x@meta.data$orig.ident, x@meta.data$integrated_snn_res.0.2), margin = 2)
        nc <- length(levels(as.factor(x@meta.data$integrated_snn_res.0.2)))}
    nr <- length(levels(as.factor(x@meta.data$orig.ident)))
    hm1 <- ComplexHeatmap::Heatmap(tab, cluster_columns = FALSE,
                            show_column_dend = TRUE, column_names_rot = 90, row_names_gp = grid::gpar(fontsize = 12), column_names_gp = grid::gpar(fontsize=12),
                            cluster_rows = FALSE, show_row_dend = FALSE, col = heatmap.cols, border_gp = gpar(col = "black", lty = 1),
                            rect_gp = gpar(col = "white", lwd = 2), column_title = c_title, height = unit(5, "mm")*nr, width = unit(5, "mm")*nc, name="Proportion of cells")
    hm1 = grid.grabExpr(draw(hm1))
    return(hm1)
}
bdim = function(x){
  if (x[["pca"]]@assay.used == "integrated"){
    DefaultAssay(x) = "integrated"
    pct=x[['pca']]@stdev/sum(x[['pca']]@stdev)*100 #Make PCT var for PCA dim selection
    csum = cumsum(pct)
    c1 = which(pct < 5 & csum > 90)[1] #PCA dim must contribute <5 variance and > 90 cumulative variance
    c2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1 #Find dim with less than .1 difference in step
    if (is.na(c1) | is.na(c2)){
      z = c(c1,c2)
      bestdims = z[!is.na(z)]
    }else{
    bestdims = min(c1,c2) #We take the lowest dimension of the two as the range of dims that capture a majority of variance
    }
    return(bestdims)
  }
  pct=x[['pca']]@stdev/sum(x[['pca']]@stdev)*100 #Make PCT var for PCA dim selection
  csum = cumsum(pct)
  c1 = which(pct < 5 & csum > 90)[1] #PCA dim must contribute <5 variance and > 90 cumulative variance
  c2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1 #Find dim with less than .1 difference in step
  if (is.na(c1) | is.na(c2)){
    z = c(c1,c2)
    bestdims = z[!is.na(z)]
  }else{
  bestdims = min(c1,c2) #We take the lowest dimension of the two as the range of dims that capture a majority of variance
  }
  return(bestdims)
}

Louvain_Annotation_HM = function(x,assay_full,file_name,hs_or_mm){
  if (hs_or_mm=="mm"){
  mrsd_1 = prop.table(table(x@meta.data$mrsd.pred_pruned, eval(parse(text=(paste0("x@meta.data$",assay_full))))))
  immgen_2 = prop.table(table(x@meta.data$immgen.pred_pruned, eval(parse(text=(paste0("x@meta.data$",assay_full))))))
  aziann="predicted.celltype_level2"
  }
  else{
  DICE_1 = prop.table(table(x@meta.data$DICE.pred_pruned, eval(parse(text=(paste0("x@meta.data$",assay_full))))))
  HPCA_2 = prop.table(table(x@meta.data$HPCA.pred_pruned, eval(parse(text=(paste0("x@meta.data$",assay_full))))))
  aziann="predicted.ann_level_3"
  }
  if (!assay_full %like% "integrated"){
  azimuth_3 <- prop.table(table(x@meta.data[aziann][,1], eval(parse(text=(paste0("x@meta.data$",assay_full))))))
  }
  pdf(file_name)
  for (obj in objects()[objects() %like% "_[1-9]$"]){
    mhm = ComplexHeatmap::Heatmap(get(obj), cluster_columns = FALSE, show_column_dend = TRUE, column_names_rot = 90, row_names_gp = grid::gpar(fontsize = 12),
                        column_names_gp = grid::gpar(fontsize=12), cluster_rows = TRUE, show_row_dend = FALSE,  col = heatmap_ann.cols, border_gp = gpar(col = "black", lty = 1),
                        rect_gp = gpar(col = "white", lwd = 2), column_title = paste0("Louvian 0.2 clusters vs. ",substr(obj,1,nchar(obj)-2)),  name="Proportion of cells")
    plot(mhm)
  }
  dev.off()
}

theme_90 <- function() {
  theme_bw(base_size=12)+
    theme(axis.text.x=element_text(angle=90,hjust = 1,vjust = 0.5),axis.text=element_text(color="black"),panel.background=element_rect(color="black"),
          strip.text = element_text(size=12),strip.background = element_rect(fill="white"))
}
find_plot_markers = function(x,assay,lfc_pdf,annot_pdf,dotplot,hs_or_mm){
  annot_colors = generate_annot_cols(x)
  #x = NormalizeData(x,normalization.method="LogNormalize",scale=10000,assay="RNA") #"HARM_LFC","HARM_ANNOT_HM","Dotplot_Features_Harm")
  iter = 0
  res = colnames(x@meta.data[,colnames(x@meta.data) %like% "*res.[0-9].[0-9]"])
  for (i in res){
      iter = iter+1
      res_str = str_match(res[iter],".[0-9].[0-9]")[1]
      print(paste0("Working on resolution: ",res_str))
      Idents(x) = i
      markers = FindAllMarkers(x, assay=assay, only.pos = TRUE, slot="counts", min.pct = 0.3, logfc.threshold = .5, max.cells.per.ident=5000, random.seed = 8675309)
      markers = filter(markers,avg_log2FC!=Inf)
      #Organize results into table
      markers %>%
          group_by(cluster) %>%
          top_n(n = 5, wt = avg_log2FC) -> top5
      markers %>%
          group_by(cluster) %>%
          top_n(n = 10, wt = avg_log2FC) -> top10
      marker.list = list(AllMarkers=markers, Top5Markers=top5)
      marker.listn10 = list(AllMarkers=markers, Top5Markers=top10)
      write.csv(marker.list[1],file=paste0(assay,"_",res_str,"_clustersDEG.csv"))
      write.csv(marker.listn10[2],file=paste0(assay,"_",res_str,"_clustersDEGtop10.csv"))
      marker.dt = as.data.table(markers)
      smry = marker.dt[, list(N=.N, mean_log2FC_cluster=mean(avg_log2FC), mean_pct_cluster=mean(pct.1)),  by=cluster]
      write.csv(smry,file=paste0(assay,"_",res_str,"_clustersDEGsummary.csv"))
      p1 = ggplot(data=marker.dt, aes(x=avg_log2FC, y=cluster, fill=cluster, color=cluster))+
                  ggridges::geom_density_ridges(jittered_points = TRUE, scale = .95, rel_min_height = .01,
                                              point_shape = "|", point_size = 2, size = 0.25, position = position_points_jitter(height = 0), alpha=0.2)+
                  geom_vline(xintercept=0.5, color="grey20", linetype="dashed", lwd=0.2)+
                  labs(y=paste0("Louvain clusters, resolution ",res_str), x="Average log2 fold change of cluster markers")
      p2 = ggplot(data=marker.dt, aes(x=pct.1, y=cluster, fill=cluster, color=cluster))+
          ggridges::geom_density_ridges(jittered_points = TRUE, scale = .95, rel_min_height = .01,
                                      point_shape = "|", point_size = 2, size = 0.25, position = position_points_jitter(height = 0), alpha=0.2)+
          labs(y=paste0("Louvain clusters, resolution ",res_str), x="Percentage of cells expressing marker in cluster")+
          geom_vline(xintercept=0.3, color="grey20", linetype="dashed", lwd=.2)
      pw = p1 / p2 + patchwork::plot_layout(guides="collect") & scale_fill_manual(values=cluster.cols) & scale_color_manual(values=cluster.cols) & theme_linedraw()
      ggsave(paste0(lfc_pdf,res_str,".pdf"),height=10,width=10,plot=pw,device='pdf')

      heatmap_data.cols = circlize::colorRamp2(breaks=c(0:5),hcl_palette="Inferno")
      heatmap_scale.cols = circlize::colorRamp2(breaks=c(-3.5,0,3.5),hcl_palette = "Lisbon")
      heatmap_marker.cols = circlize::colorRamp2(breaks=c(0:5), hcl_palette="Inferno")

      marker.genes = top5$gene
      x.hm = subset(x, downsample=100)

      mat = x.hm[[assay]]@data[marker.genes, ] %>% as.matrix()

      ra = ComplexHeatmap::rowAnnotation(`Cluster`=as.character(top5$cluster), col=list(`Cluster`=cluster.cols), show_legend=FALSE)
      if (hs_or_mm == "mm"){
      ta_immgen = ComplexHeatmap::HeatmapAnnotation(`Cluster`=x.hm@meta.data[[res[iter]]], `ImmGen`=x.hm@meta.data$immgen.pred_pruned, col=list(`ImmGen`=annot_colors[["immgen.pred_pruned"]],`Cluster`=cluster.cols))
      ta_mrsd = ComplexHeatmap::HeatmapAnnotation(`Cluster`=x.hm@meta.data[[res[iter]]], `MRSD`=x.hm@meta.data$mrsd.pred_pruned, col=list(`MRSD`=annot_colors[["mrsd.pred_pruned"]],`Cluster`=cluster.cols))
      ta_azi = ComplexHeatmap::HeatmapAnnotation(`Cluster`=x.hm@meta.data[[res[iter]]], `L2_LungmapAzi`=x.hm@meta.data$predicted.celltype_level2, col=list(`L2_LungmapAzi`=annot_colors[["predicted.celltype_level2"]],`Cluster`=cluster.cols))
      }
      else{
      ta_DICE = ComplexHeatmap::HeatmapAnnotation(`Cluster`=x.hm@meta.data[[res[iter]]], `DICE`=x.hm@meta.data$DICE.pred_pruned, col=list(`DICE`=annot_colors[["DICE.pred_pruned"]],`Cluster`=cluster.cols))
      ta_HPCA = ComplexHeatmap::HeatmapAnnotation(`Cluster`=x.hm@meta.data[[res[iter]]], `HPCA`=x.hm@meta.data$HPCA.pred_pruned, col=list(`HPCA`=annot_colors[["HPCA.pred_pruned"]],`Cluster`=cluster.cols))
      ta_azi = ComplexHeatmap::HeatmapAnnotation(`Cluster`=x.hm@meta.data[[res[iter]]], `L2_LungmapAzi`=x.hm@meta.data$predicted.ann_level_3, col=list(`L2_LungmapAzi`=annot_colors[["predicted.ann_level_3"]],`Cluster`=cluster.cols))
      }
      hm_list = list()
      for (ta in objects()[grepl("^ta_",objects())]){
          hm_list[[ta]] = ComplexHeatmap::Heatmap(mat, name = "Expression", cluster_columns = FALSE, cluster_rows=FALSE, column_split=x.hm@meta.data[[res[iter]]],
                      row_split=top5$cluster, row_names_gp = grid::gpar(fontsize = 11), column_title=character(0), column_gap = unit(0.5, "mm"),
                      col = heatmap_data.cols, top_annotation = get(ta), show_column_names = FALSE, left_annotation = ra)
      }
      nc = length(unique(markers$cluster))
      pdf(paste0(annot_pdf,res_str,".pdf"),height=round_any(nc,8,f=ceiling),width=round_any(nc,8,f=ceiling))
      lapply(hm_list,ComplexHeatmap::draw,merge_legends = TRUE, padding = unit(c(2, 2, 2, 8), "mm"))
      dev.off()

      p = DotPlot(x, assay=assay, features=unique(marker.genes), group.by=res[iter], scale=TRUE)+
          theme(axis.text.x=element_text(size=6), axis.title.x=element_blank())+
          labs(y=paste0("Louvain2 clusters - resolution ",res_str))+ theme_90()
      ggsave(paste0(dotplot,res_str,".pdf"),plot=p,height=8,width=14,device="pdf")
  }
}
feature_plots = function(x, barplot, ridgeplot,assay,hs_or_mm){ #Doesn't consider integration type i.e. Harm vs SCT
  annot_colors = generate_annot_cols(x)
  if (hs_or_mm == "mm"){gb = list("mrsd","immgen","predicted.celltype_level2")
  }else{gb=list("DICE","HPCA","predicted.ann_level_3")}
  dat <- x@meta.data
  if (assay == "RNA"){
    p1 = ggplot(data=dat, aes(x=RNA_snn_res.0.2, fill=eval(parse(text=gb[[3]]))))+
        geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.2"))+scale_fill_manual("L2_LungmapAzi", values=annot_colors[[gb[[3]]]])
    p2 = ggplot(data=dat, aes(x=RNA_snn_res.0.2, fill=get(paste0(gb[1],".pred_pruned"))))+
        geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.2"))+scale_fill_manual(gb[1], values=annot_colors[[paste0(gb[1],".pred_pruned")]])
    p3 = ggplot(data=dat, aes(x=RNA_snn_res.0.2, fill=paste0(gb[2],".pred_pruned")))+
        geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.2"))+scale_fill_manual(gb[2], values=annot_colors[[paste0(gb[2],".pred_pruned")]])
    p4 = ggplot(data=dat, aes(x=RNA_snn_res.0.5, fill=eval(parse(text=gb[[3]]))))+
        geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.5"))+scale_fill_manual("L2_LungmapAzi", values=annot_colors[[gb[[3]]]])
    p5 = ggplot(data=dat, aes(x=RNA_snn_res.0.5, fill=get(paste0(gb[1],".pred_pruned"))))+
        geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.5"))+scale_fill_manual(gb[1], values=annot_colors[[paste0(gb[1],".pred_pruned")]])
    p6 = ggplot(data=dat, aes(x=RNA_snn_res.0.5, fill=paste0(gb[2],".pred_pruned")))+
        geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.5"))+scale_fill_manual(gb[2], values=annot_colors[[paste0(gb[2],".pred_pruned")]])
    g = cowplot::plot_grid(p1,p2, p3, ncol = 1)
    g2 = cowplot::plot_grid(p4,p5, p6, ncol = 1)
    ggsave(paste0(barplot,"_resolution_0.2.pdf"),plot=g,width=10,height=18,device="pdf")
    ggsave(paste0(barplot,"_resolution_0.5.pdf"),plot=g2,width=10,height=18,device="pdf")
    p1 = RidgePlot(x, features=c("nCount_RNA", "nFeature_RNA"), group.by="RNA_snn_res.0.2", sort=TRUE) & scale_fill_manual(values=cluster.cols) & ggtitle("Resolution - 0.2")
    p2 = RidgePlot(x, features=c("nCount_RNA", "nFeature_RNA"), group.by="RNA_snn_res.0.5", sort=TRUE) & scale_fill_manual(values=cluster.cols) & ggtitle("Resolution - 0.5")
    p = cowplot::plot_grid(p1, p2, ncol = 1, nrow= 2)
    ggsave(paste0(ridgeplot,".pdf"),device='pdf',height=12,width=8,plot=p)
  }
  if (assay=="integrated"){
    #p1 = ggplot(data=dat, aes(x=RNA_snn_res.0.2, fill=eval(parse(text=gb[[3]]))))+
    #    geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.2"))+scale_fill_manual("L2_LungmapAzi", values=annot_colors[[gb[[3]]]])
    p2 = ggplot(data=dat, aes(x=integrated_snn_res.0.2, fill=get(paste0(gb[1],".pred_pruned"))))+
        geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.2"))+scale_fill_manual(gb[1], values=annot_colors[[paste0(gb[1],".pred_pruned")]])
    p3 = ggplot(data=dat, aes(x=integrated_snn_res.0.2, fill=paste0(gb[2],".pred_pruned")))+
        geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.2"))+scale_fill_manual(gb[2], values=annot_colors[[paste0(gb[2],".pred_pruned")]])
    #p4 = ggplot(data=dat, aes(x=RNA_snn_res.0.5, fill=eval(parse(text=gb[[3]]))))+
    #    geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.5"))+scale_fill_manual("L2_LungmapAzi", values=annot_colors[[gb[[3]]]])
    p5 = ggplot(data=dat, aes(x=integrated_snn_res.0.5, fill=get(paste0(gb[1],".pred_pruned"))))+
        geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.5"))+scale_fill_manual(gb[1], values=annot_colors[[paste0(gb[1],".pred_pruned")]])
    p6 = ggplot(data=dat, aes(x=integrated_snn_res.0.5, fill=paste0(gb[2],".pred_pruned")))+
        geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.5"))+scale_fill_manual(gb[2], values=annot_colors[[paste0(gb[2],".pred_pruned")]])
    g = cowplot::plot_grid(p2, p3, ncol = 1)
    g2 = cowplot::plot_grid(p5, p6, ncol = 1)
    ggsave(paste0(barplot,"_resolution_0.2.pdf"),plot=g,width=10,height=12,device="pdf")
    ggsave(paste0(barplot,"_resolution_0.5.pdf"),plot=g2,width=10,height=12,device="pdf")
    p1 = RidgePlot(x, features=c("nCount_RNA", "nFeature_RNA"), group.by="integrated_snn_res.0.2", sort=TRUE) & scale_fill_manual(values=cluster.cols) & ggtitle("Resolution - 0.2")
    p2 = RidgePlot(x, features=c("nCount_RNA", "nFeature_RNA"), group.by="integrated_snn_res.0.5", sort=TRUE) & scale_fill_manual(values=cluster.cols) & ggtitle("Resolution - 0.5")
    p = cowplot::plot_grid(p1, p2, ncol = 1, nrow= 2)
    ggsave(paste0(ridgeplot,".pdf"),device='pdf',height=12,width=8,plot=p)
  }
  if (assay=="SCT"){
    p1 = ggplot(data=dat, aes(x=SCT_snn_res.0.2, fill=eval(parse(text=gb[[3]]))))+
        geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.2"))+scale_fill_manual("L2_LungmapAzi", values=annot_colors[[gb[[3]]]])
    p2 = ggplot(data=dat, aes(x=SCT_snn_res.0.2, fill=get(paste0(gb[1],".pred_pruned"))))+
        geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.2"))+scale_fill_manual(gb[1], values=annot_colors[[paste0(gb[1],".pred_pruned")]])
    p3 = ggplot(data=dat, aes(x=SCT_snn_res.0.2, fill=paste0(gb[2],".pred_pruned")))+
        geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.2"))+scale_fill_manual(gb[2], values=annot_colors[[paste0(gb[2],".pred_pruned")]])
    p4 = ggplot(data=dat, aes(x=SCT_snn_res.0.5, fill=eval(parse(text=gb[[3]]))))+
        geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.5"))+scale_fill_manual("L2_LungmapAzi", values=annot_colors[[gb[[3]]]])
    p5 = ggplot(data=dat, aes(x=SCT_snn_res.0.5, fill=get(paste0(gb[1],".pred_pruned"))))+
        geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.5"))+scale_fill_manual(gb[1], values=annot_colors[[paste0(gb[1],".pred_pruned")]])
    p6 = ggplot(data=dat, aes(x=SCT_snn_res.0.5, fill=paste0(gb[2],".pred_pruned")))+
        geom_bar(position="fill", color="black")+labs(x=paste0("Louvain2 clusters - resolution ",".0.5"))+scale_fill_manual(gb[2], values=annot_colors[[paste0(gb[2],".pred_pruned")]])
    g = cowplot::plot_grid(p1,p2, p3, ncol = 1)
    g2 = cowplot::plot_grid(p4,p5, p6, ncol = 1)
    ggsave(paste0(barplot,"_resolution_0.2.pdf"),plot=g,width=10,height=18,device="pdf")
    ggsave(paste0(barplot,"_resolution_0.5.pdf"),plot=g2,width=10,height=18,device="pdf")
    p1 = RidgePlot(x, features=c("nCount_RNA", "nFeature_RNA"), group.by="SCT_snn_res.0.2", sort=TRUE) & scale_fill_manual(values=cluster.cols) & ggtitle("Resolution - 0.2")
    p2 = RidgePlot(x, features=c("nCount_RNA", "nFeature_RNA"), group.by="SCT_snn_res.0.5", sort=TRUE) & scale_fill_manual(values=cluster.cols) & ggtitle("Resolution - 0.5")
    p = cowplot::plot_grid(p1, p2, ncol = 1, nrow= 2)
    ggsave(paste0(ridgeplot,".pdf"),device='pdf',height=12,width=8,plot=p)
  }
}
subcluster = function(x,tx,integ,parent,hs_or_mm){ #tx should be a list() with names of CMO groups attached to CMO labels, i.e. list('model'=c("CMO_1","CMO_2"))
  annot_colors = generate_annot_cols(x)
  if (hs_or_mm == "mm"){gb=list("mrsd","immgen","predicted.celltype_level2")}
  else{gb=list("DICE","HPCA","predicted.ann_level_3")}
  q = NULL ### Assign every sample to a group name
  for (ele in tx){
      q = c(q,ele)
  }
  iter = 0
  sample2.cols = NULL ###Set sample colors
  for (ele in q){
      iter = iter + 1
      sample2.cols = c(sample2.cols,ele = sample.cols[iter])
  }
  names(sample2.cols) = q ###Ugly way to set names

  for (resolution in colnames(x@meta.data)[colnames(x@meta.data) %like% "*res.[0-9].[0-9]"]){
      for (i in unique((x@meta.data[[resolution]]))){
          eval(parse(text=paste0("sub",i,"=subset(x,subset=`",resolution,"`=='",i,"')"))) #Programatically assign subcluster objects to subclusters
      }
      for (i in objects()[objects() %like% "sub\\d{1}"]){
          iter = as.character(str_match(i,"\\d{1,}"))
          wod  = paste0(parent,"/",integ,"/",resolution,"/subcluster",iter,"/")
          print(wod)
          dir.create(file.path(wod),recursive = TRUE)
          setwd(wod)
          sub_temp = get(i)
          DefaultAssay(sub_temp) = "RNA"
          print(paste0("Starting ",i," dimensional reduction"))
          sub_temp = SCTransform(sub_temp,assay="RNA",vars.to.regress = c("percent.mt","S.Score","G2M.Score"),vst.flavor="v2")
          #sub_temp <- Seurat::NormalizeData(sub_temp, normalization.method = "LogNormalize", scale.factor = 10000, verbose=FALSE)
          #sub_temp <- Seurat::FindVariableFeatures(sub_temp, selection.method="vst", nfeatures=2000)
          #sub_temp <- Seurat::ScaleData(sub_temp, verbose=FALSE)
          sub_temp <- Seurat::RunPCA(sub_temp, assay="SCT",npcs=50)
          pdf(paste0("Cluster_",iter,"_PCA_Loadings.pdf"),width=8,height=8)
          DimHeatmap(sub_temp, dims=1:9, nfeatures=10, fast=TRUE, slot="scale.data", assays="SCT")
          dev.off()
          print(paste0("Dim Reduction Complete"))
          q = DimPlot(sub_temp, group.by="orig.ident", reduction="pca") + scale_color_manual(values=sample.cols)
          ggsave(paste0("PCA_by_CMO_harm_Cluster_",iter,".pdf"),plot=q,device='pdf')

          sub_temp_meta <- sub_temp@meta.data
          sub_temp_meta$barcode <- unlist(tstrsplit(rownames(sub_temp_meta),"_",keep=1))
          #sub_temp_meta$cell_id <- ifelse(sub_temp_meta$orig.ident=="untreated",
          #                          paste0("Untreated_", sub_temp_meta$barcode),
          #                          paste0("Treated_", sub_temp_meta$barcode))
          sub_temp_meta$rn <- rownames(sub_temp_meta)
          p1 <- DimPlot(sub_temp, group.by="orig.ident", reduction="pca") +
                  labs(title="Treatment") +
                  scale_color_manual("Treatment", values=sample.cols)

          p2 <- DimPlot(sub_temp, group.by=paste0(gb[[1]],".pred_pruned"), reduction="pca") +
                  labs(title=paste0(gb[1]," annotation")) +
                  scale_color_manual("Annotation", values=annot_colors[[paste0(gb[[1]],".pred_pruned")]])

          p3 <- DimPlot(sub_temp, group.by=paste0(gb[[2]],".pred_pruned"), reduction="pca") +
                  labs(title=paste0(gb[[2]]," annotation")) +
                  scale_color_manual("Annotation", values=annot_colors[[paste0(gb[[2]],".pred_pruned")]])

          p4 <- DimPlot(sub_temp, group.by=gb[[3]], reduction="pca") +
                  labs(title="L2_LungmapAzi") +
                  scale_color_manual("L2_LungmapAzi", values=annot_colors[[gb[[3]]]])

          p5 <- DimPlot(sub_temp, group.by="Phase", reduction="pca") +
                  labs(title="Cell cycle") +
                  scale_color_manual("Cell cycle", values=cell_cycle.cols)

          pw <- p1 + p2 + p3 + p4 + p5 + patchwork::plot_layout(ncol=2,nrow=3)
          ggsave(paste0("Sub_",iter,"_PCA_harm.pdf"),width=15,height=10,plot=pw)
          sub_temp@meta.data$orig.ident = factor(sub_temp@meta.data$orig.ident)
          if (nrow(sub_temp@meta.data) < 30){
            neighbors = round(nrow(sub_temp@meta.data)*.10,0)
          }else{
            neighbors = 30
          }
          #sub_temp <- harmony::RunHarmony(sub_temp, dims=1:bdim(sub_temp), group.by.vars = "orig.ident")
          sub_temp <- Seurat::RunUMAP(sub_temp, assay.use="SCT", dims=1:bdim(sub_temp), reduction="pca",n.neighbors=neighbors)
          sub_temp <- Seurat::FindNeighbors(sub_temp, reduction="pca", dims=1:bdim(sub_temp))
          sub_temp <- Seurat::FindClusters(sub_temp, algorithm=2, resolution=c(0.1, 0.2))

          p1 <- DimPlot(sub_temp, group.by="SCT_snn_res.0.1")+
          rcartocolor::scale_color_carto_d(palette = "Bold")

          p2 <- DimPlot(sub_temp, group.by=gb[[3]])+
          labs(title="L2_LungmapAzi") +
          scale_color_manual("L2_LungmapAzi", values=annot_colors[[gb[[3]]]])

          p3 <- DimPlot(sub_temp, group.by="SCT_snn_res.0.2")+
          rcartocolor::scale_color_carto_d(palette = "Bold")

          p4 <- DimPlot(sub_temp, group.by="orig.ident")+
          labs(title="Treatment")+
          scale_color_manual("Treatment", values=sample.cols)

          p5 <- FeaturePlot(sub_temp, features="percent.mt")+
          viridis::scale_color_viridis(option="mako")

          p6 <- FeaturePlot(sub_temp, features="Phase")+
          labs(title="Cell cycle")+
          scale_color_manual("Cell cycle", values=cell_cycle.cols)

          pw <- p1 + p2 + p3 + p4 + p5 + p6 + patchwork::plot_layout(ncol=2, nrow=3)
          ggsave(paste0("Sub_",iter,"_subclustering_All.pdf"),height=10,width=10,plot=pw)

          #Single-cell differential GEX approach (Not as good as "pseudobulk")
          Idents(sub_temp) <- sub_temp@meta.data$orig.ident
          sub_temp@meta.data$orig.ident=factor(sub_temp@meta.data$orig.ident)
          c(v1,v2) %<-% c(1,2)
          flag = FALSE
          scde_list = list()
          scde_flag = TRUE
          for (i in names(tx)){
            tx[[i]] = tx[[i]][tx[[i]] %in% levels(sub_temp@meta.data$orig.ident)]
          }
          scde_flag = tryCatch(function(){
          for (i in 1:(length(names(tx))*(length(names(tx))-1)/2)){
              if (flag == FALSE){
                  v2_og = v2
                  flag = TRUE
              }
              scde_results <- Seurat::FindMarkers(sub_temp, slot="data",ident.1=tx[[v1]],ident.2=tx[[v2]],logfc.threshold=0.25,test.use="wilcox", min.pct=0.1, min.diff.pct=-Inf, only.pos=FALSE)
              scde_results <- cbind(gene=rownames(scde_results),scde_results)
              scde_results$Sig <- data.table::fcase(scde_results$p_val_adj <= 0.05, "FDR", scde_results$p_val <= 0.05, "Nominal", default="NS")
              scde_list[[paste0(names(tx)[v1],".",names(tx)[v2])]] = scde_results
              if (v2 < length(tx)){
                  v2 = v2+1
              }else{
                  v1 = v1+1
                  v2_og = v2_og+1
                  v2 = v2_og}
              }
              scde_flag = TRUE
              #saveRDS(scde_list,"MySCDE.rds")
              return(list(scde_flag,scde_list))},
              error=function(cond){
                scde_flag = FALSE
                message(paste0("SCDE failure: Cluster ", iter, " // GROUPS: ",names(tx)[v1]," ",names(tx)[v2]))
                message(paste0(cond))
                return(scde_flag)
              },
              warning=function(cond){
                message(cond)
                return(list(scde_flag,scde_list))
              },
              finally={}
          )
          scde = scde_flag()
          try(for (i in names(scde[[2]])){
            write.csv(file=paste0(i,"_DEG.csv"),x=scde[[2]][i])
          })
          #Pseudobulk GEX approach
          ##Checking data metrics
          dat <- as.data.table(sub_temp@meta.data)
          dat[, .N, by=orig.ident]

          counts <- Seurat::AggregateExpression(sub_temp,assays="SCT",group.by="orig.ident",slot="counts")$SCT
          filt_counts <- counts[rowSums(counts>=10)>1,]

          meta <- sub_temp@meta.data
          meta <- meta[!duplicated(meta$orig.ident),]
          rownames(meta) <- meta$orig.ident
          colnames(counts) = meta$orig.ident
          colnames(filt_counts) = meta$orig.ident
          meta <- meta[match(colnames(counts), rownames(meta)),]
          # set factor levels for DE analysis
          meta$group <- factor(meta$orig.ident)#, levels=levels(sub_temp$orig.ident))
          dds <- DESeq2::DESeqDataSetFromMatrix(countData=filt_counts, colData=meta, design=~group)
          dds <- DESeq2::estimateSizeFactors(dds, type="poscounts")
          lognorm_counts <- log2(DESeq2::counts(dds, normalized=TRUE)+1)
          mlognorm_counts <- reshape2::melt(lognorm_counts)

          p = ggplot(data=mlognorm_counts, aes(x=value, y=Var2, fill=Var2, color=Var2)) + ggridges::geom_density_ridges(alpha=0.8)+
              theme_linedraw() + labs(x="log2 (normalized counts) + 1", y="") + scale_fill_manual("Condition", values=sample.cols) + scale_color_manual("Condition", values=sample.cols)
          ggsave(paste0("Sub_",iter,"_Lognorm_Counts.pdf"),plot=p)

          pdf(paste0("Sub_",iter,"_Expression_PCA.pdf"))
          pca <- FactoMineR::PCA(t(lognorm_counts), scale.unit=TRUE, ncp=10, graph=FALSE)
          plot(factoextra::fviz_screeplot(pca))
          plot(factoextra::fviz_pca_ind(pca))
          dev.off()

          if (scde[[1]] == TRUE){
              scde_list=scde[[2]]
            for (comp in names(scde_list)){
                hm_results = scde_list[[comp]][!is.na(scde_list[[comp]]$gene),]
                hm_results = hm_results[order(hm_results$p_val, decreasing=FALSE),]
                if (nrow(hm_results)>=75){
                  hm_results = hm_results[1:75,]
                }else{
                  hm_results = hm_results[1:nrow(hm_results),]
                }
                de_genes = hm_results$gene

                mat = sub_temp[["SCT"]]@data[de_genes, ] %>% as.matrix()
                mat = t(scale(t(mat)))
                top_ann = ComplexHeatmap::HeatmapAnnotation(`Sample`=sub_temp@meta.data$orig.ident, col=list(`Sample`=sample2.cols), annotation_name_gp = grid::gpar(fontsize=9, fontface="bold"))

                # row colors
                perm <- match(rownames(mat), hm_results$gene)
                hm_results <- hm_results[perm,]
                row_ann <- ComplexHeatmap::rowAnnotation(`Significance`=hm_results$Sig, col=list(`Significance`=sig.cols), annotation_name_gp = grid::gpar(fontsize=9, fontface="bold"))

                pdf(paste0("Sub_",iter,"_",comp,"_Significance_Expression_HM.pdf"),height=15,width=15)
                ch <- ComplexHeatmap::Heatmap(mat, name = "Expression", cluster_columns = TRUE, column_split=sub_temp@meta.data$orig.ident,
                                        cluster_rows=TRUE, row_split=hm_results$Sig, col=heatmap_scale.cols, row_names_gp = grid::gpar(fontsize = 11),
                                        column_title=character(0), top_annotation = top_ann, show_column_names = FALSE, left_annotation = row_ann)
                ComplexHeatmap::draw(ch, merge_legends = TRUE, padding = unit(c(2, 2, 2, 8), "mm"))
                dev.off()

                pdf(paste0("Sub_",iter,"_",comp,"_Vln_Top5_Sigdiff.pdf"),width=15)
                top5_genes <- hm_results$gene[1:6]
                q = VlnPlot(sub_temp, features=top5_genes, assay="SCT", group.by="orig.ident") &
                    scale_fill_manual(values=sample.cols) & theme_linedraw() & theme(axis.text.x=element_blank(), axis.title.x=element_blank())
                plot(q)
                dev.off()

                x_axis_lim <- abs(max(scde_list[[comp]]$avg_log2FC))+0.75
                p <- ggplot(data=scde_list[[comp]],
                            aes(x=avg_log2FC, y=-log10(p_val), color=Sig))+
                            geom_point(size=2.4, pch=19, alpha=0.7)+
                            labs(x="Average Log2 Fold Change", y="-log10(pvalue)")+
                            xlim(-x_axis_lim, x_axis_lim)+
                            scale_color_manual("Significance", values=sig.cols)


                #GSEA geneset expression analysis
                pathways <- fgsea::gmtPathways(paste0(parent,"/GMT.gmt"))
                ranks <- scde_list[[comp]]$avg_log2FC
                names(ranks) <- scde_list[[comp]]$gene
                ranks <- ranks[order(ranks, decreasing=TRUE)]
                gsea_results <- fgsea::fgsea(pathways=pathways,
                                            stats=ranks,
                                            minSize=1)
                setorder(gsea_results, pval, na.last=TRUE)
                fwrite(gsea_results,paste0("Sub_",iter,"_",comp,"_GSEA_HALLMARK.csv"))
                gsea_plot <- gsea_results
                scale_lim <- max(abs(gsea_plot$NES))
                gsea_plot$Pathway <- gsub("HALLMARK_", "", gsea_plot$pathway)
                gsea_plot$pval <- ifelse(gsea_plot$padj > 0.05, NA, gsea_plot$padj)
                p <- ggplot(data=gsea_plot, aes(x=NES, y=reorder(Pathway, NES), fill=-log10(pval)))+
                geom_bar(stat="identity", aes(color=""))+
                scale_fill_viridis_c(option="mako", na.value="#DBE2E9")+
                scale_color_manual(values="#DBE2E9")+
                guides(color=guide_legend("Not significant", override.aes=list(fill="#DBE2E9")))+
                theme_linedraw()
                pdf(paste0("Sub_",iter,"_",comp,"_Barplot_NES_Paths.pdf"))
                plot(p)
                dev.off()
              }
            }else{
              message("Skipping GSEA due to SCDE failure.")
            }
          setwd(parent)
      }
  remove(list=objects()[objects() %like% "^sub\\d{1,}"])
  }
}

doCellChat = function(seuratObject,groups,hs_mm,threads){ #May have to set .libPaths to include R_LIB like .libPaths(new=c(.libPaths(),"/work/users/c/a/came/R_LIBS"))
    cellchat = createCellChat(object=seuratObject,group.by=groups)
    if (tolower(hs_mm) == "mm"){CellChatDB = CellChatDB.mouse}else{CellChatDB = CellChatDB.human}
    cellchat@DB=CellChatDB
    cellchat = subsetData(cellchat)
    future::plan("multisession", workers = as.integer(threads))
    cellchat = identifyOverExpressedGenes(cellchat)
    cellchat = identifyOverExpressedInteractions(cellchat)
    cellchat = projectData(cellchat, PPI.human) #Uses projected data as opposed to assay data, requires raw.use=FALSE in computeCommunProb function below.
    cellchat = computeCommunProb(cellchat, type = "triMean",raw.use=FALSE)
    cellchat = filterCommunication(cellchat, min.cells = 10)
    cellchat = computeCommunProbPathway(cellchat)
    cellchat = aggregateNet(cellchat)
    return(cellchat)
}

doCellChatSummary = function(cellchatobj){
    #write communication summary
    comm_frame = subsetCommunication(cellchatobj)
    write.csv(comm_frame,"Communication_Summary.csv")
    
    #Interactions + weight chord diagram
    #groupSize = as.numeric(table(cellchatobj@idents))
    #pdf("Thisexample.pdf",width=10,height=8)
    #matr = matrix(data=c(1,2),nrow=1,ncol=2) doesnt work right currently, wont cooperate with saving an image correctly
    #layout(mat=matr)
    #netVisual_circle(cellchatobj@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    #netVisual_circle(cellchatobj@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    #dev.off()
    
    #Closely look at interactions and weights
    mat = cellchatobj@net$weight
    myl = list()
    groupSize = as.numeric(table(cellchatobj@idents))
    pdf("Interactions_Weights_Minimal.pdf",width=18,height=18)
    par(mfrow = c(5,4), xpd=TRUE)
    for (i in 1:nrow(mat)) {
        idx=paste0("P",as.character(i))
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i, ] <- mat[i, ]
        netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    }
    dev.off()
    return(comm_frame)
}