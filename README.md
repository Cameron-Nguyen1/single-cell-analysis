# single-cell-analysis

This repository is designed to facilitate the analysis of single-cell RNA-seq data that has been processed by 10x genomics count/aggr pipeline.
To run the analysis, execute each R script in the following order:
  sc_mergeqc.r (1st)
  sc_normal.r 
  sc_sct.r
  sc_harm.r
  sc_analysis.r (last)

This pipeline assumes that all work is being performed in the working directory. It assumes additional non-R files (GMT.gmt, HMD*.txt) are located within the working directory as well.
