# Data integration and UMAP Drawing
# 1. Memory Clearing
rm(
  data1.data, 
  data1.seurat, 
  data1.seurat.v3, 
  sweep.res, sweep.stats, bcmvn,
  best.pK, nExp, homotypic.prop, nExp.adj,
  df.col1
)
gc()

# 2. Data Merging
data1.singlet$Origin <- "Sample1"
