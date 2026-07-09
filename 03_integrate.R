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
...
data10.singlet$Origin <- "Sample10"
merged.seurat <- merge(x = data1.singlet, 
                       y = list(data2.singlet, ... , data10.singlet),
                       add.cell.ids = c("Sa_1", ..., "Sa_10"),
                       project = "Samples_merged"
                       )


# 3. Normalization + Batch Correction
merged.seurat <- SCTransform(
  merged.seurat, 
  assay = "RNA",
  new.assay.name = "SCT",
  method = "glmGamPoi"
)
DefaultAssay(merged.seurat) <- "SCT"
merged.seurat <- RunPCA(merged.seurat, assay = "SCT", npcs = 10)
merged.seurat <- RunHarmony(
  merged.seurat,
  group.by.vars = "Origin"
  #assay.use = "SCT"
)
merged.seurat <- RunUMAP(merged.seurat, reduction = "harmony", dims = 1:10)
merged.seurat <- FindNeighbors(merged.seurat, reduction = "harmony", dims = 1:10)
merged.seurat <- FindClusters(merged.seurat, resolution = 0.6)
DimPlot(merged.seurat, reduction = "umap", group.by = "Origin")
DimPlot(merged.seurat, reduction = "umap", label = TRUE)



