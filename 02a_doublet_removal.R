# Log Normalization
# Data 1
DefaultAssay(data1.seurat) <- "RNA"
data1.seurat <- NormalizaeData(data1.seurat, normalization.method = "LogNormalize", scale.factor = 1e4)
data1.seurat <- FindVariableFeatures(data1.seurat, selection.method = "vst", nfeatures = 2000)
data1.seurat <- ScaleData(Data1.seurat, features = rownames(data1.seurat))
data.seurat <- RunPCA(data1.seurat, features = VariableFeatures(object = data1.seurat))
ElbowPlot(data1.seurat)
data1.seurat <- FindNeighbors(data1.seurat, dims = 1:10)
data1.seurat <- FindClusters(data1.seurat, resolution = 0.6)

