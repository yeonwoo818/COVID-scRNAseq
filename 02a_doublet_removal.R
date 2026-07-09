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

# Data 10 까지 반복

# Log Normalization 하는 이유: scRNA-seq에서 cell마다 읽힌 RNA 양, UMI count, sequencing depth 가 다르기 때문.
# NormalizeData(): 각 cell의 count를 cell 전체 count 기준으로 보정. gene count / cell total count * 10000 후 로그 변환. 
# FindVariableFeatures(): cell마다 차이가 큰 gene 2000개 고름
# ScaleData(): gene별 expression 값을 평균 0, 분산 1에 가깝게 맞추는 단계.
# RunPCA(): gene을 몇개 pc로 줄이는 단계 (차원축소)
# ElbowPlot(): PCA에서 몇 개의 PC를 사용할지 판단하는 진단용 plot (PC 1~10까지 사용)
# FindNeighbors(): 각 cell이 PCA 공간에서 어떤 cell들과 가까운지 계산하는 단계. 
# FindClusters(): 비슷한 cell들을 cluster 단위로 묶는 단계.




# Doublet Deletion Using DoubletFinder
# Data 1
DefaultAssay(data1.seurat) <- "RNA"
data1.seurat.v3 <- data1.seurat
data1.seurat.v3[["RNA"]] <- as(object = data1.seurat.v3[["RNA"]], Class = "Assay")
sweep.res <- paramSweep_v3(data1.seurat.v3, PCs = 1:10)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
best.pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
nExp <- round(ncol(data1.seurat) * 0.075)
homotypic.prop <- modelHomotypic(data1.seurat.v3$seurat_clusters)
nExp.adj <- round(nExp * (1 - homotypic.prop))
data1.seurat.v3 <- doubletFinder_v3(
  data1.seurat.v3,
  pN = 0.25,
  pK = best.pK,
  nExp = nExp.adj,
  PCs = 1:10
)
head(data1.seurat.v3)
unique(data1.seurat.v3$DF.classifications_0.25_0.005_298)
df.col1 <- tail(grep("^DF\\.classifications_", colnames(data1.seurat.v3@meta.data), value = TRUE), 1)
data1.seurat$doublet_class_log <- data1.seurat.v3@meta.data[colnames(data1.seurat), df.col1]
head(data1.seurat$doublet_class_log)
data1.singlet <- subset(data1.seurat, subset = doublet_class_log == "Singlet")

# Data 10 까지 반복


