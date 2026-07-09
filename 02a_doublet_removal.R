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


# Log Normalization 하는 이유: scRNA-seq에서 cell마다 읽힌 RNA 양, UMI count, sequencing depth 가 다르기 때문.
# NormalizeData(): 각 cell의 count를 cell 전체 count 기준으로 보정. gene count / cell total count * 10000 후 로그 변환. 
# FindVariableFeatures(): cell마다 차이가 큰 gene 2000개 고름
# ScaleData(): gene별 expression 값을 평균 0, 분산 1에 가깝게 맞추는 단계.
# RunPCA(): gene을 몇개 pc로 줄이는 단계 (차원축소)
# ElbowPlot(): PCA에서 몇 개의 PC를 사용할지 판단하는 진단용 plot (PC 1~10까지 사용)
# FindNeighbors(): 각 cell이 PCA 공간에서 어떤 cell들과 가까운지 계산하는 단계. 
# FindClusters(): 비슷한 cell들을 cluster 단위로 묶는 단계.
