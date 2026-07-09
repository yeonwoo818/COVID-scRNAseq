library(Matrix)
library(Seurat)
library(harmony)
library(dplyr)
library(tidyr)
library(SingleR)
library(ggplot2)
library(DoubletFinder)
library(celldex)
library(BiocParallel)
library(spam)

# Reading Preprocessed Datasets
data1.data <- Read10X('~/COVID_preprocessing_Sample1/Solo_Sa1.out/Gene/filtered')
data2.data <- Read10X('~/COVID_preprocessing_Sample1/Solo_Sa2.out/Gene/filtered')
data3.data <- Read10X('~/COVID_preprocessing_Sample2/Solo_Sa1.out/Gene/filtered')
...
data10.data <- Read10X('~/COVID_preprocessing_Sample3/Solo_Sa4.out/Gene/filtered')


data1.seurat <- CreateSeuratObject(data1.data, project = 'Test', min.cells=3, min.features=200)
data2.seurat <- CreateSeuratObject(data2.data, project = 'Test', min.cells=3, min.features=200)
data3.seurat <- CreateSeuratObject(data3.data, project = 'Test', min.cells=3, min.features=200)
data4.seurat <- CreateSeuratObject(data4.data, project = 'Test', min.cells=3, min.features=200)
data5.seurat <- CreateSeuratObject(data5.data, project = 'Test', min.cells=3, min.features=200)
data6.seurat <- CreateSeuratObject(data6.data, project = 'Test', min.cells=3, min.features=200)
data7.seurat <- CreateSeuratObject(data7.data, project = 'Test', min.cells=3, min.features=200)
data8.seurat <- CreateSeuratObject(data8.data, project = 'Test', min.cells=3, min.features=200)
data9.seurat <- CreateSeuratObject(data9.data, project = 'Test', min.cells=3, min.features=200)
data10.seurat <- CreateSeuratObject(data10.data, project = 'Test', min.cells=3, min.features=200)


# Preprocessing Seurat Objects
VlnPlot(data1.seurat, features = c('nCount_RNA', 'nFeature_RNA'))
VlnPlot(data2.seurat, features = c('nCount_RNA', 'nFeature_RNA'))
VlnPlot(data3.seurat, features = c('nCount_RNA', 'nFeature_RNA'))
VlnPlot(data4.seurat, features = c('nCount_RNA', 'nFeature_RNA'))
VlnPlot(data5.seurat, features = c('nCount_RNA', 'nFeature_RNA'))
VlnPlot(data6.seurat, features = c('nCount_RNA', 'nFeature_RNA'))
VlnPlot(data7.seurat, features = c('nCount_RNA', 'nFeature_RNA'))
VlnPlot(data8.seurat, features = c('nCount_RNA', 'nFeature_RNA'))
VlnPlot(data9.seurat, features = c('nCount_RNA', 'nFeature_RNA'))
VlnPlot(data10.seurat, features = c('nCount_RNA', 'nFeature_RNA'))


data1.seurat <- subset(data1.seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200)
data2.seurat <- subset(data2.seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200)
data3.seurat <- subset(data3.seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200)
data4.seurat <- subset(data4.seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200)
data5.seurat <- subset(data5.seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200)
data6.seurat <- subset(data6.seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200)
data7.seurat <- subset(data7.seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200)
data8.seurat <- subset(data8.seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200)
data9.seurat <- subset(data9.seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200)
data10.seurat <- subset(data10.seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200)

