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
...
data10.seurat <- CreateSeuratObject(data10.data, project = 'Test', min.cells=3, min.features=200)


# Preprocessing Seurat Objects
VlnPlot(data1.seurat, features = c('nCount_RNA', 'nFeature_RNA'))
VlnPlot(data2.seurat, features = c('nCount_RNA', 'nFeature_RNA'))
...
VlnPlot(data10.seurat, features = c('nCount_RNA', 'nFeature_RNA'))


data1.seurat <- subset(data1.seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200)
data2.seurat <- subset(data2.seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200)
...
data10.seurat <- subset(data10.seurat, subset = nCount_RNA > 500 & nFeature_RNA > 200)

