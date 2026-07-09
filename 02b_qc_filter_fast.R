# Calculating Mitochondrial Percentage
data1.seurat[["percent.mt"]] <- PercentageFeatureSet(data1.seurat, pattern = "^mt-")
data2.seurat[["percent.mt"]] <- PercentageFeatureSet(data2.seurat, pattern = "^mt-")
...
data10.seurat[["percent.mt"]] <- PercentageFeatureSet(data10.seurat, pattern = "^mt-")


# Upper Bound: doublet(nFeature) + damaged cell(mito)
# nFeature UB
data1.singlet <- subset(data1.seurat, subset = nFeature_RNA < 5000 & percent.mt < 10)
data2.singlet <- subset(data2.seurat, subset = nFeature_RNA < 5000 & percent.mt < 10)
...
data10.singlet <- subset(data10.seurat, subset = nFeature_RNA < 5000 & percent.mt < 10)


# # nFeature UB: each sample median+3*MAD (log10)
# nf.ub <- function(obj) 10^(median(log10(obj$nFeature_RNA)) + 3 * mad(log10(obj$nFeature_RNA)))
# 
# data1.singlet <- subset(data1.seurat, subset = nFeature_RNA < nf.ub(data1.seurat) & percent.mt < 10)
# data2.singlet <- subset(data2.seurat, subset = nFeature_RNA < nf.ub(data2.seurat) & percent.mt < 10)
# data3.singlet <- subset(data3.seurat, subset = nFeature_RNA < nf.ub(data3.seurat) & percent.mt < 10)
# data4.singlet <- subset(data4.seurat, subset = nFeature_RNA < nf.ub(data4.seurat) & percent.mt < 10)
# data5.singlet <- subset(data5.seurat, subset = nFeature_RNA < nf.ub(data5.seurat) & percent.mt < 10)
# data6.singlet <- subset(data6.seurat, subset = nFeature_RNA < nf.ub(data6.seurat) & percent.mt < 10)
# data7.singlet <- subset(data7.seurat, subset = nFeature_RNA < nf.ub(data7.seurat) & percent.mt < 10)
# data8.singlet <- subset(data8.seurat, subset = nFeature_RNA < nf.ub(data8.seurat) & percent.mt < 10)
# data9.singlet <- subset(data9.seurat, subset = nFeature_RNA < nf.ub(data9.seurat) & percent.mt < 10)
# data10.singlet <- subset(data10.seurat, subset = nFeature_RNA < nf.ub(data10.seurat) & percent.mt < 10)
# 
# rm(nf.ub)


