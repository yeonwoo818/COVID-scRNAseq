# Cell Type Annotation (two references)
DefaultAssay(merged.seurat) <- "RNA"
merged.seurat.v3 <- merged.seurat
merged.seurat.v3[["RNA"]] <- as(object = merged.seurat.v3[["RNA"]], Class = "Assay")
expr <- GetAssayData(merged.seurat.v3, layer = "data")

# ref 1: MouseRNAseqData
ref.mouse <- celldex::MouseRNAseqData()
pred.mouse <- SingleR(
  test = expr,
  ref = ref.mouse,
  labels = ref.mouse$label.main,
  BPPARAM = MulticoreParam(60)
)
merged.seurat$SingleR_MouseRNAseq <- pred.mouse$labels

# ref 2: ImmGenData
ref.immgen <- celldex::ImmGenData()
pred.immgen <- SingleR(
  test = expr,
  ref = ref.immgen,
  labels = ref.immgen$label.main,
  BPPARAM = MulticoreParam(60)
)
merged.seurat$SingleR_ImmGen <- pred.immgen$labels

table(merged.seurat$SingleR_MouseRNAseq)
table(merged.seurat$SingleR_ImmGen)

p_mouse <- DimPlot(merged.seurat, reduction = "umap", group.by = "SingleR_MouseRNAseq", label = TRUE)
p_immgen <- DimPlot(merged.seurat, reduction = "umap", group.by = "SingleR_ImmGen", label = TRUE)
p_origin2 <- DimPlot(merged.seurat, reduction = "umap", group.by = "Origin")
