set.seed(1234)

### Load the libraries
library(Seurat)
library(tidyverse)
library(SingleR)
library(celldex)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)

### Load the data
data <- read.table("data/GSE132257_GEO_processed_protocol_and_fresh_frozen_raw_UMI_count_matrix.txt" , header = TRUE , sep = "\t" , stringsAsFactors = FALSE )
data <- data %>% column_to_rownames("Index")

### Create the Seurat object
seu.obj <- CreateSeuratObject(counts = data, project = "GSE132257")

### filtering our data
seu.obj <- subset(seu.obj ,  orig.ident %in% c("SMC171.T.SGI" , "SMC171.T.BEL" , "SMC13T.A1" , "SMC13T.A1.F") )
seu.obj@meta.data$orig.ident <- droplevels(seu.obj@meta.data$orig.ident)

### preprocessing
seu.obj$mt.percent <- PercentageFeatureSet(seu.obj, pattern = "^MT-")
FeatureScatter(seu.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")
VlnPlot(seu.obj, features = c("nFeature_RNA", "nCount_RNA", "mt.percent"), ncol = 3)
seu.obj <- subset(seu.obj , subset = mt.percent < 20 & nFeature_RNA > 1000 & nFeature_RNA < 5000 & nCount_RNA > 1000 )

### QC
seu.obj <- NormalizeData(seu.obj)
seu.obj <- FindVariableFeatures(seu.obj)
seu.obj <- ScaleData(seu.obj)
seu.obj <- RunPCA(seu.obj, npcs = 30)
ElbowPlot(seu.obj)
seu.obj <- FindNeighbors(seu.obj, dims = 1:20)
seu.obj <- FindClusters(seu.obj, resolution = .3)
seu.obj <- RunUMAP(seu.obj, dims = 1:20)
DimPlot(seu.obj, label = TRUE , reduction = "umap" , group.by = "RNA_snn_res.0.3" , label.size = 5)
pca_embeddings <- Embeddings(seu.obj, reduction = "pca")
duplicate_cells <- which(duplicated(pca_embeddings))
if (length(duplicate_cells) > 0) {
  seu.obj <- subset(seu.obj, cells = setdiff(Cells(seu.obj), Cells(seu.obj)[duplicate_cells]))
}
seu.obj <- RunTSNE(seu.obj, dims = 1:20 )
DimPlot(seu.obj, label = TRUE , reduction = "tsne" , group.by = "RNA_snn_res.0.3" , label.size = 5)
markers <- FindAllMarkers(seu.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1 , test.use = "wilcox")
topmarkers <- markers %>% group_by(cluster) %>% top_n(1, avg_log2FC)
DotPlot(seu.obj , features = topmarkers$gene , group.by = "RNA_snn_res.0.3" , dot.scale = 1.5) +
  RotatedAxis()

### SingleR
hpca <- celldex::HumanPrimaryCellAtlasData()
countsData <- GetAssayData(seu.obj , slot = "counts")
singleR.results <- SingleR(test = countsData, ref = hpca, labels = hpca$label.main)
table <- table(assigned = singleR.results$labels ,clusters = seu.obj$seurat_clusters)
pheatmap(log10(10 + table) , color = colorRampPalette(c("white" , "blue"))(10))

### Assign the labels
seu.obj$SingleR.labels <- singleR.results$labels[match(rownames(seu.obj@meta.data) , rownames(singleR.results))]
Idents(seu.obj) <- seu.obj@meta.data$SingleR.labels

