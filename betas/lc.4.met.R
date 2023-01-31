lc.4.met.srr <- c("SRR8325963", "SRR8325964", "SRR8325965")

seu1 <- readRDS("~/Documents/data/PASseqData/SRR832_backup/SRR8325963Solo.out/Gene/filtered/SRR8325963_seurat.rds")
seu.list <- c()
for(i in lc.4.met.srr[-1]){
  print(i)
  seu <- readRDS(paste(baseDir, i, "Solo.out/Gene/filtered/", i, "_seurat.rds", sep=""))
  seu.list <- c(seu.list, seu) 
}



merged.seu4 <- merge(seu1, y = seu.list, add.cell.ids = lc.4.met.srr, project = "LC.4")

merged.seu4 <- NormalizeData(merged.seu4, normalization.method = "LogNormalize", scale.factor = 10000)

merged.seu4 <- FindVariableFeatures(merged.seu4, selection.method = "vst", nfeatures = 2000)
#Scale data
merged.seu4 <- ScaleData(merged.seu4, features = rownames(merged.seu4))
#Linear dimension reduction - PCA
merged.seu4 <- RunPCA(merged.seu4, features = VariableFeatures(object = merged.seu4))

ElbowPlot(merged.seu4)

#Finding the neighboring structure
merged.seu4 <- FindNeighbors(merged.seu4, dims = 1:10)
#Cluster cells
merged.seu4 <- FindClusters(merged.seu4, resolution = .4)
#Non linear dimension reduction
merged.seu4 <- RunUMAP(merged.seu4, dims = 1:5)

DimPlot(merged.seu4, reduction = "umap", label=T)
DimPlot(merged.seu4, reduction = "umap",group.by = "orig.ident", label=T)

VlnPlot(merged.seu4, features = tcellmarkers)
