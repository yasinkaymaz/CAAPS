suppressPackageStartupMessages(library(Seurat,verbose=F))
args = commandArgs(T)
RunDir = args[1]
RunName = args[2]
#Load count data
scdata <- Read10X(data.dir=RunDir)
seu <- CreateSeuratObject(counts=scdata, project=RunName, min.cells = 3, min.features = 200)
#pre-filtration
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
#Normalization
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
#Variable features
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
#Scale data
seu <- ScaleData(seu, features = rownames(seu))
#Linear dimension reduction - PCA
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
#Finding the neighboring structure
seu <- FindNeighbors(seu, dims = 1:15)
#Cluster cells
seu <- FindClusters(seu, resolution = 0.5)
#Non linear dimension reduction
seu <- RunUMAP(seu, dims = 1:15)
#save the seurat object as an rds file.
saveRDS(seu, file = paste(RunDir, "/", RunName, "_seurat.rds", sep=""))

write.table(data.frame(cb=rownames(seu@meta.data), cls=seu@meta.data[,'seurat_clusters']), file=paste(RunDir, "/", RunName, "_seurat_clusters.txt", sep=""), sep="\t", quote = F, col.names = F, row.names = F)



