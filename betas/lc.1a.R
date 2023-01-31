baseDir <- "~/Documents/data/PASseqData/SRR832_backup/"

lc.1a.srr <- c("SRR8325947", "SRR8325948", "SRR8325949")

sratable <- readxl::read_excel("~/Dropbox/EGE_University/Grants/Tubitak-1002-2021/SonuçRaporu/SraRunTable.xlsx")

rownames(sratable) <- sratable$Run

sratable[lc.1a.srr,]

seu1 <- readRDS("~/Documents/data/PASseqData/SRR832_backup/SRR8325947Solo.out/Gene/filtered/SRR8325947_seurat.rds")
seu.list <- c()
for(i in lc.1a.srr[-1]){
  print(i)
  seu <- readRDS(paste(baseDir, i, "Solo.out/Gene/filtered/", i, "_seurat.rds", sep=""))
  seu.list <- c(seu.list, seu) 
}



library(Seurat)
merged.seu <- merge(seu1, y = seu.list, add.cell.ids = lc.1a.srr, project = "LC.1a")

merged.seu <- NormalizeData(merged.seu, normalization.method = "LogNormalize", scale.factor = 10000)

merged.seu <- FindVariableFeatures(merged.seu, selection.method = "vst", nfeatures = 2000)
#Scale data
merged.seu <- ScaleData(merged.seu, features = rownames(merged.seu))
#Linear dimension reduction - PCA
merged.seu <- RunPCA(merged.seu, features = VariableFeatures(object = merged.seu))

ElbowPlot(merged.seu)

#Finding the neighboring structure
merged.seu <- FindNeighbors(merged.seu, dims = 1:10)
#Cluster cells
merged.seu <- FindClusters(merged.seu, resolution = .4)
#Non linear dimension reduction
merged.seu <- RunUMAP(merged.seu, dims = 1:5)

DimPlot(merged.seu, reduction = "umap", label=T)
DimPlot(merged.seu, reduction = "umap",group.by = "orig.ident", label=T)



lcmarkers <- c("ARAF","BIRC6","EED","CD274","GPC5","MALAT1","RBM10","KRT19","CEACAM5","ENO2","CA12", "KISS1R", "LYPD3", "SLC7A11","TMPRSS3")

tcellmarkers <- c("PTPRC","CD4","CD8A","CXCR3","CCR4","CCR7")

FeaturePlot(merged.seu, features = tcellmarkers)
VlnPlot(merged.seu, features = tcellmarkers)
merged.seu[["percent.lcmarkers"]] <- PercentageFeatureSet(merged.seu, features = lcmarkers)
merged.seu[["percent.tcellmarkers"]] <- PercentageFeatureSet(merged.seu, features = tcellmarkers)
FeaturePlot(merged.seu, features = "percent.lcmarkers")

VlnPlot(merged.seu, features = "percent.lcmarkers")
FeaturePlot(merged.seu, features = c("ABHD14B","PRDX5", "BCL2L11", "PLAT", "ATP6V1D", "RUVBL1"))

VlnPlot(merged.seu, features = c("ABHD14B","PRDX5", "BCL2L11", "PLAT", "ATP6V1D", "RUVBL1"))
