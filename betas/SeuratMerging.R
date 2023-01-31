

srr.list <- c("SRR8325947", "SRR8325948", "SRR8325949", "SRR8325963", "SRR8325964", "SRR8325965", "SRR8325979", "SRR8325980", "SRR8325981", "SRR8326004", "SRR8326005", "SRR8325995", "SRR8325996", "SRR8325997")

baseDir <- "~/Documents/data/PASseqData/SRR832_backup/"

seu.list <- c()
for(i in srr.list){
  print(i)
  seu <- readRDS(paste(baseDir, i, "Solo.out/Gene/filtered/", i, "_seurat.rds", sep=""))
  seu.list <- c(seu.list, seu) 
}


seu <- readRDS("~/Documents/data/PASseqData/SRR832_backup/SRR8326006Solo.out/Gene/filtered/SRR8326006_seurat.rds")

srr.list <- c("SRR8326006", srr.list)

library(Seurat)
merged.seu <- merge(seu, y = seu.list, add.cell.ids = srr.list, project = "LCdata")

merged.seu <- NormalizeData(merged.seu, normalization.method = "LogNormalize", scale.factor = 10000)

merged.seu <- FindVariableFeatures(merged.seu, selection.method = "vst", nfeatures = 2000)
#Scale data
merged.seu <- ScaleData(merged.seu, features = rownames(merged.seu))
#Linear dimension reduction - PCA
merged.seu <- RunPCA(merged.seu, features = VariableFeatures(object = merged.seu))

ElbowPlot(merged.seu)

#Finding the neighboring structure
merged.seu <- FindNeighbors(merged.seu, dims = 1:15)
#Cluster cells
merged.seu <- FindClusters(merged.seu, resolution = 0.5)
#Non linear dimension reduction
merged.seu <- RunUMAP(merged.seu, dims = 1:10)

DimPlot(merged.seu, reduction = "umap", label=T)
DimPlot(merged.seu, reduction = "umap",group.by = "orig.ident", label=T)



lcmarkers <- c("ARAF","BIRC6","EED","CD274","GPC5","MALAT1","RBM10","KRT19","CEACAM5","ENO2","CA12", "KISS1R", "LYPD3", "SLC7A11","TMPRSS3")

tcellmarkers <- c("PTPRC","CD4","CD8A","CXCR3","CCR4","CCR7")

