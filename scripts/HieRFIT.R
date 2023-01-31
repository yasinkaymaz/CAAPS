suppressPackageStartupMessages(library(HieRFIT))
suppressPackageStartupMessages(library(tidyverse))
args <- commandArgs(TRUE)
hiermod_rds <- args[1]
seu_rds <- args[2]
ct_out <- args[3]

hiermod <- readRDS(hiermod_rds)
seu <- readRDS(seu_rds)

#Project the cell class labels on the new dataset:
hier_obj <- HieRFIT(Query = seu[["RNA"]]@data, refMod = hiermod)

write.table(data.frame(cb=rownames(hier_obj@Evaluation), cls=hier_obj@Evaluation[, 'Projection']),  # nolint
            file = ct_out,
            sep = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)
