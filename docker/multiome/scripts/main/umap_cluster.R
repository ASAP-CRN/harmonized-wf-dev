working_dir <- '/data/CARD_singlecell/scvi-rna/'; setwd(working_dir)

source('scripts/main/load_packages.r')

object <- readRDS(snakemake@input[['seurat_object']])

DefaultAssay(object) <- 'integrated'

object <- object %>% 
    RunUMAP(dims=1:30, reduction='pca') %>%
    FindNeighbors(dims=1:30, reduction='pca') %>% 
    FindClusters(resolution=0.8, algorithm=3)

saveRDS(object, snakemake@output[['seurat_object']])