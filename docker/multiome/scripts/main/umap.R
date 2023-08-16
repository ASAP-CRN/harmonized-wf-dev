working_dir <- '/data/CARD_singlecell/harmony-rna/'; setwd(working_dir)

source('scripts/main/load_packages.r')

object <- readRDS(snakemake@input[['seurat_object']]) %>% 
            RunUMAP(dims=1:50, reduction='harmony')

saveRDS(object, snakemake@output[['seurat_object']])

