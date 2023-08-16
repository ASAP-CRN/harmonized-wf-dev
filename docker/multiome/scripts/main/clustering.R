working_dir <- '/data/CARD_singlecell/harmony-rna/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=ngbs * 1000 * 1024^2)


object <- readRDS(snakemake@input[['seurat_object']]) %>% 
            FindClusters(resolution=snakemake@params[['resolution']], 
                            algorithm=snakemake@params[['algorithm']]) %>% 
            
            AnnotateClusters(genes=readRDS(snakemake@input[['markers']]), 
                                plot.out=snakemake@output[['plot']], cell.size=1.0)

saveRDS(object, snakemake@output[['seurat_object']])


