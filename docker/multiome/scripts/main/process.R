working_dir <- '/data/CARD_singlecell/harmony-rna/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=ngbs * 1000 * 1024^2)

assay <- 'RNA'

object <- readRDS(snakemake@input[['seurat_object']])


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


all.genes <- rownames(object)

object <- object %>% 

    NormalizeData() %>% 
    CleanVarGenes(nHVG=3000, use.sct=TRUE) %>% 
    ScaleData(features=all.genes, verbose=FALSE) %>% 
    CellCycleScoring(s.features=s.genes, g2m.features=g2m.genes) %>% 

    DietSeurat(assays=assay)


saveRDS(object, snakemake@output[['seurat_object']])