working_dir <- '/data/CARD_singlecell/harmony-rna/'; setwd(working_dir)

source('scripts/main/load_packages.r')
source('https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R')

future::plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=ngbs * 1000 * 1024^2)


markers <- readRDS(snakemake@input[['markers']])
object <- readRDS(snakemake@input[['seurat_object']])


all.genes <- rownames(object)

scores <- object %>% 
            ScaleData(verbose=FALSE, features=all.genes) %>% 
            AnnotateSubtypes(genes.list=markers)

umap <- copy(as.data.frame(object[['umap']]@cell.embeddings))
setDT(umap, keep.rownames='cells')

metadata <- copy(object@meta.data)
setDT(metadata, keep.rownames='cells')

metadata <- metadata[scores[, .(seurat_clusters, type)], on='seurat_clusters'][umap, on='cells']

fwrite(metadata, snakemake@output[['metadata']])
