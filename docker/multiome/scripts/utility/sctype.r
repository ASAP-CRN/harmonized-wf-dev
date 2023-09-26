
AnnotateSubtypes <- function(object, genes.list, assay='RNA', is.scaled=TRUE) {
    
    if (class(object) != 'Seurat') stop('error: please use a seurat object')

    # get cell-type by cell matrix
    counts <- object %>% GetAssayData(assay=assay, slot=ifelse(is.scaled, 'scale.data', 'data'))    
    es.max <- sctype_score(scRNAseqData=counts, scaled=is.scaled, gs=genes.list) 

    # merge by cluster
    cL_resutls <- do.call('rbind', lapply(unique(object@meta.data$seurat_clusters), function(cl) {
        
        es.max.cl <- sort(rowSums(es.max[, rownames(object@meta.data[object@meta.data$seurat_clusters == cl, ])]), decreasing=!0)
        head(data.frame(cluster=cl, type=names(es.max.cl), scores=es.max.cl, ncells=sum(object@meta.data$seurat_clusters == cl)), 10)

        }
    ))

    sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n=1, wt=scores)  

    # set low-confident (low ScType score) clusters to 'unknown'
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells / 4] <- 'Unknown'

    setDT(sctype_scores)
    setnames(sctype_scores, old='cluster', new='seurat_clusters')
    
    return(sctype_scores)
}
