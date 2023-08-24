
scrublet <- function(object, min_counts=2, min_cells=3, expected_doublet_rate=0.06, min_gene_variability_pctl=85, n_prin_comps=50, sim_doublet_ratio=2, n_neighbors=NULL) {

    if (class(object) != 'Seurat') stop('object is not of type "Seurat"') 

    reticulate::source_python('/data/ShernData/scripts_NF1/scrublet_py.py')

    X <- as(t(as.matrix(GetAssayData(object=object, slot='counts'))), 'TsparseMatrix')
    
    val <- X@x
    i <- as.integer(X@i)
    j <- as.integer(X@j)
    dim <- as.integer(X@Dim)

    if (is.null(n_neighbors)) n_neighbors <- round(0.5 * sqrt(nrow(X)))

    scrublet_py_args <- c(
        list(
                i=i, j=j, val=val, dim=dim,
                expected_doublet_rate=expected_doublet_rate, min_counts=min_counts, 
                min_cells=min_cells, min_gene_variability_pctl=min_gene_variability_pctl, 
                n_prin_comps=n_prin_comps, sim_doublet_ratio=sim_doublet_ratio, n_neighbors=n_neighbors
            )
        )

    scrublet_res <- do.call(scrublet_py, scrublet_py_args)
    names(scrublet_res) <- c('doublet_scores', 'predicted_doublets')

    d <- data.table(cells=rownames(X), doublet_scores=as.numeric(scrublet_res$doublet_scores))

    doublet_scores <- d[, doublet_scores]
    names(doublet_scores) <- d[, cells]


    object <- object %>% AddMetaData(metadata=doublet_scores, col.name='doublet_scores')
                    
    return(object)
  
}

