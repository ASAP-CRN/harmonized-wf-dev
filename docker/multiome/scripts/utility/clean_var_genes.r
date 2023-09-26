
CleanVarGenes <- function(object, nHVG=3000, assay='RNA', use.sct=TRUE) {
    
    if (class(object) != 'Seurat') stop('please provide a Seurat object...')

    DefaultAssay(object) <- assay
    
    g <- data.table(genes=rownames(object))
    
    if (use.sct) {
        genes <- object %>% 
            subset(features=g[!(genes %like% '^MT' | genes %like% '^RP[LS]'), genes]) %>% 
            SCTransform(variable.features.n=nHVG, verbose=FALSE) %>% VariableFeatures()
    } else {
        genes <- object %>% 
            subset(features=g[!(genes %like% '^MT' | genes %like% '^RP[LS]'), genes]) %>% 
            NormalizeData() %>% FindVariableFeatures(nfeatures=nHVG) %>% VariableFeatures()
    }       
    
    VariableFeatures(object, assay=assay) <- genes

    DefaultAssay(object) <- assay

    return(object)
}
