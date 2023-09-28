

FilterCells <- function(object, min.genes=300, max.genes=7500, max.mito=70, min.umi=500, max.umi=50000, filter_doublets=TRUE) {

    metadata <- copy(object@meta.data)
    setDT(metadata, keep.rownames='cells')

    new_names <- c('cells', 'laneID', 'nUMI', 'nGene')
    setnames(metadata, old=colnames(metadata)[1:4], new=new_names)

    metadata <- na.omit(metadata)


    if (filter_doublets) {
        mix <- mixtools::normalmixEM(x=metadata[, doublet_scores], k=2)
        metadata <- metadata[doublet_scores < round(NormalMixCutoff(mix), 2)]   
    }

    mito <- metadata[, .(cellNo=.I[percent.mt < max.mito])]
    UMIs <- metadata[, .(cellNo=.I[between(nUMI, min.umi, max.umi)])]
    genes <- metadata[, .(cellNo=.I[between(nGene, min.genes, max.genes)])]

    
    cells2keep <- metadata[fintersect(fintersect(UMIs, genes), mito)[, cellNo], cells]

    object <- object %>% subset(cells=cells2keep)

    return(object)

}
