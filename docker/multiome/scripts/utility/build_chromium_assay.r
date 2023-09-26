

# hub <- AnnotationHub()
# query(hub, c('EnsDb', 'sapiens', '98'))
# enDB <- hub[['AH75011', force=TRUE]]

AddChromiumAssay <- function(object, input.data, enDB, frag.file) {
    
    atac_counts <- input.data[['Peaks']]
    grange.counts <- StringToGRanges(rownames(atac_counts), sep=c(':', '-'))
    grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac_counts <- atac_counts[as.vector(grange.use), ]
    
    suppressWarnings(annotations <- GetGRangesFromEnsDb(ensdb=enDB))
    
    genome(annotations) <- 'hg38'
    annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), 'UCSC'))
    

    chrom_assay <- CreateChromatinAssay(
        counts=atac_counts,
        sep=c(':', '-'),
        genome='hg38',
        fragments=frag.file,
        min.cells=0,
        annotation=annotations
    )
    
    object[['ATAC']] <- chrom_assay

    return(object)
}

