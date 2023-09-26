

get_annotations <- function(genome='hg38', hub.id=NULL) {

    require(AnnotationHub)
    
    if (is.null(hub.id)) hub.id <- 'AH75011'

    hub <- AnnotationHub::AnnotationHub()
    reference <- hub[[hub.id]]

    annotations <- suppressWarnings(GetGRangesFromEnsDb(ensdb=reference))

    genome(annotations) <- genome
    annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), 'UCSC'))

    return(annotations)

}