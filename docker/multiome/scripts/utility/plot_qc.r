

PlotQC <- function(object.list, project.name='', outpath='', filter.doublets=TRUE) {

    m <- rbindlist(lapply(object.list, function(object) {
        m <- copy(object@meta.data); setDT(m, keep.rownames='cells')
    }))

    cutoff <- NormalMixCutoff(mixtools::normalmixEM(m[, doublet_scores], k=2))

    m[, `:=` (
        project=rep(project.name, nrow(m)),
        predicted_gmm_doublets=fifelse(doublet_scores < cutoff, 'singlet', 'doublet')
    )]

    fwrite(x=m, file='output/unfiltered_metadata.csv')
    

    object.list <- sapply(object.list, function(object) {
        object <- object %>% subset(
            cells=
                m[orig.ident %in% Project(object) & 
                predicted_gmm_doublets %chin% 'singlet', cells]
            )
        
        project <- m[, project]; names(project) <- m[, cells]
        object %>% AddMetaData(metadata=factor(project), col.name='project')

    }, simplify=FALSE)



    m <- rbindlist(lapply(object.list, function(object) {
        m <- copy(object@meta.data); setDT(m, keep.rownames='cells')
    }))

    noise <- c('nCount_ATAC', 'nCount_RNA', 'nFeature_ATAC', 'nFeature_RNA', 'percent.mt', 'percent.rb', 'doublet_scores', 'nucleosome_signal', 'TSS.enrichment')


    plot.list <- lapply(noise, function(feature) {
        m %>% ggplot(aes(x=project, y=get(feature))) + geom_violin(fill='steelblue', color='black') + 
        theme_bw() + theme(axis.text.x=element_blank(), axis.title.x=element_blank()) + ylab(feature)
    })

    p <- ggpubr::ggarrange(plotlist=plot.list, legend='none', align='hv', ncol=3, nrow=3)

    ggsave(plot=p, width=15, height=8, filename=paste0(outpath, 'qc_plot1.pdf'))

    
    plot1 <- m %>% ggplot(aes(x=get(noise[2]), y=get(noise[4]))) + geom_point(alpha=0.3) 
        theme_bw() + xlab('Number of UMIs - RNA') + ylab('Number of Genes - RNA')
    
    plot2 <- m %>% ggplot(aes(x=get(noise[1]), y=get(noise[3]))) + geom_point(alpha=0.3) 
        theme_bw() + xlab('Number of UMIs - ATAC') + ylab('Number of Genes - ATAC')
    

    ggsave(plot=plot1 + plot2, width=18, height=9, filename=paste0(outpath, 'qc_plot2.pdf'))
    

    return(object.list)

}



