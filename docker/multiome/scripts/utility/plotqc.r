

PlotQC <- function(object.list, outpath='plots/', project.name='') {

    m <- rbindlist(lapply(object.list, function(object) {
        m <- copy(object@meta.data); setDT(m, keep.rownames='cells')
    }))

    m[, project := rep(project.name, nrow(m))]

    noise <- c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.rb', 'doublet_scores')


    plot.list <- lapply(noise, function(feature) {
        m %>% ggplot(aes(x=project, y=get(feature))) + geom_violin(fill='steelblue', color='black') + 
        theme_bw() + theme(axis.text.x=element_blank(), axis.title.x=element_blank()) + ylab(feature)
    })

    p <- ggpubr::ggarrange(plotlist=plot.list, legend='none', align='hv', ncol=3, nrow=2)

    ggsave(plot=p, width=14, height=7, filename=paste0(outpath, 'qc_plot1.pdf'))

    
    plot1 <- m %>% ggplot(aes(x=get(noise[1]), y=get(noise[2]))) + geom_point(alpha=0.3) +
        theme_bw() + xlab('Number of UMIs - RNA') + ylab('Number of Genes - RNA')

    ggsave(plot=plot1, width=12, height=9, filename=paste0(outpath, 'qc_plot2.pdf'))
    

    return(object.list)

}



