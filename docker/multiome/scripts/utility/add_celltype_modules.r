
    
AnnotateClusters <- function(object, genes, plot.out=getwd(), h=9, w=13, assay=DefaultAssay(object), cell.size=0.1) {

    colors <- c('lightgrey', 'navy')
    object <- object %>% AddModuleScore(features=genes, assay=assay, name=names(genes))

    p <- object %>% 
        FeaturePlot(features=paste0(names(genes), seq(length(genes))), pt.size=cell.size, 
        min.cutoff='q10', max.cutoff='q90', ncol=3, cols=colors, reduction='umap')
    
    ggsave(plot=p, height=h, width=w, filename=plot.out)
    
    return(object)
    
}
    
    