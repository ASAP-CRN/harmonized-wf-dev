AnnotateClusters <- function(object, genes, plot.out=getwd(), h=9, w=13, assay=DefaultAssay(object), cell.size=0.1) {

    colors <- c('lightgrey', 'navy')
    object <- object %>% AddModuleScore(features=genes, assay=assay, name=names(genes))

    p <- object %>%
        FeaturePlot(features=paste0(names(genes), seq(length(genes))), pt.size=cell.size,
        min.cutoff='q10', max.cutoff='q90', ncol=3, cols=colors, reduction='umap')

    # Save plot as PDF
    ggsave(plot=p, height=h, width=w, device="pdf", filename=paste0(plot.out, ".pdf"))

    # Save plot as PNG
    ggsave(plot=p, height=h, width=w, device="png", filename=paste0(plot.out, ".png"))

    return(object)
}
