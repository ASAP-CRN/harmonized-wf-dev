working_dir <- '/data/CARD_singlecell/harmony-rna/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=ngbs * 1000 * 1024^2)


m <- fread(snakemake@input[['metadata']])

noise <- c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.rb', 'doublet_scores')


plot.list <- future.apply::future_lapply(noise, function(feature) {
    
    m %>%

        ggplot(aes(x=snakemake@params[['project_name']], y=get(feature))) + 
        
        geom_violin(fill='steelblue', color='black') + theme_bw() + 
        
        theme(axis.text.x=element_blank(), axis.title.x=element_blank()) + ylab(feature)
        
})

p <- ggpubr::ggarrange(plotlist=plot.list, legend='none', align='hv', ncol=3, nrow=3)

ggsave(plot=p, width=15, height=8, filename=snakemake@output[['plot_1']])


p <- m %>% ggplot(aes(x=get(noise[1]), y=get(noise[2]))) + geom_point(alpha=0.3) +
    theme_bw() + xlab('Number of UMIs') + ylab('Number of Genes')


ggsave(plot=p, width=18, height=9, filename=snakemake@output[['plot_2']])







