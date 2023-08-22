working_dir <- '/data/CARD_singlecell/harmony-rna/'; setwd(working_dir)

invisible(lapply(c('data.table', 'ggplot2', 'dplyr'), require, character.only=TRUE))


group <- snakemake@params[['groups']]
m <- fread(snakemake@input[['metadata']])

setkey(m, cells)

colors <- sample(colors(), nrow(m[, .N, get(group)]))


g <- m %>% ggplot(aes(x=UMAP_1, y=UMAP_2)) + theme_classic() +
        
        geom_point(aes(color=m[, get(group)]), alpha=0.8, size=0.1) + 
        scale_color_manual(values=colors) + labs(color=m[, get(group)]) +
        guides(color=guide_legend(override.aes=list(size=5)))


ggsave(plot=g, width=13, height=9, filename=snakemake@output[['plot']])

