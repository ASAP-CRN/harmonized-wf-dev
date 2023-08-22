working_dir <- '/data/CARD_singlecell/harmony-rna/'; setwd(working_dir)

invisible(lapply(c('data.table', 'ggplot2', 'dplyr'), require, character.only=TRUE))


feature <- snakemake@params[['features']]
m <- fread(snakemake@input[['metadata']])

setkey(m, cells)

middle <- m[, mean(c(max(get(feature)), min(get(feature))))]
pal <- wesanderson::wes_palette('Zissou1', type='discrete')[c(1, 3, 5)]

g <- m %>% ggplot(aes(x=UMAP_1, y=UMAP_2)) + 

        geom_point(color=m[, get(feature)], size=0.5, alpha=0.8) +
        scale_color_gradient2(low=pal[1], mid=pal[2], high=pal[3], midpoint=middle) +
        labs(color=m[, get(feature)]) + theme_classic() + ggtitle('')


ggsave(plot=g, width=13, height=9, filename=snakemake@output[['plot']])

