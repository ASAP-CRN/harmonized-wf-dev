# Parse command-line arguments
library('argparse')
parser <- ArgumentParser(description='Plot features')
parser$add_argument('--working-dir', dest='working_dir', type='character', help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser$add_argument('--metadata', dest='metadata', type='character', help='Metadata file output by sctype (annotate_clusters)')
parser$add_argument('--feature', dest='feature', type='character', help='Feature to plot umaps for')
parser$add_argument('--output-feature-umap-plot', dest='output_feature_umap_plot', type='character', help='Output file to write the feature umap plot to')
args <- parser$parse_args()

# Set working directory and load packages
setwd(args$working_dir)
invisible(lapply(c('data.table', 'ggplot2', 'dplyr'), require, character.only=TRUE))

# Set variables from args or snakemake parameters
feature <- if (is.null(args$feature)) snakemake@params[['features']] else args$feature
metadata <- fread(if (is.null(args$metadata)) snakemake@input[['metadata']] else args$metadata)
output_feature_umap_plot <- if (is.null(args$output_feature_umap_plot)) snakemake@output[['plot']] else args$output_feature_umap_plot

# Main
setkey(metadata, cells)

middle <- metadata[, mean(c(max(get(feature)), min(get(feature))))]
pal <- wesanderson::wes_palette('Zissou1', type='discrete')[c(1, 3, 5)]

g <- metadata %>% ggplot(aes(x=UMAP_1, y=UMAP_2)) +

        geom_point(color=metadata[, get(feature)], size=0.5, alpha=0.8) +
        scale_color_gradient2(low=pal[1], mid=pal[2], high=pal[3], midpoint=middle) +
        labs(color=metadata[, get(feature)]) + theme_classic() + ggtitle('')

ggsave(plot=g, width=13, height=9, filename=output_feature_umap_plot)
