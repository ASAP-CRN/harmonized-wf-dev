# Parse command-line arguments
library('argparse')
parser <- ArgumentParser(description='Plot groups')
parser$add_argument('--working-dir', dest='working_dir', type='character', help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser$add_argument('--metadata', dest='metadata', type='character', help='Metadata file output by sctype (annotate_clusters)')
parser$add_argument('--group', dest='group', type='character', help='Group to plot umaps for')
parser$add_argument('--output-group-umap-plot-prefix', dest='output_group_umap_plot_prefix', type='character', help='Output file prefix to write the group umap plot to; will write both a PDF and a PNG.')
args <- parser$parse_args()

# Set working directory and load packages
setwd(args$working_dir)
invisible(lapply(c('data.table', 'ggplot2', 'dplyr'), require, character.only=TRUE))

# Set variables from args or snakemake parameters
group <- if (is.null(args$group)) snakemake@params[['groups']] else args$group
metadata <- fread(if (is.null(args$metadata)) snakemake@input[['metadata']] else args$metadata)
output_group_umap_plot_prefix <- if (is.null(args$output_group_umap_plot_prefix)) gsub("\\.pdf$", "", snakemake@output[['plot']]) else args$output_group_umap_plot_prefix

# Main
setkey(metadata, cells)

colors <- sample(colors(), nrow(metadata[, .N, get(group)]))

g <- metadata %>% ggplot(aes(x=UMAP_1, y=UMAP_2)) + theme_classic() +

        geom_point(aes(color=as.factor(metadata[, get(group)])), alpha=0.8, size=0.1) +
        scale_color_manual(values=colors) + labs(color=metadata[, get(group)]) +
        guides(color=guide_legend(override.aes=list(size=5)))

# Save plot as PDF
ggsave(plot=g, width=13, height=9, device="pdf", filename=paste0(output_group_umap_plot_prefix, ".pdf"))

# Save plot as PNG
ggsave(plot=g, width=13, height=9, device="png", filename=paste0(output_group_umap_plot_prefix, ".png"))
