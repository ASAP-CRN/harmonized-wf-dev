# Parse command-line arguments
library('argparse')
parser <- ArgumentParser(description='Run clustering')
parser$add_argument('--working-dir', dest='working_dir', type='character', help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser$add_argument('--script-dir', dest='script_dir', type='character', help='Directory containing workflow scripts', default='scripts')
parser$add_argument('--threads', dest='threads', type='integer', help='Number of threads to use for processing')
parser$add_argument('--seurat-object', dest='seurat_object', type='character', help='Seurat object for a dataset')
parser$add_argument('--clustering-algorithm', dest='clustering_algorithm', type='integer', help='Clustering algorithm to use')
parser$add_argument('--clustering-resolution', dest='clustering_resolution', type='numeric', help='Clustering resolution')
parser$add_argument('--cell-type-markers-list', dest='cell_type_markers_list', type='character', help='Seurat object containing a list of major cell type markers')
parser$add_argument('--output-cell-type-plot', dest='output_cell_type_plot', type='character', help='Output file to save major cell type umap plot in PDF format')
parser$add_argument('--output-seurat-object', dest='output_seurat_object', type='character', help='Output file to save Seurat object to')
args <- parser$parse_args()

# Set working directory and load packages
setwd(args$working_dir)
source(paste0(args$script_dir, '/main/load_packages.r'))

# Set variables from args or snakemake parameters
threads <- if (is.null(args$threads)) snakemake@threads else args$threads
seurat_object <- if (is.null(args$seurat_object)) snakemake@input[['seurat_object']] else args$seurat_object
clustering_algorithm <- if (is.null(args$clustering_algorithm)) snakemake@params[['algorithm']] else args$clustering_algorithm
clustering_resolution <- if (is.null(args$clustering_resolution)) snakemake@params[['resolution']] else args$clustering_resolution
cell_type_markers_list <- if (is.null(args$cell_type_markers_list)) snakemake@input[['markers']] else args$cell_type_markers_list
output_cell_type_plot <- if (is.null(args$output_cell_type_plot)) snakemake@output[['plot']] else args$output_cell_type_plot
output_seurat_object <- if (is.null(args$output_seurat_object)) snakemake@output[['seurat_object']] else args$output_seurat_object

# Main
future::plan('multicore', workers=threads)
options(future.globals.maxSize=ngbs * 1000 * 1024^2)

object <- readRDS(seurat_object) %>%
            FindClusters(resolution=clustering_resolution,
                            algorithm=clustering_algorithm) %>%

            AnnotateClusters(genes=readRDS(cell_type_markers_list),
                                plot.out=output_cell_type_plot, cell.size=1.0)

saveRDS(object, output_seurat_object)
