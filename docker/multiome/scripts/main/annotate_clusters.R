# Parse command-line arguments
library('argparse')
parser <- ArgumentParser(description='Annotate clusters')
parser$add_argument('--working-dir', dest='working_dir', type='character', help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser$add_argument('--script-dir', dest='script_dir', type='character', help='Directory containing workflow scripts', default='scripts')
parser$add_argument('--threads', dest='threads', type='integer', help='Number of threads to use for processing')
parser$add_argument('--seurat-object', dest='seurat_object', type='character', help='Seurat object for a dataset')
parser$add_argument('--cell-type-markers-list', dest='cell_type_markers_list', type='character', help='Seurat object containing a list of major cell type markers')
parser$add_argument('--output-metadata-file', dest='output_metadata_file', type='character', help='Output file to write metadata to')
args <- parser$parse_args()
script_dir <- args$script_dir

# Set working directory and load packages
setwd(args$working_dir)
source(paste0(script_dir, '/main/load_packages.r'))
source('https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R')

# Set variables from args or snakemake parameters
threads <- if (is.null(args$threads)) snakemake@threads else args$threads
seurat_object <- if (is.null(args$seurat_object)) snakemake@input[['seurat_object']] else args$seurat_object
cell_type_markers_list <- if (is.null(args$cell_type_markers_list)) snakemake@input[['markers']] else args$cell_type_markers_list
output_metadata_file <- if (is.null(args$output_metadata_file)) snakemake@output[['metadata']] else args$output_metadata_file

# Main
# Using future when running with 75 samples was causing the error:
## Error: Failed to retrieve the result of MulticoreFuture (future_lapply-6) from the forked worker (on localhost; PID 201). Post-mortem diagnostic: No process exists with this PID, i.e. the forked localhost worker is no longer alive. The total size of the 15 globals exported is 2.31 GiB. The three largest globals are ‘object’ (2.30 GiB of class ‘S4’), ‘split.cells’ (7.23 MiB of class ‘list’) and ‘features’ (2.37 MiB of class ‘character’)
# future::plan('multicore', workers=threads)
# options(future.globals.maxSize=ngbs * 1000 * 1024^2)

markers <- readRDS(cell_type_markers_list)
object <- readRDS(seurat_object)

all.genes <- rownames(object)

scores <- object %>%
            ScaleData(verbose=FALSE, features=unlist(markers)) %>%
            AnnotateSubtypes(genes.list=markers)

umap <- copy(as.data.frame(object[['umap']]@cell.embeddings))
setDT(umap, keep.rownames='cells')

metadata <- copy(object@meta.data)
setDT(metadata, keep.rownames='cells')

metadata <- metadata[scores[, .(seurat_clusters, type)], on='seurat_clusters'][umap, on='cells']

fwrite(metadata, output_metadata_file)
