# Parse command-line arguments
library('argparse')
parser <- ArgumentParser(description='Call doublets')
parser$add_argument('--working-dir', dest='working_dir', type='character', help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser$add_argument('--script-dir', dest='script_dir', type='character', help='Directory containing workflow scripts', default='scripts')
parser$add_argument('--threads', dest='threads', type='integer', help='Number of threads to use for processing')
parser$add_argument('--seurat-objects', dest='seurat_objects', type='character', nargs='+', help='Set of input seurat objects for datasets')
parser$add_argument('--project-name', dest='project_name', type='character', help='Project name')
parser$add_argument('--output-metadata-file', dest='output_metadata_file', type='character', help='Output file to write metadata to')
args <- parser$parse_args()
script_dir <- args$script_dir

# Set working directory and load packages
setwd(args$working_dir)
source(paste0(script_dir, '/main/load_packages.r'))

# Set variables from args or snakemake parameters
threads <- if (is.null(args$threads)) snakemake@threads else args$threads
seurat_objects <- if (is.null(args$seurat_objects)) snakemake@input[['seurat_object']] else args$seurat_objects
project_name <- if (is.null(args$project_name)) snakemake@params[['project_name']] else args$project_name
output_metadata_file <- if (is.null(args$output_metadata_file)) snakemake@output[['metadata']] else args$output_metadata_file

# Main
future::plan('multicore', workers=threads)
options(future.globals.maxSize=ngbs * 1000 * 1024^2)

object.list <- future.apply::future_lapply(seurat_objects, readRDS)

m <- rbindlist(lapply(object.list, function(object) {
    m <- copy(object@meta.data); setDT(m, keep.rownames='cells')
}))

# cutoff <- NormalMixCutoff(mixtools::normalmixEM(m[, doublet_scores], k=2))
cutoff <- 0.06

m[, `:=` (
    project=rep(project_name, nrow(m)),
    predicted_gmm_doublets=fifelse(doublet_scores < cutoff, 'singlet', 'doublet')
)]

fwrite(x=m, file=output_metadata_file)
