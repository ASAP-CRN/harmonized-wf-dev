# Parse command-line arguments
library('argparse')
parser <- ArgumentParser(description='Call doublets')
parser$add_argument('--working-dir', dest='working_dir', type='character', help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser$add_argument('--script-dir', dest='script_dir', type='character', help='Directory containing workflow scripts', default='scripts')
parser$add_argument('--threads', dest='threads', type='integer', help='Number of threads to use for processing')
parser$add_argument('--seurat-objects-fofn', dest='seurat_objects_fofn', type='character', help='Newline-delimited paths to the set of input seurat objects (file-of-filenames)')
parser$add_argument('--project-name', dest='project_name', type='character', help='Project name')
parser$add_argument('--output-metadata-file', dest='output_metadata_file', type='character', help='Output file to write metadata to')
args <- parser$parse_args()
script_dir <- args$script_dir

# Set working directory and load packages
setwd(args$working_dir)
source(paste0(script_dir, '/main/load_packages.r'))

# Set variables from args or snakemake parameters
threads <- if (is.null(args$threads)) snakemake@threads else args$threads
seurat_objects <- if (is.null(args$seurat_objects_fofn)) snakemake@input[['seurat_object']] else scan(args$seurat_objects_fofn, what='character', sep='\n')
project_name <- if (is.null(args$project_name)) snakemake@params[['project_name']] else args$project_name
output_metadata_file <- if (is.null(args$output_metadata_file)) snakemake@output[['metadata']] else args$output_metadata_file

# Main
# Using future when running with 75 samples was causing the error:
## Error: Failed to retrieve the result of MulticoreFuture (future_lapply-1) from the forked worker (on localhost; PID 496). Post-mortem diagnostic: No process exists with this PID, i.e. the forked localhost worker is no longer alive. The total size of the 5 globals exported is 30.95 KiB. The three largest globals are ‘...future.FUN’ (22.34 KiB of class ‘function’), ‘...future.elements_ii’ (8.61 KiB of class ‘list’) and ‘future.call.arguments’ (0 bytes of class ‘list’)
# future::plan('multicore', workers=threads)
# options(future.globals.maxSize=ngbs * 1000 * 1024^2)
# object.list <- future.apply::future_lapply(seurat_objects, readRDS)

object.list <- list()
for (seurat_object in seurat_objects) {
    object.list <- append(object.list, readRDS(seurat_object))
}

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
