# Parse command-line arguments
library('argparse')
parser <- ArgumentParser(description='Run harmony analysis')
parser$add_argument('--working-dir', dest='working_dir', type='character', help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser$add_argument('--script-dir', dest='script_dir', type='character', help='Directory containing workflow scripts', default='scripts')
parser$add_argument('--threads', dest='threads', type='integer', help='Number of threads to use for processing')
parser$add_argument('--seurat-objects', dest='seurat_objects', type='character', nargs='+', help='Set of input seurat objects for datasets')
parser$add_argument('--output-seurat-object', dest='output_seurat_object', type='character', help='Output file to save Seurat object to')
# Accept positional arguments for compatibility with snakemake workflow
parser$add_argument('positional_arguments', nargs='*', type='character', help='Optionally proivde arguments positionally. Format: `harmony.R seurat_objects output_seurat_object threads`. seurat_objects is a space-separated set of input seurat objects.')
args <- parser$parse_args()

## If positional arguments are used, ensure there are at least 3 (seurat_objects, output_seurat_object, threads)
if (length(args$positional_arguments) > 0 && length(args$positional_arguments) < 3) {
    cat("Unexpected number of positional arguments provided; positional arguments should be provided in the format:\n\tseurat_objectA seurat_objectB ... seurat_objectZ output_seurat_object threads\nArgs provided:\n\t")
    cat(args$positional_arguments)
    cat("\n")
    stop()
}

# Set working directory and load packages
setwd(args$working_dir)
source(paste0(args$script_dir, '/main/load_packages.r'))

# Set variables from args
threads <- if (length(args$positional_arguments) > 0) as.numeric(tail(args$positional_arguments, n=1)) else args$threads
seurat_objects <- if (length(args$positional_arguments) > 0) head(args$positional_arguments, -2) else args$seurat_objects
output_seurat_object <- if (length(args$positional_arguments) > 0) tail(args$positional_arguments, n=2) else args$output_seurat_object

# Main
future::plan('multicore', workers=threads)
options(future.globals.maxSize=ngbs * 1000 * 1024^2)

object.list <- future.apply::future_lapply(seurat_objects, readRDS)

top.genes <- object.list %>% SelectIntegrationFeatures(nfeatures=5000)

message('merging samples...')
object <- merge(x=object.list[[1]], y=object.list[-1])
message('Done.')

batch <- 'sample'; noise <- c('percent.mt', 'percent.rb', 'nFeature_RNA', 'nCount_RNA', 'doublet_scores', 'G2M.Score', 'S.Score')

object <- object %>%

    ScaleData(vars.to.regress=noise, features=top.genes, verbose=FALSE) %>%
    RunPCA(features=top.genes, npcs=50, verbose=FALSE) %>%

    harmony::RunHarmony(reduction='pca', group.by.vars=batch, max.iter.harmony=50, reduction.save='harmony')

saveRDS(object, output_seurat_object)[1])
