# Parse command-line arguments
library('argparse')
library('pryr')
parser <- ArgumentParser(description='Run harmony analysis')
parser$add_argument('--working-dir', dest='working_dir', type='character', help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser$add_argument('--script-dir', dest='script_dir', type='character', help='Directory containing workflow scripts', default='scripts')
parser$add_argument('--threads', dest='threads', type='integer', help='Number of threads to use for processing')
parser$add_argument('--group-by-vars', dest='group_by_vars', nargs='+', type='character', help='Name of the column(s) to be used for batch correction', default='sample')
parser$add_argument('--seurat-objects-fofn', dest='seurat_objects_fofn', type='character', help='Newline-delimited paths to the set of input seurat objects (file-of-filenames)')
parser$add_argument('--output-seurat-object', dest='output_seurat_object', type='character', help='Output file to save Seurat object to')
# Accept positional arguments for compatibility with snakemake workflow
parser$add_argument('positional_arguments', nargs='*', type='character', help='Optionally proivde arguments positionally. Format: `harmony.R seurat_objects output_seurat_object threads`. seurat_objects is a space-separated set of input seurat objects.')
args <- parser$parse_args()
script_dir <- args$script_dir

## If positional arguments are used, ensure there are at least 3 (seurat_objects, output_seurat_object, threads)
if (length(args$positional_arguments) > 0 && length(args$positional_arguments) < 3) {
    cat("Unexpected number of positional arguments provided; positional arguments should be provided in the format:\n\tseurat_objectA seurat_objectB ... seurat_objectZ output_seurat_object threads\nArgs provided:\n\t")
    cat(args$positional_arguments)
    cat("\n")
    stop()
}

# Set working directory and load packages
setwd(args$working_dir)
source(paste0(script_dir, '/main/load_packages.r'))

# Set variables from args
threads <- if (length(args$positional_arguments) > 0) as.numeric(tail(args$positional_arguments, n=1)) else args$threads
group_by_vars <- if (length(args$positional_arguments) > 0) 'sample' else args$group_by_vars
seurat_objects <- if (length(args$positional_arguments) > 0) head(args$positional_arguments, -2) else scan(args$seurat_objects_fofn, what='character', sep='\n')
output_seurat_object <- if (length(args$positional_arguments) > 0) tail(args$positional_arguments, n=2) else args$output_seurat_object

# Main
# future::plan('multicore', workers=threads)
# options(future.globals.maxSize=ngbs * 1000 * 1024^2)

message("Loading object list")
print(mem_change(object.list <- future.apply::future_lapply(seurat_objects, readRDS)))
message("Total mem used: ", utils:::format.object_size(mem_used(), "auto"))

top.genes <- object.list %>% SelectIntegrationFeatures(nfeatures=5000)

message('Merging samples')
print(mem_change(object <- merge(x=object.list[[1]], y=object.list[-1])))
message("Total mem used: ", utils:::format.object_size(mem_used(), "auto"))

message("Deleting object list")
print(mem_change(rm(object.list)))
message("Total mem used: ", utils:::format.object_size(mem_used(), "auto"))

noise <- c('percent.mt', 'percent.rb', 'nFeature_RNA', 'nCount_RNA', 'doublet_scores', 'G2M.Score', 'S.Score')

message("Running data scaling, PCA, and harmony")
print(mem_change(
    object <- object %>%

        ScaleData(vars.to.regress=noise, features=top.genes, verbose=FALSE) %>%
        RunPCA(features=top.genes, npcs=50, verbose=FALSE) %>%

        harmony::RunHarmony(reduction='pca', group.by.vars=group_by_vars, max.iter=50, reduction.save='harmony')
))
message("Total mem used: ", utils:::format.object_size(mem_used(), "auto"))

saveRDS(object, output_seurat_object)
