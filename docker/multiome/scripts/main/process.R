# Parse command-line arguments
library('argparse')
parser <- ArgumentParser(description='Normalize seurat objects')
parser$add_argument('--working-dir', dest='working_dir', type='character', help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser$add_argument('--script-dir', dest='script_dir', type='character', help='Directory containing workflow scripts', default='scripts')
parser$add_argument('--threads', dest='threads', type='integer', help='Number of threads to use for processing')
parser$add_argument('--seurat-object', dest='seurat_object', type='character', help='Seurat object for a dataset')
parser$add_argument('--output-seurat-object', dest='output_seurat_object', type='character', help='Output file to save Seurat object to')
args <- parser$parse_args()
script_dir <- args$script_dir

# Set working directory and load packages
setwd(args$working_dir)
source(paste0(script_dir, '/main/load_packages.r'))

# Set variables from args or snakemake parameters
threads <- if (is.null(args$threads)) snakemake@threads else args$threads
seurat_object <- if (is.null(args$seurat_object)) snakemake@input[['seurat_object']] else args$seurat_object
output_seurat_object <- if (is.null(args$output_seurat_object)) snakemake@output[['seurat_object']] else args$output_seurat_object

# Main
future::plan('multicore', workers=threads)
options(future.globals.maxSize=ngbs * 1000 * 1024^2)

assay <- 'RNA'

object <- readRDS(seurat_object)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

all.genes <- rownames(object)

object <- object %>%

    NormalizeData() %>%
    CleanVarGenes(nHVG=3000, use.sct=TRUE) %>%
    ScaleData(features=all.genes, verbose=FALSE) %>%
    CellCycleScoring(s.features=s.genes, g2m.features=g2m.genes) %>%

    DietSeurat(assays=assay)

saveRDS(object, output_seurat_object)
