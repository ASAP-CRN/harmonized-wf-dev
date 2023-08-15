# Parse command-line arguments
library('argparse')
parser <- ArgumentParser(description='Filter')
parser$add_argument('--working-dir', dest='working_dir', type='character', help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser$add_argument('--script-dir', dest='script_dir', type='character', help='Directory containing workflow scripts', default='scripts')
parser$add_argument('--seurat-object', dest='seurat_object', type='character', help='Seurat object for a dataset')
parser$add_argument('--metadata', dest='metadata', type='character', help='Metadata file output by gmm_doublet_calling')
parser$add_argument('--output-seurat-object', dest='output_seurat_object', type='character', help='Output file to save Seurat object to')
args <- parser$parse_args()

# Set working directory and load packages
setwd(args$working_dir)
source(paste0(args$script_dir, '/main/load_packages.r'))

# Set variables from args or snakemake parameters
seurat_object <- if (is.null(args$seurat_object)) snakemake@input[['seurat_object']] else args$seurat_object
metadata <- if (is.null(args$metadata)) snakemake@input[['metadata']] else args$metadata
output_seurat_object <- if (is.null(args$output_seurat_object)) snakemake@output[['seurat_object']] else args$output_seurat_object

# Main
object <- readRDS(seurat_object)

object <- object %>% subset(
    subset=
            nCount_RNA %between% c(500, 2e5) &
            nFeature_RNA %between% c(500, 12500) &
            percent.mt < 0.5 & percent.rb < 0.5,
    cells=
            fread(metadata)[
                sample %chin% Project(object) &
                predicted_gmm_doublets %chin% 'singlet', cells
            ]
    )

saveRDS(object, output_seurat_object)
