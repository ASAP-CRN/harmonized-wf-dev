# Parse command-line arguments
library('argparse')
parser <- ArgumentParser(description='Preprocess')
parser$add_argument('--working-dir', dest='working_dir', type='character', help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser$add_argument('--script-dir', dest='script_dir', type='character', help='Directory containing workflow scripts', default='scripts')
parser$add_argument('--sample-id', dest='sample_id', type='character', help='Sample/dataset ID')
parser$add_argument('--batch', dest='batch', type='character', help="Batch from which the sample/dataset originated")
parser$add_argument('--project', dest='project', type='character', help="Project ID")
parser$add_argument('--raw-counts', dest='raw_counts', type='character', help='Unfiltered feature-barcode matrices HDF5 output by cellranger')
parser$add_argument('--filtered-counts', dest='filtered_counts', type='character', help='Filtered feature-barcode matrices HDF5 output by cellranger')
parser$add_argument('--soup-rate', dest='soup_rate', type='numeric', help='Dataset contamination rate fraction; used to remove mRNA contamination from the RNAseq data')
parser$add_argument('--output-seurat-object', dest='output_seurat_object', type='character', help='Output file to save Seurat object to')
args <- parser$parse_args()
script_dir <- args$script_dir

# Set working directory and load packages
setwd(args$working_dir)
source(paste0(script_dir, '/main/load_packages.r'))

# Set variables from args or snakemake parameters
dataset <- if (is.null(args$sample_id)) snakemake@params[['dataset']] else args$sample_id
batch <- if (is.null(args$batch)) fread(snakemake@input[['datasets']])[sample %chin% dataset, batch] else args$batch
project <- if (is.null(args$project)) NULL else gsub("-", "_", args$project)
batch_id <- gsub("-", "_", paste(project, batch, sep="_"))
raw_counts_file <- if (is.null(args$raw_counts)) paste0(data_path, batch, '/Multiome/', dataset, '/outs/raw_feature_bc_matrix.h5') else args$raw_counts
filtered_counts_file <- if (is.null(args$filtered_counts)) paste0(data_path, batch, '/Multiome/', dataset, '/outs/filtered_feature_bc_matrix.h5') else args$filtered_counts

soup_rate <- if (is.null(args$soup_rate)) snakemake@params[['soup_rate']] else args$soup_rate
output_seurat_object <- if (is.null(args$output_seurat_object)) snakemake@output[['seurat_object']] else args$output_seurat_object

# Main
raw.counts <- Read10X_h5(raw_counts_file)
filtered.counts <- Read10X_h5(filtered_counts_file)

adj.matrix <- suppressWarnings(SoupCorrect(raw.counts, filtered.counts, contamination_rate=soup_rate))
object <- CreateSeuratObject(adj.matrix, min.cells=0, min.features=0, project=dataset)

object[['percent.mt']] <- PercentageFeatureSet(object, pattern='^MT-')
object[['percent.rb']] <- PercentageFeatureSet(object, pattern='^RP[SL]')

m <- copy(object@meta.data)
setDT(m, keep.rownames='cells')

m[,
    `:=` (
        sample=dataset,
        batch=batch,
        project=project,
        batch_id=batch_id
        )
]

batch <- m[, batch]
sample <- m[, sample]
project <- m[, project]
batch_id <- m[, batch_id]

names(sample) <- names(batch) <- names(project) <- names(batch_id) <- m[, cells]

doublet_rate <- (ncol(object) / 1000) * 0.008

object <- object %>%

    AddMetaData(metadata=factor(batch), col.name='batch') %>%
    AddMetaData(metadata=factor(sample), col.name='sample') %>%
    AddMetaData(metadata=factor(project), col.name='project') %>%
    AddMetaData(metadata=factor(batch_id), col.name='batch_id') %>%

    scrublet(n_prin_comps=30, expected_doublet_rate=doublet_rate)

saveRDS(object, output_seurat_object)
