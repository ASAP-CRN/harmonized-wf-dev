# Parse command-line arguments
library('argparse')
parser <- ArgumentParser(description='Preprocess')
parser$add_argument('--working-dir', dest='working_dir', type='character', help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser$add_argument('--script-dir', dest='script_dir', type='character', help='Directory containing workflow scripts', default='scripts')
parser$add_argument('--dataset', dest='dataset', type='character', help="Name of dataset to process")
parser$add_argument('--batch-file', dest='batch_file', type='character')
parser$add_argument('--raw-counts', dest='raw_counts', type='character')
parser$add_argument('--filtered-counts', dest='filtered_counts', type='character')
parser$add_argument('--soup-rate', dest='soup_rate', type='character')
parser$add_argument('--seurat-object', dest='seurat_object', type='character', help='Output file to save Seurat object to')
args <- parser$parse_args()

setwd(args$working_dir)

source(paste0(args$script_dir, '/main/load_packages.r'))

reticulate::source_python(paste0(args$script_dir, '/utility/scrublet_py.py'))

dataset <- if (is.null(args$dataset)) snakemake@params[['dataset']] else args$dataset
batch_file <- if (is.null(args$batch_file)) snakemake@input[['datasets']] else args$batch_file

batch <- fread(batch_file)[sample %chin% dataset, batch]

raw_counts_file <- if (is.null(args$raw_counts)) paste0(data_path, batch, '/Multiome/', dataset, '/outs/raw_feature_bc_matrix.h5') else args$raw_counts
filtered_counts_file <- if (is.null(args$filtered_counts)) paste0(data_path, batch, '/Multiome/', dataset, '/outs/filtered_feature_bc_matrix.h5') else args$filtered_counts

raw.counts <- Read10X_h5(raw_counts_file)
filtered.counts <- Read10X_h5(filtered_counts_file)

soup_rate <- if (is.null(args$soup_rate)) snakemake@params[['soup_rate']] else args$soup_rate
adj.matrix <- suppressWarnings(SoupCorrect(raw.counts, filtered.counts, contamination_rate=soup_rate))
object <- CreateSeuratObject(adj.matrix, min.cells=0, min.features=0, project=dataset)


object[['percent.mt']] <- PercentageFeatureSet(object, pattern='^MT-')
object[['percent.rb']] <- PercentageFeatureSet(object, pattern='^RP[SL]')

m <- copy(object@meta.data)
setDT(m, keep.rownames='cells')

m[,
    `:=` (
        sample=dataset,
        batch=batch
        )
]

batch <- m[, batch]
sample <- m[, sample]

names(sample) <- names(batch) <- m[, cells]


doublet_rate <- (ncol(object) / 1000) * 0.008

object <- object %>% 
    
    AddMetaData(metadata=factor(batch), col.name='batch') %>% 
    AddMetaData(metadata=factor(sample), col.name='sample') %>% 
    
    scrublet(n_prin_comps=30, expected_doublet_rate=doublet_rate) 


seurat_object <- if (is.null(args$seurat_object)) snakemake@output[['seurat_object']] else args$seurat_object
saveRDS(object, seurat_object)
