import argparse
import scanpy
import sys
import pandas as pd

sys.path.append("/opt/scripts/utility")
from helpers import score_cell_cycle, update_validation_metrics


parser = argparse.ArgumentParser(description="Normalize seurat objects")

parser.add_argument(
    "--adata-input", dest="adata_input", type=str, help="AnnData object for a dataset"
)
parser.add_argument(
    "--batch-key",
    dest="batch_key",
    type=str,
    help="Key in AnnData object for batch information",
)
parser.add_argument(
    "--adata-output",
    dest="adata_output",
    type=str,
    help="Output file to save AnnData object to",
)
parser.add_argument(
    "--n-top-genes",
    dest="n_top_genes",
    type=int,
    help="number of HVG genes to keep",
    default=3000,
)
parser.add_argument(
    "--marker-genes",
    dest="marker_genes",
    type=str,
    default="resources/celltype_marker_table.csv",
    help="Path to marker_genes .csv file",
)
parser.add_argument(
    "--output-validation-file",
    dest="output_validation_file",
    type=str,
    help="Output file to write validation metrics to",
)


args = parser.parse_args()


# Set CPUs to use for parallel computing
scanpy._settings.ScanpyConfig.n_jobs = -1

adata = scanpy.read_h5ad(args.adata_input)  # type: ignore


# does this work with sparse uint8?
adata.layers["counts"] = adata.X.copy()  # type: ignore

scanpy.pp.normalize_total(adata, target_sum=1e4)
scanpy.pp.log1p(adata)

# get cell cycle scores
score_cell_cycle(adata, organism="human")

# TODO: load the top_genes from the qc plotting step and subset...
# top_genes = 'top_genes.csv'
#  if we have memory issues, consider subsetting to a gene_list

### MAKE SURE MARKER GENES ARE KEPT
# #  1. load marker_genes
#  1. load marker_genes
markers = pd.read_csv(args.marker_genes, index_col=0)
# defensive
markers = markers[~markers.index.duplicated(keep="first")].rename_axis(index=None)

batch_key = args.batch_key
# WARNING: using 'sample' can cause loess to fail in the highly_variable_genes function
# HACK: using 'batch_id' instead of 'sample' for now
if batch_key == "sample":
    print("WARNING: using 'batch_id' instead of 'sample' for now")
    batch_key = "batch_id"

hvgs_full = scanpy.pp.highly_variable_genes(
    adata,
    n_top_genes=20_000,
    batch_key=batch_key,
    flavor="seurat_v3",
    check_values=True,
    layer="counts",
    subset=False,
    inplace=False,
)

# TODO: double check that this works...
# hvgs_full.loc[markers.index, 'highly_variable_rank'] = 1.
hvgs_full.loc[markers.index, "highly_variable_nbatches"] = (
    max(hvgs_full["highly_variable_nbatches"]) + 1
)

# Sort genes by how often they selected as hvg within each batch and
# break ties with median rank of residual variance across batches
hvgs_full.sort_values(
    ["highly_variable_nbatches", "highly_variable_rank"],
    ascending=[False, True],
    na_position="last",
    inplace=True,
)


hvgs_full = hvgs_full.iloc[: args.n_top_genes].index.to_list()

# # sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
# scanpy.pp.highly_variable_genes(
#     adata,
#     batch_key=args.batch_key,
#     subset=True,
#     flavor='seurat_v3',
#     layer='counts',
#     n_top_genes=args.n_top_genes
# )

# load the raw test_ad
adata = adata[:, hvgs_full]

scanpy.pp.pca(adata, n_comps=30)

adata.write_h5ad(filename=args.adata_output, compression="gzip")

#######  validation metrics
val_metrics = pd.read_csv(args.output_validation_file)

output_metrics = update_validation_metrics(adata, "feature_selection", val_metrics)

# log the validation metrics
output_metrics.to_csv(args.output_validation_file, index=True)
