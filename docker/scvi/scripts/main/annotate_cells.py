# TODO:  implement a function that takes in a scvi model and adata and annotates the clusters
#     refer to utily/sctype.r
import os
import argparse
import pandas as pd
import numpy as np
import scanpy
import scvi


parser = argparse.ArgumentParser(description="Annotate clusters")
parser.add_argument(
    "--adata-input", dest="adata_input", type=str, help="AnnData object for a dataset"
)
parser.add_argument(
    "--marker-genes",
    dest="marker_genes",
    type=str,
    default="resources/celltype_marker_table.csv",
    help="Path to marker_genes .csv file",
)
parser.add_argument(
    "--batch-key",
    dest="batch_key",
    type=str,
    help="Key in AnnData object for batch information",
)
parser.add_argument(
    "--output-cell-types-file",
    dest="cell_type_output",
    type=str,
    help="Output file to write celltypes to",
)
parser.add_argument(
    "--adata-output",
    dest="adata_output",
    type=str,
    help="Output file to save AnnData object to",
)
parser.add_argument(
    "--output-metadata-file",
    dest="output_metadata_file",
    type=str,
    help="Output file to write metadata to",
)


args = parser.parse_args()


# Set CPUs to use for parallel computing
scanpy._settings.ScanpyConfig.n_jobs = -1

# 0. load adata
adata = scanpy.read_h5ad(args.adata_input)  # type: ignore

#  1. load marker_genes
markers = pd.read_csv(args.marker_genes, index_col=0)
# defensive
markers = markers[~markers.index.duplicated(keep="first")].rename_axis(index=None)

# TODO: add a check to make sure the marker genes are in the adata
# output now many are missing...

#  2. copy for cellassign
# bdata = adata[:, markers.index].copy() #
bdata = adata[:, adata.var.index.isin(markers.index)].copy()

#  3. get size_factor and noise
lib_size = bdata.X.sum(1)  # type: ignore
bdata.obs["size_factor"] = lib_size / np.mean(lib_size)
# noise = ['doublet_score', 'pct_counts_mt', 'pct_counts_rb'] #, 'S.Score', 'G2M.Score']

#  4. model = CellAssign(bdata, marker_genes)
scvi.external.CellAssign.setup_anndata(
    bdata,
    size_factor_key="size_factor",
    batch_key=args.batch_key,
    layer="counts",
    # continuous_covariate_keys=noise
)

#  5. model.train()
model = scvi.external.CellAssign(bdata, markers)
plan_args = {"lr_factor": 0.05, "lr_patience": 20, "reduce_lr_on_plateau": True}
model.train(
    max_epochs=1000,
    accelerator="gpu",
    early_stopping=True,
    plan_kwargs=plan_args,
    early_stopping_patience=40,
)

#  6. model.predict()
bdata.obs["cellassign_types"] = model.predict().idxmax(axis=1).values

# 7. transfer cell_type to adata
adata.obs["cell_type"] = bdata.obs["cellassign_types"]

#  8. save model & artificts
predictions = (
    bdata.obs[["sample", "cellassign_types"]]
    .reset_index()
    .rename(columns={"index": "cells"})
)
predictions.to_csv(
    args.cell_type_output, index=False
)  # # pred_file = "cellassign_predictions.csv"

# 9. write_h5ad
adata.write_h5ad(filename=args.adata_output, compression="gzip")

# 10. save metadata
adata.obs.to_csv(
    args.output_metadata_file, index=True
)  # metadata_file = "metadata.csv"
