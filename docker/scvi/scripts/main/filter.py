# import muon.pp.filter_obs as filter_obs
import muon
import scanpy
import argparse
import pandas as pd
import pandas as pd
import sys

sys.path.append("/opt/scripts/utility")
from helpers import update_validation_metrics

parser = argparse.ArgumentParser(description="Filter")
parser.add_argument(
    "--adata-input", dest="adata_input", type=str, help="AnnData object for a dataset"
)
parser.add_argument(
    "--adata-output",
    dest="adata_output",
    type=str,
    help="Output file to save AnnData object to",
)
# TODO: add filter parameters as arguments
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

# TODO: make these cutoffs arguments...
# muon api is better than the scanpyt api for this...
muon.pp.filter_obs(adata, "pct_counts_mt", lambda x: x <= 10)
muon.pp.filter_obs(adata, "doublet_score", lambda x: x < 0.2)
muon.pp.filter_obs(adata, "total_counts", lambda x: (x >= 500) & (x <= 100000))
muon.pp.filter_obs(adata, "n_genes_by_counts", lambda x: (x >= 300) & (x <= 10000))

# save the filtered adata
adata.write_h5ad(filename=args.adata_output, compression="gzip")


#######  validation metrics
val_metrics = pd.read_csv(args.output_validation_file)

output_metrics = update_validation_metrics(adata, "filter", val_metrics)

# log the validation metrics
output_metrics.to_csv(args.output_validation_file, index=True)
