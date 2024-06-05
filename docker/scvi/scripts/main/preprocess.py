import argparse
import sys
import scanpy


parser = argparse.ArgumentParser(description="Preprocess")
parser.add_argument(
    "--adata-input", dest="adata_input", type=str, help="AnnData object for a dataset"
)
# toplevel metadata to add to the adata
parser.add_argument("--sample-id", dest="sample_id", type=str, help="Sample/dataset ID")
parser.add_argument(
    "--batch",
    dest="batch",
    type=str,
    help="Batch from which the sample/dataset originated",
)
parser.add_argument("--project", dest="project", type=str, help="Project ID")
parser.add_argument(
    "--adata-output",
    dest="adata_output",
    type=str,
    help="Output file to save AnnData object to",
)

args = parser.parse_args()


sys.path.append("/opt/scripts/utility")
from helpers import anndata_from_h5

# Set CPUs to use for parallel computing
scanpy._settings.ScanpyConfig.n_jobs = -1

# load the data from cellbender output
adata = anndata_from_h5(args.adata_input)
# adata = scanpy.read_10x_h5(args.adata_input)

adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["rb"] = adata.var_names.str.startswith(("RPL", "RPS"))
scanpy.pp.calculate_qc_metrics(
    adata, qc_vars=["rb", "mt"], percent_top=None, log1p=False, inplace=True
)

# add doublet score
scanpy.external.pp.scrublet(adata)
# keep the metadata clean ad drop predicted_doublet... rely on the probability of dublet i.e. doublet_score
adata.obs.drop("predicted_doublet", axis=1, inplace=True)

# add metadata
adata.obs["sample"] = args.sample_id
adata.obs["batch"] = args.batch
adata.obs["project"] = args.project
# adata.obs['batch_id'] = args.project+args.batch
adata.obs["batch_id"] = f"{args.project}_{args.batch}"  #

# drop the celbender obs? 'background_fraction', 'cell_probability', 'cell_size', 'droplet_efficiency',

adata.write_h5ad(filename=args.adata_output)
