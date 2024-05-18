import argparse
from scvi.external import SCAR
from scanpy import read_10x_h5


parser = argparse.ArgumentParser(description="Preprocess")
parser.add_argument(
    "--working-dir",
    dest="working_dir",
    type=str,
    help="Working directory",
    default="/data/CARD_singlecell/harmony-rna/",
)
parser.add_argument(
    "--script-dir",
    dest="script_dir",
    type=str,
    help="Directory containing workflow scripts",
    default="scripts",
)
parser.add_argument("--sample-id", dest="sample_id", type=str, help="Sample/dataset ID")
parser.add_argument(
    "--batch",
    dest="batch",
    type=str,
    help="Batch from which the sample/dataset originated",
)
parser.add_argument("--project", dest="project", type=str, help="Project ID")
parser.add_argument(
    "--raw-counts",
    dest="raw_counts",
    type=str,
    help="Unfiltered feature-barcode matrices HDF5 output by cellranger",
)
parser.add_argument(
    "--filtered-counts",
    dest="filtered_counts",
    type=str,
    help="Filtered feature-barcode matrices HDF5 output by cellranger",
)
parser.add_argument(
    "--ambient-p",
    dest="soup_rate",
    type=float,
    help="Dataset contamination rate fraction; used to remove mRNA contamination from the RNAseq data",
)
parser.add_argument(
    "--adata-output",
    dest="adata_output",
    type=str,
    help="Output file to save AnnData object to",
)

args = parser.parse_args()


# use scar instead of cellbender
raw_adata = read_10x_h5(args.raw_counts)
raw_adata.var_names_make_unique()

filt_adata = read_10x_h5(args.filtered_anndata)
filt_adata.var_names_make_unique()

SCAR.setup_anndata(filt_adata)
SCAR.get_ambient_profile(adata=filt_adata, raw_adata=raw_adata, prob=0.995)
scar_vae = SCAR(filt_adata)
scar_vae.train()
filt_adata.obsm["X_scAR"] = scar_vae.get_latent_representation()
filt_adata.layers["denoised"] = scar_vae.get_denoised_counts()

# add metadata
adata.obs["sample"] = args.sample_id
adata.obs["batch"] = args.batch
adata.obs["project"] = args.project
adata.obs["batch_id"] = args.project + args.batch  #

filt_adata.write_h5ad(filename=snakemake.output.anndata, compression="gzip")  # type: ignore
