import argparse
import scanpy
import anndata
import pandas as pd
import sys

sys.path.append("/opt/scripts/utility")
from helpers import get_validation_metrics

parser = argparse.ArgumentParser(description="Call doublets")
parser.add_argument(
    "--adata-objects-fofn",
    dest="adata_objects_fofn",
    type=str,
    help="Newline-delimited TSV of sample names and paths to the set of input adata objects (file-of-filenames)",
)
parser.add_argument(
    "--adata-output",
    dest="adata_output",
    type=str,
    help="Output file to save AnnData object to",
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

scanpy.settings.verbosity = 1
scanpy.settings.figdir = "plots/"
scanpy.settings.set_figure_params(
    dpi=100, fontsize=10, dpi_save=300, format="png", figsize=("12", "8")
)  # type: ignore

metrics = [
    "n_genes_by_counts",
    "total_counts",
    "pct_counts_mt",
    "pct_counts_rb",
    "doublet_score",
]
scanpy.settings.verbosity = 1
scanpy.settings.figdir = "plots/"
scanpy.settings.set_figure_params(
    dpi=100, fontsize=10, dpi_save=300, format="png", figsize=("12", "8")
)  # type: ignore

adatas = {}
top_genes = {}


#  note that the sample id should be the official ASAP_samples
with open(args.adata_objects_fofn, "r") as file:
    for sample in file:
        # Check that sample line is not empty
        if sample.strip():
            columns = sample.strip().split("\t")
            sample_name, file_path = columns
            raw = scanpy.read_h5ad(file_path)
            # code below if memory issues with concatenating all the ge
            # adata = raw.copy()
            # scanpy.pp.normalize_total(adata, target_sum=1e4)
            # scanpy.pp.log1p(adata)
            # scanpy.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=8000)

            # # ranked_genes = adata.var[adata.var.highly_variable].dispersions_norm.argsort().to_dict()
            # for k,v in ranked_genes.items():
            #     if k in top_genes:
            #         top_genes[k] += v
            #     else:
            #         top_genes[k] = v

            # raw = minify_adata(raw)
            adatas[sample_name] = raw

# we could subset to the top_genes here before concat if we have memory issues (e.g. whole dataset harmonization.)
adata = anndata.concat(merge="same", uns_merge="same", index_unique="_", adatas=adatas)


for metric in metrics:  # type: ignore
    scanpy.pl.violin(adata, keys=metric, size=0, save="".join("_" + metric))

# top_genes = pd.DataFrame(index=top_genes.keys(), columns=['rank'], data=top_genes.values())
# TODO: export top_genes? and plots

# export concatenated data.
adata.write_h5ad(filename=args.adata_output, compression="gzip")


#######  validation metrics
val_metrics = get_validation_metrics(adata, "concatenation")
# log the validation metrics

val_metrics.to_csv(args.output_validation_file, index=True)
