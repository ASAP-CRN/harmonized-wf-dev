import argparse
import scanpy
import leidenalg


parser = argparse.ArgumentParser(description="Annotate clusters")
parser.add_argument(
    "--adata-input", dest="adata_input", type=str, help="AnnData object for a dataset"
)
parser.add_argument(
    "--adata-output",
    dest="adata_output",
    type=str,
    help="Output file to save AnnData object to",
)
parser.add_argument(
    "--latent-key",
    dest="latent_key",
    type=str,
    default="X_scvi",
    help="Latent key to the scvi latent",
)

args = parser.parse_args()


# Set CPUs to use for parallel computing
scanpy._settings.ScanpyConfig.n_jobs = -1

adata = scanpy.read_h5ad(args.adata_input)  # type: ignore

# calculate neighbor graph on scVI latent
scanpy.pp.neighbors(adata, n_neighbors=30, use_rep=args.latent_key)

# do leiden
scanpy.tl.leiden(adata, resolution=0.4)

scanpy.tl.umap(adata)

# Export figures
# scanpy.pl.umap(adata, color=, size=0, save='_major_cell_type.png')
# scanpy.pl.umap(adata, color=, size=0, save='_major_cell_type.pdf')


adata.write_h5ad(filename=args.adata_output, compression="gzip")
