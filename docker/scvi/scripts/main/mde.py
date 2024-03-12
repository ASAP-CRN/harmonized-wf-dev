import argparse
import scanpy
from math import sqrt
from statistics import median
from scvi.model.utils import mde


parser = argparse.ArgumentParser(
    description='Minimum-Distortion Embedding (MDE)'
)
parser.add_argument(
    '--latent-key',
    dest='latent_key',
    type=str,
    default='X_scvi',
    help='Latent key to save the scvi latent to'
)
parser.add_argument(
    '--adata-input',
    dest='adata_input',
    type=str,
    help='AnnData object for a dataset'
)
parser.add_argument(
    '--adata-output',
    dest='adata_output',
    type=str,
    help='Output file to save AnnData object to'
)

args = parser.parse_args()


adata = scanpy.read_h5ad(args.adata_input)


# k = round(median(sqrt(adata.obs.sample.value_counts() * 0.5)))

# scanpy.pp.neighbors(adata, n_neighbors=50, use_rep=snakemake.params.latent_key)
# scanpy.tl.leiden(adata, resolution=0.4)
# scanpy.tl.umap(adata)
mde(adata.obsm[args.latent_key])

adata.write_h5ad(filename=args.adata_output, compression='gzip')
