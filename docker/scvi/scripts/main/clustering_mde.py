import argparse
import scanpy
from scvi.model.utils import mde


parser = argparse.ArgumentParser(
    description='Run clustering'
)
parser.add_argument(
    '--working-dir',
    dest='working_dir',
    type=str,
    help='Working directory',
    default='/data/CARD_singlecell/harmony-rna/'
)
parser.add_argument(
    '--script-dir',
    dest='script_dir',
    type=str,
    help='Directory containing workflow scripts',
    default='scripts'
)
parser.add_argument(
    '--threads',
    dest='threads',
    type=int,
    help='Number of threads to use for processing'
)
parser.add_argument(
    '--adata-object',
    dest='adata_object',
    type=str,
    help='AnnData object for a dataset'
)
parser.add_argument(
    '--clustering-algorithm',
    dest='clustering_algorithm',
    type=int,
    help='Clustering algorithm to use'
)
parser.add_argument(
    '--clustering-resolution',
    dest='clustering_resolution',
    type=float,
    help='Clustering resolution'
)
parser.add_argument(
    '--cell-type-markers-list',
    dest='cell_type_markers_list',
    type=str,
    help='Seurat object containing a list of major cell type markers'
)
parser.add_argument(
    '--output-cell-type-plot-prefix',
    dest='output_cell_type_plot_prefix',
    type=str,
    help='Prefix for output file to save major cell type umap plot in PDF and PNG formats'
)
parser.add_argument(
    '--adata-output',
    dest='adata_output',
    type=str,
    help='Output file to save AnnData object to'
)
parser.add_argument(
    '--latent-key',
    dest='latent_key',
    type=str,
    default='X_scvi',
    help='Latent key to save the scvi latent to'
)

args = parser.parse_args()

adata = scanpy.read_h5ad(args.adata_input) # type: ignore

# TODO: impliment clustering


# TODO - write_h5ad option compression='gzip' is giving an error
#adata.write_h5ad(filename=args.adata_output, compression='gzip')
adata.write_h5ad(filename=args.adata_output)
