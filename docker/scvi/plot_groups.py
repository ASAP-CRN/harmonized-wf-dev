import scanpy
import anndata
import argparse
# Parse command-line arguments
parser = argparse.ArgumentParser(description='Plot groups')
parser.add_argument('--working-dir', dest='working_dir', type=str, help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser.add_argument('--metadata', dest='metadata', type=str, help='Metadata file output by sctype (annotate_clusters)')
parser.add_argument('--group', dest='group', type=str, help='Group to plot umaps for')
parser.add_argument('--output-group-umap-plot-prefix', dest='output_group_umap_plot_prefix', type=str, help='Output file prefix to write the group umap plot to; will write both a PDF and a PNG.')
args = parser.parse_args()

# Set working directory and load packages
scanpy.settings.verbosity = 1
scanpy.settings.figdir = 'plots/'
scanpy.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, format='png', figsize=('12', '8')) # type: ignore

# Parse the arguments
args = parser.parse_args()

adata = scanpy.read_h5ad(args.adata_input) # type: ignore
