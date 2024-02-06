import scanpy as sc
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Plot groups')
parser.add_argument('--working-dir', dest='working_dir', type=str, 
                    help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser.add_argument('--metadata', dest='metadata', type=str, 
                    help='Metadata file output by sctype (annotate_clusters)')
parser.add_argument('--group', dest='group', type=str, 
                    help='Group to plot umaps for')
parser.add_argument('--output-group-umap-plot-prefix', dest='output_group_umap_plot_prefix', type=str, 
                    help='Output file prefix to write the group umap plot to; will write both a PDF and a PNG.')
parser.add_argument('--feature', dest='feature', type=str, 
                    help='Feature to plot umaps for')
parser.add_argument('--output-feature-umap-plot-prefix', dest='output_feature_umap_plot_prefix', type=str, 
                    help='Output file to write the feature umap plot to')
args = parser.parse_args()

# Set working directory and load packages
sc.settings.verbosity = 1
sc.settings.figdir = 'plots/'
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, format='png', figsize=('12', '8')) # type: ignore

# Parse the arguments
args = parser.parse_args()

adata = sc.read_h5ad(args.adata_input) # type: ignore


plot_featues = [x for x in args.feature if x in adata.obs.columns]
file_name = args.output_feature_umap_plot_prefix + '_features_umap.png'
sc.pl.embedding(adata, basis="umap", color=plot_featues, frameon=False, show=False, ncols=1, save=file_name)

plot_groups = [x for x in args.group if x in adata.obs.columns]
file_name = args.output_group_umap_plot_prefix + '_groups_umap.png'
sc.pl.embedding(adata, basis="umap", color=plot_groups, frameon=False, show=False, ncols=1, save=file_name)


# TODO: impliment plotting code
