import scanpy
import anndata
import argparse

# Create the parser
parser = argparse.ArgumentParser(description='Plot features')

# Add arguments
parser.add_argument('--working-dir', dest='working_dir', type=str, 
                    help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser.add_argument('--metadata', dest='metadata', type=str, 
                    help='Metadata file output by sctype (annotate_clusters)')
parser.add_argument('--feature', dest='feature', type=str, 
                    help='Feature to plot umaps for')
parser.add_argument('--output-feature-umap-plot-prefix', dest='output_feature_umap_plot_prefix', type=str, 
                    help='Output file to write the feature umap plot to')

# Parse the arguments
args = parser.parse_args()


scanpy.settings.verbosity = 1
scanpy.settings.figdir = 'plots/'
scanpy.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, format='png', figsize=('12', '8')) # type: ignore
