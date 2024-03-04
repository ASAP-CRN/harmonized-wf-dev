import argparse
import scanpy as sc


parser = argparse.ArgumentParser(
    description='Plot groups'
)
parser.add_argument(
    '--adata-input',
    dest='adata_input',
    type=str,
    help='AnnData object for a dataset'
)
parser.add_argument(
    '--group',
    dest='group',
    type=str,
    help='Group to plot umaps for'
)
parser.add_argument(
    '--output-group-umap-plot-prefix',
    dest='output_group_umap_plot_prefix',
    type=str,
    help='Output file prefix to write the group umap plot to; will write both a PDF and a PNG.'
)
parser.add_argument(
    '--feature',
    dest='feature',
    type=str,
    help='Feature to plot umaps for'
)
parser.add_argument(
    '--output-feature-umap-plot-prefix',
    dest='output_feature_umap_plot_prefix',
    type=str,
    help='Output file to write the feature umap plot to'
)

args = parser.parse_args()


# Set CPUs to use for parallel computing
sc._settings.ScanpyConfig.n_jobs = -1
scvi.settings.dl_num_workers = -1

# Set working directory and load packages
sc.settings.verbosity = 1
sc.settings.figdir = 'plots/'
sc.settings.set_figure_params(
    dpi=100,
    fontsize=10,
    dpi_save=300,
    format='png',
    figsize=('12', '8')
) # type: ignore

adata = sc.read_h5ad(args.adata_input) # type: ignore

features_list = args.feature.split(',')
plot_features = [x for x in features_list if x in adata.obs.columns]
file_name = args.output_feature_umap_plot_prefix + '_features_umap.png'
sc.pl.embedding(
    adata,
    basis="umap",
    color=plot_features,
    frameon=False,
    show=False,
    ncols=1,
    save=file_name
)

groups_list = args.group.split(',')
plot_groups = [x for x in groups_list if x in adata.obs.columns]
file_name = args.output_group_umap_plot_prefix + '_groups_umap.png'
sc.pl.embedding(
    adata,
    basis="umap",
    color=plot_groups,
    frameon=False,
    show=False,
    ncols=1,
    save=file_name
)


# TODO: impliment plotting code
