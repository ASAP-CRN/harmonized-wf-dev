import scanpy
from scvi.model.utils import mde
import argparse
import scanpy
parser = argparse.ArgumentParser(description='Run clustering')
parser.add_argument('--working-dir', dest='working_dir', type=str, help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser.add_argument('--script-dir', dest='script_dir', type=str, help='Directory containing workflow scripts', default='scripts')
parser.add_argument('--threads', dest='threads', type=int, help='Number of threads to use for processing')
parser.add_argument('--seurat-object', dest='seurat_object', type=str, help='Seurat object for a dataset')
parser.add_argument('--clustering-algorithm', dest='clustering_algorithm', type=int, help='Clustering algorithm to use')
parser.add_argument('--clustering-resolution', dest='clustering_resolution', type=float, help='Clustering resolution')
parser.add_argument('--cell-type-markers-list', dest='cell_type_markers_list', type=str, help='Seurat object containing a list of major cell type markers')
parser.add_argument('--output-cell-type-plot-prefix', dest='output_cell_type_plot_prefix', type=str, help='Prefix for output file to save major cell type umap plot in PDF and PNG formats')
parser.add_argument('--output-seurat-object', dest='output_seurat_object', type=str, help='Output file to save Seurat object to')
args = parser.parse_args()

adata = scanpy.read_h5ad(snakemake.input.anndata)

# calculate PCS for reference
scanpy.pp.pca(adata, n_comps=50)

# calculate neighbor graph on scVI latent
scanpy.pp.neighbors(adata, n_neighbors=50, use_rep=snakemake.params.latent_key)

# do leiden
scanpy.tl.leiden(adata, resolution=0.4)

# mde is faster than umap, but needs a GPU to really shine
mde(adata.obsm[snakemake.params.latent_key]) 


adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip')
