# TODO:  implement a function that takes in a scvi model and adata and annotates the clusters
#     refer to utily/sctype.r
import argparse
import scanpy


parser = argparse.ArgumentParser(
    description='Annotate clusters'
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
    '--seurat-object',
    dest='seurat_object',
    type=str,
    help='Seurat object for a dataset'
)
parser.add_argument(
    '--cell-type-markers-list',
    dest='cell_type_markers_list',
    type=str,
    help='Seurat object containing a list of major cell type markers'
)
# parser.add_argument(
#     '--output-metadata-file',
#     dest='output_metadata_file',
#     type=str,
#     help='Output file to write metadata to'
# )
parser.add_argument(
    '--adata-output',
    dest='adata_output',
    type=str,
    help='Output file to save AnnData object to'
)

args = parser.parse_args()

adata = scanpy.read_h5ad(args.adata_input) # type: ignore

# TODO: write annotation code to clusters
# cell assign??
#  1. load marker_genes 
#  2. bdata = adata[:, marker_genes].copy()
#  3. adata.obs["size_factor"] = adata.obs["size_factor"].sum(1)/(adata.obs["size_factor"].sum(1).mean())
#  4. model = CellAssign(bdata, marker_genes)
#  5. model.train()
#  6. model.predict()


# TODO - write_h5ad option compression='gzip' is giving an error
#adata.write_h5ad(filename=args.adata_output, compression='gzip')
adata.write_h5ad(filename=args.adata_output)
