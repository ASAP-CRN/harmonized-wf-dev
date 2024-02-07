import argparse
import scanpy

parser = argparse.ArgumentParser(description='Annotate clusters')
parser.add_argument('--working-dir', dest='working_dir', type=str, 
                    help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser.add_argument('--script-dir', dest='script_dir', type=str, 
                    help='Directory containing workflow scripts', default='scripts')
parser.add_argument('--threads', dest='threads', type=int, 
                    help='Number of threads to use for processing')
parser.add_argument('--adata-input', dest='adata_input', type=str, 
                    help='AnnData object for a dataset')
parser.add_argument('--cell-type-markers-list', dest='cell_type_markers_list', type=str, 
                    help='Seurat object containing a list of major cell type markers')
parser.add_argument('--output-metadata-file', dest='output_metadata_file', type=str, 
                    help='Output file to write metadata to')
parser.add_argument('--adata-output', dest='adata_output', type=str, 
                    help='Output file to save AnnData object to')
parser.add_argument('--latent-key', dest='latent_key', type=str, default='X_scvi',
                    help='latent key to the scvi latent')
# Parse the arguments
args = parser.parse_args()

adata = scanpy.read_h5ad(args.adata_input) # type: ignore

# calculate neighbor graph on scVI latent
scanpy.pp.neighbors(adata, n_neighbors=50, use_rep=args.latent_key)

# do leiden
scanpy.tl.leiden(adata, resolution=0.4)

scanpy.tl.umap(adata)


adata.write_h5ad(filename=args.adata_output, compression='gzip') 
