# TODO:  implement a function that takes in a scvi model and adata and annotates the clusters
#     refer to utily/sctype.r
import argparse


parser = argparse.ArgumentParser(description='Annotate clusters')
parser.add_argument('--working-dir', dest='working_dir', type=str, help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser.add_argument('--script-dir', dest='script_dir', type=str, help='Directory containing workflow scripts', default='scripts')
parser.add_argument('--threads', dest='threads', type=int, help='Number of threads to use for processing')
parser.add_argument('--seurat-object', dest='seurat_object', type=str, help='Seurat object for a dataset')
parser.add_argument('--cell-type-markers-list', dest='cell_type_markers_list', type=str, help='Seurat object containing a list of major cell type markers')
parser.add_argument('--output-metadata-file', dest='output_metadata_file', type=str, help='Output file to write metadata to')
args = parser.parse_args()


# TODO: write annotation code to clusters