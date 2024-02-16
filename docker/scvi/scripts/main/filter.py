# import muon.pp.filter_obs as filter_obs
import muon
import scanpy
import argparse


parser = argparse.ArgumentParser(
    description='Filter'
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
# TODO: add filter parameters as arguments

args = parser.parse_args()

adata = scanpy.read_h5ad(args.adata_input) # type: ignore

# TODO: make these cutoffs arguments...
# muon api is better than the scanpyt api for this...
muon.pp.filter_obs(adata, 'pct_counts_mt', lambda x: x <= 10)
muon.pp.filter_obs(adata, 'doublet_score', lambda x: x < 0.2)                    
muon.pp.filter_obs(adata, 'total_counts', lambda x: (x >= 500) & (x <= 100000))
muon.pp.filter_obs(adata, 'n_genes_by_counts', lambda x: (x >= 300) & (x <= 10000))

# save the filtered adata
adata.write_h5ad(filename=args.adata_output)
