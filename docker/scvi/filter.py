# import muon.pp.filter_obs as filter_obs
import muon
import scanpy
import argparse

# Create the parser
parser = argparse.ArgumentParser(description='Filter')

# Add arguments
parser.add_argument('--working-dir', dest='working_dir', type=str, 
                    help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser.add_argument('--script-dir', dest='script_dir', type=str, 
                    help='Directory containing workflow scripts', default='scripts')
parser.add_argument('--metadata', dest='metadata', type=str, 
                    help='Metadata file output by gmm_doublet_calling')
parser.add_argument('--adata-input', dest='adata_input', type=str, 
                    help='AnnData object for a dataset')
parser.add_argument('--output-adata', dest='output_adatat', type=str, 
                    help='Output file to save AnnData object to')

# backed adata hack https://discourse.scverse.org/t/concat-anndata-objects-on-disk/400/2
# Parse the arguments
args = parser.parse_args()

adata = scanpy.read_h5ad(args.adata_input) # type: ignore

# muon api is better than the scanpyt api for this...
muon.pp.filter_obs(adata, 'pct_counts_mt', lambda x: x <= 10)
muon.pp.filter_obs(adata, 'doublet_score', lambda x: x < 0.2)                    
muon.pp.filter_obs(adata, 'total_counts', lambda x: (x >= 500) & (x <= 100000))
muon.pp.filter_obs(adata, 'n_genes_by_counts', lambda x: (x >= 300) & (x <= 10000))

adata.write_h5ad(filename=args.output_adata, compression='gzip') 
