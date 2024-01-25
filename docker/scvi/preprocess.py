import scanpy
import argparse
import os

# Create the parser
parser = argparse.ArgumentParser(description='Preprocess')

# Add arguments
parser.add_argument('--working-dir', dest='working_dir', type=str, 
                    help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser.add_argument('--script-dir', dest='script_dir', type=str, 
                    help='Directory containing workflow scripts', default='scripts')
parser.add_argument('--adata-input', dest='adata_input', type=str, 
                    help='AnnData object for a dataset')
# toplevel metadata to add to the adata
parser.add_argument('--sample-id', dest='sample_id', type=str, 
                    help='Sample/dataset ID')
parser.add_argument('--batch', dest='batch', type=str, 
                    help='Batch from which the sample/dataset originated')
parser.add_argument('--project', dest='project', type=str, 
                    help='Project ID')
                    
parser.add_argument('--output-adata', dest='output_adatat', type=str, 
                    help='Output file to save AnnData object to')


# Parse the arguments
args = parser.parse_args()

os.setwd(args.working_dir)
from util.helpers import anndata_from_h5, get_solo_results

# load the data from cellbender output
adata = anndata_from_h5(args.adata_input)
# adata = scanpy.read_10x_h5(args.adata_input)

adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['rb'] = adata.var_names.str.startswith(('RPL', 'RPS'))
scanpy.pp.calculate_qc_metrics(adata, qc_vars=['rb', 'mt'], percent_top=None, log1p=False, inplace=True)

# add doublet score
scanpy.external.pp.scrublet(adata)
# keep the metadata clean ad drop predicted_doublet... rely on the probability of dublet i.e. doublet_score
adata.obs.drop('predicted_doublet', axis=1, inplace=True)

# add metadata
adata.obs['sample'] = args.sample_id
adata.obs['batch'] = args.batch
adata.obs['project'] = args.project
adata.obs['batch_id'] = args.project+args.batch #


adata.write_h5ad(filename=args.output_adata, compression='gzip') 
