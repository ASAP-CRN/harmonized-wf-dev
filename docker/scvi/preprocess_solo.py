import scvi
import scanpy
from cellbender.remove_background.downstream import anndata_from_h5
import argparse

# Create the parser
parser = argparse.ArgumentParser(description='Preprocess')

# Add arguments
parser.add_argument('--working-dir', dest='working_dir', type=str, 
                    help='Working directory', default='/team-X/xxx')
parser.add_argument('--script-dir', dest='script_dir', type=str, 
                    help='Directory containing workflow scripts', default='scvi')
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

# load the data from cellbender output
adata = anndata_from_h5(args.inpadata_inputut)
# adata = scanpy.read_10x_h5(args.input)

adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['rb'] = adata.var_names.str.startswith(('RPL', 'RPS'))
scanpy.pp.calculate_qc_metrics(adata, qc_vars=['rb', 'mt'], percent_top=None, log1p=False, inplace=True)

# fit the scvi and solo models to create doublet scores
scvi.model.SCVI.setup_anndata(adata)
vae = scvi.model.SCVI(adata)
vae.train()

solo = scvi.external.SOLO.from_scvi_model(vae)
solo.train()

scores = solo.predict(soft=True)

# add doublet score
adata.obs['doublet_score'] = scores['doublet']
# adata.obs['singlet_score'] = scores['singlet']

# add metadata
adata.obs['sample'] = args.sample_id
adata.obs['batch'] = args.batch
adata.obs['project'] = args.project
adata.obs['batch_id'] = args.project+args.batch #


adata.write_h5ad(filename=args.output_adata, compression='gzip') 

