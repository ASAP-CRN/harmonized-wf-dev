import argparse
import sys
import scvi
import scanpy


parser = argparse.ArgumentParser(
    description='Preprocess'
)
parser.add_argument(
	'--working-dir',
	dest='working_dir',
	type=str,
    help='Working directory',
	default='/team-X/xxx'
)
parser.add_argument(
	'--script-dir',
	dest='script_dir',
	type=str,
    help='Directory containing workflow scripts',
	default='scvi'
)
parser.add_argument(
	'--adata-input',
	dest='adata_input',
	type=str,
    help='AnnData object for a dataset'
)
# toplevel metadata to add to the adata
parser.add_argument(
	'--sample-id',
	dest='sample_id',
	type=str,
    help='Sample/dataset ID'
)
parser.add_argument(
	'--batch',
	dest='batch',
	type=str,
    help='Batch from which the sample/dataset originated'
)
parser.add_argument(
	'--project',
	dest='project',
	type=str,
    help='Project ID'
)
parser.add_argument(
	'--adata-output',
	dest='adata_output',
	type=str,
    help='Output file to save AnnData object to'
)

args = parser.parse_args()


os.chdir(args.working_dir)

sys.path.append('/opt/scripts/utility')
from helpers import anndata_from_h5, get_solo_results

# load the data from cellbender output
adata = anndata_from_h5(args.adata_input)
# adata = scanpy.read_10x_h5(args.adata_input)

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

# I'm not sure that this is returning reliable results... use scrublet for now...
# scores = solo.predict(adata, return_solo=False)
scores = get_solo_results(solo, adata, vae, gen_report=False, expected_doublet_rate=None) # 5% doublet rate matches scrublet default...

# add doublet score
adata.obs['doublet_score'] = scores['softmax_scores']

# add metadata
adata.obs['sample'] = args.sample_id
adata.obs['batch'] = args.batch
adata.obs['project'] = args.project
adata.obs['batch_id'] = args.project+args.batch #

# drop the celbender obs? 'background_fraction', 'cell_probability', 'cell_size', 'droplet_efficiency',

adata.write_h5ad(filename=args.adata_output, compression='gzip')
