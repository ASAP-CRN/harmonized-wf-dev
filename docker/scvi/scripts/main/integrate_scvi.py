import os
import argparse
import scvi
import anndata as ad
import scanpy


parser = argparse.ArgumentParser(
    description='Run scVI integration'
)
parser.add_argument(
    '--latent-key',
    dest='latent_key',
    type=str,
    default='X_scvi',
    help='Latent key to save the scvi latent to'
)
parser.add_argument(
    '--batch-key',
    dest='batch_key',
    type=str,
    help='Key in AnnData object for batch information'
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
parser.add_argument(
    '--output-scvi-dir',
    dest='output_scvi_dir',
    type=str,
    help='Output folder to save `scvi` model'
)

# TODO: optional scvi arguments

args = parser.parse_args()


# Set CPUs to use for parallel computing
scanpy._settings.ScanpyConfig.n_jobs = -1

## parameters
n_latent = 30
n_layers = 2
train_size = 0.85
scvi_epochs = 200
batch_size = 512
# accelerator = 'gpu'
# devide = "cuda:0"
dispersion = 'gene'   
plan_kwargs = {'lr_factor': 0.1, 'lr_patience': 20, 'reduce_lr_on_plateau': True}
gene_likelihood = 'zinb'
latent_distribution = 'normal'
early_stopping = True
early_stopping_patience = 40
# latent_key='X_scvi'

adata = scanpy.read_h5ad(args.adata_input) # type: ignore

# just use the counts layer. does this delete the "raw"
# adata = adata.layers['counts']

## integrate the data with `scVI`
# noise = ['doublet_score', 'pct_counts_mt', 'pct_counts_rb']
# scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key=args.batch_key, continuous_covariate_keys=noise)
scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key=args.batch_key)

model = scvi.model.SCVI(
    adata, 
    n_layers=n_layers, 
    n_latent=n_latent, 
    dispersion=dispersion,
    gene_likelihood=gene_likelihood,
)

model.train(
    train_size=train_size,
    max_epochs=scvi_epochs,
    early_stopping=early_stopping,
    early_stopping_patience=early_stopping_patience,
    plan_kwargs=plan_kwargs,
)

adata.obsm[args.latent_key] = model.get_latent_representation() # type: ignore


# make minified adata
# TODO: impliment

# artifacts
model.save(args.output_scvi_dir, overwrite=True)

adata.write_h5ad(filename=args.adata_output, compression='gzip')
