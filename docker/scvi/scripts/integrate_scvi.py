import scvi
import anndata as ad
import argparse

# Create the parser
parser = argparse.ArgumentParser(description='Run harmony analysis')

# Add arguments
parser.add_argument('--latent-key', dest='latent_key', type=str, default='X_scvi',
                    help='latent key to save the scvi latent to')
parser.add_argument('--adata-input', dest='adata_input', type=str, 
                    help='AnnData object for a dataset')
parser.add_argument('--adata-output', dest='adata_output', type=str, 
                    help='Output file to save AnnData object to')
parser.add_argument('--output-scvi', dest='output_scvi', type=str, 
                    help='Output file to save `scvi` model')
# TODO: optional scvi arguments


# Parse the arguments
args = parser.parse_args()

## parameters
n_latent = 10
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
# scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key='sample', continuous_covariate_keys=noise)
scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key='sample')

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
model.save(args.output_scvi)

adata.write_h5ad(filename=args.adata_output, compression='gzip') 
