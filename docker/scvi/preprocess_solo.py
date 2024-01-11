import scvi
import scanpy


adata = scanpy.read_10x_h5(snakemake.input.cellbender) # type: ignore


adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['rb'] = adata.var_names.str.startswith(('RPL', 'RPS'))

# scvi to normalize and prep for solo
noise = ['doublet_score', 'pct_counts_mt', 'pct_counts_rb']
scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key='sample', continuous_covariate_keys=noise)
model = scvi.model.SCVI(
    adata, 
    n_layers=1, 
    n_latent=50, 
    dispersion='gene',
    gene_likelihood='zinb'
)

model.train(
    train_size=0.7,
    max_epochs=1000,
    accelerator='gpu',  
    early_stopping=True,
    early_stopping_patience=40,
    plan_kwargs={'lr_factor': 0.1, 'lr_patience': 20, 'reduce_lr_on_plateau': True}
)


adata.obsm[snakemake.params.latent_key] = model.get_latent_representation() # type: ignore


##  use solo instead of scrublet
solo = scvi.external.SOLO.from_scvi_model(model) # adata is already registered
solo.train()
doublet_score = solo.predict(soft=True)
adata.obs['doublet_score'] = doublet_score

# artifact 1
# minify the adata 

# artifact 2
# save the models? (could save model with minified adata together)
scanpy.pp.calculate_qc_metrics(adata, qc_vars=['rb', 'mt'], percent_top=None, log1p=False, inplace=True)

adata.X = model.get_normalized_expression()


adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip') # type: ignore

