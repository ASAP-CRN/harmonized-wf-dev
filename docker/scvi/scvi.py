import scvi
import anndata as ad




# https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt

# cell_cycle_genes = [x.strip() for x in open('./data/regev_lab_cell_cycle_genes.txt')]
# # Split into 2 lists
# s_genes = cell_cycle_genes[:43]
# g2m_genes = cell_cycle_genes[43:]

# cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]


# parameters
n_latent = 10
n_layers = 2
train_size = 0.85
latent_key='X_scvi'

# extract sample name from the file name.

adata = ad.concat(
    merge='same', uns_merge='same', index_unique='_',
    adatas={[item for item in dataset.split('_') if 'ARC' in item][0]:   
            ad.read_h5ad(dataset) for dataset in snakemake.input.objects} # type: ignore
)

noise = ['doublet_score', 'pct_counts_mt', 'pct_counts_rb']

scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key='sample', continuous_covariate_keys=noise)

model = scvi.model.SCVI(
    adata, 
    n_layers=n_layers, 
    n_latent=n_latent, 
    dispersion='gene',
    gene_likelihood='zinb'
)


model.train(
    train_size=0.85,
    max_epochs=1000,
    accelerator='gpu',  
    early_stopping=True,
    early_stopping_patience=40,
    plan_kwargs={'lr_factor': 0.1, 'lr_patience': 20, 'reduce_lr_on_plateau': True}
)


adata.obsm[snakemake.params.latent_key] = model.get_latent_representation() # type: ignore

# artifact 1
# minify the adata 

# artifact 2
# save the model? (could save model with minified adata together)

adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip') # type: ignore