import scanpy

adata = scanpy.read_10x_h5(snakemake.input.cellbender) # type: ignore

adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['rb'] = adata.var_names.str.startswith(('RPL', 'RPS'))

scanpy.pp.calculate_qc_metrics(adata, qc_vars=['rb', 'mt'], percent_top=None, log1p=False, inplace=True)

adata.obs['sample'] = snakemake.params.sample # type: ignore


adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip') # type: ignore

##  use solo instead of scrublet
solo = scvi.external.SOLO.from_scvi_model(model) # adata is already registered
solo.train()
doublet_score = solo.predict(soft=True)
adata.obs['doublet_score'] = doublet_score

# artifact 1
# minify the adata 

# artifact 2
# save the models? (could save model with minified adata together)

adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip') # type: ignore
