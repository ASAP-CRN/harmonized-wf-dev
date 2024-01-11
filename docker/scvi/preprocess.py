import scanpy


adata = scanpy.read_10x_h5(snakemake.input.cellbender) # type: ignore


adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['rb'] = adata.var_names.str.startswith(('RPL', 'RPS'))

scanpy.pp.calculate_qc_metrics(adata, qc_vars=['rb', 'mt'], percent_top=None, log1p=False, inplace=True)

scanpy.external.pp.scrublet(adata, expected_doublet_rate=(adata.n_obs / 1000) * 0.008)
adata.obs.drop('predicted_doublet', axis=1, inplace=True)

adata.obs['sample'] = snakemake.params.sample # type: ignore


adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip') # type: ignore