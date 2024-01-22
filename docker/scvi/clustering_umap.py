import scanpy

adata = scanpy.read_h5ad(snakemake.input.anndata)


# k = round(median(sqrt(adata.obs.sample.value_counts() * 0.5)))
# calculate PCS for reference
scanpy.pp.pca(adata, n_comps=50)

# calculate neighbor graph on scVI latent
scanpy.pp.neighbors(adata, n_neighbors=50, use_rep=snakemake.params.latent_key)

# do leiden
scanpy.tl.leiden(adata, resolution=0.4)

scanpy.tl.umap(adata)


adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip')
