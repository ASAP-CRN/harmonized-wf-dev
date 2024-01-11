import scanpy
from math import sqrt
from statistics import median

adata = scanpy.read_h5ad(snakemake.input.anndata)


# k = round(median(sqrt(adata.obs.sample.value_counts() * 0.5)))

scanpy.pp.neighbors(adata, n_neighbors=50, use_rep=snakemake.params.latent_key)

scanpy.tl.leiden(adata, resolution=0.4)

scanpy.tl.umap(adata)


adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip')
