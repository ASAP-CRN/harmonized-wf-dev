import scanpy
from scvi.model.utils import mde

adata = scanpy.read_h5ad(snakemake.input.anndata)

# calculate PCS for reference
scanpy.pp.pca(adata, n_comps=50)

# calculate neighbor graph on scVI latent
scanpy.pp.neighbors(adata, n_neighbors=50, use_rep=snakemake.params.latent_key)

# do leiden
scanpy.tl.leiden(adata, resolution=0.4)

# mde is faster than umap, but needs a GPU to really shine
mde(adata.obsm[snakemake.params.latent_key]) 


adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip')
