
from scvi.external import SCAR
import scvi
import anndata as ad


# use scar instead of cellbender

adata = ad.read_h5ad(snakemake.input.filtered_anndata) 
raw_adata = ad.read_h5ad(snakemake.input.raw_anndata)

SCAR.setup_anndata(adata, batch_key="batch")
SCAR.get_ambient_profile(adata=adata, raw_adata=raw_adata, prob=0.995)
vae = SCAR(adata)
vae.train()
adata.obsm["X_scAR"] = vae.get_latent_representation()
adata.layers["denoised"] = vae.get_denoised_counts()

adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip') # type: ignore
