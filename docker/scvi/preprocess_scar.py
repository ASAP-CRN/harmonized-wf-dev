
from scvi.external import SCAR
from scanpy import read_10x_h5 

# use scar instead of cellbender
raw_adata = read_10x_h5(snakemake.input.raw_anndata)
raw_adata.var_names_make_unique()

filt_adata = read_10x_h5(snakemake.input.filtered_anndata)
filt_adata.var_names_make_unique()

SCAR.setup_anndata(filt_adata)
SCAR.get_ambient_profile(adata=filt_adata, raw_adata=raw_adata, prob=0.995)
scar_vae = SCAR(filt_adata)
scar_vae.train()
filt_adata.obsm["X_scAR"] = scar_vae.get_latent_representation()
filt_adata.layers["denoised"] = scar_vae.get_denoised_counts()

filt_adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip') # type: ignore
