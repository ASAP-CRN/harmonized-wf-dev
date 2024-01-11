import muon
import scanpy


adata = scanpy.read_h5ad(snakemake.input.anndata) # type: ignore


muon.pp.filter_obs(adata, 'pct_counts_mt', lambda x: x <= 10)
muon.pp.filter_obs(adata, 'doublet_score', lambda x: x < 0.2)                    
muon.pp.filter_obs(adata, 'total_counts', lambda x: (x >= 500) & (x <= 100000))
muon.pp.filter_obs(adata, 'n_genes_by_counts', lambda x: (x >= 300) & (x <= 10000))

adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip') # type: ignore