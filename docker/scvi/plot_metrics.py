import scanpy
import anndata

scanpy.settings.verbosity = 1
scanpy.settings.figdir = 'plots/'
scanpy.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, format='png', figsize=('12', '8')) # type: ignore


adata = anndata.concat(
    merge='same', uns_merge='same', index_unique='_',
    adatas={[item for item in dataset.split('_') if 'ARC' in item][0]: 
            scanpy.read_h5ad(dataset) for dataset in snakemake.input.objects} # type: ignore
)


for metric in snakemake.params.metrics: # type: ignore
    scanpy.pl.violin(adata, keys=metric, size=0, save=''.join('_' + metric))


adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip') # type: ignore

