import scanpy
import anndata

scanpy.settings.verbosity = 1
scanpy.settings.figdir = 'plots/'
scanpy.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, format='png', figsize=('12', '8')) # type: ignore
