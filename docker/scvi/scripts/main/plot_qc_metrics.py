import argparse
import scanpy 
import anndata


parser = argparse.ArgumentParser(
    description='Call doublets'
)
# parser.add_argument(
#     '--working-dir',
#     dest='working_dir',
#     type=str,
#     help='Working directory',
#     default='/data/CARD_singlecell/harmony-rna/'
# )
# parser.add_argument(
#     '--script-dir',
#     dest='script_dir',
#     type=str,
#     help='Directory containing workflow scripts',
#     default='scripts'
# )
# parser.add_argument(
#     '--threads',
#     dest='threads',
#     type=int,
#     help='Number of threads to use for processing'
# )
parser.add_argument(
    '--adata-objects-fofn',
	dest='adata_objects_fofn',
	type=str,
    help='Newline-delimited paths to the set of input adata objects (file-of-filenames)'
)
parser.add_argument(
    '--project-name',
	dest='project_name',
	type=str,
    help='Project name'
)

parser.add_argument(
    '--adata-output',
	dest='adata_output',
	type=str,
    help='Output file to save AnnData object to'
)

args = parser.parse_args()


scanpy.settings.verbosity = 1
scanpy.settings.figdir = 'plots/'
scanpy.settings.set_figure_params(
    dpi=100,
    fontsize=10,
    dpi_save=300,
    format='png',
    figsize=('12', '8')
) # type: ignore

metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rb', 'doublet_score']
scanpy.settings.verbosity = 1
scanpy.settings.figdir = 'plots/'
scanpy.settings.set_figure_params(
    dpi=100,
    fontsize=10,
    dpi_save=300,
    format='png',
    figsize=('12', '8')
) # type: ignore

adatas = {}
top_genes = {}


# TODO:  change how we get sample names / file names to aggregate?
#  note that the sample id should be the official ASAP_samples 
with open(args.adata_objects_fofn, 'r') as file:
    file_contents = file.read()

# Remove empty elements in list, this will cause an error when reading in
samples = list(filter(None, file_contents.split('\n')))

for sample in samples:
    raw = scanpy.read_h5ad(sample)
    # code below if memory issues with concatenating all the ge
    # adata = raw.copy()
    # scanpy.pp.normalize_total(adata, target_sum=1e4)
    # scanpy.pp.log1p(adata)
    # scanpy.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=8000)

    # # ranked_genes = adata.var[adata.var.highly_variable].dispersions_norm.argsort().to_dict()   
    # for k,v in ranked_genes.items():
    #     if k in top_genes:
    #         top_genes[k] += v
    #     else:
    #         top_genes[k] = v
    
    # raw = minify_adata(raw)
    sample_name = sample.split(".")[0]
    adatas[sample_name] = raw

# we could subset to the top_genes here before concat if we have memory issues (e.g. whole dataset harmonization.)
adata = anndata.concat(
            merge='same',
            uns_merge='same',
            index_unique='_',
            adatas=adatas
        )


for metric in metrics: # type: ignore
    scanpy.pl.violin(adata, keys=metric, size=0, save=''.join('_' + metric))

# top_genes = pd.DataFrame(index=top_genes.keys(), columns=['rank'], data=top_genes.values())
#TODO: export top_genes? and plots

# export concatenated data.
adata.write_h5ad(filename=args.adata_output)
