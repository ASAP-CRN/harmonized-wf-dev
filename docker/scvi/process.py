import scanpy
import argparse

# Create the parser
parser = argparse.ArgumentParser(description='Normalize seurat objects')

# Add arguments
parser.add_argument('--working-dir', dest='working_dir', type=str, 
                    help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser.add_argument('--script-dir', dest='script_dir', type=str, 
                    help='Directory containing workflow scripts', default='scripts')
parser.add_argument('--threads', dest='threads', type=int, 
                    help='Number of threads to use for processing')
parser.add_argument('--adata-input', dest='adata_input', type=str, 
                    help='AnnData object for a dataset')
parser.add_argument('--output-adata', dest='output_adatat', type=str, 
                    help='Output file to save AnnData object to')

# Parse the arguments
args = parser.parse_args()

# TODO: impliment cell cycle scoring 
# https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt

# cell_cycle_genes = [x.strip() for x in open('./data/regev_lab_cell_cycle_genes.txt')]
# # Split into 2 lists
# s_genes = cell_cycle_genes[:43]
# g2m_genes = cell_cycle_genes[43:]

# cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]


adata = scanpy.read_h5ad(snakemake.input.anndata) # type: ignore

# does this work with sparse uint8?
adata.layers['counts'] = adata.X.copy() # type: ignore

scanpy.pp.normalize_total(adata, target_sum=1e4)

scanpy.pp.log1p(adata)

adata.raw = adata

# sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

scanpy.pp.highly_variable_genes(
    adata, 
    batch_key='sample', 
    subset=True, 
    flavor='seurat_v3', 
    layer='counts', 
    n_top_genes=3000
)


adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip') # type: ignore