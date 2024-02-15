import argparse
import scanpy
import sys
sys.path.append('/opt/scripts/utility')
from helpers import score_cell_cycle


parser = argparse.ArgumentParser(
    description='Normalize seurat objects'
)

parser.add_argument(
	'--adata-input',
	dest='adata_input',
	type=str,
    help='AnnData object for a dataset'
)
parser.add_argument(
	'--adata-output',
	dest='adata_output',
	type=str,
    help='Output file to save AnnData object to'
)
parser.add_argument(
    '--n-top-genes',
    dest='n_top_genes',
    type=str,
    help='number of HVG genes to keep',
    default=8000
)

args = parser.parse_args()

adata = scanpy.read_h5ad(args.adata_input) # type: ignore

# TODO: load the top_genes from the qc plotting step and subset... to the top 8k genes
#  if we have memory issues, consider subsetting to a gene_list
# gene_list = pd.read_csv(top_genes, header=None)[0].tolist()
# adata = adata[:, gene_list.index]

# does this work with sparse uint8?
adata.layers['counts'] = adata.X.copy() # type: ignore

scanpy.pp.normalize_total(adata, target_sum=1e4)
scanpy.pp.log1p(adata)

# get cell cycle scores
score_cell_cycle(adata, organism="human")

# sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
scanpy.pp.highly_variable_genes(
    adata, 
    batch_key='sample', 
    subset=True, 
    flavor='seurat_v3', 
    layer='counts', 
    n_top_genes=args.n_top_genes
)

#  write_h5ad option compression='gzip' is giving an error
adata.write_h5ad(filename=args.adata_output)
