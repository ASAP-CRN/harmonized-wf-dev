from pathlib import Path
import os
import anndata as ad
import scanpy as sc

DATA_PATH = Path("/media/ergonyc/DATA/scdata/ASAP")

HAFLER_DATA = DATA_PATH / "team-hafler/raw/1.0.0"

LEE_DATA = DATA_PATH / "team-lee/raw/1.0.0"

OUTPUT_PATH = DATA_PATH / "artifacts"

# test with HAFLER
RAW = "raw_feature_bc_matrix.h5"
FILT = "filtered_feature_bc_matrix.h5"

raw = sorted([f for f in HAFLER_DATA.glob(f"*{RAW}")])
filtered = sorted([f for f in HAFLER_DATA.glob(f"*{FILT}")])

samples = [ dict(raw=r, filt=f) for r, f in zip(raw, filtered ) if r.stem.rstrip(RAW) == f.stem.rstrip(FILT)]
samples

samples[0]['raw']


os.chdir(sample)

Path.cwd()
# args = parser.parse_args()

# os.chdir(os.path.join(os.getcwd(), 'data', args.sample))
import subprocess

# input_path = Path('/media/ergonyc/DATA/scdata/ASAP/team-hafler/raw/1.0.0/HSDG07HC.raw_feature_bc_matrix.h5')
# input = input_path.name
# sample = input_path.parent

output = 'HSDG078C_output.h5'
fpr = '0.0' 
cellbender_args = [
    '--input', input, 
    '--output', output,
    '--fpr', fpr,        
    '--cuda'
]


subprocess.run(['cellbender', 'remove-background'] + cellbender_args)
# import function
from cellbender.remove_background.downstream import anndata_from_h5

# load the data
out_adata = anndata_from_h5(output)
out_adata

# from cellbender.remove_background.downstream import load_anndata_from_input_and_output

# adata = load_anndata_from_input_and_output(
#     input_file='my_raw_10x_file.h5',
#     output_file='my_cellbender_output_file.h5',
#     input_layer_key='raw',  # this will be the raw data layer
# )
# adata