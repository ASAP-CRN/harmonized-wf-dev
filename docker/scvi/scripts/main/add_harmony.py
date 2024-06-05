import argparse
import scanpy

parser = argparse.ArgumentParser(description="Add PCA and Harmony integration")
parser.add_argument(
    "--batch-key",
    dest="batch_key",
    type=str,
    help="Key in AnnData object for batch information",
)
parser.add_argument(
    "--adata-input", dest="adata_input", type=str, help="AnnData object for a dataset"
)

parser.add_argument(
    "--adata-output",
    dest="adata_output",
    type=str,
    help="Output file to save AnnData object to",
)

args = parser.parse_args()


adata = scanpy.read_h5ad(args.adata_input)  # type: ignore

# make sure we have PCA
if "X_pca" not in adata.obsm:
    scanpy.pp.pca(adata, n_comps=30)

# add harmony
if "X_pca_harmony" not in adata.obsm:
    scanpy.external.pp.harmony_integrate(adata, args.batch_key)
    # 9. write_h5ad

adata.write_h5ad(filename=args.adata_output, compression="gzip")
