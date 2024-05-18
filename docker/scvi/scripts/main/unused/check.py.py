


# gs://asap-dev-data-cohort-feb12test 

# gs://asap-raw-data-team-lee/workflow_execution/preprocess/cellranger/
# gs://asap-raw-data-team-lee/workflow_execution/preprocess/counts_to_adata/
# gs://asap-raw-data-team-lee/workflow_execution/preprocess/remove_technical_artifacts/

# gs://asap-raw-data-team-hafler/workflow_execution/preprocess/cellranger/
# gs://asap-raw-data-team-hafler/workflow_execution/preprocess/counts_to_adata/
# gs://asap-raw-data-team-hafler/workflow_execution/preprocess/remove_technical_artifacts/

# gs://asap-raw-data-team-jakobsson/workflow_execution/preprocess/cellranger/
# gs://asap-raw-data-team-jakobsson/workflow_execution/preprocess/counts_to_adata/1.0.0
# gs://asap-raw-data-team-jakobsson/workflow_execution/preprocess/remove_technical_artifacts/1.0.0

# %%

import pandas as pd
from pathlib import Path
import anndata as ad
import scanpy as sc
import scvi

# %%

cohort_path = Path.home() / "Projects/ASAP/harmonized-wf-dev/check"



# artifacts
cell_types = pd.read_csv(cohort_path / "cohort_analysis/asap-cohort.cell_types.csv")
meta_data = pd.read_csv(cohort_path / "cohort_analysis/asap-cohort.annotate_cells.metadata.csv")

# %%
cell_types.head()
# %%
meta_data.head()




# %%
adata = ad.read_h5ad(cohort_path / "asap-cohort.merged_adata_object.h5ad")
adata

# %%
adata = ad.read_h5ad(cohort_path / "asap-cohort.merged_adata_object_filtered.h5ad")
adata

# %%
adata = ad.read_h5ad(cohort_path / "asap-cohort.merged_adata_object_filtered_normalized.h5ad")
adata

# %%
adata = ad.read_h5ad(cohort_path / "asap-cohort.adata_object.scvi_integrated.h5ad" )
adata

# %%
adata = ad.read_h5ad(cohort_path / "asap-cohort.adata_object.scvi_integrated.umap_cluster.h5ad" )
adata

# %%
adata = ad.read_h5ad(cohort_path / "cohort_analysis/asap-cohort.adata_object.scvi_integrated.umap_cluster.annotate_cells.h5ad")
adata


# %%
vae = scvi.model.SCVI.load(cohort_path / "asap-cohort.scvi_model", adata.copy())

# %%

ca_model = scvi.external.CellAssign.load(cohort_path / "asap-cohort.cellassign_model", adata.copy())

# %%
