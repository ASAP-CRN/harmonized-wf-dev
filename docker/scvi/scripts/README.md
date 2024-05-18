# scVI analysis Pipeline python scripts


### _PREPROCESS_
- _pre-preprocessing_: [`cellbender`](../../cellbender/scripts/cellbender.py)
- _preprocessing_: [`scrublet`]('preprocess.py') 

### _QC: metrics_
- 'plot_qc_metrics.py'
    - merges adatas
    - NOTE:  merging here could be inefficient.  stubs of code to create gene_list to pre-filter 
        - consider filtering cells first... 'processing' filters genes to "highly variable"

### _QC: FILTER_ 
- _filtering_: 'filter.py'

### _FEATURE SELECTION_ 
- _processing_: 'process.py'.  normalize + feature selection (i.e. identification of highly variable genes.)

### _INTEGRATION_
- 'integrate_scvi.py'.  
    - `scVI` integration 

### _CLUSTERING_
- _clustering_: `umap` ('clustering_umap.py') 
    - in the future we may choose `mde` ('clustering_mde.py') over `umap`, as it is super fast and efficient on a GPU, and the embeddings are only useful for visualization so the choice is semi-arbitrary.


### _ANNOTATION_
- NOTE:  considering that this might be an "additional curation" item for cross-team harmonized artifacts.  This annotation should probably be done separately on different "atlas" level subsets of the data.
- 'annotate_cells.py'.  Use cellassign and a list of marker genes.  Need to curate a list of marker genes apropriate for all the brain regions we are processing.  Currently using CARD cortical list of genes.  NOTE: this is not annotating the "clusters" but the cells based on marker gene expression.

### _PLOTTING_
- 'plot_feats_and_groups.py'.  
- features: 'n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rb', 'doublet_score','S_score','G2M_score'
- groups: "sample", "batch", "cell_type"

### CONSIDERATIONS.
The "best" tools are probably `scAR`+`SOLO`+`scVI`.   Tradeoffs are probably some computation.  However, withouth further testing
we will default to `cellbender`+`scrublet` for removing ambient and creating the doublet metrics. `scVI` to be used for integration.

