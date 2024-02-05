

#### _PREPROCESS_
- _pre-preprocessing_: `cellbender` ('cellbender.py') ~~or `scAR`('preprocess_scar.py').~~  (replacing `soupx` in preprocess.R)
- _preprocessing_: `scrublet` ('preprocess.py') ~~or`SOLO` ('preprocess_solo.py')~~ 

#### _QC_
- 'plot_qc_metrics.py'
    - merges adatas
    - NOTE:  merging here could be inefficient.  stubs of code to create gene_list to pre-filter 
        - consider filtering cells first... 'processing' filters genes to "highly variable"

#### _FILTER and PROCESS_ 
- _filtering_: 'filter.py'
- _processing_: 'process.py'.  normalize + identify variable genes

#### _INTEGRATION_
- 'integrate_scvi.py'.  
    - `scVI` integration 

#### _CLUSTERING_
- _clustering_: `umap` ('clustering_umap.py') or `mde` ('clustering_mde.py')
- NOTE:  mde is super fast and efficient on a GPU

#### _ANNOTATION_
- 'annotate_clusters.py'.  To implement .. Use cellassign and a list of marker genes.  Need to curate a list of marker genes apropriate for all the brain regions we are processing.

#### _PLOTTING_
- 'plot_feats_and_groups.py'.  


#### CONSIDERATIONS.
The "best" tools are probably `scAR`+`SOLO`+`scVI`.   Tradeoffs are probably some computation.  However, withouth further testing
we will default to `cellbender`+`scrublet` for removing ambient and creating the doublet metrics. `scVI` to be used for integration.

