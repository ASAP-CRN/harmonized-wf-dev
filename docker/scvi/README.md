

## _PREPROCESS_
- _pre-preprocessing_: `cellbender` ('cellbender.py') or `scAR`('prepreprocess.py').  (replacing `soupx` in preprocess.R)
- _preprocessing_: `scrublet` ('preprocess.py') or `SOLO` ('preprocess_solo.py'). 

## _QC_
- 'plot_qc_metrics.py'
- NOTE:  currently concatenating all of the adatas.  This is inefficient

## _FILTER and PROCESS_
- _filtering_: 'filter.py'
- _processing_: 'process.py'.  normalize + identify variable genes

## _INTEGRATION_
- 'scvi.py'.  
    - merge adatas
    - `scVI` integration 

## _CLUSTERING_
- _clustering_: `umap` ('clustering_umap.py') or `mde` ('clustering_mde.py')
- NOTE:  mde is super fast and efficient on a GPU

## _ANNOTATION_
- 'annotate_clusters.py'.  To implement ..

## _PLOTTING_
- 'plot_featueres.py'.  To implement ..
- 'plot_groups.py'.  To implement ..




## CONSIDERATIONS.
The "best" tools are probably `scAR`+`SOLO`+`scVI`.  Tradeoffs are probably some computation.

- filtering and processing on concatenated datasets... highly variable gene selection might be optimized
Tradeoffs between generating extra artifacts 

