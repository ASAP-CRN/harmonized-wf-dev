

There are two variations here:

1) `cellbender`->`preprocess`-> `plot_metrics` (merge) -> `filter` -> `scvi` -> `clustering_mde` (or `clustering_umap`)

and 

2) `prepreprocess`(scar)-> `preprocess_solo` -> (`plot_metrics`) (merge) -> `filter` -> `scvi` ->`clustering_mde` (or `clustering_umap`)


For #2 its still not clear what the optimal flow is.  `scVI` VAEs are trained on individual samples AND are used to harmonize teh data downstream. 