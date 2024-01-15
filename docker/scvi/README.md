

There are two variations sketched here:


1) A more generic variant: `cellbender`->`preprocess`-> `plot_metrics` (merge) -> `filter` -> `scvi` -> `mde` (visualizing only, or `clustering_mde` to cluster on the neighbor graph (slow))

and 
2) a full scvi-tools variant: `prepreprocess`(scar)-> `preprocess_solo` -> (`plot_metrics`) (merge) -> `filter` -> `scvi` ->`clustering_mde` 


For #2 its still not clear what the optimal flow is.  `scVI` VAEs are trained on individual samples AND are used to harmonize teh data downstream. 

Also, its unclear how well normalizing with scvi (in pre-process) and then using those normalized "expression" values works compared with batch correcting all at once.
