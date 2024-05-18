





python3 /opt/scripts/main/preprocess.py \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--adata-input ~{cellbender_counts} \
			--sample-id ~{sample_id} \
			--batch ~{batch} \
			--project ~{project_id} \
			--adata-output ~{sample_id}.adata_object.h5ad


python3 /opt/scripts/main/plot_qc_metrics.py \
			--working-dir "$(pwd)" \
			--script-dir "/opt/scripts" \
			--threads 8 \
			--adata-objects-fofn "adata_objects_paths.txt" \
			--project-name "test" \
			--adata-output "testing/02_test_merged_adata.h5ad"


python3 /opt/scripts/main/filter.py \
    --adata-input "testing/02_test_merged_adata.h5ad"  \
    --adata-output "testing/03_test_merged_adata_filtered.h5ad"



python3 /opt/scripts/main/process.py \
				--adata-input "testing/03_test_merged_adata_filtered.h5ad" \
				--adata-output "testing/04_test_merged_adata_filtered_normalized.h5ad" \
				--n-top-genes 3000


python3 /opt/scripts/main/integrate_scvi.py \
			--latent-key "X_scvi" \
			--batch-key "sample" \
			--adata-input "testing/04_test_merged_adata_filtered_normalized.h5ad" \
			--adata-output "testing/05_test_merged_adata_filtered_normalized_integrated.h5ad" \
			--output-scvi "test_scvi_model.pkl"


python3 /opt/scripts/main/clustering_umap.py \
				--working-dir "$(pwd)" \
				--script-dir "/opt/scripts" \
				--threads 8 \
				--adata-input "testing/05_test_merged_adata_filtered_normalized_integrated.h5ad" \
				--adata-output "testing/06_test_merged_adata_filtered_normalized_integrated_clustered.h5ad" \
				--latent-key "X_scvi" 


python3 /opt/scripts/main/annotate_cells.py \
			--adata-input "testing/06_test_merged_adata_filtered_normalized_integrated_clustered.h5ad" \
			--marker-genes "testing/celltype_marker_table.csv" \
			--output-cellassign "testing/cellassign_model.pkl" \
			--output-cell-types-file "testing/cell_types.csv" \
			--adata-output "testing/07_test_merged_adata_filtered_normalized_integrated_clustered_annotated.h5ad"


python3 /opt/scripts/main/plot_feats_and_groups.py \
			--working-dir "$(pwd)" \
			--adata-input "testing/07_test_merged_adata_filtered_normalized_integrated_clustered_annotated.h5ad" \
			--group ~{sep=',' groups} \
			--output-group-umap-plot-prefix "~{cohort_id}" \
			--feature ~{sep=',' features} \
			--output-feature-umap-plot-prefix "~{cohort_id}"
