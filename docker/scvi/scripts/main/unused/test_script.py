
# python3 /opt/scripts/main/preprocess.py \
# 			--working-dir "$(pwd)" \
# 			--script-dir /opt/scripts \
# 			--adata-input ~{cellbender_counts} \
# 			--sample-id ~{sample_id} \
# 			--batch ~{batch} \
# 			--project ~{project_id} \
# 			--adata-output ~{sample_id}.adata_object.h5ad


# python3 /opt/scripts/main/plot_qc_metrics.py \
# 			--adata-objects-fofn "adata_objects_paths.txt" \
# 			--adata-output "artifacts/asap-cohort.merged_adata_object.h5ad"
			# --output-validation-file "testing/output_validation.csv" 


python3 /opt/scripts/main/filter.py \
    --adata-input "artifacts/asap-cohort.merged_adata_object.h5ad" \
    --adata-output "artifacts/asap-cohort.merged_adata_object_filtered.h5ad" \
	--output-validation-file "artifacts/asap-cohort.output_validation.csv" 


python3 /opt/scripts/main/process.py \
				--adata-input "artifacts/asap-cohort.merged_adata_object_filtered.h5ad" \
				--adata-output "artifacts/asap-cohort.merged_adata_object_filtered_normalized.h5ad" \
                --marker-genes "testing/celltype_marker_table.csv" \
				--n-top-genes 3000 \
                --batch-key "sample" \
                --output-validation-file "artifacts/asap-cohort.output_validation.csv" 


python3 /opt/scripts/main/integrate_scvi.py \
			--latent-key "X_scvi" \
			--batch-key "sample" \
			--adata-input "artifacts/asap-cohort.merged_adata_object_filtered_normalized.h5ad"\
			--adata-output "artifacts/asap-cohort.merged_adata_object_filtered_normalized_integrated.h5ad" \
			--output-scvi "artifacts/test_scvi_model.pkl"



python3 /opt/scripts/main/clustering_umap.py \
				--adata-input "artifacts/asap-cohort.merged_adata_object_filtered_normalized_integrated.h5ad" \
				--adata-output "artifacts/asap-cohort.merged_adata_filtered_normalized_integrated_clustered.h5ad" \
				--latent-key "X_scvi" 


python3 /opt/scripts/main/annotate_cells.py \
			--adata-input "artifacts/asap-cohort.merged_adata_filtered_normalized_integrated_clustered.h5ad"  \
			--marker-genes "testing/celltype_marker_table.csv" \
			--output-cell-types-file "artifacts/cell_types.csv" \
			--adata-output "artifacts/asap-cohort.merged_adata_filtered_normalized_integrated_clustered_annotated.h5ad"

python3 /opt/scripts/main/add_harmony.py \
			--batch-key "sample" \
			--adata-input "artifacts/asap-cohort.merged_adata_filtered_normalized_integrated_clustered_annotated.h5ad"  \
			--adata-output "artifacts/asap-cohort.merged_adata_filtered_normalized_integrated_clustered_annotated_harmony.h5ad"

python3 /opt/scripts/main/artifact_metrics.py \
            --latent-key "X_scvi" \
			--batch-key "sample" \
			--adata-input "artifacts/asap-cohort.merged_adata_filtered_normalized_integrated_clustered_annotated_harmony.h5ad"  \
			--output-report-dir "artifacts/scib_report/" \





python3 /opt/scripts/main/integrate_scvi.py \
			--latent-key "X_scvi" \
			--batch-key "batch_id" \
			--adata-input "artifacts/asap-cohort.merged_adata_object_filtered_normalized.h5ad"\
			--adata-output "artifacts/asap-cohort2.merged_adata_object_filtered_normalized_integrated.h5ad" \
			--output-scvi "artifacts/test_scvi_model2.pkl"


python3 /opt/scripts/main/clustering_umap.py \
				--adata-input "artifacts/asap-cohort2.merged_adata_object_filtered_normalized_integrated.h5ad" \
				--adata-output "artifacts/asap-cohort2.merged_adata_filtered_normalized_integrated_clustered.h5ad" \
				--latent-key "X_scvi" 


python3 /opt/scripts/main/annotate_cells.py \
			--adata-input "artifacts/asap-cohort2.merged_adata_filtered_normalized_integrated_clustered.h5ad"  \
			--marker-genes "testing/celltype_marker_table.csv" \
			--output-cell-types-file "artifacts/cell_types2.csv" \
			--adata-output "artifacts/asap-cohort2.merged_adata_filtered_normalized_integrated_clustered_annotated.h5ad"

python3 /opt/scripts/main/add_harmony.py \
			--batch-key "batch_id" \
			--adata-input "artifacts/asap-cohort2.merged_adata_filtered_normalized_integrated_clustered_annotated.h5ad"  \
			--adata-output "artifacts/asap-cohort2.merged_adata_filtered_normalized_integrated_clustered_annotated_harmony.h5ad"



python3 /opt/scripts/main/artifact_metrics.py \
            --latent-key "X_scvi" \
			--batch-key "batch_id" \
			--adata-input "artifacts/asap-cohort2.merged_adata_filtered_normalized_integrated_clustered_annotated_harmony.h5ad"  \
			--output-report-dir "artifacts/scib_report2/" \




# python3 /opt/scripts/main/plot_feats_and_groups.py \
# 			--working-dir "$(pwd)" \
# 			--adata-input "testing/07_test_merged_adata_filtered_normalized_integrated_clustered_annotated.h5ad" \
# 			--group ~{sep=',' groups} \
# 			--output-group-umap-plot-prefix "~{cohort_id}" \
# 			--feature ~{sep=',' features} \
# 			--output-feature-umap-plot-prefix "~{cohort_id}"
