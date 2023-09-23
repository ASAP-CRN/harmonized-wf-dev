#lee_testing

library("tidyverse")
#---------- Sample Info
#----- Covar

covar_csv <- "team-Lee-all_samples.csv"
covar <- read_csv(covar_csv) %>%
  # add CaseID column for merging with Banner demographic info
  mutate(CaseID = gsub(ID, pattern = "(.*)_([0-9][0-9]).*([0-9][0-9])", replacement = "\\2-\\3")) %>%
  # remove PMI for now -- we'll add this next w/ the other covars
  dplyr::select(-PMI) %>%
  mutate(redo = case_when(SEQ_ID %in% c("MFGHC1225", "MFGPD1441", "MFGPD1504", 
                                        "MFGPD1921", "MFGPD2005") ~ "redo",
                          TRUE ~ "no_redo"))

# # merge with extra demographic covar information from Banner (PMI, etc.)
# covar_extra <- rbind(read_csv("PD_ASAP_Sample_batch_information_banner_cases.csv") %>% dplyr::mutate(neurological_dx_years = NA),
#                      read_csv("PD_ASAP_Sample_batch_information_banner_controls.csv")) %>%
#   dplyr::filter(!is.na(CaseID))
# covar <- left_join(covar, covar_extra, by = "CaseID")
# 
# knitr::kable(covar) %>% 
#   kableExtra::kable_styling(bootstrap_options = c("striped", 
#                                                   "hover", 
#                                                   "condensed", 
#                                                   "responsive"))

NameProject <- "PD_MFG_snRNAseq_with_redos"

#----- Count Names
NameCount <- "SEQ_ID"
NamesCount <- covar$SEQ_ID #gsub(covar$SEQ_ID, pattern = "(HC|PD)", replacement = "")
#----- Sample Names
NamesSample <- covar$SAMPLE #gsub(covar$SAMPLE, pattern = "_(HC|PD)", replacement = "")
# #---- Gene Names
# #(must be run after cellranger)
# Gene_Emsembl_Symbol <- read_tsv(
#   paste0(DirCounts,
#          NamesCount[1], 
#          "_count/outs/filtered_feature_bc_matrix/features.tsv.gz"),
#   col_names = c("Gene_Emsembl","Gene_Symbol","Assay_type"), 
#   col_types = "ccc")[,c(1:2)] %>% 
#   unique()

NameCount <- "SEQ_ID"
subdir = covar$subdir
DirCounts <- "cellranger_mRNA_premRNA/"

full_path <- covar$fullpath <- paste0("../data/team-lee/cellranger_counts/",covar$subdir,"/",NamesCount,"_count_filtered_feature_bc_matrix.h5")

#---------- Seurat_Object
#----- Count Matrix
seurat_counts <- lapply(full_path, function(path_i){
  print(path_i)
  count <- Read10X_h5(path_i)
                   # gene.column=2, # col1 Ensemble, col2 Symbol
                   # unique.features=TRUE,
                   # strip.suffix=TRUE) 
  seurat_obj <- CreateSeuratObject(counts = count, 
                                   project = path_i,
                                   min.cells = 10,
                                   min.genes = 200)
})

#----- Merge Seurat Objects
seurat_combined <- merge(x = seurat_counts[[1]], 
                         y = seurat_counts[2:length(seurat_counts)],
                         add.cell.ids = unlist(NamesCount), 
                         project = NameProject)

rm(seurat_counts)


Is_Seurat <- function(
    seurat_object
) {
  if (!inherits(what = "Seurat", x = seurat_object)) {
    cli_abort(message = "{.code seurat_object} provided is not an object of class: Seurat.")
  }
}

Fetch_Meta <- function(
    object
) {
  # Pull meta data
  object_meta <- object_meta <- slot(object = object, name = "meta.data")
  
  return(object_meta)
}


Add_Mito_Ribo_Seurat <- function(
    seurat_object,
    species,
    mito_name = "percent_mito",
    ribo_name = "percent_ribo",
    mito_ribo_name = "percent_mito_ribo",
    mito_pattern = NULL,
    ribo_pattern = NULL,
    mito_features = NULL,
    ribo_features = NULL,
    ensembl_ids = FALSE,
    assay = NULL,
    overwrite = FALSE
) {
  # Accepted species names
  accepted_names <- data.frame(
    Mouse_Options = c("Mouse", "mouse", "Ms", "ms", "Mm", "mm"),
    Human_Options = c("Human", "human", "Hu", "hu", "Hs", "hs"),
    Marmoset_Options = c("Marmoset", "marmoset", "CJ", "Cj", "cj", NA),
    Zebrafish_Options = c("Zebrafish", "zebrafish", "DR", "Dr", "dr", NA),
    Rat_Options = c("Rat", "rat", "RN", "Rn", "rn", NA),
    Drosophila_Options = c("Drosophila", "drosophila", "DM", "Dm", "dm", NA),
    Macaque_Options = c("Macaque", "macaque", "Rhesus", "macaca", "mmulatta", NA)
  )
  
  
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)
  
  # Check name collision
  if (any(duplicated(x = c(mito_name, ribo_name, mito_ribo_name)))) {
    cli_abort(message = "One or more of values provided to {.code mito_name}, {.code ribo_name}, {.code mito_ribo_name} are identical.")
  }
  
  # Overwrite check
  if (mito_name %in% colnames(x = seurat_object@meta.data) || ribo_name %in% colnames(x = seurat_object@meta.data) || mito_ribo_name %in% colnames(x = seurat_object@meta.data)) {
    if (!overwrite) {
      cli_abort(message = c("Columns with {.val {mito_name}} and/or {.val {ribo_name}} already present in meta.data slot.",
                            "i" = "*To run function and overwrite columns set parameter {.code overwrite = TRUE} or change respective {.code mito_name}, {.code ribo_name}, and/or {.code mito_ribo_name}*")
      )
    }
    cli_inform(message = c("Columns with {.val {mito_name}} and/or {.val {ribo_name}} already present in meta.data slot.",
                           "i" = "Overwriting those columns as .code {overwrite = TRUE.}")
    )
  }
  
  # Checks species
  if (is.null(x = species)) {
    cli_abort(message = c("No species name or abbreivation was provided to {.code species} parameter.",
                          "i" = "If not using default species please set {.code species = other}.")
    )
  }
  
  # Set default assay
  assay <- assay %||% DefaultAssay(object = seurat_object)
  
  # Species Spelling Options
  mouse_options <- accepted_names$Mouse_Options
  human_options <- accepted_names$Human_Options
  marmoset_options <- accepted_names$Marmoset_Options
  zebrafish_options <- accepted_names$Zebrafish_Options
  rat_options <- accepted_names$Rat_Options
  drosophila_options <- accepted_names$Drosophila_Options
  macaque_options <- accepted_names$Macaque_Options
  
  # Check ensembl vs patterns
  if (ensembl_ids && species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options) && any(!is.null(x = mito_pattern), !is.null(x = ribo_pattern), !is.null(x = mito_features), !is.null(x = ribo_features))) {
    cli_warn(message = c("When using a default species and setting {.code ensembl_ids = TRUE} provided patterns or features are ignored.",
                         "*" = "Supplied {.code mito_pattern}, {.code ribo_pattern}, {.code mito_features}, {.code ribo_features} will be disregarded.")
    )
  }
  
  # Assign mito/ribo pattern to stored species
  if (species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options) && any(!is.null(x = mito_pattern), !is.null(x = ribo_pattern))) {
    cli_warn(message = c("Pattern expressions for included species are set by default.",
                         "*" = "Supplied {.code mito_pattern} and {.code ribo_pattern} will be disregarded.",
                         "i" = "To override defaults please supply a feature list for mito and/or ribo genes.")
    )
  }
  
  if (species %in% mouse_options) {
    mito_pattern <- "^mt-"
    ribo_pattern <- "^Rp[sl]"
  }
  if (species %in% human_options) {
    mito_pattern <- "^MT-"
    ribo_pattern <- "^RP[SL]"

  }
  
  # Check that values are provided for mito and ribo
  if (is.null(x = mito_pattern) && is.null(x = mito_features) && is.null(x = ribo_pattern) && is.null(x = ribo_pattern)) {
    cli_abort(message = c("No features or patterns provided for mito/ribo genes.",
                          "i" = "Please provide a default species name or pattern/features."))
  }
  
  # Retrieve ensembl ids if TRUE
  if (ensembl_ids) {
    mito_features <- Retrieve_Ensembl_Mito(species = species)
    ribo_features <- Retrieve_Ensembl_Ribo(species = species)
  }
  
  mito_features <- mito_features %||% grep(pattern = mito_pattern, x = rownames(x = seurat_object[[assay]]), value = TRUE)
  
  ribo_features <- ribo_features %||% grep(pattern = ribo_pattern, x = rownames(x = seurat_object[[assay]]), value = TRUE)
  
  # Check features are present in object
  length_mito_features <- length(x = intersect(x = mito_features, y = rownames(x = seurat_object[[assay]])))
  
  length_ribo_features <- length(x = intersect(x = ribo_features, y = rownames(x = seurat_object[[assay]])))
  
  # Check length of mito and ribo features found in object
  if (length_mito_features < 1 && length_ribo_features < 1) {
    cli_abort(message = c("No Mito or Ribo features found in object using patterns/feature list provided.",
                          "i" = "Please check pattern/feature list and/or gene names in object.")
    )
  }
  if (length_mito_features < 1) {
    cli_warn(message = c("No Mito features found in object using pattern/feature list provided.",
                         "i" = "No column will be added to meta.data.")
    )
  }
  if (length_ribo_features < 1) {
    cli_warn(message = c("No Ribo features found in object using pattern/feature list provided.",
                         "i" = "No column will be added to meta.data.")
    )
  }
  
  # Add mito and ribo columns
  if (length_mito_features > 0) {
    good_mito <- mito_features[mito_features %in% rownames(x = seurat_object)]
    seurat_object[[mito_name]] <- PercentageFeatureSet(object = seurat_object, features = good_mito, assay = assay)
  }
  if (length_ribo_features > 0) {
    good_ribo <- ribo_features[ribo_features %in% rownames(x = seurat_object)]
    seurat_object[[ribo_name]] <- PercentageFeatureSet(object = seurat_object, features = good_ribo, assay = assay)
  }
  
  # Create combined mito ribo column if both present
  if (length_mito_features > 0 && length_ribo_features > 0) {
    object_meta <- Fetch_Meta(object = seurat_object) %>%
      rownames_to_column("barcodes")
    
    object_meta <- object_meta %>%
      mutate({{mito_ribo_name}} := .data[[mito_name]] + .data[[ribo_name]])
    
    object_meta <- object_meta %>%
      select(all_of(c("barcodes", mito_ribo_name))) %>%
      column_to_rownames("barcodes")
    
    seurat_object <- AddMetaData(object = seurat_object, metadata = object_meta)
  }
  
  # return final object
  return(seurat_object)
}


#---------- Add Metadata
#----- Add Percentage of Mitocondrial and Ribosomal Counts
seurat_combined <- Add_Mito_Ribo_Seurat(seurat_object = seurat_combined, 
                                        species = "Human")
# View(seurat_combined@meta.data)
# orig.ident: contains the sample identity if known
# nCount_RNA: number of UMIs per cell
# nFeature_RNA: number of genes detected per cell
seurat_meta <- function(seurat_obj, name_count) {
  #----- Create metadata dataframe
  metadata <- seurat_obj@meta.data
  #----- Rename 
  names(metadata)[1:3] <- c(name_count,"nUMI","nGene")
  #----- Add Cell IDs
  metadata$cells <- rownames(metadata)
  #----- Add number of genes per UMI
  metadata$log10GenesPerUMI <- 
    log10(metadata$nGene) / log10(metadata$nUMI)
  #----- Compute percent mito ratio
  # metadata <- cbind(metadata,
  #                   PercentageFeatureSet(object = seurat_combined, 
  #                                        pattern = "^MT-") / 100)
  # names(metadata)[6] <- "mitoRatio"
  #----- Merge Covar
  metadata <- metadata %>% 
    left_join(covar, by = name_count)
  rownames(metadata) <- metadata$cells
  #----- Add metadata back to Seurat object
  seurat_obj@meta.data <- metadata
  return(seurat_obj)
}
seurat_combined <- seurat_meta(seurat_obj = seurat_combined,
                               name_count = NameCount)


#---------- Senescence
#----- Senescence Markers, p16 = CDKN2A, p21 = CDKN1A
#----- + extra senescence-relevant genes of interest: p19 = CDKN2D, uPAR = PLAUR
seurat_combined$CDKN2A_Count <- 
  seurat_combined@assays$RNA@counts[rownames(seurat_combined) == "CDKN2A",]
seurat_combined$CDKN1A_Count <- 
  seurat_combined@assays$RNA@counts[rownames(seurat_combined) == "CDKN1A",]
seurat_combined$CDKN2D_Count <- 
  seurat_combined@assays$RNA@counts[rownames(seurat_combined) == "CDKN2D",]
seurat_combined$GPNMB_Count <- 
  seurat_combined@assays$RNA@counts[rownames(seurat_combined) == "GPNMB",]
seurat_combined$PLAUR_Count <- 
  seurat_combined@assays$RNA@counts[rownames(seurat_combined) == "PLAUR",]

#---------- Save Data
#----- Remove Data
#----- Save Data
saveRDS(seurat_combined, file="seurat_combined.rds")

#-----------------Part 2
#---------- QC Filtering
#----- Filtering
nUMI_low = 1000
nUMI_high = 20000
nGene_low = 700
nGene_high = 6000
percent_mito_high = 5
percent_ribo_high = 30
log10GenesPerUMI_low = 0.85
seurat_filtered <- subset(x = seurat_combined, 
                          subset = (nUMI >= nUMI_low) & 
                            (nUMI <= nUMI_high) & 
                            (nGene >= nGene_low) & 
                            (nGene <= nGene_high) &
                            (percent_mito < percent_mito_high) & 
                            (percent_ribo < percent_ribo_high) &
                            (log10GenesPerUMI > log10GenesPerUMI_low))

#---- save rds
saveRDS(seurat_filtered, file="seurat_filtered.rds")

# ------------- Part 3

#---------- Harmony
#----- Batch Correction
seurat_batch <- function(seurat_obj, number_features, sample, dimensions) {
  #----- Normalize the counts
  harmony <- NormalizeData(seurat_obj)
  #----- Identify the most variable genes (THIS IS FAILING)
  harmony <- FindVariableFeatures(object = harmony, 
                                  selection.method = "vst",
                                  nfeatures = number_features, 
                                  verbose = FALSE)
  #----- Scale the counts
  harmony <- ScaleData(object = harmony)
  #----- Perform PCA
  harmony <- RunPCA(object = harmony,
                    verbose = FALSE)
  # #----- Score cells for cell cycle
  # harmony <- CellCycleScoring(object = harmony, 
  #                             g2m.features = g2m_genes, 
  #                             s.features = s_genes)
  #---- Harmony Integration
  harmony <- RunHarmony(object = harmony, 
                        group.by.vars = sample)
  #---- Perform UMAP
  harmony <- RunUMAP(object = harmony, 
                     reduction = "harmony", 
                     dims = dimensions)
  return(harmony)
} 
seurat_harmony <- seurat_batch(seurat_obj = seurat_filtered, 
                               number_features = 2000,
                               sample = "SAMPLE",
                               dimensions = 1:30)  

number_features = 2000
sample = "SAMPLE"
dimensions = 1:30
harmony <- NormalizeData(seurat_filtered)
#----- Identify the most variable genes (THIS IS FAILING)
harmony <- FindVariableFeatures(object = harmony, 
                                selection.method = "vst",
                                nfeatures = number_features, 
                                verbose = FALSE)
#----- Scale the counts
harmony <- ScaleData(object = harmony)
#----- Perform PCA
harmony <- RunPCA(object = harmony,
                  verbose = FALSE)
# #----- Score cells for cell cycle
# harmony <- CellCycleScoring(object = harmony, 
#                             g2m.features = g2m_genes, 
#                             s.features = s_genes)
#---- Harmony Integration
harmony <- RunHarmony(object = harmony, 
                      group.by.vars = sample)
#---- Perform UMAP
harmony <- RunUMAP(object = harmony, 
                   reduction = "harmony", 
                   dims = dimensions)


# ElbowPlot(object = seurat_harmony, ndims = 50)
# DimPlot(seurat_harmony, split.by = "sample", ncol = 6) + NoLegend()

#---------- Data
#----- Save Data
saveRDS(seurat_harmony, file="seurat_harmony.rds")
#----- Remove Data
rm(seurat_filtered)
