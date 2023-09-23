# 
# 
# 
# ototypes
# 
# 
# filtered counts: gs://asap-raw-data-team-lee/cellranger_counts/SN/SN_0009_PD_count_filtered_feature_bc_matrix.h5
# 
# 
# raw counts: gs://asap-raw-data-team-lee/cellranger_counts/SN/SN_0009_PD_count_raw_feature_bc_matrix.h5

reticulate::conda_create(envname = "ret_py", python_version = 3.9)

# install SciPy
reticulate::conda_install("ret_py", "scrublet") # import SciPy (it will be automatically discovered in "r-reticulate")
install.packages('SoupX')


library(Seurat)
library(SeuratData)
library(ggplot2)


harmony_repo <- '~/Projects/SingleCell/Harmony-RNA-Workflow/'


reticulate::use_condaenv("ret_py")

source(paste0(harmony_repo, 'scripts/main/load_packages.r'))
reticulate::source_python(paste0(harmony_repo,'scripts/utility/scrublet_py.py'))

source(paste0(harmony_repo, "scripts/utility/soupx.r"))

#------------ part 1

ndataset <- 'team-Lee-samples.csv'

setwd('~/Projects/ASAP/harmonized-wf-dev/examples')

input_table = read.csv('team-Lee-samples.csv')


datasets = input_table$sample
batch = input_table$batch
dataset = datasets[47]

# dataset = "SN/SN_0009_PD"


raw_fn <- '../data/team-lee/cellranger_counts/SN/SN_0009_PD_count_raw_feature_bc_matrix.h5'
filt_fn <- '../data/team-lee/cellranger_counts/SN/SN_0009_PD_count_filtered_feature_bc_matrix.h5'


harmony_rna_wf_preprocess <- function(input_row) {
 
  sample_root <-  paste0("../data/team-lee/cellranger_counts/", input_row$subdir, "/", input_row$sample)
  
  raw_fname <-  paste0(sample_root,"_count_raw_feature_bc_matrix.h5")
  filt_fname <- paste0(sample_root,"_count_filtered_feature_bc_matrix.h5")
  
  
  raw.counts <- Read10X_h5(raw_fname)
  filtered.counts <- Read10X_h5(filt_fname)
  
  soup_rate=0.10
  
  adj.matrix <- suppressWarnings(SoupCorrect(raw.counts, filtered.counts, contamination_rate=soup_rate))
  
  object <- CreateSeuratObject(adj.matrix, min.cells=0, min.features=0, project=input_row$sample)
  
  
  object[['percent.mt']] <- PercentageFeatureSet(object, pattern='^MT-')
  object[['percent.rb']] <- PercentageFeatureSet(object, pattern='^RP[SL]')
  
  m <- copy(object@meta.data)
  setDT(m, keep.rownames='cells')
  
  m[,
    `:=` (
      sample=input_row$sample,
      batch=input_row$batch
    )
  ]
  
  batch <- m[, batch]
  sample <- m[, sample]
  
  names(sample) <- names(batch) <- m[, cells]
  
  
  doublet_rate <- (ncol(object) / 1000) * 0.008
  
  object <- object %>% 
    AddMetaData(metadata=factor(batch), col.name='batch') %>% 
    AddMetaData(metadata=factor(sample), col.name='sample') %>% 
    scrublet(n_prin_comps=30, expected_doublet_rate=doublet_rate) 
  
  
  outfnm = paste0("../objects/seurat_object_",input_row$sample,"_preprocessed_01.rds")
  saveRDS(object, outfnm)
  return(outfnm)
  
}

# simulate snakemake with lapply
object.list <- lapply(split(input_table,1:nrow(input_table)), harmony_rna_wf_preprocess)



# get_file_names <- function(input_row) {
#   outfnm = paste0("../objects/seurat_object_",input_row$sample,"_preprocessed_01.rds")
#   return(outfnm)
#   
# }
# 
# object.list <- lapply(split(input_table,1:nrow(input_table)), get_file_names)



#------------ part 2
#----- Count Matrix
# seurat_counts <- lapply(full_path, function(path_i){
#   print(path_i)
#   count <- Read10X_h5(path_i)
#   # gene.column=2, # col1 Ensemble, col2 Symbol
#   # unique.features=TRUE,
#   # strip.suffix=TRUE) 
#   seurat_obj <- CreateSeuratObject(counts = count, 
#                                    project = path_i,
#                                    min.cells = 10,
#                                    min.genes = 200)
# })


objects <- future.apply::future_lapply(object.list, readRDS)


m <- rbindlist(lapply(objects, function(object) {
  m <- copy(object@meta.data); setDT(m, keep.rownames='cells')
}))

# cutoff <- NormalMixCutoff(mixtools::normalmixEM(m[, doublet_scores], k=2))
cutoff <- 0.06
metadata_out='unfiltered_metadata.csv'
project_name = 'testing00'
m[, `:=` (
  project=rep(project_name, nrow(m)),
  predicted_gmm_doublets=fifelse(doublet_scores < cutoff, 'singlet', 'doublet')
)]

fwrite(x=m, file=metadata_out)
