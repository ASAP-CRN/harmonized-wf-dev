
working_dir <- '/data/CARD_singlecell/harmony-rna/'; setwd(working_dir)

source('scripts/main/load_packages.r')


object <- readRDS(snakemake@input[['seurat_object']])
        
object <- object %>% subset(
    
    subset=
            nCount_RNA %between% c(500, 2e5) &
            nFeature_RNA %between% c(500, 12500) &
            percent.mt < 0.5 & percent.rb < 0.5,
    cells=
            fread(snakemake@input[['metadata']])[
                sample %chin% Project(object) & 
                predicted_gmm_doublets %chin% 'singlet', cells
            ]  
    )

saveRDS(object, snakemake@output[['seurat_object']])
