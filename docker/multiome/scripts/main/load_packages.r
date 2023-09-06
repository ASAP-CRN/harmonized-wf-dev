script_dir <- if (exists("script_dir")) script_dir else "scripts"

PKGS <- c('Seurat', 'dplyr', 'ggplot2', 'data.table')

invisible(sapply(PKGS, require, character.only=TRUE))
invisible(lapply(list.files(paste0(script_dir, '/utility'), full.names=TRUE, pattern='\\.r$'), source))

data_path <- '/data/CARD_singlecell/Brain_atlas/NABEC_multiome/'

setDTthreads(threads=1)

ngbs <- 50
