working_dir <- '/data/CARD_singlecell/harmony-rna/'; setwd(working_dir)

source('scripts/main/load_packages.r')
arguments <- commandArgs(trailingOnly=TRUE)

future::plan('multicore', workers=as.numeric(tail(arguments, n=1)))
options(future.globals.maxSize=ngbs * 1000 * 1024^2)


input.files <- Filter(Negate(is.null), lapply(arguments, function(argument) {
    if (argument %like% '_03.rds') argument
}))


object.list <- future.apply::future_lapply(input.files, readRDS)

top.genes <- object.list %>% SelectIntegrationFeatures(nfeatures=5000)

message('merging samples...')
object <- merge(x=object.list[[1]], y=object.list[-1])
message('Done.')

batch <- 'sample'; noise <- c('percent.mt', 'percent.rb', 'nFeature_RNA', 'nCount_RNA', 'doublet_scores', 'G2M.Score', 'S.Score')


object <- object %>% 

            ScaleData(vars.to.regress=noise, features=top.genes, verbose=FALSE) %>% 
            RunPCA(features=top.genes, npcs=50, verbose=FALSE) %>% 
            
            harmony::RunHarmony(reduction='pca', group.by.vars=batch, max.iter.harmony=50, reduction.save='harmony')


saveRDS(object, tail(arguments, n=2)[1])


