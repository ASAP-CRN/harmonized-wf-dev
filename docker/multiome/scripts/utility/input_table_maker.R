


require(data.table)

batch_list <- list.files()

batches <- sapply(batch_list, function(batch) {
    list.files(paste0(batch, '/Multiome'))
}, simplify=FALSE)


batches <- sapply(batches, data.table, simplify=FALSE)

batches <- sapply(seq_along(batches), function(batch_num) {
    batches[[batch_num]][, batch := names(batches)[batch_num]]
}, simplify=FALSE)

samples <- rbindlist(batches)
setnames(samples, old='V1', new='sample')





fwrite(samples, '/data/CARD_singlecell/snakemake_multiome/input/samples.csv')


fwrite(x, '/data/CARD_singlecell/multiome-test/input/random_samples.csv')
