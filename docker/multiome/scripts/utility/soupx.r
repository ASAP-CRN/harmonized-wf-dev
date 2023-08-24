

SoupCorrect <- function(raw.matrix, filt.matrix, contamination_rate=NULL) {
  
  object  <- CreateSeuratObject(counts=filt.matrix[['Gene Expression']])
  soup.channel  <- SoupX::SoupChannel(raw.matrix[['Gene Expression']], filt.matrix[['Gene Expression']])
  
  object <- object %>% 
    SCTransform(verbose=FALSE) %>% 
    RunPCA(verbose=FALSE) %>% 
    RunUMAP(dims=1:30, verbose=FALSE) %>% 
    FindNeighbors(dims=1:30, verbose=FALSE) %>% 
    FindClusters(verbose=FALSE, algorithm=3, resolution=0.8)
  
  meta <- object@meta.data
  umap <- object@reductions$umap@cell.embeddings
  soup.channel <- SoupX::setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel <- SoupX::setDR(soup.channel, umap)
  
  if (is.null(contamination_rate)) {
    soup.channel  <- SoupX::autoEstCont(soup.channel, forceAccept=TRUE, doPlot=FALSE)
  } else {
    soup.channel <- SoupX::setContaminationFraction(soup.channel, contamination_rate, forceAccept=TRUE)
  }

  adj.matrix  <- SoupX::adjustCounts(soup.channel)

  return(adj.matrix)
  
}

  