# All integration methods to compare, containing class I: Raw, fastMNN, Seurat, Harmony, LIGER, BBKNN, Conos, and scVI;
# and class II: scANVI, scPoli. 
library(reticulate)
use_condaenv("/home/server/anaconda3/envs/zy/bin/python")
anndata = import("anndata",convert=FALSE)
bbknn = import("bbknn", convert=FALSE)
sc = import("scanpy",convert=FALSE)
library(batchelor)
library(conos)
library(harmony)
library(rliger)
library(Seurat)
library(scater)
library(irlba)
library(BiocNeighbors)
library(SingleCellExperiment)
library(Matrix)
library(umap)
library(dplyr)
library(Rcpp)

# input data: (1) list of filtered gene by cell expression matrix; (2) data frame of meta information.

# batches
# meta
# vargenes

###################################################################################
### Raw ###
###################################################################################
run_Raw = function(batches, meta, is.normalize = TRUE, vargs = NULL, out.npcs = 30) {
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  # modified rownames
  if (!is.null(vargs)) vargs = rownames(batch_seurat)[which(rownames(batches[[1]]) %in% vargs)]
  if (is.normalize == TRUE) {
    batch_seurat <- NormalizeData(batch_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  VariableFeatures(batch_seurat) = rownames(batch_seurat)
  if (!is.null(vargs)) VariableFeatures(batch_seurat) = vargs
  batch_seurat = ScaleData(batch_seurat)
  batch_seurat = RunPCA(object = batch_seurat, features = VariableFeatures(batch_seurat), npcs = out.npcs, verbose = FALSE)
  raw_res = t(as.data.frame(batch_seurat@reductions$pca@cell.embeddings))
  colnames(raw_res) = rownames(meta)
  return(raw_res)
}

###################################################################################
### fastMNN ###
###################################################################################
# k (20 by default)
run_fastMNN = function(batches, meta, is.normalize = TRUE, nfeatures = 2000, vargs = NULL,
                    k = 20, out.npcs = 30)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  # modified rownames
  if (!is.null(vargs)) vargs = rownames(batch_seurat)[which(rownames(batches[[1]]) %in% vargs)]
  if (is.normalize == TRUE) {
    batch_seurat <- NormalizeData(batch_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  batch_seurat <- FindVariableFeatures(object = batch_seurat, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  if (!is.null(vargs)) VariableFeatures(batch_seurat) = vargs
  vargenes <- batch_seurat@assays[["RNA"]]@var.features

  # log-normalized data matrices as input
  data_list = lapply(1:length(batches), function(i) batch_seurat@assays[["RNA"]]@data[vargenes, colnames(batches[[i]])])
  ##########################################################
  # run fastMNN
  t1 = Sys.time()
  out_mnn_total = do.call(batchelor::fastMNN, c(data_list, k = k, d = out.npcs))
  t2 = Sys.time()
  print(t2-t1)
  
  fastmnn_res = t(out_mnn_total@assays@data@listData[["reconstructed"]]@seed@components)
  colnames(fastmnn_res) = rownames(meta)
  return(fastmnn_res)
}

###################################################################################
### Seurat ###
###################################################################################
# k.filter (200 by default)
run_Seurat = function(batches, meta, is.normalize = TRUE, nfeatures = 2000, vargs = NULL, reduction = c("cca", "rpca", "rlsi"),
                      k.filter = 200, out.npcs = 30)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  # modified rownames
  if (!is.null(vargs)) vargs = rownames(batch_seurat)[which(rownames(batches[[1]]) %in% vargs)]
  batch_list = SplitObject(batch_seurat, split.by = "Batch")
  for(i in 1:length(batch_list)) {
    if(is.normalize == TRUE) {
      batch_list[[i]] <- NormalizeData(batch_list[[i]],normalization.method = "LogNormalize", scale.factor = 10000)
    }
    batch_list[[i]] <- FindVariableFeatures(batch_list[[i]], selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  }
  
  # run Seurat 
  t1 = Sys.time()
  if (is.null(vargs)) {
    features <- SelectIntegrationFeatures(batch_list)
    if (reduction == "rpca") {
      for(i in 1:length(batch_list)) {
        batch_list[[i]] <- ScaleData(batch_list[[i]], features = features, verbose = FALSE)
        batch_list[[i]] <- RunPCA(batch_list[[i]], features = features, verbose = FALSE)
      }
    }
    cell_anchors = FindIntegrationAnchors(object.list = batch_list, k.filter = k.filter, reduction = reduction)
  } else {
    if (reduction == "rpca") {
      for(i in 1:length(batch_list)) {
        batch_list[[i]] <- ScaleData(batch_list[[i]], features = vargs, verbose = FALSE)
        batch_list[[i]] <- RunPCA(batch_list[[i]], features = vargs, verbose = FALSE)
      }
    }
    cell_anchors = FindIntegrationAnchors(object.list = batch_list, k.filter = k.filter, reduction = reduction, anchor.features = vargs)
  }
  rm(batch_list, batch_seurat)
  gc()
  batch_correct = IntegrateData(anchorset = cell_anchors)
  t2 = Sys.time()
  print(t2-t1)
  
  DefaultAssay(batch_correct) = "integrated"
  batch_correct = ScaleData(object = batch_correct)
  batch_correct = RunPCA(object = batch_correct, npcs = out.npcs, verbose = FALSE)
  seurat_res = t(as.data.frame(batch_correct@reductions$pca@cell.embeddings))
  colnames(seurat_res) = rownames(meta)
  return(seurat_res)
}

###################################################################################
### Harmony ###
###################################################################################
# key parameters: group.by.vars ("Batch" by default), theta (2 by default)
run_Harmony = function(batches, meta, is.normalize = TRUE, nfeatures = 2000, group.by.vars = "Batch", theta = 2, vargs = NULL, 
                        out.npcs = 30)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  # modified rownames
  if (!is.null(vargs)) vargs = rownames(batch_seurat)[which(rownames(batches[[1]]) %in% vargs)]
  if (is.normalize == TRUE) {
    batch_seurat <- NormalizeData(object = batch_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  batch_seurat <- FindVariableFeatures(object = batch_seurat, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  if (!is.null(vargs)) VariableFeatures(batch_seurat) = vargs
  batch_seurat <- ScaleData(object = batch_seurat)
  batch_seurat = RunPCA(object = batch_seurat, npcs = out.npcs, features = VariableFeatures(object = batch_seurat))
  
  #run Harmony
  t1 = Sys.time()
  batch_seurat = RunHarmony(object = batch_seurat, group.by.vars = group.by.vars, theta = theta, plot_convergence = TRUE, 
                             nclust = 50, max.iter.cluster = 100)
  t2 = Sys.time()
  print(t2-t1)
  
  harmony_res = t(as.data.frame(batch_seurat@reductions[["harmony"]]@cell.embeddings))
  colnames(harmony_res) = rownames(meta)
  return(harmony_res)
}

###################################################################################
### LIGER ###
###################################################################################
# k (20 by default), lambda (5 by default)
run_LIGER = function(batches, meta, is.normalize = TRUE,
                      k = 20, lambda = 5, nrep = 1) 
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  if (is.normalize == TRUE) {
    batch_seurat = NormalizeData(batch_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  batch_list = SplitObject(batch_seurat, split.by = "Batch")

  liger_object = createLiger(batches, remove.missing = F)
  liger_object = rliger::normalize(liger_object)
  if (is.normalize == TRUE) {
    liger_object@norm.data = lapply(batch_list, function(x) x@assays[["RNA"]]@data)
  } else {
    liger_object@norm.data = liger_object@raw.data
  }
  liger_object = rliger::selectGenes(liger_object, var.thresh = 0.1)
  liger_object = rliger::scaleNotCenter(liger_object, remove.missing = F)
  
  ##########################################################
  # run LIGER
  t1 = Sys.time()
  liger_object = rliger::optimizeALS(liger_object, k = k, lambda = lambda, nrep = nrep)
  liger_object = rliger::quantile_norm(liger_object)
  t2 = Sys.time()
  print(t2-t1)
  
  liger_res = t(liger_object@H.norm)
  colnames(liger_res) = rownames(meta)
  return(liger_res)
}

###################################################################################
### Conos ###
###################################################################################
run_Conos = function(batches, meta, is.normalize = TRUE, nfeatures = 2000, vargs = NULL,
                     k = 20, out.npcs = 30)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  # modified rownames
  if (!is.null(vargs)) vargs = rownames(batch_seurat)[which(rownames(batches[[1]]) %in% vargs)]
  batch_list = SplitObject(batch_seurat, split.by = "Batch")
  for(i in 1:length(batch_list)) {
    if(is.normalize == TRUE) {
      batch_list[[i]] <- NormalizeData(batch_list[[i]],normalization.method = "LogNormalize", scale.factor = 10000)
    }
    batch_list[[i]] <- FindVariableFeatures(batch_list[[i]], selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    if (!is.null(vargs)) VariableFeatures(batch_list[[i]]) = vargs
  }
  batch_list = lapply(batch_list, function(x) ScaleData(x) %>% RunPCA())
  
  # Construct Conos object
  con = Conos$new(batch_list, n.cores=1)
  
  # run Conos
  t1 = Sys.time()
  # Build joint graph
  con$buildGraph(k = k)
  t2 = Sys.time()
  print(t2-t1)
  
  # Find communities
  con$findCommunities()
  # Generate embedding
  conos_res = t(con$embedGraph(target.dims = out.npcs, method = "largeVis", verbose = FALSE)[Reduce('c', lapply(batch_list, colnames)),])
  colnames(conos_res) = rownames(meta)
  return(conos_res)
}

###################################################################################
### BBKNN ###
###################################################################################
# return UMAP embedding
run_BBKNN = function(batches, meta, is.normalize = TRUE, nfeatures = 2000, vargs = NULL)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  # modified rownames
  if (!is.null(vargs)) vargs = rownames(batch_seurat)[which(rownames(batches[[1]]) %in% vargs)]
  if (is.normalize == TRUE) {
    batch_seurat <- NormalizeData(batch_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  batch_seurat <- FindVariableFeatures(object = batch_seurat, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  if (!is.null(vargs)) VariableFeatures(batch_seurat) = vargs
  vargenes <- batch_seurat@assays[["RNA"]]@var.features
  batch_seurat <- ScaleData(object = batch_seurat)
  batch_seurat = RunPCA(object = batch_seurat, features = VariableFeatures(object = batch_seurat))
  data_pca = batch_seurat@reductions$pca@cell.embeddings
  rownames(data_pca) = rownames(meta)
  adata = anndata$AnnData(X=data_pca, obs=meta)
  sc$tl$pca(adata)
  adata$obsm$X_pca = data_pca

  # run BBKNN
  t1 = Sys.time()
  bbknn$bbknn(adata,batch_key='Batch')
  t2 = Sys.time()
  print(t2-t1)
  
  sc$tl$umap(adata)
  bbknn_res = t(py_to_r(adata$obsm[["X_umap"]]))
  rownames(bbknn_res) = c("UMAP_1", "UMAP_2")
  colnames(bbknn_res) = rownames(meta)
  return(bbknn_res)
}