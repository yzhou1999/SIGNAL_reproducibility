library(Seurat)
library(symphony)
library(harmony)
library(tibble)
library(matrixStats)
scale_data = SIGNAL:::scale_data

# vargenes must in the rownames of seurat object.
run_Seurat_mapping = function(reference.dataset, reference.meta, query.dataset, query.meta, cts = "CellType",
                               vars = "batchlb", nfeatures = 2000, vargenes = NULL, is.normalize = T, npcs = 20) {
  message('Build Seurat reference.')
  t1 = Sys.time()
  ref_seurat = CreateSeuratObject(do.call(cbind, reference.dataset), meta.data = reference.meta, 
                                     min.cells = 0, min.features = 0, project = "reference")
  ref_seurat = SplitObject(ref_seurat, split.by = vars)
  for(i in 1:length(ref_seurat)) {
    # normalization
    if(is.normalize == TRUE) {
      ref_seurat[[i]] = NormalizeData(ref_seurat[[i]],normalization.method = "LogNormalize", scale.factor = 10000)
    }
    # select hvgs
    if(!is.null(vargenes)) {
      VariableFeatures(ref_seurat[[i]]) = vargenes
    } else {
      ref_seurat[[i]] = FindVariableFeatures(ref_seurat[[i]], selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    }
    ref_seurat[[i]] = ScaleData(ref_seurat[[i]], verbose = FALSE)
    ref_seurat[[i]] = RunPCA(ref_seurat[[i]], verbose = FALSE)
  }
  if(!is.null(vargenes)) {
    vargs = vargenes
    cell_anchors = FindIntegrationAnchors(object.list = ref_seurat, anchor.features = vargenes)
  } else {
    vargs = SelectIntegrationFeatures(ref_seurat)
    cell_anchors = FindIntegrationAnchors(object.list = ref_seurat, anchor.features = nfeatures)
  }
  
  ref_seurat = IntegrateData(anchorset = cell_anchors)
  t2 = Sys.time()
  print(t2-t1)
  
  DefaultAssay(ref_seurat) = "integrated"
  ref_seurat = ScaleData(ref_seurat, verbose = FALSE)
  ref_seurat = RunPCA(ref_seurat, npcs = npcs, verbose = FALSE)
  ref_seurat = RunUMAP(ref_seurat, reduction = "pca", dims = 1:20, verbose = FALSE)
  
  message('Do Seurat mapping.')
  t1 = Sys.time()
  query.dataset = lapply(query.dataset, function(x) x[intersect(vargs, rownames(x)), ])
  query_seurat = CreateSeuratObject(do.call(cbind, query.dataset), meta.data = query.meta, 
                                   min.cells = 0, min.features = 0, project = "query")
  # set hvgs
  VariableFeatures(query_seurat) = rownames(query_seurat)
  query_seurat = SplitObject(query_seurat, split.by = vars)
  
  # mapping and label transfer
  anchors = list()
  pred_list = list()
  for(i in 1:length(query_seurat)) {
    if(is.normalize == TRUE) {
      query_seurat[[i]] = NormalizeData(query_seurat[[i]],normalization.method = "LogNormalize", scale.factor = 10000)
    }
    
    anchors[[i]] = FindTransferAnchors(
      reference = ref_seurat,
      query = query_seurat[[i]],
      k.filter = NA,
      reduction = "pcaproject",
      reference.reduction = "pca",
      dims = 1:npcs,
      verbose = FALSE
    )
    
    pred_list[[i]] = TransferData(anchorset = anchors[[i]], 
                                  refdata = ref_seurat@meta.data[[cts]],
                                  dims = 1:npcs,
                                  verbose = FALSE)[["predicted.id"]]
    
    query_seurat[[i]] = MapQuery(
      anchorset = anchors[[i]], 
      query = query_seurat[[i]],
      reference = ref_seurat, 
      refdata = list(cell_type = cts),
      reference.reduction = "pca",
      verbose = FALSE
    )
  }
  
  all_query = merge(query_seurat[[1]], query_seurat[2:length(query_seurat)], merge.dr = "ref.pca")
  all_query = all_query@reductions[["ref.pca"]]@cell.embeddings
  
  pred_labels = Reduce('c', pred_list)
  t2 = Sys.time()
  print(t2-t1)
  
  mapped_res = rbind(ref_seurat@reductions[["pca"]]@cell.embeddings, all_query) %>% t()
  return(list("pred_labels" = pred_labels, "mapped_res" = mapped_res))
}

run_Symphony_mapping = function(reference.dataset, reference.meta, query.dataset, query.meta, cts = "CellType",
                                 nfeatures = 2000, vargenes = NULL, is.normalize = T, npcs = 20, pred_k = 5,
                                 ref.vars = c('donor', 'batchlb'), theta = c(2, 4), query.vars = 'batchlb') {
  message('Build Symphony reference.')
  t1 = Sys.time()
  reference.dataset = Reduce(cbind2, reference.dataset)
  if (!is.null(vargenes)) {
    reference.dataset = reference.dataset[vargenes, ]
  } else {
    vargenes = vargenes_vst(reference.dataset, groups = as.character(reference.meta$batchlb), topn = 1000)
    reference.dataset = reference.dataset[vargenes, ]
  }
  vargenes_means_sds = tibble(symbol = vargenes, mean = Matrix::rowMeans(reference.dataset))
  vargenes_means_sds$stddev = rowSds(as.matrix(reference.dataset), center = vargenes_means_sds$mean)
  exp_scaled = scale_data(reference.dataset, row.means = vargenes_means_sds$mean, row.sds = vargenes_means_sds$stddev)
  s = irlba::irlba(exp_scaled, nv = npcs)
  Z_pca = diag(s$d) %*% t(s$v)
  loadings = s$u
  
  message('Harmonize reference')
  # run Harmony
  set.seed(0)
  harmObj = harmony::HarmonyMatrix(
    data_mat = t(Z_pca),
    meta_data = reference.meta,
    theta = theta,
    vars_use = ref.vars,
    nclust = 100,
    max.iter.harmony = 10,
    tau = 5,
    return_object = TRUE,
    do_pca = FALSE
  )
  reference = symphony::buildReferenceFromHarmonyObj(
    harmObj,
    reference.meta,
    vargenes_means_sds, 
    loadings,
    verbose = TRUE,
    do_umap = FALSE)
  t2 = Sys.time()
  print(t2-t1)
  
  message('Do Symphony mapping.')
  t1 = Sys.time()
  query.dataset = Reduce(cbind2, query.dataset)
  query_combined = mapQuery(exp_query = query.dataset, 
                            metadata_query = query.meta,
                            ref_obj = reference,
                            vars = query.vars, 
                            do_normalize = FALSE, do_umap = FALSE)
  query_combined = knnPredict(query_combined, reference, reference[["meta_data"]][[cts]], k = pred_k)
  
  pred_labels = as.character(query_combined[["meta_data"]][["cell_type_pred_knn"]])
  all_Z = rbind(t(reference$Z_corr), t(query_combined$Z))
  t2 = Sys.time()
  print(t2-t1)
  
  mapped_res = t(all_Z)
  return(list("pred_labels" = pred_labels, "mapped_res" = mapped_res))
}