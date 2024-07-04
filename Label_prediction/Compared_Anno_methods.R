# PGTV, SciBet, scmap-cluster, scmap-cell, SingleR, scPred
library(Seurat)
library(scibet)
library(scmap)
library(SingleR)
library(scPred)
library(SIGNAL)

pgtv_pred <- function(ref_data, ref_group, query_data, query_group, do.cosine = T, do.scale = F, 
                      conf.score = F, proj.dims = 30, eta = 0.25, smooth.k = 5) {
  t1 = Sys.time()
  pred_lab <- Run.LabelTransfer.Single(ref_data, ref_group, query_data, do.cosine, do.scale, conf.score, 
                                       proj.dims, eta, smooth.k)
  t2 = Sys.time()
  print(t2 - t1)
  if (!conf.score) {
    true_lab <- query_group
    res <- evaluate(true_lab, pred_lab)
  } else {
    true_lab <- query_group
    res <- list("preds" = evaluate(true_lab, pred_lab$Prediction),
                "conf.scores" = pred_lab$Confidence)
  }
  return(res)
}

scibet_pred <- function(ref_data, ref_group, query_data, query_group) {
  t1 = Sys.time()
  true_lab <- query_group
  labels <- ref_group
  ref_expr <- as.data.frame(t(ref_data))
  ref_expr$label <- labels
  query_expr <- t(query_data)
  pred_lab <- SciBet(ref_expr, query_expr)
  t2 = Sys.time()
  print(t2 - t1)
  res <- evaluate(true_lab, pred_lab)
  return(res)
}

scmapcluster_pred <- function(ref_data, ref_group, query_data, query_group) {
  t1 = Sys.time()
  true_lab <- query_group
  celltype <- as.data.frame(as.factor(ref_group))
  colnames(celltype) <- "cell_type1"
  ref_expr <- SingleCellExperiment(assays = list(normcounts = as.matrix(ref_data)), colData = celltype)
  logcounts(ref_expr) <- as.matrix(ref_data)
  rowData(ref_expr)$feature_symbol <- rownames(ref_expr)
  ref_expr <- ref_expr[!duplicated(rownames(ref_expr)), ]
  ref_expr <- selectFeatures(ref_expr, suppress_plot = FALSE)
  ref_expr <- indexCluster(ref_expr)
  
  query_expr <- SingleCellExperiment(assays = list(normcounts = as.matrix(query_data)))
  logcounts(query_expr) <- as.matrix(query_data)
  rowData(query_expr)$feature_symbol <- rownames(query_expr)
  
  scmap_pre <- scmapCluster(
    query_expr,
    list(
      ref_expr@metadata$scmap_cluster_index
    )
  )
  t2 = Sys.time()
  print(t2 - t1)
  scmap_results <- scmap_pre$scmap_cluster_labs
  res <- evaluate(true_lab, scmap_results)
  return(res)
}

scmapcell_pred <- function(ref_data, ref_group, query_data, query_group) {
  t1 = Sys.time()
  true_lab <- query_group
  celltype <- as.data.frame(as.factor(ref_group))
  colnames(celltype) <- "celltype1"
  ref_expr <- SingleCellExperiment(assays = list(normcounts = as.matrix(ref_data)), colData = celltype)
  logcounts(ref_expr) <- as.matrix(ref_data)
  rowData(ref_expr)$feature_symbol <- rownames(ref_expr)
  ref_expr <- ref_expr[!duplicated(rownames(ref_expr)), ]
  ref_expr <- selectFeatures(ref_expr, suppress_plot = FALSE)
  ref_expr <- indexCell(ref_expr)
  
  query_expr <- SingleCellExperiment(assays = list(normcounts = as.matrix(query_data)))
  logcounts(query_expr) <- as.matrix(query_data)
  rowData(query_expr)$feature_symbol <- rownames(query_expr)
  
  scmap_pre <- scmapCell(
    query_expr,
    list(
      ref_expr@metadata$scmap_cell_index
    )
  )
  t2 = Sys.time()
  print(t2 - t1)
  cells <- as.data.frame(scmap_pre[[1]][[1]])
  scores <- as.data.frame(scmap_pre[[1]][[2]])
  cells <- t(cells)
  scores <- t(scores)
  max_scores <- max.col(scores)
  pre_cells <- c()
  for (i in 1:nrow(cells)) {
    pre_cells[i] <- cells[i, max_scores[i]]
  }
  pre_cells <- as.numeric(pre_cells)
  scmap_results <- list()
  for (j in 1:length(pre_cells)) {
    scmap_results[[j]] <- celltype[pre_cells[j], 1]
  }
  scmap_results <- as.character(unlist(scmap_results))
  res <- evaluate(true_lab, scmap_results)
  return(res)
}

singler_pred <- function(ref_data, ref_group, query_data, query_group) {
  t1 = Sys.time()
  cts <- unique(ref_group)
  ref.data <- sapply(cts, function(ct) {
    if (length(which(ref_group == ct)) > 1) {
      rowMeans(ref_data[, which(ref_group == ct)])
    } else {
      ref_data[, which(ref_group == ct)]
    }
  })
  pred_lab <- SingleR(test = query_data, ref = ref.data, labels = cts)$labels
  t2 = Sys.time()
  print(t2 - t1)
  true_lab <- query_group
  res <- evaluate(true_lab, pred_lab)
  return(res)
}

scpred_pred <- function(ref_data, ref_group, query_data, query_group) {
  t1 = Sys.time()
  true_lab <- query_group
  label <- ref_group
  ref_expr <- CreateSeuratObject(counts = ref_data, min.cells = 0, min.features = 0) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
    ScaleData(verbose = F) %>%
    RunPCA(verbose = F, npcs = 50) %>%
    RunUMAP(dims = 1:30)
  ref_expr$cell_type <- label
  ref_expr <- getFeatureSpace(ref_expr, "cell_type")
  ref_expr <- trainModel(ref_expr)
  query_expr <- CreateSeuratObject(counts = query_data, min.cells = 0, min.features = 0) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
    ScaleData(verbose = F) %>%
    RunPCA(verbose = F, npcs = 50) %>%
    RunUMAP(dims = 1:30)
  query_expr <- scPredict(query_expr, ref_expr)
  t2 = Sys.time()
  print(t2 - t1)
  pred_lab = query_expr$scpred_prediction
  res <- evaluate(true_lab, pred_lab)
  return(res)
}

# evaluate
evaluate <- function(true_lab, pred_lab) {
  unique_true <- unlist(unique(true_lab))
  unique_pred <- unlist(unique(pred_lab))
  
  unique_all <- unique(c(unique_true, unique_pred))
  conf <- table(true_lab, pred_lab)
  pop_size <- rowSums(conf)
  
  conf_F1 <- table(true_lab, pred_lab)
  F1 <- vector()
  sum_acc <- 0
  
  for (i in c(1:length(unique_true))) {
    findLabel <- colnames(conf_F1) == row.names(conf_F1)[i]
    if (sum(findLabel)) {
      prec <- conf_F1[i, findLabel] / colSums(conf_F1)[findLabel]
      rec <- conf_F1[i, findLabel] / rowSums(conf_F1)[i]
      if (prec == 0 || rec == 0) {
        F1[i] <- 0
      } else {
        F1[i] <- (2 * prec * rec) / (prec + rec)
      }
      sum_acc <- sum_acc + conf_F1[i, findLabel]
    } else {
      F1[i] <- 0
    }
  }
  
  pop_size <- pop_size[pop_size > 0]
  names(F1) <- names(pop_size)
  macro_F1 <- mean(F1)
  total <- length(pred_lab)
  num_unlab <- sum(pred_lab == "unassigned") + sum(pred_lab == "Unassigned") + sum(pred_lab == "rand") + sum(pred_lab == "Unknown") + sum(pred_lab == "unknown") + sum(pred_lab == "Node") + sum(pred_lab == "ambiguous") + sum(pred_lab == "unassign")
  per_unlab <- num_unlab / total
  acc <- sum_acc / sum(conf_F1)
  
  result <- list(Acc = acc, Macro_F1 = macro_F1, Conf = conf, F1 = F1, True_labels = true_lab, Pred_labels = pred_lab)
  return(result)
}