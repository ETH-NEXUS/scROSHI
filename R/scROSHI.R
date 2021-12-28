scROSHI <- function(sce_data,
                    celltype_lists,
                    type_config,
                    min_genes=5,
                    min_var=1.5,
                    n_top_genes=2000,
                    n_nn=5,
                    thresh_unknown=0.05,
                    thresh_uncert=0.1,
                    thresh_uncert_second=0.8,
                    verbose=0){
  ## load celltype gene list -> all types in one file, subtype distinction via config file
  if (is.list(celltype_lists)){
    cell.type <- celltype_lists
  } else if (endsWith(celltype_lists, ".gmt")) {
    tmp <- readLines(celltype_lists)
    tmp <- lapply(tmp, function(x) strsplit(x, "\\\t")[[1]])
    names(tmp) <- sapply(tmp, function(x) x[1])
    cell.type <- sapply(tmp, function(x) x[-1])
  } else if (endsWith(celltype_lists, ".gmx")) {
    cell.type <- utils::read.table(celltype_lists, sep = "\t", head = T, stringsAsFactors = F)
    cell.type <- apply(cell.type[-1, ], 2, function(x) x[x != ""])
  } else {
    stop("Error. File type of cell type list must be either .gmt or .gmx.")
  }

  all.ct.genes <- as.character(unlist(cell.type))
  nr.genes <- sapply(cell.type, length)
  if(verbose == 1){
    cat("\n\nNumber of genes on each cell type list:\n\n")
    print(nr.genes)
    cat("\n\n")
  }

  # celltype config file, sub types are individual for each major type
  major_types <- type_config$Major
  if(verbose == 1){
    print("str(major_types):")
    print(utils::str(major_types))
  }
  minor_types <- lapply(type_config$Subtype, function(x) strsplit(x, ",")[[1]])
  names(minor_types) <- type_config$Major
  # remove "none"
  minor_types <- lapply(minor_types, function(x) x[x != "none"])
  # remove NA
  idx <- sapply(minor_types, function(x) length(which(is.na(x))))
  minor_types <- minor_types[which(idx == 0)]
  # remove empty list items
  idx <- sapply(minor_types, length)
  minor_types <- minor_types[idx > 0]
  # final list
  if(verbose == 1){
    print("str(minor_types):")
    print(utils::str(minor_types))
  }
  # Keep list of all possible final cell types
  all_celltypes_full_ct_name <- apply(type_config, 1, function(x) paste(x, collapse = ","))
  all_celltypes_full_ct_name <- paste(all_celltypes_full_ct_name, collapse = ",")
  all_celltypes_full_ct_name <- gsub("^NA,|,NA", "", all_celltypes_full_ct_name)
  all_celltypes_full_ct_name <- strsplit(all_celltypes_full_ct_name, ",")[[1]]
  idx <- which(all_celltypes_full_ct_name != "none")
  all_celltypes_full_ct_name <- unique(all_celltypes_full_ct_name[idx])

  if (utils::hasName(S4Vectors::metadata(sce_data), "all_celltypes_full_ct_name")) {
    S4Vectors::metadata(sce_data)$all_celltypes_full_ct_name <- all_celltypes_full_ct_name
  } else {
    S4Vectors::metadata(sce_data) <- c(S4Vectors::metadata(sce_data), list(all_celltypes_full_ct_name = all_celltypes_full_ct_name))
  }
  xx <- gsub(pattern = "([^_]+)_.*", "\\1", all_celltypes_full_ct_name)
  if (utils::hasName(S4Vectors::metadata(sce_data), "all_celltypes")) {
    S4Vectors::metadata(sce_data)$all_celltypes <- xx
  } else {
    S4Vectors::metadata(sce_data) <- c(S4Vectors::metadata(sce_data), list(all_celltypes = xx))
  }
  # Keep list of all possible major cell types
  if (utils::hasName(S4Vectors::metadata(sce_data), "all_major_types_full_ct_name")) {
    S4Vectors::metadata(sce_data)$all_major_types_full_ct_name <- major_types
  } else {
    S4Vectors::metadata(sce_data) <- c(S4Vectors::metadata(sce_data), list(all_major_types_full_ct_name = major_types))
  }
  xx <- gsub(pattern = "([^_]+)_.*", "\\1", major_types)
  if (utils::hasName(S4Vectors::metadata(sce_data), "all_major_types")) {
    S4Vectors::metadata(sce_data)$all_major_types <- xx
  } else {
    S4Vectors::metadata(sce_data) <- c(S4Vectors::metadata(sce_data), list(all_major_types = xx))
  }

  # filter out non-unique genes as long as more than min_genes are left
  genes_dupl <- which(duplicated(unlist(cell.type[major_types])))
  genes_dupl <- as.character(unlist(cell.type[major_types])[genes_dupl])
  # prelim new celltype-specific gene list
  test_major_types <- lapply(cell.type[major_types], function(x) setdiff(x, genes_dupl))
  # check if length > min_genes
  tmp <- which(sapply(test_major_types, length) >= min_genes)
  tmp <- setdiff(names(test_major_types), names(tmp))
  if (length(tmp) > 0) {
    # final new celltype specific gene list
    # if there is a cell type in temp that has less than min_genes genes
    # replace it with the cell type list BEFORE filtering for unique genes
    for (ii in seq(length(tmp))) {
      test_major_types[[tmp[ii]]] <- as.character(cell.type[[tmp[[ii]]]])
    }
  }
  # perform the first celltyping
  first.score <- f_score_ctgenes_U(sce_data, cell.type[major_types], min_genes)
  first.class <- f_annot_ctgenes(first.score, thresh_unknown, thresh_uncert)
  first.class$cell.type <- factor(first.class$cell.type, levels = c(major_types, "unknown", "uncertain"))
  # attach major celltype to SCE object
  SummarizedExperiment::colData(sce_data)$celltype_major_full_ct_name <- first.class$cell.type
  SummarizedExperiment::colData(sce_data)$celltype_major <- as.factor(gsub(pattern = "([^_]+)_.*", "\\1", first.class$cell.type))
  all.major.types <- sort(c(S4Vectors::metadata(sce_data)$all_major_types, "uncertain", "unknown"))
  if(verbose == 1){
    cat("\n\nall.major.types in alphabetical order:\n\n")
    print(all.major.types)
  }
  sce_data$celltype_major <- factor(sce_data$celltype_major, levels = all.major.types)
  #major score
  bad <- apply(first.score, 2, function(x) length(unique(x)))
  notbad <- which(bad != 1)
  #first_celltype_scores
  first_celltype_scores <- t(first.score[, notbad])
  # Rphenograph of celltype-specific genes.
  tmp <- SummarizedExperiment::assay(sce_data, "pearson_resid")
  idx <- match(all.ct.genes, rownames(tmp))
  idx <- idx[!is.na(idx)]
  tmp <- tmp[idx, ]
  r.umap <- uwot::umap(t(tmp), n_neighbors = n_nn, spread = 1, min_dist = 0.01,
                 #y = sce_data$celltype.major,
                 ret_nn = T)
  SingleCellExperiment::reducedDim(sce_data, "umap_ct") <- r.umap$embedding
  if (utils::hasName(S4Vectors::metadata(sce_data), "umap_ct")) {
    S4Vectors::metadata(sce_data)$umap_ct <- r.umap$nn$euclidean
  } else {
    S4Vectors::metadata(sce_data) <- c(S4Vectors::metadata(sce_data), list(umap_ct = r.umap$nn$euclidean))
  }

  #########################################################
  ## assign unknown or uncertain cells to the majority celltype of their nearest neighbors
  prelim.celltype <- sce_data$celltype_major_full_ct_name
  idx <- which(prelim.celltype %in% c("unknown", "uncertain"))
  uu.nn <- r.umap$nn$euclidean$idx[idx, , drop = F]
  # get preliminary cell types of the nearest neighbors of each uncertain/unknown cell
  # disregard the cell itself
  uu.nn.ct <- prelim.celltype[as.numeric(uu.nn[, -1])]
  # have preliminary cell types shaped as matrix
  # for each uncertain/unknown cell there is a row with (n_nn - 1) columns.
  # the columns contain the cell types of the nearest neighbours of the cell.
  uu.nn.ct <- matrix(uu.nn.ct, ncol = ncol(uu.nn) - 1, byrow = F)
  # deal with definitely unknown/uncertain cases:
  id_clear <- apply(uu.nn.ct, 1, function(x) all(x %in% c("unknown", "uncertain")))
  id_unclear <- which(!id_clear)
  # deal with rest:
  f.vote <- function(x) names(sort(table(x[!x %in% c("unknown", "uncertain")]), decreasing = T)[1])
  # f.vote(uu.nn.ct[1,])
  vote.ct <- apply(uu.nn.ct[id_unclear, , drop = F], 1, f.vote)
  # replace "unknown" and "uncertain" with the majority celltype of their nearest neighbors
  prelim.celltype[idx[id_unclear]] <- vote.ct
  if(verbose == 1){
    cat("\nResults of first cell type classification:\n\n")
    as.matrix(table(sce_data$celltype_major_full_ct_name))
  }

  # perform second celltyping
  # prefill final ct column with results of major cell type classification
  SummarizedExperiment::colData(sce_data)$celltype_final_full_ct_name <- as.character(sce_data$celltype_major_full_ct_name)
  # iterate over all major cell types (ii is index of major cell type in list `minor_types`)
  for (ii in seq_len(length(minor_types))) {
    # for each major cell type with minor defined, get all cells with this major cell type as prelim.celltype
    idy <- which(prelim.celltype %in% names(minor_types)[ii])
    if(verbose == 1){
      cat("\n\nPerform second cell typing for", length(idy), "cells of major type:   ", names(minor_types)[ii], "\n")
    }
    if (length(idy) > 0) {
      these.cells <- sce_data[, idy]
      these_major <- these.cells$celltype_major_full_ct_name
      these.ct <- match(minor_types[[ii]], names(cell.type))
      these.ct <- these.ct[!is.na(these.ct)]
      second.score <- f_score_ctgenes_U(these.cells, cell.type[these.ct], min_genes)
      second.class <- f_annot_ctgenes(second.score, thresh_unknown, thresh_uncert_second)
      if(verbose == 1){
        table(second.class$cell.type)
      }
      # cells with (second) label "unknown", "uncertain" keep the original first.class as final cell type
      idx <- which(second.class$cell.type %in% c("unknown", "uncertain"))
      if (length(idx) > 0) second.class$cell.type[idx] <- as.character(these_major[idx])
      sce_data$celltype_final_full_ct_name[idy] <- second.class$cell.type
      # for all successful second classifications, replace major celltype with parent of minor type
      idx <- which(!second.class$cell.type %in% c("unknown", "uncertain"))
      if (length(idx) > 0) sce_data$celltype_major_full_ct_name[idy][idx] <- names(minor_types)[ii]
    }
  }

  SummarizedExperiment::colData(sce_data)$celltype_final <- gsub(pattern = "([^_]+)_.*", "\\1", SummarizedExperiment::colData(sce_data)$celltype_final_full_ct_name)
  all.cell.types <- sort(c(S4Vectors::metadata(sce_data)$all_celltypes, "uncertain", "unknown"))
  if(verbose ==1){
    print("all.cell.types in alphabetical order:")
    print(all.cell.types)
    sce_data$celltype_final <- factor(sce_data$celltype_final, levels = all.cell.types)
    cat("\n\nsce_data$celltype_final after second cell typing:\n\n")
    as.matrix(table(sce_data$celltype_final))
    cat("\nAdapted first cell type classification:\n\n")
    as.matrix(table(sce_data$celltype_major_full_ct_name))
  }
  # have colData(sce_data)$celltype_final_full_ct_name as factors
  all.cell.types.full <- sort(c(S4Vectors::metadata(sce_data)$all_celltypes_full_ct_name, "uncertain", "unknown"))
  # print("all.cell.types.full in alphabetical order:")
  # print(all.cell.types.full)
  SummarizedExperiment::colData(sce_data)$celltype_final_full_ct_name <- factor(SummarizedExperiment::colData(sce_data)$celltype_final_full_ct_name, levels = all.cell.types.full)


  # keep sce object including cluster 0
  sce_data_incl_cluster0 <- sce_data

  # remove cells that are in cluster 0
  mask_keep_cells <- SummarizedExperiment::colData(sce_data)$phenograph_clusters != 0
  if(verbose == 1){
    cat("\n\nNumber of cells that remain after filtering out cluster 0:\n\n")
    print(sum(mask_keep_cells))
    cat("\n\nsce object before filtering out cells that are in cluster 0:\n\n")
    print(sce_data)
  }
  sce_data <- sce_data[, mask_keep_cells]
  SingleCellExperiment::reducedDims(sce_data)$phenodist <- SingleCellExperiment::reducedDims(sce_data)$phenodist[, mask_keep_cells]
  if(verbose == 1){
    cat("\n\nsce object after filtering out cells that are in cluster 0:\n\n")
    print(sce_data)
  }
  ## write final celltype to disk
  res <- SummarizedExperiment::colData(sce_data)[, c("barcodes", "celltype_final")]
  return(res)
}
