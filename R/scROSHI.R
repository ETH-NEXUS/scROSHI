#' Robust Supervised Hierarchical Identification of Single Cells
#'
#'@description
#' scROSHI identifies cell types based on expression profiles of single cell analysis by
#' utilizing previously obtained cell type specific gene sets. It takes into account the
#' hierarchical nature of cell type relationship and does not require training or
#' annotated data.
#'
#' @param sce_data A SingleCellExperiment object containing the expression
#' profiles of the single cell analysis.
#' @param celltype_lists Marker gene list for all cell types. It can be provided
#' as a list of genes with cell types as names or as a path to a file containing the
#' marker genes. Supported file formats are .gmt or .gmx files.
#' @param type_config Config file to define major cell types and hierarchical subtypes.
#' It should be provided as a two-column data.frame where the first column are the
#' major cell types and the second column are the subtypes. If several subtypes exists,
#' they should be separated by comma.
#' @param gene_symbol Variable name in the row data of the SingleCellExperiment object containing the gene names.
#' @param count_data Assay name in the SingleCellExperiment object containing the count data.
#' @param cell_scores Boolean value determining if the scores should be saved.
#' @param min_genes scROSHI filters out non-unique genes as long as more than min_genes
#' are left. If there is a cell type that has less than min_genes genes, it will be
#' replaced with the cell type list BEFORE filtering for unique genes (default 5).
#' @param min_var Minimum variance for highly variable genes (default 1.5).
#' @param n_top_genes Maximum number of highly variable genes (default 2000).
#' @param n_nn Number of nearest neighbors for umap for assignment of cell types
#' (default 5).
#' @param thresh_unknown If none of the probabilities is above this threshold,
#' the cell type label is assigned to the class unknown (default 0.05).
#' @param thresh_uncert If the ratio between the largest and the second largest
#' probability is below this threshold, the cell type label is assigned to the
#' class uncertain for the major cell types (default 0.1).
#' @param thresh_uncert_second If the ratio between the largest and the second largest
#' probability is below this threshold, the cell type label is assigned to the
#' class uncertain for the subtypes (default 0.8).
#' @param verbose Level of verbosity. Zero means silent, one makes a verbose output.
#' @param output Defines the output. sce: The output is a SingleCellExperiment
#' object with the cell types appended to the meta data. df: The output id a
#' data.frame with two columns. The first column contains the barcode of the cell
#' and the second column contains the cell type labels.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data("test_sce_data")
#' data("config")
#' data("marker_list")
#'
#' results <- scROSHI(sce_data = test_sce_data,
#'                    celltype_lists = marker_list,
#'                    type_config = config)
#' table(results$celltype_final)
#' }

scROSHI <- function(sce_data,
                    celltype_lists,
                    type_config,
                    count_data = "normcounts",
                    gene_symbol = "SYMBOL",
                    cell_scores = FALSE,
                    min_genes=5,
                    min_var=1.5,
                    n_top_genes=2000,
                    n_nn=5,
                    thresh_unknown=0.05,
                    thresh_uncert=0.1,
                    thresh_uncert_second=0.8,
                    verbose=0,
                    output="sce"){
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
  first.score <- f_score_ctgenes_U(sce = sce_data,
                                   gset = cell.type[major_types],
                                   count_data = count_data,
                                   gene_symbol = gene_symbol,
                                   min_genes = min_genes,
                                   verbose = verbose)
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
  # attach major celltype scores to SCE object
  if (cell_scores==TRUE){
    first_ct_scores <- first_celltype_scores
    colnames(first_ct_scores) <- paste("score_major",colnames(first_ct_scores),sep="_")
    SummarizedExperiment::colData(sce_data) <- cbind(SummarizedExperiment::colData(sce_data),first_ct_scores)
  }
  # Rphenograph of celltype-specific genes.
  tmp <- SummarizedExperiment::assay(sce_data, count_data)
  idx <- match(all.ct.genes, rownames(tmp))
  idx <- idx[!is.na(idx)]
  tmp <- tmp[idx, ]

  # Fix sparse matrix
  if (is(tmp, "sparseMatrix")) {
    tmp <- Matrix::as.matrix(tmp)
  }

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
      second.score <- f_score_ctgenes_U(sce = these.cells,
                                        gset = cell.type[these.ct],
                                        count_data = count_data,
                                        gene_symbol = gene_symbol,
                                        min_genes = min_genes,
                                        verbose = verbose)
      # attach second celltype scores to SCE object
      if (cell_scores==TRUE){
        second_ct_scores <- t(second.score)
        second_ct_scores_2 <- matrix(NA,ncol=ncol(second_ct_scores),nrow = nrow(SummarizedExperiment::colData(sce_data)))
        second_ct_scores_2[idy,] <- second_ct_scores
        colnames(second_ct_scores_2) <- paste("score_subtype",colnames(second_ct_scores),sep="_")
        SummarizedExperiment::colData(sce_data) <- cbind(SummarizedExperiment::colData(sce_data),second_ct_scores_2)
      }
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

  ## write final celltype to disk
  if(output == "sce"){
    res <- sce_data
  }
  if(output == "df"&cell_scores==FALSE){
    res <- SummarizedExperiment::colData(sce_data)[, c("barcodes", "celltype_major", "celltype_final")]
    res <- as.data.frame(res)
    rownames(res) <- NULL
  }
  if(output == "df"&cell_scores==TRUE){
    s_major <- colnames(SummarizedExperiment::colData(sce_data))[grep("score_major",colnames(SummarizedExperiment::colData(sce_data)))]
    s_sub <- colnames(SummarizedExperiment::colData(sce_data))[grep("score_sub",colnames(SummarizedExperiment::colData(sce_data)))]
    res <- SummarizedExperiment::colData(sce_data)[, c("barcodes", "celltype_major", s_major,  "celltype_final", s_sub)]
    res <- as.data.frame(res)
    rownames(res) <- NULL
  }
  return(res)
}
