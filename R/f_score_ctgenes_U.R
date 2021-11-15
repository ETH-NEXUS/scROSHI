#' Calculate celltype score: U-test
#'
#' @param sce sce object
#' @param gset gene set
#' @param min_genes minimum number of genes
#'
#' @return
#' @export
#'
#' @examples
f_score_ctgenes_U <- function(sce, gset, min_genes = 5) {
  cat("\nCalculate scores for cell type classification.")
  # gset = cell.type[these.ct]
  # sce = these.cells
  # min_genes = 5
  # select count table with all.ct.genes
  all.genes <- unique(as.character(unlist(gset)))
  cat("\n\nNumber of genes on cell type specific lists:", length(all.genes))
  # count matrix including all cells and only the cell type specific genes
  mat <- SummarizedExperiment::assay(sce, "normcounts")[SingleCellExperiment::rowData(sce)$SYMBOL %in% all.genes, , drop = F]
  cat("\n\nNumber of genes included in matrix:", dim(mat)[1])
  # keep genes only if they have counts in at least one cell
  keep <- rowSums(mat)
  mat <- mat[keep > 0, , drop = F]
  cat("\n\nNumber of genes with with a sum of counts > 0:", dim(mat)[1], "\n\n\n")
  colnames(mat) <- sce$barcodes
  # get row index of dd for each gene in a list of celltypes
  idxs <- limma::ids2indices(gene.sets = gset, rownames(mat), remove.empty = F)
  # generate gene x celltype matrix
  # all values are 0
  ds <- matrix(0, nrow = nrow(mat), ncol = length(idxs))
  rownames(ds) <- rownames(mat)
  colnames(ds) <- names(idxs)
  # fill matrix with 1 where a gene is specific for a cell type
  for (cell_type in seq(length(idxs))) {
    ds[idxs[[cell_type]], cell_type] <- 1
  }
  # perform wilcox test using gene x sample and gene x celltype matrices.
  # wilcox.test(mat[ds[,1]==1,1],mat[ds[,1]==0,1])$p.value
  m.cts <- apply(mat, 2, function(x) {
    # on each column of mat (each cell)
    # perform wilcox test with each list of cell type genes against all other cell type genes
    apply(ds, 2, function(y) f_my_wilcox_test(x[y == 1], x[y == 0], min_genes))
  })
  m.cts[is.na(m.cts)] <- 1
  return(m.cts)
}
