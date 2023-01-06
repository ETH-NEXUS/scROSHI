#' Choose best matching cell type
#'
#' @param m.cts Matrix containing the cell type scores. The rows represent the cell types, whereas the columns represent the samples.
#' @param unknown If none of the probabilities is above this threshold,
#' the cell type label is assigned to the class unknown.
#' @param uncertain If the ratio between the largest and the second largest
#' probability is below this threshold, the cell type label is assigned to the
#' class uncertain for the major cell types.
#' @return Data frame containing the best matching cell type for each sample.
#'
#' @export
#'
#' @examples
#' m.cts <- matrix(c(0.2, 0.001, 0.002, 0.1), nrow=2)
#' colnames(m.cts) <- c("sample1", "sample2")
#' rownames(m.cts) <- c("cell_type1", "cell_type2")
#' f_annot_ctgenes(m.cts, 0.05, 0.1)

f_annot_ctgenes <- function(m.cts, unknown, uncertain) {
  ct.annot <- apply(m.cts, 2, function(x) names(which.min(x)))
  ct.annot <- data.frame(barcodes = names(ct.annot),
                         cell.type = as.character(ct.annot),
                         stringsAsFactors = F)
  best <- apply(m.cts, 2, min)
  ct.annot$celltype.score <- -log10(best)
  # deal with low certainty annotations:
  secondbest <- apply(m.cts, 2, function(x) sort(x, decreasing = F)[2])
  ct.annot$cts.2 <- -log10(secondbest)
  ct.annot$reclass <- NA
  # if best and second best are too similar & secondbest is at all relevant
  # cell gets reclassified as "uncertain"
  ct.annot$reclass[best / secondbest > uncertain & secondbest < unknown] <- "uncertain"
  # deal with low confidence annotations:
  # if total cell type score is too low cell gets reclassified as "unknown"
  ct.annot$reclass[ct.annot$celltype.score < -log10(unknown)] <- "unknown"
  # cell.type col is overwritten with `reclass` col where necessary
  ct.annot$cell.type[!is.na(ct.annot$reclass)] <- ct.annot$reclass[!is.na(ct.annot$reclass)]
  return(ct.annot[, 1:2])
}
