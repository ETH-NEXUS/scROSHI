#' Choose best matching cell type
#'
#' @param m.cts Cell type score
#' @param unknown Threshold unknown
#' @param uncertain Threshold uncertain
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- f_annot_ctgenes(m.cts, unknown, uncertain)
#' }
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
