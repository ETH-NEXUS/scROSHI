#' Wilcox test
#'
#' Calculate spearman correlation p-value or return preset values
#'
#' @param x value 1
#' @param y value 2
#' @param min_genes minimum genes
#'
#' @return
#' @export
#'
#' @examples
#' f_my_wilcox_test(rnorm(10,1,2),rnorm(10,5,2))
f_my_correlation_test <- function(x, y, min_genes = 5) {
  if (length(x) >= min_genes & length(y) >= min_genes) {
    return(stats::cor.test(x, y, alternative = "greater", method = "spearman", exact=F)$p.value)
  } else {
    return(1)
  }
}
