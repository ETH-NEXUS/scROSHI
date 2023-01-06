#' Wilcox test
#'
#' Calculate Wilcox p-value or return preset values
#'
#' @param x A numeric vector of values.
#' @param y A numeric vector of values.
#' @param min_genes A numeric value defining the threshold for the minimum number of genes.
#' @return A numeric value representing the p value of the Wilcox test.
#'
#' @export
#'
#' @examples
#' f_my_wilcox_test(rnorm(10, 1, 2), rnorm(10, 5, 2))

f_my_wilcox_test <- function(x, y, min_genes = 5) {
  if (length(x) >= min_genes & length(y) >= min_genes) {
    return(stats::wilcox.test(x, y, alternative = "greater")$p.value)
  } else {
    return(1)
  }
}
