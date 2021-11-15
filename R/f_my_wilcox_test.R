f_my_wilcox_test <- function(x, y, min_genes = 5) {
  if (length(x) >= min_genes & length(y) >= min_genes) {
    return(wilcox.test(x, y, alternative = "greater")$p.value)
  } else {
    return(1)
  }
}
