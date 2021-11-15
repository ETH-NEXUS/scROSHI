#' Score heatmap plot function
#'
#' @param final.score final score
#' @param major.type major type
#'
#' @return
#' @export
#'
#' @examples
f_plot_final_score <- function(final.score, major.type) {
  bad <- apply(final.score, 2, function(x) length(unique(x)))
  notbad <- which(bad != 1)
  if (length(notbad) > 1) {
    hm1 <- pheatmap(-log10(final.score[, notbad, drop = F]), scale = "none", color = col.pal,
                    silent = T, show_rownames = T, show_colnames = F,
                    clustering_method = "ward.D2", treeheight_row = 0, treeheight_col = 20)
    hm2 <- pheatmap(-log10(final.score[, notbad, drop = F]), scale = "column", color = col.pal,
                    silent = T, show_rownames = T, show_colnames = F,
                    clustering_method = "ward.D2", treeheight_row = 0, treeheight_col = 20)
    pp <- plot_grid(hm1$gtable, hm2$gtable, nrow = 2)
    ggsave(path %&% "." %&% major.type %&% "_final_celltype_scores.png", pp,
           width = 30, height = 13, units = "cm", dpi = 300)
  }
  write.table(t(final.score), path %&% "." %&% major.type %&% "_final_celltype_scores.tsv",
              sep = "\t", row.names = T, col.names = NA)
}
