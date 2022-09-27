#' Calculate celltype score: U-test
#'
#' @param sce sce object
#' @param lprof list of profiles (= named expression vectors)
#' @param min_genes minimum number of genes
#'
#' @export
#'
#' @examples
#' #' \dontrun{
#' result <- f_score_profile_cor(sce, lprof, min_genes = 5,verbose = 0)
#' }

f_score_profile_cor = function(sce, lprof, min_genes = 5, verbose = 0) {
  if(verbose == 1){
    cat("\nCalculate scores for cell type classification.")
  }
  # lprof = l_prof
  # sce = these.cells
  # min_genes = 5
  
  all_genes = unique(as.character(sapply(lprof, names)))
  if(verbose == 1){
    cat("\n\nNumber of genes on cell type specific lists: ", length(all_genes))
  }
  
  # count matrix including all cells and only the cell type specific genes
  m_cnt = SummarizedExperiment::assay(sce, "normcounts")[rownames(sce) %in% all_genes, , drop = F]
  if(verbose == 1){
    cat("\n\nNumber of genes included in matrix: ", dim(m_cnt)[1])
  }
  # keep genes only if they have counts in at least one cell
  keep = rowSums(m_cnt)
  m_cnt = m_cnt[keep > 0, , drop = F]
  if(verbose == 1){
    cat("\n\nNumber of genes with a sum of counts over all cells > 0: ", dim(m_cnt)[1], "\n\n\n")
  }
  

  m_cts = matrix(NA, nrow = length(lprof), ncol = ncol(sce))
  dimnames(m_cts) = list(names(lprof), colnames(m_cnt))
  for(jj in seq_along(lprof)){
    # jj = 1
    this_profile = lprof[[jj]]
    this_genes = names(this_profile)
    for(kk in seq_len(ncol(m_cts))){
      # kk = 1
      this_cell = m_cnt[, kk]

      # perform correlation test using gene x cell matrix <m_cnt> and 
      # celltype list of named (gene) expression vectors lprof.
      id_match = match(names(this_cell), names(this_prof))
      cor.test(this_cell, this_profile[id_match], alternative="greater", method="spearman", exact=F)$p.value
      m_cts[jj, kk] = f_my_correlation_test(this_cell, this_profile[id_match], min_genes = 5)
    }
  
  }
  m_cts[is.na(m_cts)] = 1
  return(m_cts)
}
