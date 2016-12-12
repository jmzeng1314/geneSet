#' Hypergeometric Tests for a single geneSet
#'
#' It'll do Hypergeometric Tests for a single geneSet, such as DNA replication. And we need 3 parameters, first is this_geneSet which contains all of the genes in this geneSet.
#' Then is choose_geneIds,which probably is the significantly differently expression genes. Last is the backgroud_geneIds, which should be all of the genes which have been detected by microarray or RNA-seq, and so on .
#' The N is the size of backgroud_geneIds
#' The n is the size of choose_geneIds
#' The M is the size of intersect(this_geneSet,backgroud_geneIds)
#' The K is the size of intersect(this_geneSet,choose_geneIds)
#' The exp_count is n*M/N
#' The OddsRatio is k/exp_count
#' The p is the statistic P value for this Hypergeometric Tests
#'
#'
#' @param this_geneSet  A vector which contains all of the genes in this geneSet.
#' @param choose_geneIds A vector which contains all of gene which defined by user, probably significantly differently expression genes.
#' @param backgroud_geneIds A vector which contains all of the genes have been detected by microarray or RNA-seq, and so on
#' @return a vector include(N,n,M,k,exp_count,OddsRatio,p)
#' @export
#' @keywords hyperGtest
#' @examples
#' hyperGTest_single()

hyperGTest_single <- function(this_geneSet=sample(1:100,10),
                              choose_geneIds= sample(1:100,20),
                              backgroud_geneIds=1:100) {
  if(!all(choose_geneIds %in% backgroud_geneIds)){
    stop('There are some genes in choose_geneIds are not belong to  backgroud_geneIds!!!')
  }
  N <- length(backgroud_geneIds)
  n <- length(choose_geneIds)
  M <- length(intersect(this_geneSet,backgroud_geneIds))
  exp_count=n*M/N
  k <- length(intersect(this_geneSet,choose_geneIds))
  OddsRatio=k/exp_count
  p=phyper(k-1,M, N-M, n, lower.tail=F)
  return(c(N,n,M,k,exp_count,OddsRatio,p))
}

