#' Hypergeometric Tests for a single geneSet
#'
#' It'll  calculate signal to noise value for each gene in the expression matrix according to the proper grouplist.
#'
#'
#' @param exprSet   Matrix for microarray expression values,rownames must be genes, colnames must be samples
#' @param group_list Factors for two groups that are tested for differential expression.(first group is control)
#' @return a vector of signal to noise value for each gene .
#' @export
#' @keywords signal2noise
#' @examples
#' calc_s2n(exprSet,group_list )

calc_s2n <- function( exprSet,group_list ){
  dat=exprSet
  group1= which(group_list==levels(group_list)[1])
  group2= which( group_list==levels(group_list)[2])
  dat1=dat[,group1];dat2=dat[,group2]
  dat=cbind(dat1,dat2)
  if (ncol(dat1)<3 | ncol(dat2)<3  ){
    stop('both of the group shoud have more than 2 samples for the calculation of signal 2 noise value !!!')
  }
  s2n <- unlist(lapply(1:nrow(dat), function(i){
    (mean(dat2[i,])-mean(dat1[i,]))/(sd(dat1[2,])+sd(dat1[i,]))
  }))
  s2n <- scale(s2n)[,1]
  names(s2n) <- rownames(exprSet)

  s2n <- sort(s2n,decreasing = T)

  return( s2n )
}

