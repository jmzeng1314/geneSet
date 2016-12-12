#' Do PCA for a single geneSet based on  bootstraps by genes
#'
#' It'll random choose the same number of gene with the size of this geneSet to extract the expression matrix and then do PCA for this small expression matrix .
#'
#'
#'
#' @param prefix    The prefix for all of the output files.( we don't need it now actually,just placeholder)
#' @param exprSet   Matrix for microarray expression values,rownames must be genes, colnames must be samples
#' @param group_list Factors for two groups that are tested for differential expression.
#' @param this_geneSet   A vector which contains all of the genes in this geneSet.
#' @param n Times for random permuted or bootstraps by genes,default:1000
#' @return The p vaule for this PCA test
#' @export
#' @keywords PCA
#' @examples
#' PCA_gene_single(exprSet=exprSet,group_list=group_list, this_geneSet=this_geneSet)

PCA_gene_single <- function(prefix='test',exprSet,group_list, this_geneSet, n=1000) {
 if( F ){
   library(CLL)
   data(sCLLex)
   suppressMessages(library(limma))
   exprSet = exprs(sCLLex)
   pdata=pData(sCLLex)
   group_list = pdata$Disease
   this_geneSet <- sample(rownames(exprSet),200)
 }
  size <- length(this_geneSet)
  library(gmodels)
  backgroud_PCA=sapply(1:n,function(y) {
    dat=t(exprSet[sample(row.names(exprSet), size, replace=TRUE), ]);
    round(100*summary(fast.prcomp(dat))$importance[2,1],2)
  })
  this_geneSet_gene=intersect(rownames(exprSet),as.character( this_geneSet ))
  dat=t(exprSet[this_geneSet_gene,]);
  this_PCA <- round(100*summary(fast.prcomp(dat))$importance[2,1],2)

  p <- 1-pnorm((this_PCA-mean(backgroud_PCA))/sd(backgroud_PCA))
  return( p )
}





