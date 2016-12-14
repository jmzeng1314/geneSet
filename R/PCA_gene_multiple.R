#' Do PCA for a single geneSet based on  bootstraps by genes
#'
#' It'll random choose the same number of gene with the size of this geneSet to extract the expression matrix and then do PCA for this small expression matrix .
#'
#'
#'
#' @param prefix    The prefix for all of the output files.( we don't need it now actually,just placeholder)
#' @param exprSet   Matrix for microarray expression values,rownames must be genes, colnames must be samples
#' @param geneSet_list  A list which contains all of the vector for each geneSet.
#' @param n Times for random permuted or bootstraps by genes,default:1000
#' @return A numeric vector of P value for each PCA test about all of the geneSet .
#' @export
#' @keywords PCA
#' @examples
#' PCA_gene_multiple(exprSet=exprSet, geneSet_list=geneSet_list)

PCA_gene_multiple <- function(prefix='test',exprSet, geneSet_list, n=1000) {
  if( F ){
    library(CLL)
    data(sCLLex)
    suppressMessages(library(limma))
    exprSet = exprs(sCLLex)
    pdata=pData(sCLLex)
    group_list = pdata$Disease
    geneSet_list <- list(set1 = sample(rownames(exprSet),50),
                         set2 = sample(rownames(exprSet),100),
                         set3 = sample(rownames(exprSet),150)
                           )
  }
  p <- unlist(lapply(geneSet_list, function(this_geneSet){
    PCA_gene_single(exprSet=exprSet, this_geneSet=this_geneSet)
  }))
  size <- unlist(lapply(geneSet_list, function(this_geneSet){
    length(this_geneSet)
  }))

  return( data.frame(geneSet_name = names(geneSet_list),
                     p = p,
                     size = size
                     )
          )
}





