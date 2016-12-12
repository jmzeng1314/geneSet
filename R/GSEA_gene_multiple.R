#' Do GSEA for a list of geneSets based on random permuted or bootstraps by genes
#'
#' It'll do GSEA Tests for a list of geneSets according to the signal 2 noise values for all of the genes or other values,such as foldchange.
#' It'll return ES ,the running ES score for all of the genes, a numeric vector .
#' rand_es ,the maximal running ES score for each permutation, a numeric vector.
#' p which calculated by formule:p <- 1-pnorm((max(ES)-mean(rand_es))/sd(rand_es))
#'
#'
#' @param gene_values  A vector of values for each gene,most of time it's  signal 2 noise.
#' @param gene_names   A vector of names  for each gene
#' @param geneSet_list  A list which contains all of the vector for each geneSet.
#' @param n Times for random permuted or bootstraps by genes,default:1000
#' @return a list contain 4 elements: ES,rand_es, p and ismem
#' @export
#' @keywords GSEA
#' @examples
#' GSEA_gene_multiple()

GSEA_gene_multiple <- function(gene_values=rnorm(100),
                             gene_names=1:100,
                             geneSet_list=list(set1=sample(1:100,10),set2=sample(1:100,10)),
                             n=1000) {
  lapply(geneSet_list, function(this_geneSet) GSEA_gene_single(gene_values,gene_names,this_geneSet,n))
}
