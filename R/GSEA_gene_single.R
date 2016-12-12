#' Do GSEA for a single geneSet based on random permuted or bootstraps by genes
#'
#' It'll do GSEA Tests for a single geneSet according to the signal 2 noise values for all of the genes or other values,such as foldchange.
#' It'll return ES ,the running ES score for all of the genes, a numeric vector .
#' rand_es ,the maximal running ES score for each permutation, a numeric vector.
#' p which calculated by formule:p <- 1-pnorm((max(ES)-mean(rand_es))/sd(rand_es))
#' ismem is a logical vector shows whether the gene belong to the current geneSet or not.
#'
#'
#'
#' @param gene_values  A vector of values for each gene,most of time it's  signal 2 noise.
#' @param gene_names   A vector of names  for each gene
#' @param this_geneSet   A vector which contains all of the genes in this geneSet.
#' @param n Times for random permuted or bootstraps by genes,default:1000
#' @return a list contain 4 elements: ES,rand_es, p and ismem
#' @export
#' @keywords GSEA
#' @examples
#' GSEA_gene_single()

GSEA_gene_single <- function(gene_values=rnorm(100),
                             gene_names=1:100,
                             this_geneSet=sample(1:100,10),
                             n=1000) {
  if (any(duplicated( gene_names  ))){
    stop("The gene names for genes should not be duplicated!!!")
  }
  names(gene_values) <- gene_names
  gene_values <- sort(gene_values,decreasing = TRUE)

  n1 <- length(gene_values)
  n2 <- length(this_geneSet)

  ismem<-is.element(names(gene_values),this_geneSet)
  ES<-c(0)
  for (i in 1:n1) {
    if (ismem[i]) x <- sqrt((n1-n2)/n2) else x<- -sqrt(n2/(n1-n2)) # 70 For brca
    ES[i+1] <- ES[i]+x
  }
  ES<-ES[-1]  ## the running ES score for all of the genes, a numeric vector

  ## plot(ES,type='l',xlim=c(0,n1))
  ## random permutation
  rand_es<-c()
  max.rES <- c()
  for (j in 1:n) {
    #myData<-cbind(myData[,1:2],is.element(c(1:n1),sample(1:n1,n2)))
    ismem <- sample(ismem)
    rES<-c(0)
    for (i in 1:n1) {
      if (ismem[i]) x <- sqrt((n1-n2)/n2) else x<- -sqrt(n2/(n1-n2))
      rES[i+1] <- rES[i]+x
    }
    rand_es <- rbind(rand_es,max(rES)) ## the maximal running ES score for each permutation, a numeric vector.
    #max.rES <- c(max.rES,max(rES))
  }
  #b<-density(rand_es, adjust=1.5)
  #p = sum(max.rES>max(ES))/n
  p <- 1-pnorm((max(ES)-mean(rand_es))/sd(rand_es))
  ismem<-is.element(names(gene_values),this_geneSet)
  return(list(ES=ES,rand_es=rand_es,p=p,ismem=ismem))
}





