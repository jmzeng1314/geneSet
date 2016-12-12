#' draw 3 figures for the results of GSEA_gene_single function
#'
#' Firstly you need to do GSEA Tests for a list of geneSets,and get the results, then just use this function to draw some pictures.
#'
#'
#' @param GSEA_single_results  A list contain 4 elements: ES,rand_es, p and ismem
#' @param  gene_values  A vector of values for each gene,most of time it's  signal 2 noise.
#' @param prefix   The prefix for the out images,default:test
#' @param imageType choose png,pdf,emf,tiff,jpeg,bmp
#' @return 3 figures:
#' @export
#' @keywords GSEA,plot
#' @examples
#' plot_GSEA_single();lapply(GSEA_gene_multiple(),plot_GSEA_single)

plot_GSEA_single <- function(GSEA_single_results,gene_values=rnorm(100),prefix ='test',imageType='png'){
  if(missing(GSEA_single_results)){
    GSEA_single_results<- GSEA_gene_single()
  }
  gene_values <- sort(gene_values,decreasing = TRUE)
  ES <- GSEA_single_results$ES
  rand_es <- GSEA_single_results$rand_es
  p <- GSEA_single_results$p
  ismem <- GSEA_single_results$ismem
  n1 <- length(ES)
  #print(str(GSEA_single_results))
  # figure 1:
  tmp_create_figures = get(imageType) ## A function
  if (  imageType %in% c('png','tiff','jpeg','bmp') ) {
    tmp_create_figures( paste0(prefix,'_Enrichment_score','.',imageType))
  } else if (  imageType %in% c('pdf','emf')  ){
    library(devEMF)
    tmp_create_figures( paste0(prefix,'_Enrichment_score','.',imageType) )
  }else{
    stop(' we just accept png,pdf,emf,tiff,jpeg,bmp !!!')
  }
  plot(ES,type='l',xlim=c(0,n1))
  abline(h=0,col='red')
  dev.off()


  # figure 2:
  cols <- rep("grey",n1)
  cols[ismem] <- "red" # samples of interest are colored red
  tmp_create_figures = get(imageType) ## A function
  if (  imageType %in% c('png','tiff','jpeg','bmp') ) {
    tmp_create_figures( paste0(prefix,'_values_of_samples','.',imageType))
  } else if (  imageType %in% c('pdf','emf')  ){
    library(devEMF)
    tmp_create_figures( paste0(prefix,'_values_of_samples','.',imageType) )
  }else{
    stop(' we just accept png,pdf,emf,tiff,jpeg,bmp !!!')
  }
  #pdf(file = "values of samples.pdf",width = 8,height = 6)
  barplot(as.numeric(gene_values),border = NA,space=0.5,col=cols,main=paste0("The size of this geneSet is :",length(cols[cols=='red'])))
  dev.off()

  # figures 3:
  tmp_create_figures = get(imageType) ## A function
  if (  imageType %in% c('png','tiff','jpeg','bmp') ) {
    tmp_create_figures( paste0(prefix,'_Distribution_of_max_score_in_permutations','.',imageType))
  } else if (  imageType %in% c('pdf','emf')  ){
    library(devEMF)
    tmp_create_figures( paste0(prefix,'_Distribution_of_max_score_in_permutations','.',imageType) )
  }else{
    stop(' we just accept png,pdf,emf,tiff,jpeg,bmp !!!')
  }
  #pdf(file = "Distribution of max score in permutations.pdf",width = 8,height = 6)
  b<-density(rand_es, adjust=1.5)
  plot(b,xlim=c(min(rand_es),max(c(rand_es,ES))),main=paste0("The p value is : ",p))
  par(new=TRUE) ## red line will indicate the max value of Enrichment score as a comparision with random distribution
  plot(c(max(ES),max(ES)),c(0,0.1),"l",col="red",axes=FALSE,ann=FALSE,xlim=c(min(rand_es),max(c(rand_es,ES))))
  dev.off()

}

