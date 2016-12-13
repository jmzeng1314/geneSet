#' draw 3 figures for the results of GSEA_gene_single function
#'
#' Firstly you need to do GSEA Tests for a list of geneSets,and get the results, then just use this function to draw some pictures.
#'
#'
#' @param GSEA_single_results  A list contain 4 elements: ES,rand_es, p and ismem
#' @param  gene_values  A vector of values for each gene,most of time it's  signal 2 noise.
#' @param this_geneSet   The name of current geneSet,default:geneset1
#' @param imageType choose png,pdf,emf,tiff,jpeg,bmp
#' @param multiple  plot 3 ugly figures or just 1(pretty),default:T
#' @return image file:
#' @export
#' @keywords GSEA,plot
#' @examples
#' plot_GSEA_single();lapply(GSEA_gene_multiple(),plot_GSEA_single)

plot_GSEA_single <- function(GSEA_single_results,
                             gene_values=rnorm(100),
                             geneSet_name='geneset1',
                             phen1 = 'case',phen2 = 'control',
                             imageType='png',
                             multiple=T){
  library(devEMF)
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
  if(multiple){
    # figure 1:
    tmp_create_figures = get(imageType) ## A function
    if (  imageType %in% c('png','tiff','jpeg','bmp') ) {
      tmp_create_figures( paste0(geneSet_name,'_Enrichment_score','.',imageType))
    } else if (  imageType %in% c('pdf','emf')  ){
      library(devEMF)
      tmp_create_figures( paste0(geneSet_name,'_Enrichment_score','.',imageType) )
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
      tmp_create_figures( paste0(geneSet_name,'_values_of_samples','.',imageType))
    } else if (  imageType %in% c('pdf','emf')  ){
      library(devEMF)
      tmp_create_figures( paste0(geneSet_name,'_values_of_samples','.',imageType) )
    }else{
      stop(' we just accept png,pdf,emf,tiff,jpeg,bmp !!!')
    }
    #pdf(file = "values of samples.pdf",width = 8,height = 6)
    barplot(as.numeric(gene_values),border = NA,space=0.5,col=cols,main=paste0("The size of this geneSet is :",length(cols[cols=='red'])))
    dev.off()

    # figures 3:
    tmp_create_figures = get(imageType) ## A function
    if (  imageType %in% c('png','tiff','jpeg','bmp') ) {
      tmp_create_figures( paste0(geneSet_name,'_Distribution_of_max_score_in_permutations','.',imageType))
    } else if (  imageType %in% c('pdf','emf')  ){
      library(devEMF)
      tmp_create_figures( paste0(geneSet_name,'_Distribution_of_max_score_in_permutations','.',imageType) )
    }else{
      stop(' we just accept png,pdf,emf,tiff,jpeg,bmp !!!')
    }
    #pdf(file = "Distribution of max score in permutations.pdf",width = 8,height = 6)
    b<-density(rand_es, adjust=1.5)
    plot(b,xlim=c(min(rand_es),max(c(rand_es,ES))),main=paste0("The p value is : ",p))
    par(new=TRUE) ## red line will indicate the max value of Enrichment score as a comparision with random distribution
    plot(c(max(ES),max(ES)),c(0,0.1),"l",col="red",axes=FALSE,ann=FALSE,xlim=c(min(rand_es),max(c(rand_es,ES))))
    dev.off()
  }else{
    cex <- 0.8
    tmp_create_figures = get(imageType) ## A function
    if (  imageType %in% c('png','tiff','jpeg','bmp') ) {
      tmp_create_figures( paste0(geneSet_name,'_GSEA','.',imageType),width = 680,height = 680 ,res=160)
    } else if (  imageType %in% c('pdf','emf')  ){
      library(devEMF)
      tmp_create_figures( paste0(geneSet_name,'_GSEA','.',imageType) )
    }else{
      stop(' we just accept png,pdf,emf,tiff,jpeg,bmp !!!')
    }

    par( mgp=c(3, 1, 0))
    layout(matrix(c(1,2,3), 3, 1, byrow = TRUE),
           heights=c(4,1,3))
    par(mar=c(0,5,2,2)+0.1)
    max_ES_pos <- which(abs(ES) == max(abs(ES)))
    col <- ifelse(ES[max_ES_pos] > 0, 2, 4)

    ## first figure:
    plot(ES,main=geneSet_name,ylab = "Running Enrichment Score (RES)",xaxt='n',
         xlim=c(1, n1), type = "l", lwd = 2, col = col,cex=cex)
    abline(h=0,lwd = 1, lty = 2, col = 1)# zero RES line
    abline(v=max_ES_pos,lwd = 1, lty = 3, col = col) # max enrichment vertical line

    ## second figure:
    par(mar=c(0,5,0,2)+0.1)
    plot(gene_values,type = 'n',xaxt='n', ann=FALSE,yaxt='n' )
    lapply(1:length(ismem), function(i) if(ismem[i]){abline(v=i)})# enrichment tags

    ## last figure:
    par(mar=c(3,5,0,2)+0.1)
    plot(gene_values,type = 'h',ylab = "Ranked list mertix \n (signal2noise)" ,cex=cex)
    temp <- order(abs(gene_values), decreasing=T)
    arg.correl <- temp[n1]
    abline(v=arg.correl, lwd = 1, lty = 3,col = 3) # zero crossing correlation vertical line

    leg.txt <- paste("\"", phen1, "\" ","(positively correlated)", sep="", collapse="")
    text(x=1, y=max(gene_values), adj = c(0, 1), labels=leg.txt,col='red')

    leg.txt <- paste("\"", phen2, "\" ","(negatively correlated)", sep="", collapse="")
    text(x=n1, y=min(gene_values), adj = c(1, 0), labels=leg.txt, col = 'blue')

    leg.txt <- paste("Zero crossing at ", arg.correl, sep="", collapse="")
    text(x=arg.correl, y=0, adj = c(0, -0.5), labels=leg.txt  )
    dev.off()

    #par(old_par)
  }

  ## don't run the code below:
  if (F){
    for (i in 1:length(gs.names)){

      GSEA_single_results <- list(
        ES = Obs.RES[i,],  ## the running ES score for each gene, (genes have been sorted.)
        rand_es= Obs.RES[i,], ## we don't need it to draw figures.
        p = 0.1, ## we also don't need it .
        ismem =Obs.indicator[i,] ## a logical vector

      )
      plot_GSEA_single(GSEA_single_results,
                       gene_values=obs.s2n,
                       geneSet_name=gs.names[i],
                       phen1 = 'CBX6.5',phen2 = 'control',
                       imageType='emf',
                       multiple=F)
    }
  }


}


