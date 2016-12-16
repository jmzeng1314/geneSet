#' create GMT file
#'
#' create gmt format file for GSEA software by broad institute.
#'
#' @param kegg2symbol_list A list contains all of the geneSet, each geneSet has some genes.(HUGO symbol)
#' @param gmt_file The filename for gmt file.
#' @return write gmt file
#' @export
#' @keywords gmt
#' @examples
#' #' write_gmt()


write_gmt <- function(geneSet=kegg2symbol_list,gmt_file='kegg2symbol.gmt'){

  sink( gmt_file )
  for (i in 1:length(geneSet)){
    cat(names(geneSet)[i])
    cat('\tNA\t')
    cat(paste(geneSet[[i]],collapse = '\t'))
    cat('\n')

  }

  sink()

}
