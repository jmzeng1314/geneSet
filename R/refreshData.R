update_kegg <- function(refresh=F){
  if(refresh){
    dir_bin='data'
    path2gene_file=file.path(dir_bin,'KEGG_update','kegg2geneID.txt')
    path2name_file=file.path(dir_bin,'KEGG_update','kegg_hierarchical.txt')
    if (file.exists(path2gene_file)){
      keggID2geneID=read.table(path2gene_file,sep="\t",colClasses=c('character'))
      keggID2geneID_df=keggID2geneID[,c(2,1)]
      names(keggID2geneID_df)=c('gene_id','path_id')
      tmp=read.table(path2gene_file,sep="\t",colClasses=c('character'))
      #tmp=toTable(org.Hs.egPATH)
      # first column is kegg ID, second column is entrez ID
      GeneID2kegg_list<<- tapply(tmp[,1],as.factor(tmp[,2]),function(x) x)
      kegg2GeneID_list<<- tapply(tmp[,2],as.factor(tmp[,1]),function(x) x)
    }else{stop("we can not find the file:path2gene_file")}
    if (file.exists(path2name_file)){
      kegg2name<<- read.delim(path2name_file,header=F,sep="\t",colClasses=c('character'),stringsAsFactors =F)
      colnames(kegg2name)=c('parent1','parent2','pathway_id','pathway_name')
      ###kegg2name$pathway_id=as.numeric(kegg2name$pathway_id)
      rownames(kegg2name)=kegg2name$pathway_id
    }else{stop("we can not find the file:path2name_file")}
    kegg2symbol_list <- lapply(kegg2GeneID_list,function(x) as.character(geneAnno(x)$symbol))
    devtools::use_data(keggID2geneID_df,GeneID2kegg_list,kegg2GeneID_list,kegg2symbol_list, kegg2name,overwrite =T)
  }
}


