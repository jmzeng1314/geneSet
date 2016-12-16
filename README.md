# geneSet(A R package for comprehensive geneSet analysis)
------
So far, we provide 3 basic statistic analysis for the expression matrix according to user defined gene set, below are the function I defined:
> * hyperGTest_single
> * hyperGTest_multiple
> * GSEA_gene_single
> * GSEA_gene_multiple
> * plot_GSEA_single
> * PCA_gene_single
> * PCA_gene_multiple 


Please follow the link below to install this package and have a fun !

### [如何安装别人开发的未发表的包](http://www.bio-info-trainee.com/2092.html)

Also please feel free to contact with me if there's a bug or it's not clear to use.
>  I perfer email communication: jmzeng1314@163.com 

you can follow the codes below to get familar with this package:
```R
suppressMessages(library(CLL))
data(sCLLex)
suppressMessages(library(limma))
exprSet = exprs(sCLLex)
pdata=pData(sCLLex)
group_list = pdata$Disease
library(humanid) ## https://github.com/jmzeng1314/humanid
exprSet = get_symbol_exprSet(sCLLex,platformDB='hgu95av2.db')

## firstly do GSEA:
s2n <- calc_s2n(exprSet,group_list ) 
tmp <- GSEA_gene_single(s2n,names(s2n),kegg2symbol_list[[1]],n = 100)
## it's too slow, I need to modify it.
plot_GSEA_single(GSEA_single_results = tmp,gene_values = s2n,multiple = F)
## don't run this, I guess it will consume more than 2 hours
GSEA_gene_multiple(s2n,names(s2n),head(kegg2symbol_list),100)

## Then do PCA:
PCA_gene_single(exprSet=exprSet,this_geneSet=kegg2symbol_list[[1]])
PCA_gene_multiple(exprSet=exprSet,geneSet_list=kegg2symbol_list)

## lastly do hyperGTest:
design=model.matrix(~factor(sCLLex$Disease))
fit=lmFit(exprSet,design)
fit=eBayes(fit)
DEG=topTable(fit,coef=2,adjust='BH',n=Inf)
logFC_cutoff  <- mean(abs(DEG$logFC)) + 2*sd(abs(DEG$logFC))
pvalue_cutoff  <- 0.05
DEG$symbol <- rownames(DEG)
DEG$sigORnot <- ifelse(abs(DEG$logFC) > logFC_cutoff & DEG$P.Value <pvalue_cutoff ,
                ifelse(DEG$logFC >logFC_cutoff,'UP','DOWN'),'NOT')
choose_geneIds <- DEG[DEG$sigORnot != 'NOT','symbol']
backgroud_geneIds <- DEG$symbol
hyperGTest_single(kegg2symbol_list[[1]],choose_geneIds,backgroud_geneIds)
hyperGTest_multiple(kegg2symbol_list,choose_geneIds,backgroud_geneIds)

 
```



