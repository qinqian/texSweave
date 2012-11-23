library(affy)
library(affyPLM)
library(AnnotationDbi)                 # annotation db
library(hgu133plus2hsrefseqcdf)        # refseq annotation
library(hgu133plus2hsrefseq.db)        # chrom and genesname info  
library(hgu133plus2hsentrezgcdf)       # entrez
library(hgu133plus2hsentrezg.db)       # chrom and genesname info
library(limma)
library(ggplot2)

# set directory to .CEL directory
setwd("~/Desktop/2012_fall_training/3rd/GSE10437/")
rna.workflow <- function(pat, pd, method, cdf){
  # rna expression process workflow
  # pat for files suffix
  # pd for phenotype data
  # cdf for annotation files
  files <- dir(pattern=pat)
  pd <- read.AnnotatedDataFrame(pd,header=T,
                               row.names=1, sep="", as.is=T)
  data.affy <- read.affybatch(filenames=files, phenoData=pd)
  
  #data.affy <- ReadAffy()
  data.affy@cdfName <- cdf
  
  if (method=='rma') {
    data.result <- rma(data.affy)
    data.expr <- exprs(data.result)
    data.expr.PA <- ''
  }
  if (method=='mas'){
    data.result <- mas5(data.affy)
    data.expr <- exprs(data.result)   # expression value
    # use log 2 instead of e
    data.expr <- log2(data.expr)
    data.expr.PA <- exprs(mas5calls(data.affy)) # expression pattern, see help
  }
  
  # ajust ids
  genes_t = matrix(rownames(data.expr))
  genes.refseq=apply(genes_t, 1, function(x) sub("_at", "", x))
  orig.refseq=rownames(data.expr)
  
  rownames(data.expr) <- genes.refseq
  ## QC measurements, boxplot, MVA plot
  boxplot(data.frame(data.expr))
  
  if (all(data.expr.PA!='')) {
    data<-list(refseq=genes.refseq, expr=data.expr, eset=data.result, PA=data.expr.PA, orig=orig.refseq)
  }
  else{
    data<-list(refseq=genes.refseq, expr=data.expr, eset=data.result,orig=orig.refseq)
  }
  return(data)
}

rna.expression.rma <- rna.workflow('.CEL', 'GSE10437.txt',
                                   'rma', 'hgu133plus2hsentrezgcdf') # "hgu133plus2hsrefseqcdf"
print(length(rna.expression.rma$expr[,1]))

rna.expression.mas5 <- rna.workflow('.CEL', 'GSE10437.txt',
                                    'mas', 'hgu133plus2hsentrezgcdf') # "hgu133plus2hsrefseqcdf"
write.table(rna.expression.mas5$expr, file='rna_mas5.txt', sep='\t', quote=F)
entrezid <- featureNames(rna.expression.mas5$eset)
affy.control <- grep('^A', entrezid, perl=T)

# remove affymetrix control probes
rna.expression.mas5$expr <- rna.expression.mas5$expr[-affy.control,]
rna.expression.mas5$PA <- rna.expression.mas5$PA[-affy.control,]
rna.expression.mas5$orig <- rna.expression.mas5$orig[-affy.control]
# remove Absent gene and missing genes
rna.expression.mas5$expr <- rna.expression.mas5$expr[-affy.control,]
rna.absent <- apply(rna.expression.mas5$PA, 1, paste, collapse='')
rna.absent.index <- grep('AAAA', rna.absent)
rna.expression.mas5$expr <- rna.expression.mas5$expr[-rna.absent.index,]
rna.expression.mas5$orig <- rna.expression.mas5$orig[-rna.absent.index]

## rna differential expression 
rna.diff <- function(expr='', ref='', eset='', control, treat, method){
  if (method=='foldT'){
    # fold change and T-test
    foldchange=apply(expr, 1, function(x) mean(x[treat])-mean(x[control]))
    T.p.value=apply(expr, 1, 
                    function(x) t.test(x[treat], x[control], var.equal=T)$p.value)
    #fdr=p.adjust(T.p.value, method="BH") # BH, Bonferroni, fdr
    fdr = T.p.value

    # genes
    genes.up = expr[which(fdr<0.05 & foldchange>0)]

    genes.down = expr[which(fdr<0.05 & foldchange<0)]
    genes.id = c(which(fdr<0.05 & foldchange>0), which(fdr<0.05 & foldchange<0))    
    return(expr[genes.id,])
  }
  # limma group means
  if (method=='gm'){
    gm.design<-model.matrix(~ 0+factor(c(1,1,2,2)))
    colnames(gm.design) <- c("control","treat")
    #gm.design = cbind(control = control, treat = treat)
    print(gm.design)
    gm.fit = lmFit(expr, gm.design)
    print(gm.fit)
    gm.matrix = makeContrasts(CvsT = control-treat, levels=gm.design)
    gm.fit = contrasts.fit(gm.fit, gm.matrix)
    gm.fit = eBayes(gm.fit)
    # coef=2 for between-group difference for two group comparison, coef=1 for intercept
    # coef=1,2,3 for multiple group between-group differences
    diff.expr <- topTable(gm.fit, p.value=0.05, lfc=2, number=length(expr[,1]), adjust.method="BH") #coef='CvsT',
    return(diff.expr)
  }
  # limma para
}
diff.result.ft <- rna.diff(rna.expression.mas5$expr, rna.expression.mas5$refseq, 
                           rna.expression.mas5$eset,
                           c(1,2), c(3,4),'foldT')

diff.result.gm <- rna.diff(rna.expression.mas5$expr, rna.expression.mas5$refseq, 
                           rna.expression.mas5$eset,
                           c(1,1,0,0), c(0,0,1,1),'gm')

# heatmap of differential expression genes
colnames(rna.expression.mas5$expr) = c(rep('control', 2), rep('treat', 2))
heatmap(rna.expression.mas5$expr[as.numeric(rownames(diff.result.gm)),], margins=c(5,10), 
        col=c('red','blue'))



# mappedKeys from AnnotateDbi
# geneplotter cplot for chromesome mapping
# whole annotation for hgu133plus2 
myAnnot <- data.frame(
    ID=sapply(contents(hgu133plus2hsentrezgENTREZID), paste, collapse=", "), 
    SYMBOL=sapply(contents(hgu133plus2hsentrezgSYMBOL), paste, collapse=", "), 
    DESC=sapply(contents(hgu133plus2hsentrezgGENENAME), paste, collapse=", "),
    CHR=sapply(contents(hgu133plus2hsentrezgCHR), paste, collapse=", ")
)
head(myAnnot)

# chrosome location for start and end 
entrez.diff <- rna.expression.mas5$orig[as.numeric(rownames(diff.result.gm))]
chr <- toTable(hgu133plus2hsentrezgCHR[entrez.diff])
#chrloc.end <- toTable(hgu133plus2hsrefseqCHRLOCEND[entrez.diff])
#chrloc.start <- toTable(hgu133plus2hsentrezgCHRLOC[entrez.diff])
#chrloc.end <- toTable(hgu133plus2hsrefseqCHRLOCEND[entrez.diff])

# find the probe annotationsead
# for (i in 1:length(entrez.diff)){a<-c(a, grep(paste('^',entrez.diff[i],'$',sep=""), rownames(myAnnot)))}
anno <- myAnnot[entrez.diff,]

## use merge and table to do statistics on chromosome genes number
merge.diffwithanno <- merge(myAnnot[,c(1,4)],diff.result.gm[,c(1,3)],all=F)

chrNames=c("Y","X", "22","21","20","19","18","17","16","15","14","13",
           "12","11","10","9","8","7","6","5","4","3","2","1")

chr.statstable <- table(merge.diffwithanno$CHR)[chrNames]
barplot(chr.statstable)

chr.stats <- list()
for (i in chrNames){
  pat=paste('^',i,'$', sep='')
  print(pat)
  print(length(grep(pattern=pat, chr[,2], perl=T)))
  if (length(grep(pattern=pat, chr[,2], perl=T))>=1)
    {chr.stats[[i]]<-length(grep(pattern=pat, chr[,2], perl=T))}
}

# overall chromesome distribution
chr.all <- sapply(chr.stats, paste, collapse="")

#position<-barplot(as.numeric(chr.all),plot=F)
barplot(as.numeric(chr.all), names=names(chr.all),horiz=T, xlim=c(0,100), width=5, ylim=c(-5, 130),col='black',cex.axis=2)
p
for (i in 1:length(chrNames)){
  mtext(chrNames[i], at=c(-1+5.9*i,1),side=2,line=0.2, outer=F, cex=0.7,padj=0)
}

#chrlen <- hgu133plus2hsrefseqCHRLENGTHS[c(1:22,'M','X','Y')]

### barplot for exercise 4, fisher.test chromosome number
chr.number.test <- function(diffgenes,chrdata,n) {
  # chrdata = for differentiatial expressed genes on chromosomes 
  # n denotes which chromosomes

  chrn.diff.len <- length(grep(paste('^', n, '$', sep=""),chrdata[,2],perl=T))
  print(chrn.diff.len)
  CHR=sapply(contents(hgu133plus2hsentrezgCHR), paste, collapse=", ")
  chrn.len <- length(grep(paste('^', n, '$', sep=""),CHR,perl=T))
  nonchrn.diff <- length(diffgenes) - chrn.diff.len
  nonchrn.len <- length(CHR) - chrn.len  
  chrn.nondiff.len <- chrn.len - chrn.diff.len
  nonchrn.nondiff <- nonchrn.len - nonchrn.diff
  # matrix for fisher test
  chrndiff.fisher <- matrix(c(chrn.diff.len, chrn.nondiff.len, nonchrn.diff, nonchrn.nondiff),
                            nrow=2, dimnames=list(chr=c(paste('chr', n), paste('nonchr',n)), diff=c('diff','nondiff')))
  cat('greater side p.value', fisher.test(chrndiff.fisher,alternative='greater')$p.value,'\n')
  cat('less side p.value', fisher.test(chrndiff.fisher,alternative='less')$p.value,'\n')
  cat('two sided p.value',fisher.test(chrndiff.fisher,alternative='two.sided')$p.value,'\n')
  return(chrndiff.fisher)
}
diff.fisher<-chr.number.test(entrez.diff, chr, 1)
