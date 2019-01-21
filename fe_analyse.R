library(AnnotationDbi)
library(GeneAnswers)
library(org.Hs.eg.db)

args=commandArgs(trailingOnly = T)
DEGenes=read.table(args[1], header = T, sep ="\t")

DEGenes$SYMBOL=row.names(DEGenes)
ak=AnnotationDbi::select(org.Hs.eg.db, keys=DEGenes$SYMBOL, columns="ENTREZID", keytype = "SYMBOL")
annotated_DEGenes=merge(DEGenes, ak, by.x='SYMBOL',by.y='SYMBOL', all.x=T)
annotated_DEGenes=annotated_DEGenes[match(DEGenes$SYMBOL, annotated_DEGenes$SYMBOL),]
ids=annotated_DEGenes$ENTREZID
ids=unique(ids)
ids=ids[!is.na(ids)]

#GO FEA
GAinstance=geneAnswersBuilder(ids, annotationLib = 'org.Hs.eg.db', categoryType = 'GO.BP', testType = "hyperG", FDR.correction = T)
ei=GAinstance@enrichmentInfo
terms=AnnotationDbi::select(GO.db, row.names(ei), columns = c('GOID','TERM'))
enrichmentInfo=cbind(terms,ei$`genes in Category`,ei$`percent in the observed List`,ei$`percent in the genome`,ei$`fold of overrepresents`,ei$`odds ratio`,ei$`p value`, ei$`fdr p value`)
colnames(enrichmentInfo)[3:9]=colnames(ei)
write.table(enrichmentInfo, file=paste(args[2],'table', sep='.'),sep = '\t')
#represent top ten
pdf(paste(args[2], 'pdf', sep = '.'))
par(mai=c(1,4.5,2,1))
bp=barplot(-log10(enrichmentInfo$`fdr p value`[1:10]),horiz=T,col=rainbow(10),xlab = 'pFDR-p-value of enrichment',names.arg=as.vector(enrichmentInfo$TERM[1:10]),las=1, main=args[2])
dev.off()
