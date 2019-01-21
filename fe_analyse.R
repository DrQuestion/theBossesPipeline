
###MIT License
###
###Copyright (c) 2019 Alessio Albanese
###
###Permission is hereby granted, free of charge, to any person obtaining a copy
###of this software and associated documentation files (the "Software"), to deal
###in the Software without restriction, including without limitation the rights
###to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
###copies of the Software, and to permit persons to whom the Software is
###furnished to do so, subject to the following conditions:
###
###The above copyright notice and this permission notice shall be included in all
###copies or substantial portions of the Software.
###
###THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
###IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
###FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
###AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
###LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
###OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
###SOFTWARE.




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
