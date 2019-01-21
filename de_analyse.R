args = commandArgs(trailingOnly=TRUE)
if (length(args)==2) {
    args[3]="~/analysis/results"
}
system(sprintf('mkdir -p %s', args[3]))


library(DESeq2)
matx= read.table(args[1], sep="\t", header=T, row.names = 1, check.names =F)
control_string= strsplit(strsplit(args[2], split = "_")[[1]],split = ":")[[2]][1]
treatment_string=strsplit(strsplit(args[2], split = "_")[[1]],split = ":")[[2]][2]
col= colnames(matx)
conditions= ifelse(grepl(control_string, col), control_string, treatment_string)
coldata= data.frame(columns=col, conditions= ifelse(grepl(control_string, col), control_string, treatment_string))
matx[is.na(matx)]=0
dds= DESeqDataSetFromMatrix(countData = matx, colData = coldata, design = ~conditions)
dds= DESeq(dds)
res= results(dds)
validated=res[!is.na(res$padj),]
validated=validated[validated$padj < 0.05,]
upreg=validated[validated$log2FoldChange >= 1,]
downreg=validated[validated$log2FoldChange <= -1,]
write.table(res, file=paste(args[3],"results.table", sep="/"), sep = "\t" )
write.table(upreg, file = paste(args[3],"upregulated.table", sep="/"), sep = "\t")
write.table(downreg, file = paste(args[3],"downregulated.table", sep="/"), sep = "\t")