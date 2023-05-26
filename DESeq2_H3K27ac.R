countdata <- read.table("H3K27ac_counts_peaks_reps.txt", header=TRUE, row.names=1)

### remove bam
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))

###remove first 5 columns from featurecounts output - do not contain data
countdata <- countdata[ ,6:ncol(countdata)] 

###set conditions
countdata <- as.matrix(countdata)
(condition <- factor(c(rep("Fhfl", 3), rep("pFh", 3), rep("cl1", 3), rep("cl19", 3))))

library(DESeq2) #load DESeq2

(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
keep <- rowSums(counts(dds)) >= 10 ###filter genes with low counts
dds <- dds[keep,]

dds <- DESeq(dds)

sizeFactors(dds) #output scale factors
vsd <- vst(dds, blind=FALSE)
write.table(as.data.frame(assay(vsd)), file="VST_values_H3K27ac_reps.txt", sep = '\t')

#change condition here depending on samples
res <- results(dds, contrast=c("condition","cl19","Fhfl"))
resOrdered <- res[order(res$pvalue),]
write.table(as.data.frame(resOrdered), file="allresults_Fhfl_cl19.txt", sep ='\t')
resSig <- subset(resOrdered, padj < 0.05)
write.table(as.data.frame(resSig), file="sigresults_Fhfl_cl19.txt", sep ='\t')
(summary(res))

#plotMA(res, ylim=c(-10,10), cex=0.8, xlab="Mean of normalised counts", ylab="Log2(fold change)", colNonSig = "gray32", colSig = c(ifelse(res$log2FoldChange>=-2.3219, "red3" ,"gray32"))))
pdf("Fhfl_cl19_q0.05_fc2_maplot.pdf")
plotMA(res, ylim=c(-7,7), cex=0.8, xlab="Mean of normalised counts", ylab="Log2(fold change)", colNonSig = "gray32", colSig = c(ifelse(abs(res$log2FoldChange)>1, "red3" ,"gray32"))) #this highlights significant points over fold change 2
dev.off()


