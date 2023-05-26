countdata <- read.table("/rds/user/cjr78/hpc-work/ATAC/DEseq2/ATAC_counts_peaks_merged.txt", header=TRUE, row.names=1)

### remove bam
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))

###remove first 5 columns from featurecounts output - do not contain data
countdata <- countdata[ ,6:ncol(countdata)] 

###set conditions
countdata <- as.matrix(countdata)
(condition <- factor(c(rep("Fhfl", 1), rep("pFh", 1), rep("cl1", 1), rep("cl19", 1))))

library(DESeq2) #load DESeq2

(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
keep <- rowSums(counts(dds)) >= 10 ###filter genes with low counts
dds <- dds[keep,]

dds <- DESeq(dds)

sizeFactors(dds) #output scale factors
vsd <- vst(dds, blind=FALSE)
write.table(as.data.frame(assay(vsd)), file="VST_values.txt", sep = '\t')


res <- results(dds, contrast=c("condition","Fhfl","cl1"))
resOrdered <- res[order(res$pvalue),]
write.table(as.data.frame(resOrdered), file="allresults_Fhfl_cl1.txt", sep ='\t')
resSig <- subset(resOrdered, padj < 0.05)
write.table(as.data.frame(resSig), file="sigresults_Fhfl_cl1.txt", sep ='\t')
(summary(res))

#plotMA(res, ylim=c(-10,10), cex=0.8, xlab="Mean of normalised counts", ylab="Log2(fold change)", colNonSig = "gray32", colSig = c(ifelse(res$log2FoldChange>=-2.3219, "red3" ,"gray32"))))

#plotMA(res, ylim=c(-10,10), cex=0.8, xlab="Mean of normalised counts", ylab="Log2(fold change)", colNonSig = "gray32", colSig = c(ifelse(abs(res$log2FoldChange)>2.3219, "red3" ,"gray32"))) #this highlights significant points over fold change 5
