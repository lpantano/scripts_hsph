project_summary = "project-summary.csv"
counts_file = "combined.counts"
summarydata = data.frame(read.table(project_summary, header=TRUE, sep=","), row.names="Name", check.rows=FALSE)
summarydata$Name = rownames(summarydata)
summarydata = summarydata[order(summarydata$Name),]
summarydata$time = as.factor(summarydata$time)
counts = read.table(counts_file, header=TRUE, row.names="id", check.names=FALSE)
counts = counts[, order(colnames(counts))]
colnames(counts) = gsub(".counts", "", colnames(counts))

library(DESeq2)
library(reshape)
library(ggplot2)
counts <- counts[rowSums(counts>0)>3,]

# Differential expression by time only for DMSO
samples = summarydata$Name[summarydata$treatment=="DMSO"]
dds = DESeqDataSetFromMatrix(countData=counts[,samples],
                             colData=summarydata[samples,], design = ~ time)
dds = DESeq(dds, reduced = ~ 1, test = "LRT")
res = results(dds)
#to give a summary of the DE analysis
summary(res)
#plot one gene to see if it makes sense
head(res[order(res$padj),])
plotCounts(dds, row.names(res[order(res$padj),])[1], intgroup="time")


# Differential expression by time for other group not seen in DMSO
samples = summarydata$Name[(summarydata$concentration==3 & summarydata$treatment=="Sora") |
                           summarydata$treatment=="DMSO"]
coldata = droplevels(summarydata[samples,c("time", "treatment","concentration")])
full <- model.matrix(~ time + time:treatment, coldata)
reduced <- full[,c(-2,-3,-4)]

dds = DESeqDataSetFromMatrix(countData=counts[,samples],
                             colData=summarydata[samples,], design = ~ time + time:treatment)
dds = DESeq(dds, full=full, reduced = reduced, test = "LRT")

res = results(dds)
#to give a summary of the DE analysis
summary(res)
#plot one gene to see if it makes sense
head(res[order(res$padj),])
dd = plotCounts(dds, row.names(res[order(res$padj),])[1], intgroup="time", returnData = T)
dd$treatment = summarydata[row.names(dd), "treatment"]
ggplot(dd, aes(x=time,y=count,color=treatment)) +
  geom_point()