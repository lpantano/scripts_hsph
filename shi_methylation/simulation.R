# Mike Love
# May 2, 2016

# make some simulated data across samples

# let's say all the methyl probabilities are set to 30%
n <- 400
meth <- 0.4
hmeth <- 0.4
p <- meth + hmeth
p.oxy <- p-hmeth
pct <- round(n*0.3)
COV <- 50
# make a pair of total count and methylation count for 'n' sites
makeSample <- function(n,p) {
    counts <- rpois(n, lambda=COV)
    methyl <- rbinom(n, size=counts, prob=p)
    cbind(counts,methyl)
}

# our colData
sample <- factor(rep(1:3, each=4))
oxy <- factor(rep(rep(0:1,each=2),3))
methyl <- factor(rep(0:1,6))
coldata <- data.frame(sample,oxy,methyl)
coldata

# make a matrix of data
data <- do.call(cbind, lapply(1:6, function(i) makeSample(n,p)))
head(data)

# make some true positives
# data[1:10,c(4,8,12)] <- round(data[1:10,c(3,7,11)]/2)
data[1:pct,c(3,4,7,8,11,12)] <-  do.call(cbind, lapply(1:3, function(i) makeSample(pct,p.oxy)))

# useful for splitting total counts and methylation counts
evens <- 1:6 * 2

# calculate empirical proportions: methyl counts over total counts
empirical.props <- data[,evens] / data[,evens-1]
hist(as.vector(empirical.props), breaks=0:100/100, col="black")

#############################################
# a binomial GLM per row approach
# test <- function(i) {
#     fit <- glm(empirical.props[i,] ~ oxy[evens], 
#                weights=data[i,evens-1], family=binomial(logit))
#     # do not use a Wald test here! R's glm() is bad for this...
#     # use a likelihood ratio test:
#     sumfit <- summary(aov(fit)) 
#     sumfit[[1]][1,5] # the LRT pvalue
# }

# do it and plot pvalues
#pvals <- sapply(seq_len(n), function(i) test(i))
col <- ifelse(seq_len(n) <= pct, "red", "black")
#plot(-log10(pvals), col=col, pch=16)

############################################
# using a NB and interaction term approach

# the size factors should take care of sample diffs in sequencing depth
# (it doesn't really matter that it corrects each *column*,
# you could tweak to share size factors between every *pair* of columns
# but inference will look nearly the same)
library(DESeq2)
rownames(coldata) = paste(rownames(coldata), coldata$sample, coldata$oxy, coldata$methyt ,sep = ":")
colnames(data) = rownames(coldata)
dds <- DESeqDataSetFromMatrix(data, coldata, ~ oxy + methyl + oxy:methyl)
dds <- DESeq(dds, fitType="mean")
res <- results(dds, name="oxy1.methyl1")
hist(res$pvalue)
plot(-log10(res$pvalue), col=col, pch=16)
abline(h=-log10(0.05), col="green")
pass = res$padj<0.1 & !(is.na(res$padj))
if(sum(pass)>0){
    mpval=max(res$pvalue[pass])
    abline(h=-log10(mpval), col="blue")
}
