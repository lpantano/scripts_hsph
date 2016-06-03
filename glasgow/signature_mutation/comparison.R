library(pheatmap)
library(dplyr)
library(ggplot2)
library(gridExtra)
setwd("~/orch/scratch/glasgow/signature_mutation")

genome.signature = read.csv("v5/signatures_v5_none_norm.csv", row.names = 1)
results = read.table("signatures_exons_none.csv", sep=",", row.names = 1, header=1)
variants = read.table("signatures_exons_exome_cor.csv", sep=",", row.names = 1, header=1)

metadata = read.table("files_to_test_right_v5.txt") %>% 
    tidyr::separate(snp_file, into=c("sample"), sep="_",extra = "warn")
rownames(metadata) = metadata$sample

compare = function(results, genome.signature, variants, metadata){
    cor = sapply(rownames(results), function(sample){
        a = unlist(results[sample,])
        b = unlist(genome.signature[sample,])
        # keep = a>0 | b>0
        # e=(cor.test(a[keep],b[keep], method="spearman")$estimate)
        a_top = names(sort(a[a>0.099], decreasing = TRUE))[1:3]
        a_top = a_top[!is.na(a_top)]
        b_top = names(sort(b[b>0.099], decreasing = TRUE))[1:3]
        b_top = b_top[!is.na(b_top)]
        print(intersect(a_top, b_top))
        common=intersect(a_top, b_top)
        e = length(common)
        e = e/mean(c(length(a_top), length(b_top)))
        if (a_top[1]==b_top[1])
            e =1
        # plot(unlist(results[sample,]), unlist(genome.signature[sample,]),main=paste(sample,e))
        e
    })
    variants$cor = cor
    
    variants$cancer = metadata[row.names(variants), "cancer"]
    variants$max = rowMax(as.matrix(results))
    variants$max_none = rowMax(as.matrix(genome.signature))
    variants$max_name = apply(as.matrix(results),1,function(x){names(results)[which.max(x)]})
    variants$max_name_none = apply(as.matrix(genome.signature),1,function(x){names(results)[which.max(x)]})
    variants$strong = rowSums(as.matrix(results)>0)
    variants$max_none = rowMax(as.matrix(genome.signature))
    variants$strong_none = rowSums(as.matrix(genome.signature)>0)
    
    variants$cor_cat = "low"
    variants$cor_cat[variants$cor>0.9] = "high"
    
    # ggplot(variants, aes(x=max,y=log10(total), color=strong)) + geom_point() +facet_wrap(~cor_cat)
    
    
    p1 = ggplot(variants, aes(y=max,x=log10(total), color=as.factor(cor),shape=cancer)) + geom_point() 
    
    p2 = ggplot(variants, aes(x=cancer, fill=max_name)) + geom_bar() +facet_wrap(~cor)
    
    p3 = ggplot(variants, aes(x=cancer, fill=max_name_none)) + geom_bar() +facet_wrap(~cor)
    
    grid.arrange(p1,p2,p3)
    
    variants
}

var = compare(results, genome.signature, variants, metadata)
