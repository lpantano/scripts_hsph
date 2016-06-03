library(deconstructSigs)
# library(useful)
data(signatures.nature2013)
setwd("~/orch/scratch/glasgow/signature_mutation")

# VCF files needs to be transformed to TSV files with vcf2tsv command
fns = list.files("~/orch/scratch/glasgow/signature_mutation", "tsv$", full.names = TRUE)
# this file is the one oliver gives meß
metadata = read.table("files_to_test_right_v5.txt") %>% 
    tidyr::separate(snp_file, into=c("sample"), sep="_",extra = "warn")
row.names(metadata) = metadata$sample
norm_ref = function(mut.counts){
    multiplicative.ratio      <- 1/tri.counts.genome
    norm.mut.counts           <- sapply(colnames(mut.counts), function(x) {norm.it(mut.counts[,x,drop=F], trimer.ratio = multiplicative.ratio)})
    norm.mut.counts           <- data.frame(norm.mut.counts, row.names = rownames(mut.counts))
    colnames(norm.mut.counts) <- colnames(mut.counts)
    
    # make each row sum to 1
    norm.mut.counts           <- norm.mut.counts/rowSums(norm.mut.counts)
    return(norm.mut.counts)
}

signatures.nature2013.genome = norm_ref(signatures.nature2013)

results = data.frame()
variants  = data.frame()
for (fn in fns){
    print(fn)
    df = read.table(fn, header=T, sep="\t", stringsAsFactors = F)
    df$CHROM = paste0("chr", df$CHROM)
    df$sample = strsplit(basename(fn), split = "_")[[1]][1]
    sigs.input <- mut.to.sigs.input(mut.ref = df, 
                                    sample.id = "sample", 
                                    chr = "CHROM", 
                                    pos = "POS", 
                                    ref = "REF", 
                                    alt = "ALT")
    test = whichSignatures(tumor.ref = sigs.input,
                           signatures.ref = signatures.nature2013, 
                           tri.counts.method = "exome2genome",
                           sample.id = row.names(sigs.input)[1], 
                           contexts.needed = TRUE)
    results = rbind(results, data.frame(test$weights, row.names = row.names(sigs.input)[1]))
    variants = rbind(variants, data.frame(row.names = row.names(sigs.input[1]), total=nrow(df)))
}


write.csv(results, "signatures_exons_none.csv")

# correlation with genome
pheatmap(as.matrix(results), annotation_row = metadata[,c("cancer","project")], show_rownames = F,
         cluster_cols = F,  cluster_rows = F,clustering_method = "ward.D", clustering_distance_rows = "correlation")


