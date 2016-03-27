library(deconstructSigs)
# library(useful)
data(signatures.nature2013)

results = data.frame()
# VCF files needs to be transformed to TSV files with vcf2tsv command
fns = list.files("~/orch/scratch/glasgow/signature_mutation", "tsv$", full.names = TRUE)[1:2]
# this file is the one oliver gives meß
metadata = read.table("files_to_test_right_v5.txt") %>% 
    tidyr::separate(snp_file, into=c("sample"), sep="_",extra = "warn")

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
                           sample.id = row.names(sigs.input)[1], contexts.needed = TRUE)
    results = rbind(results, data.frame(test$weights, row.names = row.names(sigs.input)[1]))
    
}


# write.csv(results, "signatures_v5_none_norm.csv")
library(pheatmap)
library(dplyr)
rownames(metadata) = metadata$sample
pheatmap(results, annotation_row = metadata[,c("cancer","project")], cluster_cols = F, clustering_method = "ward.D", clustering_distance_rows = "correlation")
