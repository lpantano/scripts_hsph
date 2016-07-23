root_path = "~/orch/scratch/walsh_asd/batch_1/batch1/final/2016-04-12_batch1/"

g806 = read.csv(file.path(root_path, "autism_dataset_806.csv"), stringsAsFactors = F)
g580 = read.csv(file.path(root_path, "autistic_disorder_580.csv"), stringsAsFactors = F)
ghp = read.csv(file.path(root_path, "genes_HP-0001249.csv"), header = F, stringsAsFactors = F)

keep = unlist(c(g806$Gene.Symbol, g580$symbol, ghp[,1]))

for (f in list.files(file.path(root_path,"filtered"),full.names = T, pattern = "tsv$")){
    print(f)
    fn = read.table(f, sep="\t", header = T, stringsAsFactors = F, quote = "")
    fn_f = subset(fn, gene %in% keep)
    write.table(fn_f, file.path(root_path, "filtered", "genes", basename(f)), sep="\t", quote=F, row.names=F)
}
