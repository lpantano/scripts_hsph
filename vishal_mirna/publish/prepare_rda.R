# root = "~/orch/scratch/vishal_mirna_kidney/publish"

# FA model
root_path = file.path(root, "FA-model")
order_group=c("control", "day1", "day2", "day3", "day7", "day14", "day28")
mrna_path = "rnaseq/final/2016-02-05_mrna_bio/files_publish"
mrna_matrix = read.csv(file.path(root_path, mrna_path, "fa_model_log2_counts.tsv"), row.names = 1)
mrna_clusters = read.csv(file.path(root_path, mrna_path, "clusters_genes.tsv"), row.names = 1)
mrna_col = data.frame(row.names=colnames(mrna_matrix), samples=colnames(mrna_matrix)) %>% tidyr::separate(samples,into = c("group"), extra = "drop") %>% mutate(group=factor(group, levels=order_group))
rownames(mrna_col) = colnames(mrna_matrix)

prot_path = "protein/files_publish"
prot_matrix = read.csv(file.path(root_path, prot_path, "fa_model_log2_counts.tsv"), row.names = 1)
prot_clusters = read.csv(file.path(root_path, prot_path, "clusters_genes_cluster.tsv"), row.names = 1)
prot_col = data.frame(row.names=colnames(prot_matrix), samples=colnames(prot_matrix)) %>% tidyr::separate(samples,into = c("group"), extra = "drop") %>% mutate(group=factor(group, levels=order_group))
rownames(prot_col) = colnames(prot_matrix)

mirna_path = "srnaseq/final/files_publish"
mirna_matrix = read.csv(file.path(root_path, mirna_path, "fa_model_mirna_log2_counts.tsv"), row.names = 1)
load(file.path(root_path, mirna_path, "ma.rda"))
mirna_col = data.frame(row.names=colnames(mirna_matrix), samples=colnames(mirna_matrix)) %>% tidyr::separate(samples,into = c("id", "group"), extra = "drop") %>% mutate(group=factor(gsub("day0", "control", tolower(group)), levels=order_group))
rownames(mirna_col) = colnames(mirna_matrix)

fa_model = list(mirna = SummarizedExperiment(assays=SimpleList(norm=as.matrix(mirna_matrix)), 
                                             colData=mirna_col),
                gene = SummarizedExperiment(assays=SimpleList(norm=as.matrix(mrna_matrix)), 
                                             colData=mrna_col),
                protein = SummarizedExperiment(assays=SimpleList(norm=as.matrix(prot_matrix)), 
                                             colData=prot_col)
                )

# UUO model
root_path = file.path(root, "UUO-model")
order_group=c("normal", "day3", "day7", "day14")
mrna_path = "rnaseq/final/2016-02-03_mrna/files_publish"
mrna_matrix = read.csv(file.path(root_path, mrna_path, "uuo_model_log2_counts.tsv"), row.names = 1)
mrna_clusters = read.csv(file.path(root_path, mrna_path, "clusters_genes.tsv"), row.names = 1)
mrna_col = data.frame(row.names=colnames(mrna_matrix), samples=colnames(mrna_matrix)) %>% tidyr::separate(samples,into = c("group"), extra = "drop") %>% mutate(group=factor(group, levels=order_group))
rownames(mrna_col) = colnames(mrna_matrix)

prot_path = "protein/files_publish"
prot_matrix = read.csv(file.path(root_path, prot_path, "uuo_model_log2_counts.tsv"), row.names = 1)
prot_clusters = read.csv(file.path(root_path, prot_path, "clusters_genes_cluster.tsv"), row.names = 1)
prot_col = data.frame(row.names=colnames(prot_matrix), samples=colnames(prot_matrix)) %>% tidyr::separate(samples,into = c("group"), extra = "drop") %>% mutate(group=factor(group, levels=order_group))
rownames(prot_col) = colnames(prot_matrix)

mirna_path = "srnaseq/final/files_publish"
mirna_matrix = read.csv(file.path(root_path, mirna_path, "uuo_model_mirna_log2_counts.tsv"), row.names = 1)
mirna_col = data.frame(row.names=colnames(mirna_matrix), samples=colnames(mirna_matrix)) %>% tidyr::separate(samples,into = c("group"), extra = "drop") %>% mutate(group=factor(gsub("day0", "control", tolower(group)), levels=order_group))
rownames(mirna_col) = colnames(mirna_matrix)


uuo_model = list(mirna = SummarizedExperiment(assays=SimpleList(norm=as.matrix(mirna_matrix)), 
                                             colData=mirna_col),
                gene = SummarizedExperiment(assays=SimpleList(norm=as.matrix(mrna_matrix)), 
                                            colData=mrna_col),
                protein = SummarizedExperiment(assays=SimpleList(norm=as.matrix(prot_matrix)), 
                                               colData=prot_col)
)


obj = list(fa=fa_model, uuo=uuo_model)

save(obj, file=file.path(root, "se.rda"))
