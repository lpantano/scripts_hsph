mirna_results = "~/orch/scratch/vishal_mirna_kidney/FA_model/mirna_bio/final/files/mirna_log_matrix.txt"
mirna_design = "~/orch/scratch/vishal_mirna_kidney/FA_model/mirna_bio/final/2016-02-08_mirna_bio/report/summary_re.csv"
mirna_de = "~/orch/scratch/vishal_mirna_kidney/FA_model/mirna_bio/final/files/mirna_de.tsv"

mrna_results = "~/orch/scratch/vishal_mirna_kidney/FA_model/mrna_bio/final/2016-02-05_mrna_bio/files/rlog_counts.tsv"
mrna_design = "~/orch/scratch/vishal_mirna_kidney/FA_model/mrna_bio/final/2016-02-05_mrna_bio/project-summary.csv"
mrna_de = "~/orch/scratch/vishal_mirna_kidney/FA_model/mrna_bio/final/2016-02-05_mrna_bio/files/fa_model_simple.tsv"

mirna_ma = as.matrix(read.table(mirna_results, header=T, row.names=1))
mrna_ma = as.matrix(read.table(mrna_results, header=T, row.names=1))
mrna_ma=mrna_ma[rowSums(mrna_ma)>0,]

mirna_des = read.csv(mirna_design, row.names="sample_id")
mrna_des = read.csv(mrna_design, row.names="descritiption")

mirna_sign = read.table(mirna_de, header=T, sep="\t", row.names = 1)
mrna_sign = read.table(mrna_de, header = T, sep = "\t", row.names = 1)

library(dplyr)
mirna_keep = row.names(mirna_sign)[mirna_sign$padj<0.05]
write.table(data.frame(m=mirna_keep), "mirna_diff.txt",row.names=F, col.names=F, quote=F)
mrna_keep = row.names(mrna_sign)[mrna_sign$padj<0.05 & abs(mrna_sign$log2FoldChange)]
mrna_keep = intersect(mrna_keep, row.names(mrna_ma)) #' fix this


library(CHBUtils)
pairs=read.table("mirna_pairs/matrix.tsv") %>% filter(V3< -.1)
map = data.frame(tx=unique(pairs$V1))
map_gene = annotate_df(map, "tx", 'mmusculus_gene_ensembl', "ensembl_transcript_id", "ensembl_gene_id")
idx = match(pairs$V1, map_gene$tx)
pairs$gene = map_gene[idx,2]
pairs$value = 1
library(tidyr)
ma = pairs %>% dplyr::select(gene,V2,value) %>%
  dplyr::distinct() %>% spread(V2,value, fill=0) %>% filter(!is.na(gene))
row.names(ma) = ma$gene
save(ma, file="mirna_pairs/ma.rda")
load("mirna_pairs/ma.rda")

# get median expresion by group
.median_by_group = function(e, g){
  sapply(g, function(i){
    idx = which(g == i)
    mean(e[idx], na.rm=TRUE)
  })
}

mirna_norm = t(apply(mirna_ma, 1, function(x, g){
  .median_by_group(x, g)
}, levels(mirna_des$group)))

mrna_norm = t(apply(mrna_ma, 1, function(x, g){
  .median_by_group(x, g)
}, levels(mrna_des$group)))


# GO enrichment with clusterprofile with total mrna expression
# as universe, then # number of mirna targeting that network.
mrna_target = intersect(rownames(ma), mrna_keep)
mirna_target = intersect(colnames(ma), mirna_keep)
cor = cor(t(mirna_norm[mirna_target, -5]), t(mrna_norm[mrna_target,]), method="kendall")

#get only cor values when gene is a target
cor_target = ma[mrna_target,mirna_target] * cor

# get genes that only get a negative correlation
keep = apply(cor_target, 1, min) < -.6
is_target_and_de = rownames(cor_target)[keep]

library(clusterProfiler)
library(KEGG.db)
library(org.Mm.eg.db)
library(CHBUtils)
library(reshape)

uni = convertIDs(rownames(mrna_ma)[rowSums(mrna_ma>0)>3],
                 from = "ENSEMBL", "ENTREZID",db = org.Mm.eg.db )
genes = convertIDs(is_target_and_de,
                 from = "ENSEMBL", "ENTREZID",db = org.Mm.eg.db )

ego = enrichKEGG(genes, "mmu", use_internal_data = TRUE,universe = uni[!is.na(uni)])
ego = enrichGO(genes, "org.Mm.eg.db", ont = "BP", universe = uni[!is.na(uni)])
# ego = enrichGO(genes, "org.Mm.eg.db", ont = "BP")

tab = ego@result %>% filter(Count<30)
ego_s = ego
ego_s@result = tab
ego_s = clusterProfiler::simplify(ego_s)
res =summary(ego_s)

net = do.call(rbind, apply(res, 1, function(x){
  .genes = unlist(strsplit(x[8],split = "/"))
  .genes = convertIDs(.genes,
             "ENTREZID", "ENSEMBL", db = org.Mm.eg.db )
  .genes = .genes[!is.na(.genes)]
  .idx = which(colSums(cor_target[.genes,] < -.6) > 0)
  if (length(.idx) > 0)
    .tab = data.frame(melt(as.matrix(cor_target[.genes,.idx]))) %>%
    filter(value < -.6) %>%
    dplyr::select(gene=X1,mir=X2) %>%
    mutate(goid=x[1],go=x[2])
    return(.tab)
}))

res_by_mir = net %>% group_by(mir, go) %>% summarise() %>%
  group_by(go) %>% summarise(nmir=n())

res[match(res$Description, res_by_mir$go), "nmir"] = res_by_mir$nmir

res %>% dplyr::select(Description, Count, nmir)


lapply(res_by_mir$go[1:2], function(x){
  .m = as.character(net[net$go==x, "mir"])
  pheatmap(mirna_ma[.m,], annotation_col = mirna_des[,"group", drop=F])
})

# kegg annotation
# change to the specie you want
library(goseq) #For getting KEGG pathways easily
library(KEGG.db)

kegg = getgo(row.names(mrna_ma), fetch.cats = "KEGG", genome = "mm10", id = "ensGene")
kegg2ens = Biobase::reverseSplit(kegg)
kname = mget(names(kegg2ens), env = KEGGPATHID2NAME)
names(kegg2ens) = kname

# mirPath
cor = cor(t(mirna_norm[, -5]), t(mrna_norm), method="kendall")

target = lapply(colnames(ma)[2:ncol(ma)], function(x){
  ma$gene[ma[,x]==1]
})
names(target) = colnames(ma)[2:ncol(ma)]

Zmi = mirna_sign$stat ; names(Zmi) = row.names(mirna_sign)
Zg = mrna_sign$stat; names(Zg) = row.names(mrna_sign)

output = mirPath(DEmirs = mirna_keep, Zmi = Zmi, Zg = Zg,
                 cor = cor,
                 targets = target, pathways = kegg2ens,
                 fdr="bonferroni")
rownames(output) = output$path
output$true = 0
output[names(others),"true"] = others
output$score = output$genes/output$weight
resEnrichment = output

# network
library(igraph)
g = graph.data.frame(pairs[,c("V2", "gene")] %>% filter(gene %in% mrna_keep), directed = F)

# plot.igraph(g, vertex.color=V(g)$color, vertex.label=NA)
ebc <- edge.betweenness.community(g, directed=T)
dist = table(ebc$membership)
big = names(dist[dist>3])
del = ebc$names[!(ebc$membership %in% big)]
g = delete.vertices(g, del)
V(g)$color[V(g)$name %in%  mirna_keep]<-"grey80"
V(g)$color[V(g)$name %in%  mrna_keep]<-"orange"
g
# plot.igraph(g, vertex.color=V(g)$color, vertex.label=NA)

# custom
sum(colSums(cor[mirna_keep,mrna_keep]< -.7) > 1)

#' go through all kegg names
enrich = lapply(names(kegg2ens)[3], function(x){
   print(x)
  .genes = kegg2ens[[x]]
  if (length(.genes)<3)
    return(NULL)
  .is_regulator = rowSums(cor[mirna_keep,.genes]< -.7) > 1
  .is_target = colSums(cor[mirna_keep,.genes]< -.7) > 1
  .ts = ma[intersect(.genes, rownames(ma)), ]
  if (nrow(.ts)==0)
    return(NULL)
  .ts_target = names(rowSums(.ts[,intersect(mirna_keep, colnames(.ts))] > 0))
  .total_target = intersect(.genes[.is_target], .ts_target)
  # .total_mirs = rownames(cor)[.is_regulator]
  # .de_reg_mirna = intersect(.total_mir, mirna_keep)
  .de_target_gene = intersect(.total_target, mrna_keep)
  .de_here = intersect(.genes, mrna_keep)
  phyper(length(.de_target_gene),
         length(.total_target),
         length(.genes) - length(.total_target),
         length(.de_here)
  )
})


#
library(STATegRa)

mirna_b1 <- createOmicsExpressionSet(Data=as.matrix(mirna_ma[,row.names(mirna_de)]),
                               pData=mirna_de[,"group",drop=FALSE],
                               pDataDescr=c("classname"))

mrna_b2 <- createOmicsExpressionSet(Data=as.matrix(mrna_ma[,row.names(mrna_de)]),
                                     pData=mrna_de[,"group",drop=FALSE],
                                     pDataDescr=c("classname"))


cc <- selectCommonComps(X=mirna_ma, Y=mrna_ma, Rmax=3)
cc$common