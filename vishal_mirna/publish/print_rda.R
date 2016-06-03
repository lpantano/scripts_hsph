require(cowplot)

gene = "ENSMUSG00000000628"
mirna = "mmu-let-7g-5p"

.summarize = function(exp, group){
    v = sapply(levels(group), function(i){
        idx = which(group == i)
        median(exp[idx], na.rm=TRUE)
    })
    v[!is.na(v)]
}

# fa_g = .summarize(assay(obj$fa$gene, "norm")[gene,], colData(obj$fa$gene)$group)
# fa_p = .summarize(assay(obj$fa$prot, "norm")[gene,], colData(obj$fa$prot)$group)
# fa_m = .summarize(assay(obj$fa$mirna, "norm")[mirna,], colData(obj$fa$mirna)$group)
fa_g = assay(obj$fa$gene, "norm")[gene,]
fa_p = assay(obj$fa$prot, "norm")[gene,]
fa_m = assay(obj$fa$mirna, "norm")[mirna,]

# uuo_g = .summarize(assay(obj$uuo$gene, "norm")[gene,], colData(obj$uuo$gene)$group)
# uuo_p = .summarize(assay(obj$uuo$prot, "norm")[gene,], colData(obj$uuo$prot)$group)
# uuo_m = .summarize(assay(obj$uuo$mirna, "norm")[mirna,], colData(obj$uuo$mirna)$group)
uuo_g = assay(obj$uuo$gene, "norm")[gene,]
uuo_p = assay(obj$uuo$prot, "norm")[gene,]
uuo_m = assay(obj$uuo$mirna, "norm")[mirna,]

fa_df = rbind( data.frame( exp = as.vector(scale(fa_g)), group=colData(obj$fa$gene)$group,  molecule = "transcript", model="FA" ),
       data.frame( exp = as.vector(scale(fa_p)), group=colData(obj$fa$prot)$group, molecule = "protein", model="FA" ),
       data.frame( exp = as.vector(scale(fa_m)), group=colData(obj$fa$mirna)$group, molecule = "miRNA", model="FA" )
)
fa_df$group = factor(fa_df$group, levels=levels(colData(obj$fa$mirna)[,"group"]))

uuo_df = rbind(data.frame( exp = as.vector(scale(uuo_g)), group=colData(obj$uuo$gene)$group, molecule = "transcript", model="UUO" ),
               data.frame( exp = as.vector(scale(uuo_p)), group=colData(obj$uuo$prot)$group, molecule = "protein", model="UUO" ),
               data.frame( exp = as.vector(scale(uuo_m)), group=colData(obj$uuo$mirna)$group, molecule = "miRNA", model="UUO" ))
uuo_df$group = factor(uuo_df$group, levels=levels(colData(obj$uuo$mirna)[,"group"]))

fa_plot = ggplot(fa_df, aes(y=exp, x=group, group=molecule, colour=molecule, fill=molecule)) +
    geom_point() + geom_smooth(method = "lm",formula = y~poly(x,3), alpha=0.2) +  ggtitle("FA-model")
# print(fa_plot)
uuo_plot = ggplot(uuo_df, aes(y=exp, x=group, group=molecule, colour=molecule, fill=molecule)) +
    geom_point() + geom_smooth(method = "lm",formula = y~poly(x,2), alpha=0.2) + ggtitle("UUO-model")
# print(uuo_plot)

grid.arrange(fa_plot, uuo_plot)
