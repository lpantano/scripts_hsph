plot_cluster  = function(norm_sign, g_in_c, xs ,groups, title) {
    # g_in_c = names(groups[groups==1])
    ma = as.data.frame(norm_sign)[g_in_c,]
    ma_long = suppressMessages(melt(cbind(gene=row.names(ma), ma), variable_name = "sample"))
    ma_long$x = xs[ma_long$sample]
    ma_long$group = groups[ma_long$sample]
    ggplot(ma_long, aes(x=x, y=value, fill=group, color=group)) + 
        geom_violin(alpha=0.3) + 
        stat_smooth(aes(x=x, y=value, group=group, color=group),method = "lm",formula = y~poly(x,2)) +
        ggtitle(paste("Group:", title, "(", length(g_in_c), " genes )")) +
        theme_bw(base_size = 11) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


make_patterns = function(ma, metadata, minc=15, summarize="group", time="time", col="condition"){
    
    counts_group = t(sapply(rownames(ma), function(g){
        sapply(unique(metadata[,summarize]), function(i){
            idx = which(metadata[,summarize] == i)
            mean(ma[g,idx], na.rm=TRUE)
        })
    }))
    
    m = (1-cor(t(counts_group), method = "kendall"))
    m[m<0] = 0
    d = as.dist(m^2)
    c = diana(d, diss = TRUE, stand = FALSE)
    
    groups = cutree(as.hclust(c), h = c$dc)
    
    norm_sign = t(apply(counts_group, 1, function(e){
        (e - min(e))/(max(e) - min(e))
    }))
    metadata_groups = metadata %>% distinct(metadata[,summarize])
    rownames(metadata_groups) = metadata_groups[,summarize]

    to_plot = names(table(groups))[table(groups) > minc]
    plots = lapply(to_plot, function(x){
        plot_cluster(norm_sign, as.character(names(groups[groups==x])), 
                     metadata_groups[,time], metadata_groups[,col], x)
    })
    
    do.call(grid.arrange, plots)
}
