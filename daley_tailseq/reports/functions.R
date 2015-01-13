read_summary_d = function(dd){
    dd = dd %>% filter(V3 > 1 & !grepl("N", V2) & !is.na(V2))
    dd$type = as.character(dd$V2)
    dd$type[ dd$type!="counts" & dd$type!="polyA" ] = "Mod"
    dd_melt = dd %>% group_by(type) %>% summarise(genes=length(unique(V1)))
    p = ggplot(dd_melt, aes(y=genes, x=type)) + 
        geom_bar( stat = 'identity') +
        ggtitle("Number of genes")
    print(p)
    dd_melt = dd %>% group_by(V1,type) %>% summarise(total=sum(V3))
    p <- ggplot(as.data.frame(dd_melt) , aes(x=type, y=total, fill=type))+
        geom_boxplot() + 
        scale_fill_brewer( palette = "Set1" ) + 
        scale_y_log10() +
        ggtitle("Average expression")
    print(p)
    dd_melt = dd %>% filter(V2=="AA" | V2=="A" | V2=="C") %>% 
        dplyr::select(V1,type=V2) %>% group_by(type) %>% 
        summarise(genes=length(unique(V1)))
    p = ggplot(dd_melt, aes(y=genes, x=type)) + 
        geom_bar( stat = 'identity', positions="dodge" ) +
        ggtitle("Number of genes")
    print(p)
}

summary_all = function(dd_list){
    all = lapply(dd_list, get_mod_table)
    names(all) = names(dd_list)
    tab = data.frame()
    for (e in names(all)) {
        tab = rbind( tab, cbind( all[[e]], samples=rep(e,nrow(all[[e]])) ) )
    }
    dd_melt = tab %>% group_by(type,samples) %>% summarise(total=sum(counts)) %>% ungroup()
    dd_melt$type = factor(dd_melt$type, levels=rev(c("U>3","UU","U","polyA","exp")))
    p = ggplot(dd_melt, aes(y=total, x=samples, fill=type, order=type )) + 
        geom_bar( stat = 'identity', position = 'dodge') +
        scale_y_log10()+
        scale_fill_brewer(palette = "Set1")+
    ggtitle("Number of reads: log(10)")
    print(p)
    p = ggplot(dd_melt, aes(y=total, x=samples, fill=type, order=type )) + 
        geom_bar( stat = 'identity') +
        scale_fill_brewer(palette = "Set1")+
    ggtitle("Number of reads: ")
    print(p)
    
}

frequency_u = function(dd_list){
    all = lapply(names(dd_list) , function (name){
        #x = normalization(dd_list[[name]])
        x = (dd_list[[name]])
        rbind(get_size(x,"<15") %>% mutate(sample=name),
              get_size(x,"<25") %>% mutate(sample=name),
              get_size(x,">25") %>% mutate(sample=name) )
        
    })
    tab = do.call(rbind,all)
    tab$type = factor(tab$type, levels=rev(c("U>3","UU","U","polyA","polyA-NoU","exp")))
    ggplot(tab %>% filter(type!="exp" & type!="polyA" & type!="polyA-NoU"), aes(x=size, y=freq, fill = type, order=type )) +
        geom_bar(stat = 'identity')+
        scale_fill_brewer(palette = "Blues")+
        facet_wrap(~sample) +
        theme_bw()
    
    ggplot(tab %>% filter(type!="exp" & type!="polyA"), aes(x=size, y=freq, fill = type, order=type )) +
        geom_bar(stat = 'identity')+
        scale_fill_brewer(palette = "Blues")+
        facet_wrap(~sample) +
        theme_bw()
    
}

get_size = function(x, size){
    t = get_mod_table(x %>% filter(V2==size | V3==size | V2=="total")) %>% mutate(size=size)    
    total = t %>% group_by(type,size) %>% summarize(total=sum(counts))
    total$freq = total$total / total$total[total$type=="polyA"] * 100
    total
}

get_size_gene = function(x, size=NULL, limit=50){
    if (is.null(size)){
        t = get_u_table(x,limit=limit)    
    }else{
        t = get_u_table(x %>% filter(V2==size | V3==size | V2=="total" ),limit=limit) %>% mutate(size=size)    
    }
    t
}

normalization = function(x,L=100000){
    t = x %>% filter(V2=="total")
    th = quantile(t$V4,c(.25,.75))
    ls = sum((t %>% filter(V4>th[1] & V4<th[2]))[,4])
    x$V4=x$V4/ls * L
    x
}

deseq_factors = function(x, d){
    x = x[rowMeans(x)>5,]
    dse = DESeqDataSetFromMatrix(countData = x, colData = d,  formula("~ condition") )   
    dse = DESeq::estimateSizeFactors(dse)
    sizeFactors(dse)
}

deseq_norm = function(x, d, size){
    dse = DESeqDataSetFromMatrix(countData = x, colData = d,  formula("~ condition") )
    sizeFactors(dse) = size
    dse
}

get_set_genes = function(dd){
    poly = unique((dd %>% filter(V2=="polyA"))[,1])
    aa = as.character(unique((dd %>% filter(V3=="AA"))[,1]))
    a = as.character(unique((dd %>% filter(V3=="A"))[,1]))
    a3 = as.character(unique((dd %>% filter(V5 > 0.85 & nchar(as.character(V3))>3 & V2!="total" & V2!="polyA"))[,1]))
    genes = unique(c(a,aa,a3))
    #c = unique((dd %>% filter(V2=="C"))[,1])
    #mod = unique((dd %>% filter(V2!="counts" & V2!="polyA"))[,1])
    n10 = round(length( c(aa,a,a3) ) * 0.1)
    n5 = round(length( c(aa,a,a3) ) * 0.05)
    top10 = unique(( dd[ dd$V1 %in% as.character(genes), ] %>%
                         arrange( desc(V4) ))[1:n10,1])
    top5 = unique(( dd[ dd$V1 %in% as.character(genes), ] %>%
                        arrange( desc(V4) ))[1:n5,1])
    return( list("polyA"=poly, "a"=a, "aa"=aa, "a3"=a3, "top10"=top10, "top5"=top5) )
}

get_set_genes_public = function(dd){
    poly = unique((dd %>% filter(V2=="polyA"))[,1])
    aa = as.character(unique((dd %>% filter(V2=="UU"))[,1]))
    a = as.character(unique((dd %>% filter(V2=="U"))[,1]))
    a3 = as.character(unique((dd %>% filter(V2=="U3"))[,1]))
    genes = unique(c(a,aa,a3))
    #c = unique((dd %>% filter(V2=="C"))[,1])
    #mod = unique((dd %>% filter(V2!="counts" & V2!="polyA"))[,1])
    n10 = round(length( c(aa,a,a3) ) * 0.1)
    n5 = round(length( c(aa,a,a3) ) * 0.05)
    top10 = unique(( dd[ dd$V1 %in% as.character(genes), ] %>%
                         arrange( desc(V3) ))[1:n10,1])
    top5 = unique(( dd[ dd$V1 %in% as.character(genes), ] %>%
                        arrange( desc(V3) ))[1:n5,1])
    return( list("polyA"=poly, "a"=a, "aa"=aa, "a3"=a3, "top10"=top10, "top5"=top5) )
    
}

filter_out = function(dd, cuta=30, cutmod=1){
    dd = dd %>% filter(V4 > cutmod & !grepl("N", V3) & !is.na(V3))
    genes = (dd %>% filter(V2=="polyA") %>% group_by(V1) %>% summarise(counts=sum(V4)) %>% ungroup() %>% filter(counts>=cuta))
    dd[ dd$V1 %in%  as.character(genes$V1), ]
}

filter_out_public = function(dd, cuta=30, cutmod=1){
    dd = dd %>% filter(V3 > cutmod & !grepl("N", V3) & !is.na(V3))
    genes = (dd %>% filter(V2=="polyA") %>% group_by(V1) %>% summarise(counts=sum(V3)) %>% ungroup() %>% filter(counts>=cuta))
    dd[ dd$V1 %in%  as.character(genes$V1), ]
}

get_u_table = function(dd,by_size=NULL,limit=50){
    dd$type = as.character(dd$V2)
    uf = get_mod_table(dd)
    total = uf %>% filter(type!="polyA" & type!="exp" & type!="polyA-NoU") %>% group_by(V1) %>% summarize(total=sum(counts))
    poly = uf %>% filter(type=="polyA")
    t = merge(total,poly, by=1)
    t = t %>% filter(total >= limit)
    t$ufreq = t$total/t$counts * 100
    t[,c(1,5)]
}

get_mod_table = function(dd,by_size=NULL){
    dd$type = as.character(dd$V2)
    exp = dd %>% filter(V2=="total") %>% 
        group_by(V1) %>% summarise(counts=sum(V4)) %>% ungroup() %>% mutate(type="exp")
    poly = dd %>% filter(V2=="polyA") %>% 
        group_by(V1) %>% summarise(counts=sum(V4)) %>% ungroup() %>% mutate(type="polyA")
    #    t = merge(exp,poly, by=1,all.x=TRUE)
    #    t[is.na(t)] = 0
    #    t$counts = t[,2] - t[,4]
    #    exp = t[,c("V1","counts","type.x")]
    #    names(exp)[3] = "type"
    
    #    mod = dd %>% filter(V2!="polyA" & V2!="total") %>% 
    #                group_by(V1) %>% summarise(counts=sum(V4)) %>% ungroup() %>% mutate(type="mod")
    #    t = merge(poly,mod,by=1,all.x=TRUE)
    #    t[is.na(t)] = 0
    #    t$counts = t[,2] - t[,4]
    #    poly = t[,c("V1","counts","type.x")]
    #    names(poly)[3] = "type"
    
    aa = dd %>% filter(V3=="AA") %>% 
        group_by(V1) %>% summarise(counts=sum(V4)) %>% ungroup() %>% mutate(type="U")
    a = dd %>% filter(V3=="A") %>% 
        group_by(V1) %>% summarise(counts=sum(V4)) %>% ungroup() %>% mutate(type="UU")
    a3 = dd %>% filter(V5 > 0.85 & nchar(as.character(V3))>3 & V2!="total" & V2!="polyA") %>% 
        group_by(V1) %>% summarise(counts=sum(V4)) %>% ungroup() %>% mutate(type="U>3")
    nou = dd %>% filter(V5 < 0.11 & V2!="total" & V2!="polyA") %>% 
        group_by(V1) %>% summarise(counts=sum(V4)) %>% ungroup() %>% mutate(type="polyA-NoU")
    
    
    dd = rbind(exp,poly,a,aa,a3,nou)
    dd[dd$counts>0,]
}

poly_read_summary = function(dd){
    dd = filter_out(dd)
    dd_melt = get_mod_table(dd)
    p = ggplot(dd_melt, aes(y=counts, x=type)) + 
        geom_boxplot() +
        scale_y_log10()+
        ggtitle("Number of reads")
    print(p)
    
    dd_melt = dd_melt %>% group_by(type) %>% summarise(total=sum(counts)) %>% ungroup()
    p = ggplot(dd_melt, aes(y=total, x=type)) + 
        geom_bar( stat = 'identity') +
        geom_text(aes(label=total), position=position_dodge(width=0.9), size=6, color="orange") +
        scale_y_log10()+
        ggtitle("Number of reads")
    print(p)
    #dd_melt = dd %>% filter(V2=="AA" | V2=="A" | V2=="C") %>% 
    #          dplyr::select(mod=V2, V3) %>% group_by(mod) %>%
    #          summarise( total=sum(V3) )
    #p = ggplot(dd_melt, aes(y=total, x=mod)) + 
    #    geom_bar(stat="identity", position = "dodge") +
    #    ggtitle("Number of reads")
    #print(p)
}

rm_na = function(x){
    x[!is.na(x)]
}

get_counts = function(dd_list){
    all = lapply(names(dd_list) , function (name){
        x = dd_list[[name]]
        t = get_mod_table(x) %>% mutate(sample=name)
    })
    tab = do.call(rbind,all)
    e= reshape(as.data.frame(tab %>% filter(type=="exp") %>% dplyr::select(V1,counts,sample)), direction = "wide", idvar = c("V1"),  timevar = "sample")
    e[is.na(e)] = 0
    u = reshape(as.data.frame(tab %>% filter(type=="U") %>% dplyr::select(V1,counts,sample)), direction = "wide", idvar = c("V1"),  timevar = "sample")
    u[is.na(u)] = 0
    uu= reshape(as.data.frame(tab %>% filter(type=="UU") %>% dplyr::select(V1,counts,sample)), direction = "wide", idvar = c("V1"),  timevar = "sample")
    uu[is.na(uu)] = 0
    u3= reshape(as.data.frame(tab %>% filter(type=="U>3") %>% dplyr::select(V1,counts,sample)), direction = "wide", idvar = c("V1"),  timevar = "sample")
    u3[is.na(u3)] = 0
    p= reshape(as.data.frame(tab %>% filter(type=="polyA") %>% dplyr::select(V1,counts,sample)), direction = "wide", idvar = c("V1"),  timevar = "sample")
    p[is.na(p)] = 0
    list(exp=e,polya=p,u=u,uu=u,u3=u3)
}

get_norm_counts = function(counts, design){
    e = counts[[1]]
    library_size = deseq_factors(e[,2:ncol(e)],design)
    norm_list = lapply(counts, function(x){
        dse = deseq_norm(x[,2:ncol(x)],design,library_size)
        counts(dse, normalized=TRUE)
    })
    names(norm_list) = names(counts)
    norm_list
}

get_DEG = function(counts, design){
    e = counts[[1]]
    library_size = deseq_factors(e[,2:ncol(e)],design)
    deg_list = lapply(counts, function(x){
        row.names(x) = x[,1]
        x = x[,2:ncol(x)]
        dse = deseq_norm(x,design,library_size)
        dse = estimateDispersions(dse)
        dse = nbinomWaldTest(dse)
    })
    names(deg_list) = names(counts)
    deg_list
}


library(org.Hs.eg.db)
cols <- c("ENSEMBL", "SYMBOL")
make_table = function(dd, prefix){
    all = lapply(c("<15","<25",">25"), function(size){
        t = get_mod_table(dd %>% filter(V2==size | V3==size | V2=="total")) %>% mutate(size=size)  
        t$size= sub(size,paste0("-pA",size),t$size)
        t
    }) 
    data = do.call(rbind,all)
    raw = paste0(prefix,"_raw.txt")
    wide = paste0(prefix,"_wide.txt")
    write.table(data, raw,row.names=F, col.names=T)
    len = sapply(as.character(data$V2), nchar)
    data$type = paste0(data$type,data$size)
    data = reshape(as.data.frame(data %>% dplyr::select(V1, counts,type)),  timevar = "type", idvar = "V1", direction = "wide")
    names(data) = sub("counts.","",names(data))
    data[is.na(data)] = 0
    data = data[,!grepl("exp<",names(data))]
    names(data) = sub("exp-pA>25","exp",names(data))
    symbol = AnnotationDbi::select(org.Hs.eg.db, as.character(data$V1), cols, keytype="ENSEMBL")
    symbol = symbol %>% distinct(ENSEMBL)
    data$symbol = symbol[match(data$V1,symbol$ENSEMBL),2]
    write.table(data, wide,row.names=F, col.names=T)
    c(raw,wide)
}


plot_family = function(dd){
    top = head(dd %>% filter(V2=="total") %>% arrange(desc(V4)),1000 ) 
    top$V4 = top$V4 / sum(dd %>% filter(V2=="total") %>% dplyr::select(V4) %>% ungroup()) * 100
   
    mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
    g <- getBM( attributes=c("ensembl_gene_id", "hgnc_symbol","gene_biotype","interpro_short_description") , filters=
                    "ensembl_gene_id"    , values =as.character(top$V1) ,mart=mart)
    g$family = 
        sapply(g$interpro_short_description,function(x){
            strsplit(x,"_")[[1]][1]
        })
    g = g %>% distinct(ensembl_gene_id)
 
    dd = merge(top,g,by=1) %>% arrange(desc(V4))
    dd$hgnc_symbol = factor(dd$hgnc_symbol,levels=dd$hgnc_symbol)
    dd = dd %>% distinct(V1) %>% group_by(family) %>% summarise(counts=sum(V4), genes = n())
    dd$family[dd$counts<0.5] = "Others"
    p = ggplot(dd , aes( x=family,y=(counts),fill=family )) +
        geom_bar(stat = 'identity') +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(p)
}


get_biotype = function(dd)
{
    require(biomaRt)
    mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
    g <- getBM( attributes=c("ensembl_gene_id","gene_biotype") , filters=
                    "ensembl_gene_id"    , values =as.character(dd$V1) ,mart=mart)
    dd$biotype = g[match(dd$V1, g[,1]),2]
    dd
}

get_description = function(dd)
{
    if (nrow(dd)==0) return(dd)
    require(biomaRt)
    mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
    g <- getBM( attributes=c("ensembl_gene_id","external_gene_name","gene_biotype","description") , filters=
                    "ensembl_gene_id"    , values =as.character(dd$gene) ,mart=mart)
    dd= cbind(dd,g[match(dd$gene, g[,1]),2:4])
    dd
}


get_size_gene = function(dd,ref_fn){
    ref=read.table(ref_fn,sep="\t") %>% filter(grepl("transcript",V8)) %>% mutate(size=V3-V2+1)
    dd$size = ref[match(dd$V1, ref[,4]),"size"]
}

