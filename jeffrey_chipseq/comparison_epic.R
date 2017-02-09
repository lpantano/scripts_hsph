library(knitr)
library(ChIPseeker)
library(dplyr)
library(ggplot2)
library(UpSetR)
library(GenomicRanges)
library(readr)
library(rio)

# set to source
epic_0h = read.table("filtered/SP140_0h_peaks.tsv", sep = " ", header = T) %>%
    filter(FDR < 0.01)
epic_4h = read.table("filtered/SP140_4h_peaks.tsv", sep = " ", header = T) %>%
    filter(FDR < 0.01)
mac_0h = import("filtered/SP140_0h.txt")
mac_4h = import("filtered/SP140_4h.txt")

gepic0h = makeGRangesFromDataFrame(epic_0h, keep.extra.columns = TRUE,
                                   ignore.strand = TRUE, end.field = "End",
                                   start.field = "Start", seqnames.field = "Chromosome")
gepic4h = makeGRangesFromDataFrame(epic_4h, keep.extra.columns = TRUE,
                                   ignore.strand = TRUE, end.field = "End",
                                   start.field = "Start", seqnames.field = "Chromosome")
gmac0h = makeGRangesFromDataFrame(mac_0h %>%
                                      select(seqnames,start,end,foldchange,score,log10qvalue),
                                  keep.extra.columns = TRUE,
                                  ignore.strand = TRUE, end.field = "end",
                                  start.field = "start", seqnames.field = "seqnames")
gmac4h = makeGRangesFromDataFrame(mac_4h %>%
                                      select(seqnames,start,end,foldchange,score,log10qvalue),
                                  keep.extra.columns = TRUE,
                                  ignore.strand = TRUE, end.field = "end",
                                  start.field = "start", seqnames.field = "seqnames")

sum(countOverlaps(gmac0h, gepic0h))/length(gmac0h)
sum(countOverlaps(gmac4h, gepic4h))/length(gmac4h)

sum(width(intersect(gmac0h, gepic0h)))/sum(width(gmac0h))
sum(width(intersect(gmac4h, gepic4h)))/sum(width(gmac4h))

hist(log10(width(gepic0h)))
hist(log10(width(gmac0h)))

hist(log10(width(gepic4h)))
hist(log10(width(gmac4h)))
summary(width(gmac4h))
summary(width(gepic4h))

#  correlation of FC size - NO correlation
cor.test(log10(width(gepic0h)), log10(mcols(gepic0h)$Fold_change))
plot(log10(width(gepic0h)), log10(mcols(gepic0h)$Fold_change))
cor.test(log10(width(gmac0h)), log10(mcols(gmac0h)$foldchange))
plot(log10(width(gmac0h)), log10(mcols(gmac0h)$foldchange))

# 0h
idx0h = findOverlaps(gmac0h, gepic0h)
j0h=data.frame(im=queryHits(idx0h), sm=mcols(gmac0h[queryHits(idx0h),])$score,
               wm=width(gmac0h[queryHits(idx0h)]),
               ie=subjectHits(idx0h), se=mcols(gepic0h[subjectHits(idx0h),])$Score,
               we=width(gepic0h[subjectHits(idx0h)])) %>%
    group_by(im) %>% summarise(ie=max(ie), sm=mean(sm), se=mean(se),wm=max(wm),we=sum(we)) %>%
    group_by(ie) %>% summarise(im=max(im),sm=mean(sm), se=mean(se),wm=sum(wm),we=max(we))

ggplot(j0h,aes(x=wm,y=we)) + geom_point(alpha=0.1) +
    scale_x_log10() + scale_y_log10() +
    theme_bw() + ggtitle("SP140_0h epic peak size vs macs2 peak size") +xlab("macs") + ylab("epic")

cor.test(log10(j0h$wm),log10(j0h$we))

ggplot(j0h,aes(x=sm,y=se,color=log10(we))) + geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10()
ggplot(j0h,aes(x=sm*wm/1000,y=se,color=log10(we))) + geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10()
cor.test(j0h$sm, j0h$se)

# 4h
idx4h = findOverlaps(gmac4h, gepic4h)
j4h=data.frame(im=queryHits(idx4h), sm=mcols(gmac4h[queryHits(idx4h),])$score,
               wm=width(gmac4h[queryHits(idx4h)]),
               ie=subjectHits(idx4h), se=mcols(gepic4h[subjectHits(idx4h),])$Score,
               we=width(gepic4h[subjectHits(idx4h)])) %>%
    group_by(im) %>% summarise(ie=max(ie), sm=mean(sm), se=mean(se),wm=max(wm),we=sum(we)) %>%
    group_by(ie) %>% summarise(im=max(im), sm=mean(sm), se=mean(se),wm=sum(wm),we=max(we))

ggplot(j4h,aes(x=wm,y=we)) + geom_point(alpha=0.1) +
    scale_x_log10() + scale_y_log10() +
    theme_bw() + ggtitle("SP140_4h epic peak size vs macs2 peak size") +xlab("macs") + ylab("epic")

cor.test(log10(j4h$wm),log10(j4h$we))

ggplot(j4h,aes(x=sm,y=se)) + geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10()
ggplot(j4h,aes(x=sm*wm/1000,y=se,color=log10(we))) + geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10()
ggplot(j4h,aes(x=sm,y=se,color=log10(we))) + geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10()
