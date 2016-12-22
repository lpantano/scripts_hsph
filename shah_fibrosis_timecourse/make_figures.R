library(CHBUtils)
library(ggplot2)
library(dplyr)

late_path = "~/orch/scratch/shah-rnaseq2/late/final/2016-04-11_late/files"
early_path = "~/orch/scratch/jshah_fibrosis_rnaseq/time_course_1/final/2016-03-28_time_course_1/de_early/"

late_exp = read.table(file.path(late_path, "mice_model_jck_wt_log2_counts.csv"), sep=",", header=T, row.names=1)
early_exp = read.table(file.path(early_path, "mice_model_jck_wt_log2_counts.csv"), sep=",", header=T, row.names=1)

late_meta = read.csv(file.path(late_path, "..", "project-summary.csv"))
early_meta = read.csv(file.path(early_path, "..", "project-summary.csv"))

early_de = read.csv(file.path(early_path, "mice_model_jck_wt_de.csv"), row.names=1)
late_de = read.csv(file.path(late_path, "mice_model_jck_wt_de.csv"), row.names=1)

names = paste(gsub("-",".",as.character(early_meta$Name)), early_meta$condition, early_meta$time, sep="_")
ma = early_exp[,names]
colnames(ma) = paste(early_meta$condition, early_meta$time, 1:ncol(ma), sep="_")
early_meta$time = relevel(early_meta$time, "P5")
pdf("figures/batch1_de_rnaseq.pdf")
df = mds(ma, plot=FALSE)
df$table %>% mutate(genotype=ifelse(grepl("WT", label), "WT", "JCK")) %>%
    mutate(time = unlist(early_meta$time)) %>%
    ggplot(aes(x=one, y=two, color=time, shape=genotype)) +
    geom_point(size=4) + theme_bw() + xlab(df$labels[1]) + ylab(df$labels[2])
dev.off()

names = paste0("X", as.character(late_meta$Name), "__")
ma = late_exp[,names]
colnames(ma) = paste(late_meta$genotype, late_meta$age, 1:ncol(ma), sep="_")
pdf("figures/batch2_de_rnaseq.pdf")
df = mds(ma, plot=FALSE)
df$table %>% mutate(genotype=ifelse(grepl("WT", label), "WT", "JCK")) %>%
    mutate(time = unlist(early_meta$time)) %>%
    ggplot(aes(x=one, y=two, color=time, shape=genotype)) +
    geom_point(size=4) + theme_bw() + xlab(df$labels[1]) + ylab(df$labels[2])

dev.off()

### custom clusters

prot <- read.csv("clusters/ProtListOut.csv")
ma_long1 = reshape::melt.data.frame(prot[prot$flipped==1,c(4:13,16)])
ma_long0 = reshape::melt.data.frame(prot[prot$flipped==0,c(4:13,16)])
.clean = . %>% mutate(time=gsub("JCK", "", gsub("WTP", "", variable)))

levels = (unique((ma_long1  %>% .clean())$time))[]

pdf("figures/proteins_cluster_flipped1.pdf", width = 9)
ggplot(ma_long1 %>% .clean(), aes(x=factor(time, levels = levels), y=value, color=cluster)) +
    geom_boxplot(alpha=0.3, outlier.size = 0, outlier.shape = NA) +
    stat_smooth(aes(x=factor(time, levels = levels), y=value, group=cluster),method = "lm",formula = y~poly(x,3)) +
    theme_bw(base_size = 11) +
    ylim(-1,2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylab("log2(JCK/WTP)") + xlab("Time") +
    scale_color_brewer(palette = "Set1")
dev.off()

pdf("figures/proteins_cluster_flipped0.pdf", width = 9)
ggplot(ma_long0 %>% .clean(), aes(x=factor(time, levels = levels), y=value, color=cluster)) +
    geom_boxplot(alpha=0.3, outlier.size = 0, outlier.shape = NA) +
    stat_smooth(aes(x=factor(time, levels = levels), y=value, group=cluster),method = "lm",formula = y~poly(x,3)) +
    theme_bw(base_size = 11) +
    ylim(-1,2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylab("log2(JCK/WTP)") + xlab("Time") +
    scale_color_brewer(palette = "Set1")
dev.off()

pprot <- read.csv("clusters/PhosListOut.csv")
map_long1 = reshape::melt.data.frame(pprot[pprot$flipped==1,c(6:15,23)])
map_long0 = reshape::melt.data.frame(pprot[pprot$flipped==0,c(6:15,23)])

pdf("figures/phospho_cluster_flipped1.pdf", width = 9)
ggplot(map_long1 %>% .clean(), aes(x=factor(time, levels = levels), y=value, color=cluster)) +
    geom_boxplot(alpha=0.3, outlier.size = 0, outlier.shape = NA) +
    stat_smooth(aes(x=factor(time, levels = levels), y=value, group=cluster),method = "lm",formula = y~poly(x,3)) +
    theme_bw(base_size = 11) +
    ylim(-1,2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylab("log2(JCK/WTP)") + xlab("Time") +
    scale_color_brewer(palette = "Set1")
dev.off()

pdf("figures/phospho_cluster_flipped0.pdf", width = 9)
ggplot(map_long0 %>% .clean(), aes(x=factor(time, levels = levels), y=value, color=cluster)) +
    geom_boxplot(alpha=0.3, outlier.size = 0, outlier.shape = NA) +
    stat_smooth(aes(x=factor(time, levels = levels), y=value, group=cluster),method = "lm",formula = y~poly(x,3)) +
    theme_bw(base_size = 11) +
    ylim(-1,2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylab("log2(JCK/WTP)") + xlab("Time") +
    scale_color_brewer(palette = "Set1")
dev.off()

