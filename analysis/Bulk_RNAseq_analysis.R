## install.packages("BiocManager")
## BiocManager::install(version = "3.14")
## BiocManager::install("edgeR")
## install.packages("rlang")
## updateR()

library(biomaRt)
library(data.table)
library(limma)
library(edgeR)
library(RColorBrewer)
library(dplyr)

setwd("~/Dropbox/Luo_lab/acj6/RNAseq/input/RNA-seq_files/")
path <- "~/Dropbox/Luo_lab/acj6/RNAseq/input/RNA-seq_files/"
files <- list.files(path=path)

file.renames <- gsub("*.tab", "", files)

for (i in 1:length(files)) {
  # read data file
  df <- read.csv(files[i], header=F, sep ="")
  # find drosophila genes in dro
  df_trim <- df[df$V1 %like% "FBgn",]
  write.table(df_trim, file=paste("~/Dropbox/Luo_lab/acj6/RNAseq/input/trimmed/", file.renames[i], sep=""), row.names=F, sep="\t")
}

files_trim = list.files("~/Dropbox/Luo_lab/acj6/RNAseq/input/trimmed/")
setwd("~/Dropbox/Luo_lab/acj6/RNAseq/trimmed/")
dge = readDGE(files_trim, columns=c(1,2), header=T)
dge$samples

dge$samples$group <- factor(c(rep("Ctrl", 3),
                              rep("Acj66", 3)), levels = c("Ctrl", "Acj66"))

cpm <- cpm(dge)
# creates density plots using cpm
log2_cpm <- cpm(dge, log=TRUE)

density(log2_cpm[,1])
plot(density(log2_cpm[,1]))

# filter out genes that have cpm <=1. That is roughly 5-10 reads in each sample given the lib size.
# To achieve this, we first create index matrix (logical) that have at least three samples have cpm>1
#cpm>1 gives a logical matrix, then rowSums>=3 gives logical rows,
dge_sel_idx <- rowSums(cpm>1)>=3 
## dge_sel_idx
# select rows that have more than 3 cpm>1 [TRUE]
dge_sub <- dge[dge_sel_idx, ]

#normalize data with TMM method
dge_sub_norm <- calcNormFactors(dge_sub, method = "TMM")
dge_sub_norm$samples$norm.factors

save(dge_sub_norm, file="~/Dropbox/Luo_lab/acj6/RNAseq/output/dge_all_log2cpm_norm.RData")
#make design matrix that incorporate the grouping info for later comparision
grps <- dge_sub$samples$group
design <- model.matrix(~0 + grps)
colnames(design) <- gsub("grps", "", colnames(design))

dge_sub_norm$samples
contrast <- makeContrasts(
  Acj66_vs_Ctrl = Acj66 - Ctrl,
  levels = colnames(design))

contrast

# Transform RNA-Seq Data Ready for Linear Modelling
dge_voom <- voom(dge_sub_norm, design)
# Add a linear fit 

dge_fit <- lmFit(dge_voom, design)
grp_fit <- contrasts.fit(dge_fit, contrasts=contrast)
efit <- eBayes(grp_fit)

save(efit, file="~/Dropbox/Luo_lab/acj6/RNAseq/output/efit_all.RData")


## MAKE PCA PLOT
library(ggplot2)
## library(svglite)

cpm_sub_norm <- cpm(dge_sub_norm)
log2_cpm_sub_norm <- cpm (dge_sub_norm, log=TRUE)

pca <- prcomp(t(log2_cpm_sub_norm), scale=TRUE) 
class(pca)
head(pca)
## plot pc1 and pc2
plot(pca$x[,1], pca$x[,2])

## make a plot showing all varience of all PC components
pca.var <- pca$sdev^2
## pca.var
pca.var.per <- round(pca.var/sum(pca.var)*100, 2)
pc_label <- paste("PC", c(1:8), sep="")
pca_all <- data.frame(PC=pc_label, Variance=pca.var.per)
pca_all$PC <- factor(pca_all$PC, level=pca_all$PC)
#barplot(pca_all, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
library(RColorBrewer)
ggplot(data=pca_all, aes(x=PC, y=Variance, fill=Variance)) +
  theme_bw()+
  ylab("%Variance (Log2CPM)") + xlab(NULL) +
  geom_bar(stat="Identity", width=0.7) +
  scale_fill_gradient(low="skyblue1",high="royalblue4", guide=FALSE)  +
  theme(
    axis.text.x = element_text(angle = 90, hjust=1,
                               vjust=0.5, colour = "black", size=12),
    axis.text.y = element_text(colour = "black", size=12),
    axis.title.x = element_text(colour = "black", size=9),
    axis.title.y = element_text(colour = "black", size=14),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor.y  = element_line(size =0.2, colour="grey88"),
    panel.grid.minor = element_blank(),
    aspect.ratio = 0.75) +
  scale_y_continuous(limits = c(0, 60),
                     expand = c(0, 0))


ggsave("PC_all_var_log2CPM.png", plot = last_plot(), device = "png", path = "~/Dropbox/Luo_lab/acj6/RNAseq/fig/",
       scale = 1, width = 8, height = 4, units = c("in"),
       dpi = 600, limitsize = FALSE)

dev.off()
## End make a plot showing all varience of all PC components

## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)
groups <- factor(dge$samples$group,
                 levels=c("Ctrl", "Acj66"))

## make PCA plot, PC1 vs PC2
library(ggrepel)
pca.data <- data.frame(Sample=rownames(pca$x), group = groups,
                       X=pca$x[,1],
                       Y=pca$x[,2])
# pca.data_sep <- subset(pca.data, pca.data$group=="iMG_ctrl" |
#                          pca.data$group=="iMG_LPS" |
#                          pca.data$group=="iMG_fibril" |
#                          pca.data$group=="iMG_HD")

my_color <- brewer.pal(n=8, "Set2") #[c(3,2,4,5,6)]
ggplot(data=pca.data, aes(x=X, y=Y)) +
  geom_point(aes(fill=group), size=2, alpha=0.8, colour="black", pch=21) + #pch is the shape stype
  scale_fill_manual(values = my_color, breaks = c("Ctrl", "Acj66")) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA (Log2CPM)") +
  theme(axis.text.x = element_text(colour="black",size=15),
        axis.text.y = element_text(colour="black",size=15),
        axis.title.x = element_text(colour="black",size=15),
        axis.title.y = element_text(colour="black",size=15)) +
  geom_text_repel(aes(x=X, y=Y, label = Sample), 
                  color="black", fontface = 'bold',size = 5, box.padding = 0.4,
                  point.padding = 0.5, segment.size=0.5, segment.colour="black")+
  theme(aspect.ratio = 1)

ggsave("PCA12_plot_Log2CPM.png", plot = last_plot(), device = "png", path = "~/Dropbox/Luo_lab/acj6/RNAseq/fig/",
       scale = 0.8, width = 5, height = 5, units = c("in"),
       dpi = 600, limitsize = FALSE)

ggsave("PCA12_plot_Log2CPM.pdf", plot = last_plot(), device = "pdf", path = "~/Dropbox/Luo_lab/acj6/RNAseq/fig/",
       scale = 0.8, width = 5, height = 5, units = c("in"),
       dpi = 600, limitsize = FALSE)
dev.off()

## make PCA plot, PC2 vs PC3
pca.data <- data.frame(Sample=rownames(pca$x),group =groups,
                       X=pca$x[,2],
                       Y=pca$x[,3])
my_color <- brewer.pal(n=8, "Set2") #[c(3,2,4,5,6)]
ggplot(data=pca.data, aes(x=X, y=Y)) +
  geom_point(aes(fill=group), size=2, alpha=0.8, colour="black", pch=21) + #pch is the shape stype
  scale_fill_manual(values = my_color, breaks = c("Ctrl", "Acj66")) +
  xlab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ylab(paste("PC3 - ", pca.var.per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA (Log2CPM)") +
  theme(axis.text.x = element_text(colour="black",size=15),
        axis.text.y = element_text(colour="black",size=15),
        axis.title.x = element_text(colour="black",size=15),
        axis.title.y = element_text(colour="black",size=15)) +
  theme(aspect.ratio = 1)

ggsave("PCA23_plot_Log2CPM.png", plot = last_plot(), device = "png", path = "~/Dropbox/Luo_lab/acj6/RNAseq/fig/",
       scale = 0.8, width = 5, height = 5, units = c("in"),
       dpi = 600, limitsize = FALSE)

ggsave("PCA23_plot_Log2CPM.pdf", plot = last_plot(), device = "pdf", path = "~/Dropbox/Luo_lab/acj6/RNAseq/fig/",
       scale = 0.8, width = 5, height = 5, units = c("in"),
       dpi = 600, limitsize = FALSE)
dev.off()


# retreving differential gene expression table

n <- length(efit$coefficients)

df_genes <- topTable(efit, number=n, sort.by = c("logFC"))

df_genes

# let's annotate the gene set

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl")
head(listDatasets(ensembl))
#
mart<-useMart(dataset="dmelanogaster_gene_ensembl", biomart='ensembl')

# head(listAttributes(mart), 80)
#
biomart_ID <- getBM(attributes=c("flybase_gene_id", "external_gene_name", 
                                 "description"), mart = mart)

# give df_genes a column that matches the flybase_gene_id with biomart_ID
df_genes$flybase_gene_id <- rownames(df_genes)

df_genes_anno <- dplyr::right_join(biomart_ID, df_genes, by="flybase_gene_id")

df_genes_anno <- df_genes_anno[order(df_genes_anno$logFC, decreasing = T),]
write.csv(df_genes_anno, file=paste(path, "../../output/DE_all.csv", sep=""))
save(df_genes_anno, file=paste(path, "/Rdata/DE_all.Rdata", sep=""))

# find genes that are abs(logFC)>1, and FDR <0.05

head(df_genes_anno)
df_genes_anno_sig <- subset(df_genes_anno, logFC >= 1 | logFC <= (-1))
  
df_genes_anno_sig <- subset(df_genes_anno_sig, adj.P.Val < 0.05)


# find genes that are associated with membrane

# head(listAttributes(mart), 80)
#
# biomart_GO <- getBM(attributes=c("flybase_gene_id", "external_gene_name", 
#                                 "description", "name_1006"), mart = mart)

# membrane_genes <- biomart_GO[biomart_GO$name_1006 %like% "membrane", ]

