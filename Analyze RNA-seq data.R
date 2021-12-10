#### Process using DESEQ2 ####
data = read.table("RNAseq_cGSCs_GBMdc_no_names.csv", sep = ";", header = TRUE, row.names = 1)
#####data = data[,c(4,5,6,13,14,15)] # Maintain only cGSCs and GBMdc
meta = read.table("Metadata_geneNames.csv", sep = ";", quote = "", header = TRUE)
data = merge(meta, data, by.x = 1, by.y = 0)
#write.table(data, "RNAseq_cGSC_and_GBMdc_names.csv", sep = ";", row.names = FALSE)

meta = read.table("GeneList hg19.txt", sep = "\t", header = TRUE)

data = read.table("RNAseq_cGSC_and_GBMdc_names.csv", sep = ";", header = TRUE)
data = data[,-c(1,4)]
data = data[which(data$SYMBOL %in% meta$name2),]
data = data[which(duplicated(data$SYMBOL) == "FALSE"),]
data = data[which(duplicated(data$ENSEMBL) == "FALSE"),]
data = data[which(is.na(data$ENSEMBL) == FALSE),]

meta = meta[which(meta$name2 %in% data$SYMBOL),]

library("countToFPKM")
library("biomaRt")
library("dplyr")

#### Obtain feature length ####
## Import feature counts matrix
counts = data[,-c(1,2)]
rownames(counts) = data[,1]
#######counts = counts[,c(4,5,6,13,14,15)]
colnames(counts) = c("GBMdc1", "GBMdc2", "GBMdc3", "icGSC1", "icGSC3", "icGSC7")

## Build a biomart query 
ens_build = "sep2015"
dataset="hsapiens_gene_ensembl"
mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset, 
                        host="www.ensembl.org", path = "/biomart/martservice", archive = FALSE)

gene.annotations = biomaRt::getBM(mart = mart, 
                                   attributes=c("ensembl_gene_id", "external_gene_name", 
                                                "start_position", "end_position"))
gene.annotations = dplyr::transmute(gene.annotations, external_gene_name,  ensembl_gene_id, 
                                     length = end_position - start_position)

# Filter and re-order gene.annotations to match the order in feature counts matrix
gene.annotations = gene.annotations %>% dplyr::filter(ensembl_gene_id %in% rownames(counts))
gene.annotations = gene.annotations[order(match(gene.annotations$ensembl_gene_id, rownames(counts))),]

# Assign feature lenghts into a numeric vector.
featureLength = gene.annotations$length

# Import sample metrics.
# Assign mean fragment length into a numeric vector.
samples.metrics = read.table("MeanFragmentLength.csv", sep=";", header=TRUE, row.names = 1)
meanFragmentLength = samples.metrics$Mean.Length
meanFragmentLength = meanFragmentLength##############[c(4,5,6,13,14,15)]

# Return FPKM into a numeric matrix.
counts = counts[which(rownames(counts) %in% gene.annotations$ensembl_gene_id),]
fpkm_matrix = fpkm (counts, featureLength, meanFragmentLength)
fpkm_all = merge(gene.annotations, fpkm_matrix, by.x = 2, by.y = 0)
fpkm_all = fpkm_all[,-c(1,3)]
fpkm_all = fpkm_all[which(duplicated(fpkm_all$external_gene_name) == FALSE),]
rownames(fpkm_all) = fpkm_all$external_gene_name
fpkm_all = fpkm_all[,-1]
colnames(fpkm_all) = c("GBMdc1", "GBMdc2", "GBMdc3", "icGSC1", "icGSC3", "icGSC7")

meanExp = apply(fpkm_all, 1, mean)
fpkm_all = fpkm_all[which(!meanExp == 0),]
#write.table(fpkm_all, "RNA-seq_FPKM_icGSC.csv", sep = ";")

## Import feature counts matrix
counts = merge(gene.annotations, counts, by.x = 2, by.y = 0)
rownames(counts) = make.names(counts$external_gene_name, unique = TRUE)	
counts = counts[,-c(1,2,3)]
library(DESeq2)
coldata = as.data.frame(cbind("Sample" = colnames(counts), "batch" = 1, 
                              "condition" = c("GBMdc", "GBMdc", "GBMdc", "icGSC", "icGSC", "icGSC")))

dds = DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design= ~ condition)

# Filter genes with low counts --> At least 6 counts, so at least one count per #
keep = rowSums(counts(dds)) >= 6 # Discard 4,655 genes, maintain 17,465 genes
dds = dds[keep,]

# Differential expression analysis #
dds = DESeq(dds)
res = results(dds)
res

qv = p.adjust(res$pvalue, method = "fdr", n = length(res$pvalue))

resLFC = lfcShrink(dds, coef="condition_icGSC_vs_GBMdc", type="apeglm")
resLFC

summary(res)

# Select signficant #
res05 = results(dds, alpha=0.05)
summary(res05)

# Independent hypothesis weighting # 
library("IHW")
resIHW = results(dds, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.05, na.rm=TRUE)
metadata(resIHW)$ihwResult

# MA-plot #
plotMA(res)

# Plot Volcano #
filter_by_Zratio = abs(res$log2FoldChange) >= 1.5
filter_by_qvalue = qv <= 0.0000000001
filter_combined = filter_by_Zratio & filter_by_qvalue

plot(res$log2FoldChange, -log10(qv), main = "cGSCs vs GSC - Volcano", xlab="Zratio")
points (res$log2FoldChange[filter_combined & res$log2FoldChange > 0],
        -log10(qv[filter_combined & res$log2FoldChange > 0]),
        pch = 16, col = "red")
points (res$log2FoldChange[filter_combined & res$log2FoldChange < 0],
        -log10(qv[filter_combined & res$log2FoldChange < 0]),
        pch = 16, col = "blue")

abline(v = 1, col = "red", lwd = 3)
abline(v = -1, col = "blue", lwd = 3)
abline(h = -log10(0.05), col = "grey", lwd = 3)


# Plot counts of the most different gene#
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
d = plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

# More information on results columns #
mcols(res)$description

#write.table(as.data.frame(res), file = "DESEQ2 cGSC vs GSC.csv", sep = ";")

# Heatmap of the count matrix #
ntd = normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd)) # Transformate so the ones with top number of count do not DESTROY everything


library("pheatmap")
select = order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df = as.data.frame(colData(dds)[,c("condition")])
rownames(df) = colnames(counts)
pheatmap(assay(ntd), cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Heatmap of sample-to-sample distances
sampleDists = dist(t(assay(ntd)))
library("RColorBrewer")
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(ntd$condition, ntd$type, sep="-")
colnames(sampleDistMatrix) = NULL
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#pdf("Distance Matrix.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()
plotPCA(ntd, intgroup=c("condition"))

# Correlation matrix #
library(corrplot)
library("RColorBrewer")
library(pheatmap)
resCor = cor(assay(ntd), method = "spearman") # 5000 most variable regions
colors = colorRampPalette( (brewer.pal(9, "Blues")) )(255)

#pdf("RNA-seq Correlation Matrix.pdf")
pheatmap(as.matrix(resCor),
         col=colors)
dev.off()

# Heatmap of all DEG #
resSign = res[which(abs(res$log2FoldChange)>1),]
resSign = resSign[which(resSign$padj<1e-10),]
resSign = resSign[which(abs(resSign$log2FoldChange)>1.5),]

dim(resSign)

FPKM = read.table("RNA-seq_FPKM_icGSC.csv", sep = ";", header = TRUE, row.names = 1)
FPKMsign = FPKM[which(rownames(FPKM) %in% rownames(resSign)),]

data_meta = as.data.frame(as.matrix(t(FPKM)))
data_meta$Group = c("GBMdc", "GBMdc", "GBMdc", "icGSC", "icGSC", "icGSC")
data_meta$Color = NA
data_meta$Color[which(data_meta$Group == "GBMdc")] = "blue"
data_meta$Color[which(data_meta$Group == "icGSC")] = "red"

FPKMsign2 = as.data.frame(as.matrix(t(scale(as.data.frame(as.matrix(t(FPKMsign)))))))

colv = as.dendrogram(hclust(as.dist(1-cor(t(as.matrix(FPKMsign2))))))
rowv = as.dendrogram(hclust(as.dist(1-cor(as.matrix(FPKMsign2)))))

library(gplots)
#pdf("Heatmap All DEG Scaled icGSC vs GBMdc.pdf")
heatmap.2(as.matrix(FPKMsign2), cexCol=0.7,
          col = rev(redblue(256)), scale = "none", trace=NULL, tracecol=NULL, ColSideColors = data_meta$Color)
legend('topleft', legend = unique(data_meta$Group), fill = unique(data_meta$Color), border=T, title='Subtype')
dev.off()

#### FPKM to TPM ####
FPKM = read.table("RNA-seq_FPKM_cGSC vs GSC.csv", sep = ";", header = TRUE, row.names = 1)
colTotal = colSums(FPKM)
TPM = as.data.frame(matrix(NA, ncol = ncol(FPKM), nrow = nrow(FPKM)))
for(i in 1:ncol(FPKM)){
  TPM[,i] = (FPKM[,i]/colTotal[i])*10^6
}

TPMlog = log2(TPM+1)
#write.table(TPMlog, file = "Log2-TPM cGSC vs GSC.csv", sep = ";")