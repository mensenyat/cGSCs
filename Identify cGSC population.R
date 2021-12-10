GBM=read.table("GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv", row.names=1, header=TRUE)
GBM.mean = apply(GBM, 1, mean)
head(GBM.mean)
GBM_t = as.data.frame(as.matrix(t(GBM)))

count0 = sapply(GBM_t, function(x) sum(x==0))
GBMred = GBM[which(count0/ncol(GBM) < 0.95),]
GBMred_t = as.data.frame(as.matrix(t(GBMred)))
high_cutoff = quantile(as.matrix(GBMred), prob=0.75)

NANOG_SOX2_POU = GBMred_t[which(GBMred_t$NANOG > high_cutoff & 
                                  GBMred_t$SOX2 > high_cutoff & 
                                  GBMred_t$POU5F1 > high_cutoff),]

cGSC = GBM_t[which(rownames(GBM_t) %in% rownames(NANOG_SOX2_POU)),]
noncGSC = GBM_t[which(!rownames(GBM_t) %in% rownames(NANOG_SOX2_POU)),]
noncGSCmean = as.data.frame(colSums(noncGSC)/nrow(noncGSC))

#### Separate per patient ####
GBMpat = cbind("Patient" = substr(rownames(GBMred_t), 1, 6), GBMred_t)
cellsxpat = as.data.frame(table(GBMpat$Patient))

NANOG_SOX2_POU_pat = GBMpat[which(GBMpat$NANOG > high_cutoff & 
                                    GBMpat$SOX2 > high_cutoff & 
                                    GBMpat$POU5F1 > high_cutoff),]

cGSCxpat = as.data.frame(table(NANOG_SOX2_POU_pat$Patient))

CellsAndcGSC = merge(cellsxpat, cGSCxpat, all = TRUE, by = 1)
colnames(CellsAndcGSC) = c("Patient", "Total cells", "cGSCs")
CellsAndcGSC[is.na(CellsAndcGSC)] = 0
CellsAndcGSC$Percentage = CellsAndcGSC$`cGSCs`/CellsAndcGSC$`Total cells`*100
CellsAndcGSC = CellsAndcGSC[order(CellsAndcGSC$Percentage),]

library(ggplot2)
p = ggplot(CellsAndcGSC, aes(x = reorder(Patient, -Percentage), y = Percentage, fill = cGSCs)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1.25, vjust = 1.25))
p
#ggsave(filename = "Percentage of cGSCs per patient.pdf", plot = p)

#### Compare cGSCs and Rest of tumor ####
GSCbarcode = rownames(NANOG_SOX2_POU)
allstembarcode = rownames(GBMred_t)[which(!rownames(GBMred_t) %in% GSCbarcode)]

GSC = GBMred_t[which(rownames(GBMred_t) %in% GSCbarcode),]
GSC_t = as.data.frame(as.matrix(t(GSC)))
all = GBMred_t[which(rownames(GBMred_t) %in% allstembarcode),]
all_t = as.data.frame(as.matrix(t(all)))
GSC_and_all = cbind(GSC_t, all_t)

GSC.mean = apply(GSC_t, 1, mean)
all.mean = apply(all_t, 1, mean)

fold = (GSC.mean-all.mean)

# Compute statistical significance #
pvalue = NULL
tstat = NULL
for(i in 1 : nrow(GSC_t)) {
  x = as.numeric(GSC_t[i,])
  y = as.numeric(all_t[i,])
  
  t = wilcox.test(x, y) 
  pvalue[i] = t$p.value
  tstat[i] = t$statistic
}

qvalue = p.adjust(pvalue, method = "fdr", n = length(pvalue)) # Correction by False Discovery Rate
combined  = cbind(qvalue, fold, abs(fold), GSC_and_all)

# Filter per Z-Ratio and qvalue #
fold_cutoff = 0.5
qvalue_cutoff = 0.05

# fold-change filter for "biological" significance
filter_by_fold = abs(fold) >= fold_cutoff
dim(GSC_and_all[filter_by_fold, ]) 

# P-value filter for "statistical" significance
filter_by_qvalue = qvalue <= qvalue_cutoff
dim(GSC_and_all[filter_by_qvalue, ]) 

# Combined filter (both biological and statistical)
filter_combined = filter_by_fold & filter_by_qvalue 

filtered = GSC_and_all[filter_combined,]
dim(filtered)
filtered_combined = merge(combined, filtered, by=0)
genes_with_changes = filtered_combined[,1:4]
#write.table(genes_with_changes, file = "DEG_cGSC_vs_allPlur.csv", sep = ";", row.names = FALSE)
#write.table(rownames(GSC_t), file = "Background_genes.csv", sep = ";") 

# Create a file with a median cGSC and all DEG #
filteredcGSC = filtered[which(colnames(filtered) %in% GSCbarcode)]
GSCmean = as.data.frame(rowSums(filteredcGSC)/ncol(filteredcGSC))
#write.table(GSCmean, "Mean cGSC All DEG.csv", sep = ";")

#pdf("Volcano_Plot_cGSC_vs_allPlur.pdf")
plot(fold, -log10(qvalue), main = "cGSC vs allPlur - Volcano", xlim=c(-1.4, 1.4), ylim=c(0,50), xlab="fold")
points (fold[filter_combined & fold > 0],
        -log10(qvalue[filter_combined & fold > 0]),
        pch = 16, col = "red")
points (fold[filter_combined & fold < 0],
        -log10(qvalue[filter_combined & fold < 0]),
        pch = 16, col = "blue")

abline(v = fold_cutoff, col = "red", lwd = 3)
abline(v = -fold_cutoff, col = "blue", lwd = 3)
abline(h = -log10(qvalue_cutoff), col = "grey", lwd = 3)
dev.off()
