#### cGSC vs icGSC Pathways ####
data = read.table("Pathways by group.csv", sep = ",", header = TRUE, row.names = 1)
data = data[which(data$Number.of.0s == 1),]
data = data[which(duplicated(data$GO) == FALSE),]
data$EnrichmentCoreMinus1[which(data$EnrichmentCoreMinus1 > 50)] = 50
data$EnrichmentMinus1[which(data$EnrichmentMinus1 > 50)] = 50
data$EnrichmentCoreMinus1[which(data$EnrichmentCoreMinus1 < -50)] = -50
data$EnrichmentMinus1[which(data$EnrichmentMinus1 < -50)] = -50

data2 = as.data.frame(cbind("Enrichment_cGSC" = as.numeric(data$Log2cGSC), 
                            "Enrichment_icGSC" = as.numeric(data$Log2icGSC)))

library(ggplot2)
pd = ggplot(data, aes(x=(Log2cGSC), y=(Log2icGSC), color = Group)) + 
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "Log2 Enrichment cGSC",
       y = "Log2 Enrichment icGSC") +
  #geom_smooth(method= "lm", formula= y~x) + 
  geom_hline(yintercept = log2(1.5)) +
  geom_hline(yintercept = -log2(1.5)) +
  geom_vline(xintercept = log2(1.5)) +
  geom_vline(xintercept = -log2(1.5)) +
  scale_color_manual(values=c('green2', "deepskyblue", "darkgoldenrod2", "darkorchid1", "grey", "red"))+
  theme_minimal()
pd
#ggsave(plot=pd,height=NA,width=NA, filename="Plot Pathways Enrichment 3.pdf", useDingbats=FALSE) 

data$cGSC = NA
data$cGSC[which(data$Enrichment.cGSC > 1.5)] = "Upregulated"
data$cGSC[which(data$Enrichment.cGSC < -1.5)] = "Downregulated"
data$icGSC = "Non Significant"
data$icGSC[which(data$Enrichment.Positive.icGSC > 1.5)] = "Upregulated"
data$icGSC[which(data$Enrichment.Negative.icGSC < -1.5)] = "Downregulated"
data$icGSCUp[which(data$Enrichment.Positive.icGSC >= 1.5)] = "Upregulated"
data$icGSCUp[which(data$Enrichment.Positive.icGSC < 1.5)] = "Non-Upregulated"
data$icGSCDown[which(data$Enrichment.Negative.icGSC > -1.5)] = "Non-Downregulated"
data$icGSCDown[which(data$Enrichment.Negative.icGSC <= -1.5)] = "Downregulated"

data2 = data[which(!data$icGSC == "Non Significant"),]

library(questionr)
table(data$cGSC, data$icGSCUp)
odds.ratio(table(data$cGSC, data$icGSCUp))
table(data$cGSC, data$icGSCDown)
odds.ratio(table(data$cGSC, data$icGSCDown))
odds.ratio(table(data2$cGSC, data2$icGSC))
caret::confusionMatrix(table(data2$cGSC, data2$icGSC))

