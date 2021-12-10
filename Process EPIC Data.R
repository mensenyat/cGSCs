library(ChAMP)

# To obtain beadcount, detection P matrix, or Meth UnMeth matrix
testDir = "D:/R/cGSCs/EPIC/All IDAT"
myImport = champ.import(directory = testDir, arraytype = "EPIC")

myLoad = champ.filter(arraytype = "EPIC")
CpG.GUI(CpG=rownames(myLoad$beta),arraytype="EPIC")
myNorm = champ.norm(beta=myLoad$beta,arraytype="EPIC",cores=5, plotBMIQ = TRUE)
champ.QC(myNorm)

# Check affecting covariates
Covar = champ.SVD(beta=myNorm,pd=myLoad$pd,
                  RGEffect = TRUE) # There is not batch effect due to array or Slide, only for Cell Type and Sample Group

goodProbes = read.table("Non-repeated no SNP EPIC Annotation.txt", sep = "\t", header = TRUE)
normGood = myNorm[which(rownames(myNorm) %in% goodProbes$probeID),]

write.table(normGood, file = "EPIC Beta-values icGSCs GBM-DCs.txt", sep = "\t")
