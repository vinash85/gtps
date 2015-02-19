hcm.gene  = c("GLA", "MYL2", "MYL3", "MYLK2", "MYOZ2", "PRKAG2", "TTR", "RYR2", "CASQ2")
both.gene = c("ACTC1", "ACTN2", "ANKRD1", "BAG3", "CAV3", "CSRP3", "LAMP2", "LDB3", "MYBPC3","MYH6", "MYH7", "NEXN", "PLN", "TNCC1", "TNN13", "TNNT2", "TPM1", "VCL" , "TTN", "JUP")
dcm.gene = c("ABCC9", "CRYAB", "CTF1", " EMD", "FHL2", "GATAD1", "LAMA4", "LMNA", "NEBL", "RBM20","SCN5A", "SGCD", "TAZ", "TCAP","TMPO", "DES", "DSC2", "DSG2", "DSP", "PKP2", "TMEM43")
arvc.gene = c("ARVC")
library(data.table)
cardio.genes = data.table(geneNames = c(hcm.gene, both.gene, dcm.gene, arvc.gene), 
	type = c(
		rep("hcm", length(hcm.gene)),
		rep("dcm", length(dcm.gene)),
		rep("both", length(both.gene)),
		rep("arvc", length(arvc.gene))
		))
# find the affyid
library(hugene10sttranscriptcluster.db)
x.temp <- hugene10sttranscritpclusterALIAS2PROBE
xx <- as.list(x.temp)
all.genes  <-  unlist(xx[additional.gene])

# find differential expression 
# figure with lvef
# predict cardiomyopathy
# predict lvef 
# run with different phenotype if there is a correlation with anything
# find gwas signal of these genes
# find other gene correleated with these genes
# linc database compound exclusively targetting these subset of genes 
