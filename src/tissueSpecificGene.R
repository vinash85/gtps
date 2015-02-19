library(avinash)
setwd("/cbcb/project2-scratch/vinash85/project/1000gene/data/tissues")
load("/cbcb/project2-scratch/vinash85/project/1000gene/data/tissues/gtex.pc.RData")
load("/cbcb/project2-scratch/vinash85/project/1000gene/data/tissues/sample.info.RData")
gtex.mat = as.matrix(gtex.pc[,dt.merge$SAMPID,with=F]) 
gtex.mat =  sweep(gtex.mat, 2, colMeans(gtex.mat), FUN="/")
gtex.mat =  sweep(gtex.mat, 1, rowMeans(gtex.mat), FUN="/")
# aggregate(gtex.mat[1:4,], by=dt.merge$SMTS, FUN=rowMean)
tissues=sort(unique(dt.merge$SMTS))

tissues.type=dt.merge$SMTS
numGenes = nrow(gtex.mat)
tissue.exp.orig  =  numeric( length(tissues) *nrow(gtex.mat))

# colnames(tissue.exp) = tis
for (inx in seq(length(tissues))) {
	tissues.inx = which(tissues.type == tissues[inx])
 tissue.exp.orig[((inx -1) * numGenes ) + seq(numGenes) ]   = rowMeans(gtex.mat[,tissues.inx])
}

# null distribution by permuation test
tissues.index = list()
for (inx in seq(length(tissues))) {
 tissues.index[[inx]] =  which(tissues.type == tissues[inx])
}		
maxiter = 100000
iters = 1:maxiter
library(foreach)
library(doMC)
registerDoMC(cores=32)

numSamp = 10* length(tissues)
tissue.exp.pert = foreach (iter = iters, .inorder=F, .combine="cbind") %dopar% {
# tissue.exp.pert = foreach (iter = iters, .inorder=F) %dopar% {
	tissues.type = unlist( lapply(tissues.index, function(xx) sample(xx, 10, replace=T)))[sample(numSamp)]
 print(head(tissues.type))
	# tissues.type=dt.merge$SMTS[sample(length(dt.merge$SMTS))]
	tissue.exp = numeric( length(tissues) *nrow(gtex.mat))

	# colnames(tissue.exp) = tissues
	for (inx in seq(length(tissues))) {
		# tissues.inx = which(tissues.type == tissues[inx])
		tissue.inx = tissues.type[((inx -1) * 10 ) + seq(10)]
		tissue.exp[((inx -1) * numGenes ) + seq(numGenes) ] = rowMeans(gtex.mat[,tissue.inx])
	}
	return(tissue.exp)

}
p.val =  rowSums( sweep(tissue.exp.pert, 1, tissue.exp.orig, FUN="<"))/maxiter
save(file="tissue.Specific.RData", p.val, tissue.exp.orig, tissues )
save(file="tissue.Specific.pert.RData", tissue.exp.pert)
# 