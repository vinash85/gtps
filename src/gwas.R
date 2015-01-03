iters  <- obj$mask[[1715]]  
eQTL <- mclapply( iters,  function(xx) {
		       aa = cor.test(obj$x[,xx], obj$y[, obj$mask2[xx,1]]);
		       return(c(aa$p.value, aa$estimate))
}, mc.cores=12)
eqtl <- t(do.call(cbind, eQTL))




#### 
eeSNP <- which(gamMa > 0.5)
load("~/project/gtps/result/2014-11-05/pedtrait.info.RData")
pedtrait.info1 <- pedtrait.info[rownames(obj$y)]
pedtrait.info1$disease1 <- pedtrait.info1$disease -1
GWAS <- mclapply(iters ,  function(xx) {
		       aa = cor.test(obj$x[,xx], (pedtrait.info1$disease1));
		       return(c(aa$p.value, aa$estimate))
}, mc.cores=12)
gwas <- t(do.call(cbind, GWAS))

gwas.eeSNP <- mclapply(eeSNP ,  function(xx) {
		       aa = cor.test(obj$x[,xx], (pedtrait.info1$disease1));
		       return(c(aa$p.value, aa$estimate))
}, mc.cores=12)
gwas.eeSNP <- t(do.call(cbind, gwas.eeSNP))


iter <- ncol(obj$x)
iters <- 1:iter
#iters <- eeSNP

library(limma)
exp.mat  <- t(as.matrix(obj$y))
design <- cbind(grp1=1, disease=pedtrait.info1$disease1) 
fit  <-  lmFit(exp.mat, design)
# Fold-change thresholding
fit2 <- treat(fit,lfc=0.1)
topTreat(fit2,coef=2)
fit2 <- eBayes(fit2) 
date=Sys.Date()
outDir <- paste0("result/", date)
dir.create(outDir)

topGene <- topTreat(fit2,coef=2,200)

