iter <- ncol(obj$x)
iters <- 1:iter
#iters <- eeSNP
eQTL <- mclapply(iters ,  function(xx) {
		       aa = cor.test(obj$x[,xx], obj$y[, obj$mask2[xx,1]]);
		       return(c(aa$p.value, aa$estimate))
}, mc.cores=1)
eqtl <- t(do.call(cbind, eQTL))


