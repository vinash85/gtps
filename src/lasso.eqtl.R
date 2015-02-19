lasso.new.eqtl <- function (x, y, gene.map=NULL, nthreads = 1)
{   
    require(glmnet)
    require(doMC)
    require(foreach)
    registerDoMC(cores = nthreads)
    n <- nrow(x)
    coefs <- numeric(ncol(x))
    genes <- ncol(y)
    coefs <- foreach(gene = genes, .inorder = T) %dopar% {
        if(!is.null(gene.map)) {
            curr.inx <- which(gene.map == gene)
            }else curr.inx = seq(ncol(x))
        fit <- cv.glmnet(x = x[, curr.inx], y = y[, gene], alpha = 0.5,
            parallel = T)
        coef = coef(fit, s = "lambda.min")
        coef <- coef[-1]
        coef
    }
    coef = do.call(c, coefs)
    coefs
}
