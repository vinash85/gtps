cv.glm <- function( formula, family = gaussian, data, weights = NULL, 
         start = NULL,  model = TRUE, method = "glm.fit",
         x = FALSE, y = TRUE, contrasts = NULL, nfolds=10, parallel=T)
{
  require(glmnet)
  N <- nrow(data)
  foldid=sample(rep(seq(nfolds),length=N))
  predicted <- numeric(N)
   if(nfolds<3)stop("nfolds must be bigger than 3; nfolds=10 recommended")
   outlist=as.list(seq(nfolds))
   if (  parallel && require(foreach)) { # parallel does not work 
    outlist = foreach (i=seq(nfolds), .inorder=T, .combine=c) %dopar% {
      which=foldid==i
      bb <- glm(y~., data=data[!which,,drop=F], family,  weights=weights, 
		start=start, model=model, method=method, x=x, y=y, contrasts=contrasts)  
        list(list( which=which, predicted=predict(bb, data.frame(data[which,]), type="response")))
    }
    for(i in seq(nfolds)){
      outlist.curr <- outlist[[i]]
      predicted[outlist.curr$which]  <-  outlist.curr$predicted
    }
  }else{
    for(i in seq(nfolds)){
      which=foldid==i
      bb <- glm(y~., data=data[!which,,drop=F], family,  weights=weights, 
		start=start, model=model, method=method, x=x, y=y, contrasts=contrasts)  
      predicted[which] <- predict(bb, data.frame(data[which,]), type="response")
    }
  }
   glmobj <- glm(y~., data=data, family,  weights=weights, 
		start=start, model=model, method=method, x=x, y=y, contrasts=contrasts) 
 list(predicted= predicted, glmobj=glmobj) 
}
glm.cv.predict <- function(data,train.inx, family=gaussian)
{
   N <- nrow(data)
  predicted  <- numeric(N)  
  dd <- cv.glm(y~., data=data[train.inx,], nfolds=nrow(glm.df), family=family)
  predicted[train.inx] <- dd$predicted
  other.inx  <-  !( seq(N) %in% train.inx)
  if(sum(other.inx) > 0) 
    predicted[ other.inx] <-   predict(dd$glmobj, data.frame(data[other.inx,]), type="response")
  predicted
} 

model.gene <- function(xx, yy.all, inx.curr, family="gaussian")
{
  
  require(glmnet);
  yy <-yy.all[inx.curr] 
  inx1 <- inx.curr[ which(!is.na(yy)  )]
  inx <- inx1
  cvobj = cv.glmnet(xx[inx,], yy, parallel=T, family=family, alpha=0.5)
  aa  <-  coef(cvobj, s="lambda.min") 
  gene.selected  <- aa[which(aa !=0),,drop=F]
  gene.selected <- rownames(gene.selected)[-1]
  if(length(gene.selected > 0) ) {
     glm.df <- data.frame(xx[, gene.selected])
  glm.df$y <- yy.all
  if(family == "gaussian") predicted  <- glm.cv.predict(data=glm.df, train.inx=inx.curr, family=gaussian) 
  if(family == "binomial"){
    predicted  <- glm.cv.predict(data=glm.df, train.inx=inx.curr, family=binomial) 
    yy.all  <-  as.numeric(yy.all)
    yy <-  as.numeric(yy)
  }
  correlation <- cbind( unlist((cor.test(yy.all, predicted))[c("estimate","p.value")]),
  unlist((cor.test(yy.all[healthy.inx], predicted[healthy.inx]))[c("estimate","p.value")]),
  unlist((cor.test(yy.all[failure.inx], predicted[failure.inx]))[c("estimate","p.value")])
  )
  #browser()
  err <- (predicted- yy.all)
  r2.all  <- 1 - mean((err - mean(err))^2 )/mean((yy.all - mean(yy.all))^2)
  r2.healthy  <- 1 - mean((err[healthy.inx] - mean(err[healthy.inx]))^2 )/mean((yy.all[healthy.inx] - mean(yy.all[healthy.inx]))^2)
  r2.failure  <- 1 - mean((err[failure.inx] - mean(err[failure.inx]))^2 )/mean((yy.all[failure.inx] - mean(yy.all[failure.inx]))^2)
  r2.cv <- 1 - min(cvobj$cvm)/mean( ((yy) - mean(yy))^2)
  r2 <- c(r2.cv, r2.all, r2.healthy, r2.failure)
		}else{
		  predicted = 0; r2 = 0 ; correlation =0
		  warning("no gene seelected")
		}
  list(predicted = predicted, r2=r2, correlation=correlation, gene.selected=gene.selected) 
 }
