#find variable gene in heart failure
setkey(pedtrait, sample_name)

expression <- pedtrait[, grep(x=colnames(pedtrait), pattern="^X"), with=F]
pedtrait.hf <- pedtrait[disease==2, c(1:5, grep(x=colnames(pedtrait), pattern="^X")), with=F]
expression.hf <- pedtrait.hf[,  grep(x=colnames(pedtrait.hf), pattern="^X"), with=F]
exp.sd <- apply(expression.hf,2,sd)
expression.hf.var <- expression.hf[, order(exp.sd, decreasing =T)[1:1000], with=F]
pca.gene <- PCA(expression.hf.var,  ncp=ncol(expression.hf.var), graph=F)
gene.proj.df <- as.data.table(pca.gene$ind$coord)
eigen.gene <- colnames(pca.gene$ind$coord)
gene.proj.df$sample_name  <-  pedtrait.hf$sample_name

physio.imputed2 <- physio.imputed1
setkey(physio.imputed2, sample_name)
physio.imputed.hf <- physio.imputed2[as.character(pedtrait.hf$sample_name)]
aa  <-  sapply(gene.proj.df[, eigen.gene,with=F],  function(xx) {bb <- cor.test(xx , physio.imputed.hf$lvef.p); c(bb$estimate, bb$p.value)} )

gene.proj.df <- as.data.table(as.matrix(expression[,  order(exp.sd, decreasing =T)[1:1000], with=F] ) %*% pca.gene$var$coord)
eigen.gene <- colnames(pca.gene$ind$coord)
gene.proj.df$sample_name  <-  pedtrait$sample_name
expression.merge <- merge(x=physio.reduced, y=gene.proj.df, by = "sample_name", sort = F)
donor.eigen.gene <- colMeans( expression.merge[disease.p==1, eigen.gene, with=F])
expression.merge$tdi <- apply(expression.merge[,eigen.gene, with=F], 1 , function(tt) sum((tt - donor.eigen.gene)^2)) 

params = c( "lvef.p","lvedd.p","lvesd.p","mitral_regurgitation.p","tricuspid_regurgitation.p","cardiac_index.p","creatinine_level.p") 

for(curr.param in params){
  curr.name <- toupper(gsub(curr.param, pattern="\\.p$", replace=""))
#donor.di <- colMeans( expression.merge[disease.p==1, curr.param, with=F])
expression.merge$pdi <- eval(parse(text=paste0("expression.merge$",curr.param)))
aa <-  cor.test(expression.merge$pdi, expression.merge$tdi, ,method="spearman")
p <- ggplot(data=expression.merge, aes(x=rank(tdi), y=rank(pdi)))
p+geom_point(aes(col=factor(etiology)) ) + geom_smooth(method=lm) +
labs(y= curr.name, x = "Transcriptomic Deviation index (Rank)",   title=paste0(curr.name , "\n correlation =",  format( aa$estimate, digits=3) , " \n p =", format(aa$p.value, digits=3, scientific=T)))   
ggsave(paste0( "tdi_weighted_",curr.name, ".pdf"))
}





# in heart failure perform PCA
# plot
