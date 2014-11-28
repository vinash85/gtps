setwd("independent_test")
sample.name <- rownames(expression.z)  
train <- sample(sample.name, size=200, replace=F)  
test <- sample.name[!( sample.name %in%  train)]
expression.train  <- expression.diff.hf[train, ]; 
expression.test  <- expression.diff.hf[test, ];
pedtrait.info.train <- pedtrait.info[train]
gene.weight <- apply(expression.train ,  2,function(tt) fisher.score(tt, pedtrait.info.train$disease) )
pca.gene <- PCA(expression.train, col.w=gene.weight, ncp=ncol(expression.train), graph=F)
gene.proj.df <- as.data.table(expression.test %*% pca.gene$var$coord)
eigen.gene <- colnames(pca.gene$ind$coord)
#proj.gene  <- expression.diff.hf[test,] %*%  pca.gene$rotation
gene.proj.df$sample_name  <-  test

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
labs(y="Phenotypic deviation index (Rank) ", x = "Transcriptomic Deviation index (Rank)",   title=paste0(curr.name , "\n correlation =",  format( aa$estimate, digits=3) , " \n p =", format(aa$p.value, digits=3, scientific=T)))   
ggsave(paste0( "tdi_weighted_imputed_",curr.name, ".pdf"))
}



setwd("../")
