load("/cbcbhomes/vinash85/shortcuts/lasso/data/Data/MikeData/magnet_PedTraits.RData")
library(data.table)
expression <- pedtrait[, grep(x=colnames(pedtrait), pattern="^X"), with=F]
expression.z <- scale(expression, scale=apply(expression, 2, sd), center=T)

library(samr)
library(erer)
gst = 6
geneno = 16271
gen = gst+geneno-1
#genexp = t(pedtrait[,gst:gen])
#diseaseid = pedtrait$disease 
osam = SAM(t(expression.z), y=pedtrait$disease, resp.type=c("Two class unpaired"), genenames = colnames(expression.z))
write(t(osam$siggenes.table$genes.up),file="data_differential_gene_up.txt", ncolumns=7, sep="\t") # disease=2 has genes upregulated
write(t(osam$siggenes.table$genes.lo),file="data_differential_gene_low.txt", ncolumns=7, sep="\t") # disease=2 has genes downregulated

osam3 = SAM(genexp, y=diseaseid, resp.type=c("Two class unpaired"), genenames = gname, testStatistic=c("wilcoxon"))
write(t(osam3$siggenes.table$genes.up),file="data_differential_wilcoxon_gene_up.txt", ncolumns=7, sep="\t") # disease=2 has genes upregulated
write(t(osam3$siggenes.table$genes.lo),file="data_differential_wilcoxon_gene_low.txt", ncolumns=7, sep="\t") # disease=2 has genes downregulated



print(osam)
plot(osam)

print(osam3)
plot(osam3)

# manual checking for genes
temp = pedtrait[, list(X8160168,disease)]
library(ggplot2) 
ggplot(data=temp,aes(factor(disease), X8160168)) + geom_boxplot()  
wilcox.test(pedtrait$X8160168[pedtrait$disease==2 ], pedtrait$X8160168[pedtrait$disease==1 ], alternative = "greater")
summary(pedtrait[disease==2, X8160168])
summary(pedtrait[disease==1, X8160168])

temp = pedtrait[, list(X7909104,disease)]
library(ggplot2) 
ggplot(data=temp,aes(factor(disease), X7909104)) + geom_boxplot()  
wilcox.test(pedtrait$X7909104[pedtrait$disease==2 ], pedtrait$X7909104[pedtrait$disease==1 ], alternative = "greater")
summary(pedtrait[disease==2, X7909104])
summary(pedtrait[disease==1, X7909104])


# analysis with limma

library(limma)
exp.mat  <- t(as.matrix(expression))
design <- cbind(grp1=1, disease=pedtrait$disease) 
fit  <-  lmFit(exp.mat, design)
# Fold-change thresholding
fit2 <- treat(fit,lfc=0.1)
topTreat(fit2,coef=2)

date=Sys.Date()
outDir <- paste0("result/", date)
dir.create(outDir)
write.fit(fit2,file=paste0(outDir, "/limMaFit.MagNet.Genes.txt"))
pdf(paste0(outDir, "/limMaFit.MagNet.Genes.pdf"))
volcanoplot(fit2,coef=2,highlight=100, names=rownames(exp.mat))
dev.off()



# linear model of phenotype variance by gene expression

#gene.df  <- copy(expression)
gene.df  <- copy(expression)
gene.df$disease  <- diseaseid$disease -1 
annova.gene  <-  aov(disease ~ ., data=gene.df)
temp3 <- cbind(gene.df[,c(1:311, 313), with=F], gene.df[,disease])
lm.gene <- lm(V2 ~., data=temp3)
#explained variance
xx <- as.matrix(expression.z[, top.diff.genes]) ; yy = gene.df$disease
yy <- scale(yy, center=T,scale=sd(yy))
library(MASS)
#to big to execute
bb  <- ginv(t(xx) %*% xx)  %*% xx %*% yet 
ev.gene <- 1- var( yy - xx %*% bb)/var(yy)
#manova.gene  <-  manova(disease ~ ., data=gene.df)


# princomp 
top.genes.tab <- topTable(fit2, coef=2, number=200)
top.diff.genes <- rownames(top.genes.tab)
pca.gene <- princomp(t(exp.mat[top.diff.genes,]), cor=TRUE)

cor.test(pca.gene$loading, lm.gene$disease)

