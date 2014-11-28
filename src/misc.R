## this is code that is not removed from anaslysis ## 

# samtools 


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


