setwd('/cbcb/project2-scratch/jooslee/srescues/test')
tcga.sur = fread("Curtis_survival.txt")
tcga.sur = tcga.sur[order(V2)]
tcga.mrna = fread('Curtis_mRNA.txt')
tcga.rna1 = as.matrix(tcga.mrna[ ,2:1982,with=F]) 
rownames(tcga.rna1) = tcga.mrna$V1
mrna.header = read.table('Curtis_mRNA.txt',nrow=1 )
colnames(tcga.rna1) = as.character(unlist(mrna.header))
# tcga.rna1 = sapply(t(tcga.rna), rank)

# setkey(tcga.sur, V1)
tcga.sur = tcga.sur[V1 %in% colnames(tcga.rna1)]
tcga.rna = tcga.rna1[,tcga.sur$V1]
known.sl = as.data.frame( fread('SL_list.txt', header=F))
livnat.sl = as.data.frame(fread('Livnat_SLs.txt', header=F))

# mrna.q = ncol(tcga.rna)/2 

#### divide in  quadrant
sl.pair = known.sl
for (inx in seq(nrow(sl.pair))) {

 sl = unlist(sl.pair[inx, ])
 aa1 = tcga.rna[sl[1],] 
 aa1 = ifelse( aa1 > median(aa1) ,1,0)   
 aa2 = tcga.rna[sl[2],] 
 aa2 = ifelse( aa2 > median(aa2) ,1,0)   
 quadrant = 2*aa2 + aa1 + 1
 quadrant = quadrant[tcga.sur$V1] 
 n.temp = table(quadrant)
 n = n.temp[c("1", "2", "3", "4")]
 # aa = interactionTest(n=n,events=tcga.sur$V3,quadrant=quadrant) 
 aa = interactionTest.ks(n=n,events=tcga.sur$V3,quadrant=quadrant, plot=T)
 print(plotInteractionTest.ks(aa$cur))
 print(sl)
 print(aa[c("greater.p", "less.p", "two.sided.p")])
 browser()
}
 
####random gene
for (inx in seq(100)) {

rr = sample.int(17091, 2,replace=F)
sl = rownames(tcga.rna)[rr]
 # sl = unlist(sl.pair[inx, ])
 aa1 = tcga.rna[sl[1],] 
 aa1 = ifelse( aa1 > median(aa1) ,1,0)   
 aa2 = tcga.rna[sl[2],] 
 aa2 = ifelse( aa2 > median(aa2) ,1,0)   
 quadrant = 2*aa2 + aa1 + 1
 quadrant = quadrant[tcga.sur$V1] 
 n.temp = table(quadrant)
 n = n.temp[c("1", "2", "3", "4")]
  aa = interactionTest.ks(n=n,events=tcga.sur$V3,quadrant=quadrant, plot=T)
 print(plotInteractionTest.ks(aa$cur))
 print(sl)
 print(aa[c("greater.p", "less.p", "two.sided.p")])
 browser()
 # aa = interactionTest(n=n,events=tcga.sur$V3,quadrant=quadrant) 
 # aa = interactionTest.chisq(n=n,events=tcga.sur$V3,quadrant=quadrant)
 # print(sl)
 # print(aa[c("greater.p", "less.p", "two.sided.p")])
 # plot(aa$cur[,1], type="l")
 # lines(aa$cur[,2],col="red") 
 # lines(aa$cur[,3],col="green") 
 # browser()
}
#### all gene pairs
med.rna = apply(tcga.rna,1,median, na.rm=T)
tcga.rna.bin = sweep(tcga.rna, 1, med.rna, FUN=">" )

gene1 = 1
 aa1 = tcga.rna.bin[gene1,]
num.genes = nrow(tcga.rna.bin)
# num.genes = 10 
interaction.info = matrix(NA, ncol=6, nrow=num.genes)
 # aa1 = ifelse( aa1 > median(aa1) ,1,0)   
system.time(
	for (gene2 in seq(num.genes)) {
		if(gene1 !=gene2){
			aa2 = tcga.rna.bin[gene2,]
			inx = (!is.na(aa2)) & (!is.na(aa1))
			if(sum(inx) < num.samples){
				aa1.q = aa1[inx]
				aa2.q = aa2[inx]
				tcga.sur.q = tcga.sur[inx] 
				quadrant = 2*aa2.q + aa1.q + 1
				quadrant = quadrant[tcga.sur.q$V1] 
				n.temp = table(quadrant)
				n = n.temp[c("1", "2", "3", "4")]
				interaction.info[gene2,]= unlist(interactionTest.ks(n=n,events=tcga.sur.q$V3,quadrant=quadrant, plot=F))
			}else{
				quadrant = 2*aa2 + aa1 + 1
				quadrant = quadrant[tcga.sur$V1] 
				n.temp = table(quadrant)
				n = n.temp[c("1", "2", "3", "4")]
				interaction.info[gene2,]= unlist(interactionTest.ks(n=n,events=tcga.sur$V3,quadrant=quadrant, plot=F))

			}
		}
	}
)



#####for each gene#### 
#
#
num.samples = ncol(tcga.rna.bin)
tcga.mrna.bin = matrix(0, num.genes, num.samples)
tcga.mrna.bin[tcga.rna.bin] = 1


system.time(
	for (gene2 in seq(num.genes)) {
		print (gene2)
		if(gene1 !=gene2){
			aa2 = tcga.rna.bin[gene2,]
			inx = (!is.na(aa2)) & (!is.na(aa1))
			if(sum(inx) < num.samples){
				aa1.q = aa1[inx]
				aa2.q = aa2[inx]
				tcga.sur.q = tcga.sur[inx] 
				quadrant = 2*aa2.q + aa1.q + 1
				quadrant = quadrant[tcga.sur.q$V1] 
				n.temp = table(quadrant)
				n = n.temp[c("1", "2", "3", "4")]
				interaction.info[gene2,]= unlist(interactionTest_ks(n=n,events=tcga.sur.q$V3,quadrant=quadrant)) 

			}else{
				quadrant = 2*aa2 + aa1 + 1
				quadrant = quadrant[tcga.sur$V1] 
				n.temp = table(quadrant)
				n = n.temp[c("1", "2", "3", "4")]
				interaction.info[gene2,]= unlist(interactionTest_ks(n=n,events=tcga.sur$V3,quadrant=quadrant))

			}
		}
	}
)

tt = matrix(NA, 2,2); tt1 = rep(0,2)
xx = interactionTestks_gene(tt, tt1);

########


interactionTest_ks(n=n,events=tcga.sur$V3,quadrant=quadrant)


sl = unlist(SDRdt[1,1:2,with=F])



################### interacting partner of SLs########
#
drivers = fread("/cbcbhomes/vinash85/project/tcga/data/oncogenic_drivers.csv")



genes = drivers[gene %in% rownames(tcga.rna.bin)]$gene
library(foreach)
library(doMC)
registerDoMC(cores=32)
num.genes = nrow(tcga.rna.bin)
interaction.result = foreach (gene = genes, .inorder=T, .combine="rbind") %dopar% {
    gene1 = which( rownames(tcga.rna.bin) %in% gene) 
    print(gene)
	aa1 = tcga.rna.bin[gene1,]
	interaction.info = matrix(NA, ncol=6, nrow=num.genes)
	for (gene2 in seq(num.genes)) {
		# print (gene2)
		aa2 = tcga.rna.bin[gene2,]
		inx = (!is.na(aa2)) & (!is.na(aa1))
		if(sum(inx) < num.samples){
			aa1.q = aa1[inx]
			aa2.q = aa2[inx]
			tcga.sur.q = tcga.sur[inx] 
			quadrant = 2*aa2.q + aa1.q + 1
			quadrant = quadrant[tcga.sur.q$V1] 
			n.temp = table(quadrant)
			n = n.temp[c("1", "2", "3", "4")]
			if(all(!is.na(n)))
			interaction.info[gene2,]= unlist(interactionTest_ks(n=n,events=tcga.sur.q$V3,quadrant=quadrant)) 

		}else{
			quadrant = 2*aa2 + aa1 + 1
			quadrant = quadrant[tcga.sur$V1] 
			n.temp = table(quadrant)
			n = n.temp[c("1", "2", "3", "4")]
			if(all(!is.na(n)))
			interaction.info[gene2,]= unlist(interactionTest_ks(n=n,events=tcga.sur$V3,quadrant=quadrant))

		}
	}
	return(interaction.info)

}

setkey(drivers, gene)
setnames(drivers, 4, "tumor.suppressor")
drivers.inter=data.table(interaction.result)
setnames(drivers.inter, 1:6, 
	c("stat.g", "stat.l", "stat.t", "p.g", "p.l", "p.t"))
drivers.inter$gene1=rep(genes,each=num.genes)
drivers.inter$gene2=rep(rownames(tcga.rna.bin), length(genes))
drivers.inter= drivers.inter[!(is.na(p.t))]
require(fdrtool)
drivers.inter$qval = fdrtool(drivers.inter$p.t, statistic = "pvalue", plot = F)$qval
drivers.inter$oncogene = drivers[drivers.inter$gene1]$oncogene 
drivers.inter$tumor.suppressor = drivers[drivers.inter$gene1]$tumor.suppressor 

drivers.g = drivers.inter[order(stat.l,decreasing=T)[1:1000]]
drivers.g.onco = drivers.g[oncogene > tumor.suppressor]
write.table(file="drivers.l", unique(drivers.g$gene2), quote=F, row.names=F, col.names=F)
drivers.g.tumor.s = drivers.g[oncogene < tumor.suppressor]
write.table(file="drivers.g.tumor.s", drivers.g.tumor.s$gene2, quote=F, row.names=F, col.names=F)
write.table(file="drivers.less.txt", drivers.inter[], quote=F, row.names=F, col.names=F)

write.table(file="background.pcg", rownames(tcga.rna.bin), quote=F, row.names=F, col.names=F)

drivers.t = drivers.inter[order(stat.t,decreasing=F)[1:1000]]

write.table(file="drivers.greater.txt", drivers.inter[order(stat.g,decreasing=T)[1:1000]], quote=F, row.names=F, col.names=F)
write.table(file="drivers.less.txt", drivers.inter[order(stat.l,decreasing=T)[1:1000]], quote=F, row.names=F, col.names=F)
write.table(file="drivers.match.txt", drivers.inter[order(stat.t,decreasing=F)[1:1000]], quote=F, row.names=F, col.names=F)



