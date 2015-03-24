library(avinash)
tab = fread("../metadata/inst.info")
cid  = fread("../level3/lincs.cid", header=F)
cid[,distil_id:=paste(V1,V2,sep=":")]
cid[,myid:=1:nrow(cid)]
setkey(cid, distil_id)
setkey(tab, distil_id)
tab1 = merge(cid[,3:4,with=F], tab)
tab.skl = tab1[cell_id=="SKL"]

rid  = unlist(fread("../level3/lincs.rid", header=F))
load("~/shortcuts/Data/MikeData/MAGnet_combat_eset.Rdata") 
eset.combat
expression <- exprs(eset.combat)
physiological.tab <- fread("~/project/gtps/gtps/doc/magnet_redcap.txt")
setkey(physiological.tab, sample_name)
expression.annot = physiological.tab[colnames(expression)]
expression.annot[,disease.p:=ifelse(etiology=="donor",1,2)]
expression.annot[,history_of_afib_aflutter.p:=sapply(history_of_afib_aflutter, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
expression.annot[,history_of_diabetes.p:=sapply(history_of_diabetes, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
expression.annot[,history_of_hypertension.p:=sapply(history_of_hypertension, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
# expression = t(expression)
expression.annot[, lvef.p:=ifelse(grepl(lvef, pattern="normal", ignore.case=T), 58, ifelse(grepl(lvef, pattern="-"),
					sapply(lvef, function(tt) mean(as.numeric(strsplit(unique(tt)[1],split="\\-")[[1]]))), 
					ifelse( grepl(lvef, pattern="<25%"), 25, ifelse(grepl(lvef, pattern="<"), as.numeric(gsub(lvef, pattern="<", replacement="")),
									     as.numeric(lvef)))))]
expression.annot[,lvef.p:=ifelse(lvef.p==0,NA,lvef.p)]
expression.annot[,lvef.p:=ifelse(lvef.p<1, lvef.p * 100,lvef.p)]



 library(biomaRt)
 ensembl <- useMart("ensembl", dataset= "hsapiens_gene_ensembl")
 affy2ensembl = getBM(attributes = c("affy_hg_u133a",  "affy_hugene_1_0_st_v1"), mart=ensembl)
 affy2ensembl = getBM(attributes = c("ensembl_id",  "entrezgene"), mart=ensembl)
 affy2ensembl1 = getBM(attributes = c("hgnc_symbol",  "entrezgene"), mart=ensembl)
affy2ensembl = affy2ensembl[affy2ensembl$affy_hg_u133a %in%  unlist(rid), ]
affy2ensembl = affy2ensembl[!is.na(affy2ensembl$affy_hugene_1_0_st_v1),]
affy2ensembl = data.table(affy2ensembl)
affy2ensembl = affy2ensembl[ affy_hugene_1_0_st_v1%in% rownames(expression)]


if(F){

	library(hugene11stprobeset.db)
xx <- data.table(toTable(hugene11stprobesetALIAS2PROBE))
setkey(xx,alias_symbol)

setnames(xx, 1, "hugene11stprobeset")
xx = xx[ hugene11stprobeset%in% rownames(expression)]

library(hgu133a2.db)
 yy <- data.table(toTable(hgu133a2ALIAS2PROBE))
setkey(yy,alias_symbol)
setnames(yy, 1, "hgu133a2")
yy = yy[hgu133a2 %in% rid]
}

#### skeletol muscle tissue expression 
sklMat = fread("lincs.sklMat.csv", header=F)
sklMat <- as.matrix(sklMat)


colnames(sklMat) = tab.skl$distil_id
rownames(sklMat) =  rid
sklMat.match = sklMat[affy2ensembl$affy_hg_u133a, ]
expression.match = expression[as.character(affy2ensembl$affy_hugene_1_0_st_v1),] 
plotSKL=F
if(plotSKL){
p1 <- hist(sklMat.match[,tab.skl.ctl$distil_id[1] ], 100)                     # centered at 4
p1 <- hist(sklMat.match[,tab.skl.ctl$distil_id[1] ], 100)                     # centered at 4
p3 <- hist(sklMat.match[,1 ], 100)                     # centered at 4
p2 <- hist( expression.match[,1],100)                     # centered at 6
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,15.1))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,15.1), add=T)
plot( p3, col=rgb(0,1,0,1/4), xlim=c(0,15.1), add=T)
}
gse.exp = sklMat.match
gseConv.exp =  apply(gse.exp, 2, function(tt) mrs(nonreference=tt, reference=exp.median))
rownames(gseConv.exp) = rownames(expression.match)
aa = colnames(gseConv.exp) 
bb = gsub(aa, pattern=":", replacement=".")
colnames(gseConv.exp) = bb
### from geo database
#Download GDS file, put it in the current directory, and load it: viral cardiomyocytes data
library(Biobase)
library(GEOquery)
if(F){


gse57781 <- getGEO('GSE57781', destdir=".")
gse.exp = exprs(gse57781$GSE57781_series_matrix.txt.gz)

gpl16686 <- getGEO('GPL16686', destdir=".")
gpl16686 = Table(gpl16686)

gpl16686.dt = data.table(do.call(rbind, 
	strsplit(as.character(gpl16686$SPOT_ID), ":|-" )))
gpl16686.dt$ID = gpl16686$ID
gpl16686.dt = gpl16686.dt[!(V3 =="unknown")]


gpl11532 <- getGEO('GPL11532', destdir=".")
gpl11532 = Table(gpl11532)
gpl11532 = data.table(gpl11532)
gpl11532.dt =  gpl11532[,list(seqname, RANGE_START, RANGE_STOP, ID)]
annot.int = intersectBed(a=gpl16686.dt, b=gpl11532.dt ,opt.string="-wa -wb -f .95" )

expression.match = expression[as.character(annot.int$V8),]
gse.exp = gse.exp[as.character(annot.int$V4),]
 

}

##### if lincs 
 ## hg133 UA for LINCS data ### no coordiantes available 
# gpl96 <- getGEO('GPL96', destdir=".")
# gpl96 = Table(gpl96)



### given expression match and gse.ex
exp.median = sort(apply(expression.match,1, median))
gseConv.exp =  data.table(apply(gse.exp, 2, function(tt) mrs(nonreference=tt, reference=exp.median)))
gseConv.exp$gseid=as.character(annot.int$V4); gseConv.exp$magnetid=as.character(annot.int$V8)
aa = (phenoData(gse57781$GSE57781_series_matrix.txt.gz))@data
gseConv.exp[,del_rep1:= GSM1388431 - GSM1388430 ]
gseConv.exp[,del_rep2:= GSM1388429 - GSM1388432 ]
gseConv.exp$donor=apply(expression.match[,expression.annot$disease.p==1],1 ,mean)
gseConv.exp$hf=apply(expression.match[,expression.annot$disease.p==2],1 ,mean)
gseConv.exp$min=apply(expression.match,1 ,min)
gseConv.exp$max=apply(expression.match,1 ,max)
gseConv.exp$sd=apply(expression.match,1 ,sd)
threshold = .25
gseConv.exp[,sel:= ifelse(((del_rep2 > threshold) & (del_rep1 > threshold)) | ((del_rep2 <= -threshold)& (del_rep1 <= -threshold)),1 ,0 )]
gseConv.exp[,del:= (del_rep1 + del_rep2)/2]
gseConv.exp[,donor.trt := ifelse(sel, ifelse(del > 0, max, min), donor)]
gseConv.exp[,hf.trt := ifelse(sel, ifelse(del > 0, max, min), hf)]
gseConv.exp[,donor.trt1 := donor+del_rep1] 
gseConv.exp[,donor.trt2 := donor+del_rep2] 
gseConv.exp[,hf.trt1 := hf+del_rep1] 
gseConv.exp[,hf.trt2 := hf+del_rep2] 
# gseConv.exp[,donor.trt:=ifelse(sel donor+(del_rep1 + del_rep2)/2]
# gseConv.exp[,hf.trt:=hf+(del_rep1 + del_rep2)/2]





tab1.skl = tab.skl
tab1.skl[,distil_id:=gsub(distil_id, pattern=":", replacement=".")]

####differential analysis

library(limma)
exp.mat  <- expression.match
design <- cbind(grp1=1, disease=expression.annot$disease.p) 
fit  <-  lmFit(exp.mat, design)
fit2 <- treat(fit,lfc=0.1)
fit2 <- eBayes(fit2) 
topGene <- topTreat(fit2,coef=2,2000)
diff.hf.gene = topGene$ID
m = length(diff.hf.gene)
n = nrow(exp.mat) - m

 #### pca analysis
	source("~/project/gtps/gtps/src/fisher.score.R")
	library(FactoMineR)
	expression.diff.hf = t(expression.match)
	gene.weight <- apply(expression.diff.hf,  2,function(tt) fisher.score(tt, expression.annot$disease) )
# see the google group discussion for more information https://groups.google.com/forum/#!topic/factominer-users/N8xzRCxCfQM



##### controls
geneNames = unique(tab.skl[pert_type =="trt_sh"]$pert_desc)

	ctlid = tab1.skl[pert_type == "ctl_untrt"]$distil_id
out = matrix(0, nrow=length(geneNames), ncol=3)
aadt = list()
bbdt = list()
max.change = list()
max.del.change = list()
for (gene in seq(length(geneNames))) {
	geneName = geneNames[gene]
	caseid = tab1.skl[pert_desc == geneName ]$distil_id
	exp.mat1  <- gseConv.exp[, c(ctlid, caseid)]
	design <- cbind(grp1=1, pert=c(rep(0, length(ctlid)), rep(1, length(caseid)))) 
	fit  <-  lmFit(exp.mat1, design)
	fit2 <- treat(fit,lfc=0.1)
	fit2 <- eBayes(fit2) 
	topGene1 <- topTreat(fit2,coef=2,p.value=.001, n=10000)
	diff.gene = topGene1$ID 
	x = sum(diff.gene %in% diff.hf.gene)
	k = length(diff.gene)
	p.val = dhyper(x, m , n, k)
	out[gene,] = c(x, k, p.val)
	print(c(gene, geneName))

	# out = data.table(out)
	# out$gene = geneNames
	# out = out[order(V3)]


# rownames(expression.join) = c(rownames(expression.diff.hf),
# "donor", "hf") 
	# paste0("donor", colnames(exp.donor.trt)), 
	# paste0("hf", colnames(exp.hf.trt)) )
	exp.mat2 = data.table(exp.mat1 )
	exp.mat2$ctl.med=apply(exp.mat2[,ctlid,with=F], 1, median)
	for (cas in caseid) 
	{
		eval(parse(text=paste0("exp.mat2[,",cas,".del:=",cas," - ctl.med]") ))

	}
	exp.df = exp.mat2 
	exp.df$sel = 0
# exp.df[rownames(expression.match) %in% topGene1$ID, sel:=1]
	exp.df[as.numeric(rownames(topGene1)), sel:=1]
	exp.df$del.med = apply(exp.df[,paste0(caseid,".del"),with=F],1,median)
	exp.df$donor.med = apply(expression.match[,expression.annot$disease.p==1],1 ,median)
	exp.df$hf.med = apply(expression.match[,expression.annot$disease.p==2],1 ,median)
	exp.df$max = apply(expression.match,1,max)
	exp.df$min = apply(expression.match,1,min)
	exp.df[,del:=ifelse(sel>0, del.med,0 )]
	exp.df[, `:=`(
		donor.sel_.01=donor.med + .01 *del, 
		donor.sel_.1=donor.med + .1 * del, 
		donor.sel_.5=donor.med + .5 * del, 
		donor.sel_1=donor.med + 1 * del, 
		donor.sel_2=donor.med + 2 * del 
		)]

	exp.df[, `:=`(
		hf.sel_.01=hf.med + .01 *del, 
		hf.sel_.1=hf.med + .1 *del, 
		hf.sel_.5=hf.med + .5 *del, 
		hf.sel_1=hf.med + 1 * del, 
		hf.sel_2=hf.med + 2 * del 
		)]


	exp.df[, `:=`(
		donor.thr_.01=
		ifelse( donor.sel_.01 > max, max, 
			ifelse( donor.sel_.01 < min, min, donor.sel_.01)), 
		donor.thr_.1=
		ifelse( donor.sel_.1 > max, max, 
			ifelse( donor.sel_.1 < min, min, donor.sel_.1)), 
		donor.thr_.5=
		ifelse( donor.sel_.5 > max, max, 
			ifelse( donor.sel_.5 < min, min, donor.sel_.5)), 
		donor.thr_1=
		ifelse( donor.sel_1 > max, max, 
			ifelse( donor.sel_1 < min, min, donor.sel_1)), 
		donor.thr_2=
		ifelse( donor.sel_2 > max, max, 
			ifelse( donor.sel_2 < min, min, donor.sel_2))
		)]
	exp.df[, `:=`(
		hf.thr_.01=
		ifelse( hf.sel_.01 > max, max, 
			ifelse( hf.sel_.01 < min, min, hf.sel_.01)), 
		hf.thr_.1=
		ifelse( hf.sel_.1 > max, max, 
			ifelse( hf.sel_.1 < min, min, hf.sel_.1)), 
		hf.thr_.5=
		ifelse( hf.sel_.5 > max, max, 
			ifelse( hf.sel_.5 < min, min, hf.sel_.5)), 
		hf.thr_1=
		ifelse( hf.sel_1 > max, max, 
			ifelse( hf.sel_1 < min, min, hf.sel_1)), 
		hf.thr_2=
		ifelse( hf.sel_2 > max, max, 
			ifelse( hf.sel_2 < min, min, hf.sel_2))
		)]

	exp.df[, `:=`(
		hf_.01=hf.med + .01 * del.med, 
		hf_.1=hf.med + .1 * del.med, 
		hf_.5=hf.med + .5 * del.med, 
		hf_1=hf.med + 1 * del.med, 
		hf_2=hf.med + 2 * del.med 
		)]

	exp.df[, `:=`(
		donor_.01=donor.med + .01 * del.med, 
		donor_.1=donor.med + .1 * del.med, 
		donor_.5=donor.med + .5 * del.med, 
		donor_1=donor.med + 1 * del.med, 
		donor_2=donor.med + 2 *del.med 
		)]


	exp.df[, `:=`(
		donor.thr.val_.01=
		ifelse( donor_.01 > max, max, 
			ifelse( donor_.01 < min, min, donor_.01)), 
		donor.thr.val_.1=
		ifelse( donor_.1 > max, max, 
			ifelse( donor_.1 < min, min, donor_.1)), 
		donor.thr.val_.5=
		ifelse( donor_.5 > max, max, 
			ifelse( donor_.5 < min, min, donor_.5)), 
		donor.thr.val_1=
		ifelse( donor_1 > max, max, 
			ifelse( donor_1 < min, min, donor_1)), 
		donor.thr.val_2=
		ifelse( donor_2 > max, max, 
			ifelse( donor_2 < min, min, donor_2))
		)]

	exp.df[, `:=`(
		hf.thr.val_.01=
		ifelse( hf_.01 > max, max, 
			ifelse( hf_.01 < min, min, hf_.01)), 
		hf.thr.val_.1=
		ifelse( hf_.1 > max, max, 
			ifelse( hf_.1 < min, min, hf_.1)), 
		hf.thr.val_.5=
		ifelse( hf_.5 > max, max, 
			ifelse( hf_.5 < min, min, hf_.5)), 
		hf.thr.val_1=
		ifelse( hf_1 > max, max, 
			ifelse( hf_1 < min, min, hf_1)), 
		hf.thr.val_2=
		ifelse( hf_2 > max, max, 
			ifelse( hf_2 < min, min, hf_2))
		)]

	exp.mat3 = t(as.matrix (exp.df[, grep( colnames(exp.df), pattern="donor|hf"), with=F]))



# exp.mat3 = t(as.matrix(exp.mat2))
# exp.mat3 = rbind(exp.mat3,del=exp.df$del)
	expression.join = rbind(expression.diff.hf, exp.mat3) 
	ind.sup = (nrow(expression.diff.hf) + 1): nrow(expression.join)
	pca.gene <- PCA(expression.join, scale.unit=T, ind.sup=ind.sup , col.w=gene.weight, ncp=ncol(expression.diff.hf), graph=F)
# pca.gene1 = prinp.comp(expression.join) 
	gene.proj.df <- pca.gene$ind.sup$coord
	eigen.gene <- colnames(pca.gene$ind.sup$coord)

	magnet.gene.df = as.data.table(pca.gene$ind$coord)
	donor.eigen.gene <- colMeans( magnet.gene.df[expression.annot$disease.p==1, eigen.gene, with=F])
	hf.eigen.gene <- colMeans( magnet.gene.df[expression.annot$disease.p==2, eigen.gene, with=F])
	expression.annot$tdi = apply(magnet.gene.df[,eigen.gene, with=F], 1 , function(tt) sum((tt - donor.eigen.gene)^2)) 
	expression.annot$hftdi = apply(magnet.gene.df[,eigen.gene, with=F], 1 , function(tt) sum((tt - hf.eigen.gene)^2)) 
	tdi.donor <- apply(gene.proj.df[,eigen.gene], 1 , function(tt) sum((tt - donor.eigen.gene)^2))
	tdi.hf <- apply(gene.proj.df[,eigen.gene], 1 , function(tt) sum((tt - hf.eigen.gene)^2))
	st =  seq(3,42, 5)
	he = rep(st,each=4) + 0:3 
	tdi.donor.diff = tdi.donor[he] - tdi.donor[ he + 1 ]
	tdi.hf.diff = tdi.hf[he] - tdi.hf[ he + 1 ]
    donor.diff =c(
    	 tdi.donor[1 ] - tdi.donor[c(3:7, 13:17, 28:37)],
    	 tdi.donor[2 ] - tdi.donor[c(8:12, 18:27, 38:42)])
    hf.diff =c(
    	 tdi.hf[1 ] - tdi.hf[c(3:7, 13:17, 28:37)],
    	 tdi.hf[2 ] - tdi.hf[c(8:12, 18:27, 38:42)])

	max.change[[geneName]] = max(c(donor.diff, hf.diff))
	max.del.change[[geneName]] = max(c(tdi.hf.diff, tdi.donor.diff))
	aa = data.table(cond=names(tdi.donor), tdi.donor, tdi.hf)
	aadt[[geneName]] = aa
	bb = data.table(cond=names(tdi.donor.diff), tdi.donor.diff, tdi.hf.diff)
	bbdt[[geneName]] = bb
	print(c(gene, geneName))
	if(gene %in% c(1,10,50, 100, 200, 300, 400, 500))
	save(file="intermediate.RData", max.change, max.del.change, aadt, bbdt)

}
# tdi.donor = tdi[1:ncol(exp.donor.trt)]
# tdi.hf =  tdi[-(1:ncol(exp.donor.trt))]
# gene.pca = data.table(t(gene.proj.df) )
# gene.pca$max = apply(magnet.gene.df, 2, max) 
# gene.pca$min = apply(magnet.gene.df, 2, min)

tab.erg = tab1[(pert_desc=="ERG")&(pert_type=="trt_sh")]
write.table(file="erg.cid", tab.erg$myid, quote=F, row.names=F, col.names=F)
tab.cntl = tab1[pert_type=="ctl_untrt"]
write.table(file="control.cid", tab.cntl$myid, quote=F, row.names=F, col.names=F)
###finding correlation


erg.mat = fread("lincs.mat.ERG") 
tab.erg.96 = tab.erg[pert_time==96]
erg.mat.96=erg.mat[,tab.erg$pert_time==96,with=F]
tissues = unique(tab.erg.96$cell_id)
erg.med = matrix(0,nrow=nrow(erg.mat.96), ncol=length(tissues)) 
colnames(erg.med) = tissues
for (inx in seq(length(tissues))) {
tissue = tissues[inx]
erg.med[,inx] = apply(erg.mat.96[,tab.erg.96$cell_id == tissue, with=F], 1, median)
}

ctl.mat = fread("lincs.mat.control") 
tissues1 = unique(tab.cntl$cell_id)
ctl.med = matrix(0,nrow=nrow(ctl.mat), ncol=length(tissues1)) 
colnames(ctl.med) = tissues1
for (inx in seq(length(tissues1))) {
	tissue = tissues1[inx]
	ctl.med[,inx] = apply(ctl.mat[,tab.cntl$cell_id == tissue, with=F], 1, median)
}
control.median.exp =ctl.med

save(file="control.median.exp.RData", control.median.exp)

erg.cor = cor(erg.med)
erg.cor = erg.cor[lower.tri(erg.cor)]

ctl.med.match = ctl.med[,colnames(erg.med)]
ctl.cor = cor(ctl.med.match)
ctl.cor = ctl.cor[lower.tri(ctl.cor)]

del.cor = cor(ctl.med.match - erg.med)
del.cor = del.cor[lower.tri(del.cor)]

dt.avi = data.table(erg.cor, ctl.cor, del.cor)

p = ggplot(data=dt.avi, aes(x=ctl.cor, y=erg.cor)) + geom_point()
p = ggplot(data=dt.avi, aes(x=ctl.cor, y=del.cor)) + geom_point()


### KRAS
tab.kras = tab1[(pert_desc=="KRAS")&(pert_type=="trt_sh")]
write.table(file="kras.cid", tab.kras$myid, quote=F, row.names=F, col.names=F)
###finding correlation


kras.mat = fread("lincs.mat.kras") 
tab.kras.96 = tab.kras[pert_time==96]
kras.mat.96=kras.mat[,tab.kras$pert_time==96,with=F]
tissues = unique(tab.kras.96$cell_id)
kras.med = matrix(0,nrow=nrow(kras.mat.96), ncol=length(tissues)) 
colnames(kras.med) = tissues
for (inx in seq(length(tissues))) {
tissue = tissues[inx]
kras.med[,inx] = apply(kras.mat.96[,tab.kras.96$cell_id == tissue, with=F], 1, median)
}

kras.cor = cor(kras.med)
kras.cor = kras.cor[lower.tri(kras.cor)]

ctl.med.match = ctl.med[,colnames(kras.med)]
ctl.cor = cor(ctl.med.match)
ctl.cor = ctl.cor[lower.tri(ctl.cor)]

del.cor = cor(ctl.med.match - kras.med)
del.cor = del.cor[lower.tri(del.cor)]

dt.avi = data.table(kras.cor, ctl.cor, del.cor)

p = ggplot(data=dt.avi, aes(x=ctl.cor, y=kras.cor)) + geom_point()
ggsave("kras.cntl.pdf",p)
p = ggplot(data=dt.avi, aes(x=ctl.cor, y=del.cor)) + geom_point()
ggsave("kras.del.pdf",p)


hf.donor.d = 1924

p = ggplot(data=expression.annot, aes(x=tdi, y=hftdi, color=etiology)) + geom_point(alpha=.2)+ 
geom_segment(aes(y = 4, x = hf.donor.d, yend = 923, xend = 1390), arrow = arrow())


