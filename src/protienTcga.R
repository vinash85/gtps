library(avinash)
data.base.dir = "/cbcb/project2-scratch/vinash85/data/tcga/luad/24FebMatchUmMatch/"
output.base.dir ="~/project/tcga/result/24FebMatchUnmatch/"
dir.create(output.base.dir)
setwd(output.base.dir)
files.list = fread(paste0(data.base.dir,"file_manifest.txt"))
header = colnames(files.list)
header.new = gsub(header, pattern=" ", replacement="_")
setnames(files.list, header, header.new)
cnv.files = files.list[grep(Platform_Type, pattern="^CNV")]	
cnv.files = cnv.files[grep(File_Name, pattern="hg19")]
cnv.files = cnv.files[!grepl(File_Name, pattern="nocnv")]
cnv.files = cnv.files[!duplicated(cnv.files$Sample)]
chromHMM="/cbcb/lab/sridhar/Avinash/data/chromHMM//wgEncodeBroadHmmNhlfHMM.bed"
chromHMM.tab=fread(chromHMM, sep="\t")
chromHMM="wgEncodeBroadHmmNhlfHMM.bed"
dir.create("")
chromHMM.tab[,id:=paste(V1,V2,V3,sep=":")]
write.table(file=chromHMM, x=chromHMM.tab, sep="\t", quote=F, row.names=F, col.names=F)
cnv.mat = matrix(0, nrow=nrow(cnv.files), ncol=nrow(chromHMM.tab),
dimnames=list(cnv.files$Sample, chromHMM.tab$id))
seg.order = chromHMM.tab$id
for (inx in seq(nrow(cnv.files))) {
	file.curr = paste0(data.base.dir, "CNV_SNP_Array/BI__Genome_Wide_SNP_6/Level_3/", cnv.files[inx]$File_Name)
	aa = fread(file.curr)
	aa[, Chromosome:=paste0("chr",Chromosome)]
	file1=paste0("./chromHMM/", cnv.files[inx]$File_Name)
	write.table(file=file1, x=aa[,list(Chromosome, Start, End, Segment_Mean)], quote=F, row.names=F, col.names=F, sep="\t")
	file.out = paste0(file1, ".chromHMM")
	command = paste("bedtools map -b", file1, "-a", chromHMM, " -null 0 -c 4 -o mean  >",file.out)
    system(command)
	print(inx)
	bb = fread(file.out, sep="\t", select=c(11))
	cnv.mat[inx,] = bb$V11 
}
save(file="cnv.mat.RData", cnv.mat)

##### for RNASEQ V1 #######################################
gene.files = files.list[grep(File_Name, pattern="gene.quantification")]	
gene.files = gene.files[!duplicated(gene.files$Sample)]
geneid = unlist(fread( paste0(data.base.dir, "RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/", gene.files[1]$File_Name), select=1))
entrezid = sapply( strsplit(geneid, "\\|" ), function(tt) tt[2])
library(org.Hs.eg.db)
x <- org.Hs.egCHRLOC
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
x <- org.Hs.egCHRLOCEND
mapped_genes <- mappedkeys(x)
xx.end <- as.list(x[mapped_genes])

sub.inx = which(entrezid %in% names(xx))
entrezid.sub = entrezid[sub.inx] 
xx1 = sapply( entrezid.sub, function(tt) { t1 = (xx[[tt]])[1];  c( names(t1), abs(t1) ,  sign(t1)) }) 
gene.annot = data.table(t(xx1))
setnames(gene.annot, 1:3, c("chr", "start", "strand"))
gene.annot$entrezid = entrezid.sub 
gene.annot$end =  sapply( entrezid.sub, function(tt) { t1 = (xx.end[[tt]])[1];   abs(t1)  })  
gene.annot[, chr:=paste0("chr",chr)]
gene.annot[,  strand:=ifelse(strand==1, "+", "-")]
gene.annot$entrezid = entrezid.sub

gene.mat = matrix(0, nrow=nrow(gene.files), ncol= length(entrezid),
dimnames=list(gene.files$Sample, entrezid ))
for (inx in seq(nrow(gene.files))) {
	file.curr = paste0(data.base.dir, "RNASeq/UNC__IlluminaHiSeq_RNASeq/Level_3/", gene.files[inx]$File_Name)
	aa = fread(file.curr, sep="\t", select=4)
	gene.mat[inx,] = aa$RPKM
	print(inx) 
}


##### for RNASEQ V2 ################################
gene.files = files.list[grep(File_Name, pattern="genes.normalized_results")]	
gene.files = gene.files[!duplicated(gene.files$Sample)]
geneid = unlist(fread( paste0(data.base.dir, "RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/", gene.files[1]$File_Name), select=1))
entrezid = sapply( strsplit(geneid, "\\|" ), function(tt) tt[2])
hgncid = sapply( strsplit(geneid, "\\|" ), function(tt) tt[1])
library(org.Hs.eg.db)
x <- org.Hs.egCHRLOC
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
x <- org.Hs.egCHRLOCEND
mapped_genes <- mappedkeys(x)
xx.end <- as.list(x[mapped_genes])

sub.inx = which(entrezid %in% names(xx))
entrezid.sub = entrezid[sub.inx] 
xx1 = sapply( entrezid.sub, function(tt) { t1 = (xx[[tt]])[1];  c( names(t1), abs(t1) ,  sign(t1)) }) 
gene.annot = data.table(t(xx1))
setnames(gene.annot, 1:3, c("chr", "start", "strand"))
gene.annot$entrezid = entrezid.sub 
gene.annot$end =  sapply( entrezid.sub, function(tt) { t1 = (xx.end[[tt]])[1];   abs(t1)  })  
gene.annot[, chr:=paste0("chr",chr)]
gene.annot[,  strand:=ifelse(strand==1, "+", "-")]
gene.annot$entrezid = entrezid.sub

gene.mat = matrix(0, nrow=nrow(gene.files), ncol= length(entrezid),
dimnames=list(gene.files$Sample, entrezid ))

for (inx in seq(nrow(gene.files))) {
	file.curr = paste0(data.base.dir, "RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/", gene.files[inx]$File_Name)
	aa = fread(file.curr, sep="\t", select=2)
	gene.mat[inx,] = aa$normalized_count
	print(inx) 
}
save(file="gene.mat.RData", gene.mat)






#############matrix-eqtl########

gene.mat = gene.mat[ , apply(gene.mat,2,sd)  >0]
chromHMM.sd = apply(cnv.mat,2 ,sd)
cnv.mat = cnv.mat[, chromHMM.sd >0]
chromHMM.tab.sub = chromHMM.tab[ chromHMM.sd >0]

dir.create("MATRIXEQTL")
snp= cbind( segid=colnames(cnv.mat), t(cnv.mat))
write.table(file=paste0("MATRIXEQTL/snp"), x = snp , row.names = F, col.names =T,  sep="\t", quote=F )
snp.loc = chromHMM.tab.sub[,list(snp=id, chr=V1, pos=floor((V2+V3)/2))]
write.table(file=paste0("MATRIXEQTL/snp.loc"), x = snp.loc, row.names = F, col.names =T,  sep="\t", quote=F )

 ####################entrezid###################
gene.mat.sub = t(gene.mat[rownames(cnv.mat), entrezid.sub])


####################entrezid$$$$$$$$$$$$$
# gene.mat.sub = t(gene.mat)
stopifnot(identical(colnames(gene.mat.sub), rownames(cnv.mat))) 

write.table(file="MATRIXEQTL/exp.loc", x=gene.annot[,list(entrezid, chr, start, end)], quote=F, 
		row.names=F, col.names=F, sep="\t")
exp = cbind( entrezid=entrezid.sub, gene.mat.sub)
source("~/project/tcga/result/24febluad/MATRIXEQTL/sample.cis.r", echo=T) 

#####regressing out copy number variation of gene from gene-expression#####
luad.gene.scna = fread("~/project/tcga/data/luad/LUAD_SCNA.dat")
luad.gene.scna.mat =  as.matrix(luad.gene.scna[,-1,with=F])
rownames(luad.gene.scna.mat) = luad.gene.scna$genes
 


 hgncid.sub = hgncid[sub.inx]
 gene.mat = gene.mat[, sub.inx]
 hgncid.sub1 = hgncid.sub[!(hgncid.sub=="?")]
gene.mat = gene.mat[,!(hgncid.sub=="?")]
colnames(gene.mat) =hgncid.sub1
gene.mat = gene.mat[,1:60]
# luad.gene.scna.mat = luad.gene.scna.mat[rownames(luad.gene.scna.mat) %in% colnames(gene.mat), ]

#### use this if genes need to be removed not needed usually######
gene.mat = gene.mat[ , apply(gene.mat,2,sd)  >0]
gene.mat.sub = gene.mat[,colnames(gene.mat) %in% rownames(luad.gene.scna.mat)]

luad.gene.t = t(luad.gene.scna.mat)
sampid = gsub(x=rownames(gene.mat.sub), pattern="-0[1,2]$", replacement="",  perl=T)
aa = t(luad.gene.t)

samp.inx = (!duplicated(sampid) & (sampid %in% colnames(aa)) )
bb = gene.mat.sub[samp.inx,]
xx= t(aa)
xx[is.na(xx)] <- 0
rownames(bb) = sampid[samp.inx]
bb = bb[colnames(aa),]
yy = log(bb+.001)
# yy = bb

yy= scale( bb,center=T, scale=apply(bb,2,sd))

library(glmnet)
library(doMC)
registerDoMC(cores=8)
yy.residue = yy
for (gene in seq(ncol(yy))) {
	yyi = yy[,gene]
	cvobj = cv.glmnet(x=xx, y=yyi, parallel=T, alpha=0.5, nfold=5 )
	# dt = data.table(y=bb[,gene], x=aa[,gene])
	# coef = lm(y~x, data=dt)
	# fit(dt$x)
	yy.residue[,gene] <-  yy[,gene] - predict(cvobj,newx=xx, s="lambda.min")
	print(gene)
}
option(error=recover)
gene.mat.back = gene.mat
gene.mat = yy.residue 

load("cnv.mat.RData")


sampid = gsub(x=rownames(cnv.mat), pattern="-0[1,2]$", replacement="",  perl=T)

samp.inx = (!(duplicated(sampid)) & (sampid %in% rownames(gene.mat)) )
cnv.mat = cnv.mat[samp.inx,]
 sampid.sub = sampid[samp.inx]
 rownames(cnv.mat) =sampid.sub
# rownames(cnv.mat) = rownames(gene.mat)
yy1 = gene.mat[rownames(cnv.mat),]


yy1 = gene.mat1[,1:10]
yy1.residue = yy1
xx1 = cnv.mat
for (gene in seq(ncol(yy1))) {
	yyi = yy1[,gene]
	cvobj = cv.glmnet(x=xx1, y=yyi, parallel=T, alpha=0.5, nfold=5 )
	# dt = data.table(y=bb[,gene], x=aa[,gene])
	# coef = lm(y~x, data=dt)
	# fit(dt$x)
	yy1.residue[,gene] <-  yy1[,gene] - predict(cvobj,newx=xx1, s="lambda.min")
	print(gene)
}



#####top associated genes ##########
load("../24febluad/MATRIXEQTL/matrixeqtl.RData")
aa = me.linear$trans$eqtls
temp = unique(aa$gene[1:10000])

gene.files = files.list[grep(File_Name, pattern="genes.normalized_results")]	
gene.files = gene.files[!duplicated(gene.files$Sample)]
geneid = unlist(fread( paste0(data.base.dir, "RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/", gene.files[1]$File_Name), select=1))
entrezid = sapply( strsplit(geneid, "\\|" ), function(tt) tt[2])
hgncid = sapply( strsplit(geneid, "\\|" ), function(tt) tt[1])
gene.sel = which(entrezid %in% temp)
hgncid1 = hgncid[gene.sel]
load("gene.mat.RData")
gene.mat = gene.mat[, as.character(temp)]
colnames(gene.mat) =hgncid1
gene.sd = apply(gene.mat,2,sd)





#############matrix-eqtl V3 ########

gene.mat = gene.mat[ , apply(gene.mat,2,sd)  >0]
chromHMM.sd = apply(cnv.mat,2 ,sd)
cnv.mat = cnv.mat[, chromHMM.sd >0]
chromHMM.tab.sub = chromHMM.tab[ chromHMM.sd >0]

dir.create("MATRIXEQTL")
snp= cbind( segid=colnames(cnv.mat), t(cnv.mat))
write.table(file=paste0("MATRIXEQTL/snp"), x = snp , row.names = F, col.names =T,  sep="\t", quote=F )
snp.loc = chromHMM.tab.sub[,list(snp=id, chr=V1, pos=floor((V2+V3)/2))]
write.table(file=paste0("MATRIXEQTL/snp.loc"), x = snp.loc, row.names = F, col.names =T,  sep="\t", quote=F )

gene.mat.sub = t(gene.mat)
stopifnot(identical(colnames(gene.mat.sub), rownames(cnv.mat))) 
id.mat = data.table(entrezid = entrezid, hgncid=hgncid)
setkey(id.mat, hgncid)
id.mat = id.mat[as.character(colnames(gene.mat))]
setkey(gene.annot, entrezid)
gene.annot1 = gene.annot[id.mat$entrezid,list(entrezid, chr, start, end)]
gene.annot1$hgncid = id.mat$hgncid
write.table(file="MATRIXEQTL/exp.loc", x=gene.annot1[,list(hgncid, chr, start, end)], quote=F, 
		row.names=F, col.names=F, sep="\t")
exp = cbind( entrezid=rownames(gene.mat.sub), gene.mat.sub)
write.table(file="MATRIXEQTL/exp", x=exp,
 row.names = F, col.names =T,  sep="\t", quote=F )
source("~/project/tcga/result/24febluad/MATRIXEQTL/sample.cis.r", echo=T) 

#####regressing out copy number variation of gene from gene-expression#####
luad.gene.scna = fread("~/project/tcga/data/luad/LUAD_SCNA.dat")
luad.gene.scna.mat =  as.matrix(luad.gene.scna[,-1,with=F])
rownames(luad.gene.scna.mat) = luad.gene.scna$genes




#####oncogene and tumor supressor gene#######
oncogene = fread("~/project/tcga/data/oncogenic_drivers.csv")
load("gene.mat.RData")
gene.mat = log(gene.mat + 1)
output.base.dir ="~/project/tcga/result/24FebMatchUnmatch/oncogene"
gene.sd = apply(gene.mat,2,sd)

oncogene = oncogene[gene %in% hgncid]
setnames(oncogene, 4, "tumor_suppressor")
oncogene[,label:=ifelse(oncogene > tumor_suppressor, "oncogene", "tumor_suppressor")]
gene.sel = which(hgncid !="?")

hgncid1 = hgncid[gene.sel]
load("gene.mat.RData")
gene.mat = gene.mat[,gene.sel]
colnames(gene.mat) =hgncid1
gene.sd = apply(gene.mat,2,sd)
notonco = which(!(hgncid1 %in% oncogene$gene))
control.inx = avinash::sample.control(target=gene.sd[oncogene$gene], data= gene.sd[!(hgncid1 %in% oncogene$gene)], size=500)

genetab = data.table( gene=c(oncogene$gene, hgncid1[c(notonco[control.inx], sample(notonco, 500))]),  
	label = c(oncogene$label, rep("control",500), rep("random", 500)))
gene.mat = gene.mat[,genetab$gene]
