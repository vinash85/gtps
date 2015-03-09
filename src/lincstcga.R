library(avinash)
lincRNA = fread(file="~/project/tcga/data/gencode.v21.lincRNA.gtf")
lincRNA[,id:=seq(nrow(lincRNA))]
write.table(file="~/project/tcga/data/gencode.v21.lincRNA.bed", lincRNA[, list(V1, V4, V5, V7, id)], quote=F, sep="\t", row.names=F, col.names=F)
data.base.dir = "/cbcb/project2-scratch/vinash85/data/tcga/luad/24FebMatchUmMatch/"
output.base.dir ="~/project/tcga/result//26FebLincRNA/"
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
lRNA="~/project/tcga/data/gencode.v21.lincRNA.bed"
lRNA.tab=fread(lRNA, sep="\t")

dir.create("lRNA")
lRNA.tab[,id:=paste(V1,V2,V3,sep=":")]
lRNA = "gencode.v21.lincRNA.bed"
write.table(file=lRNA, x=lRNA.tab[,list(V1,V2,V3,V4,V5, id)], sep="\t", quote=F, row.names=F, col.names=F)
lRNA = "gencode.v21.lincRNA1.bed"
system("bedSort gencode.v21.lincRNA.bed gencode.v21.lincRNA1.bed")
cnv.mat = matrix(0, nrow=nrow(cnv.files), ncol=nrow(lRNA.tab),
dimnames=list(cnv.files$Sample, lRNA.tab$id))
seg.order = lRNA.tab$id
for (inx in seq(nrow(cnv.files))) {
	file.curr = paste0(data.base.dir, "CNV_SNP_Array/BI__Genome_Wide_SNP_6/Level_3/", cnv.files[inx]$File_Name)
	aa = fread(file.curr)
	aa[, Chromosome:=paste0("chr",Chromosome)]
	file1=paste0("./lRNA/", cnv.files[inx]$File_Name)
	write.table(file=file1, x=aa[,list(Chromosome, Start, End, Segment_Mean)], quote=F, row.names=F, col.names=F, sep="\t")
	file.out = paste0(file1, ".lRNA")
	command = paste("bedtools map -f .4 -b", file1, "-a", lRNA, " -null 0 -c 4 -o mean  >",file.out)
    system(command)
	print(inx)
	bb = fread(file.out, sep="\t", select=c(7))
	cnv.mat[inx,] = bb$V7 
}
save(file="cnv.mat.RData", cnv.mat)

write.table(file="lincRNA.cnv.mat", x=cnv.mat, sep="\t", quote=F, row.names=T, col.names=T)

load("../24FebMatchUnmatch//gene.mat.RData")


##### for RNASEQ V2 ################################
gene.files = files.list[grep(File_Name, pattern="genes.normalized_results")]	
gene.files = gene.files[!duplicated(gene.files$Sample)]
geneid = unlist(fread( paste0(data.base.dir, "RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/", gene.files[1]$File_Name), select=1))
entrezid = sapply( strsplit(geneid, "\\|" ), function(tt) tt[2])
hgncid = sapply( strsplit(geneid, "\\|" ), function(tt) tt[1])
library(org.Hs.eg.db)
x <- org.Hs.egCHRLOC
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes)]
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




#############matrix-eqtl########

gene.mat = gene.mat[ , apply(gene.mat,2,sd)  >0]
lRNA.sd = apply(cnv.mat,2 ,sd)
cnv.mat = cnv.mat[, lRNA.sd >0]
lRNA.tab.sub = lRNA.tab[ lRNA.sd >0]

dir.create("MATRIXEQTL")
snp= cbind( segid=colnames(cnv.mat), t(cnv.mat))
write.table(file=paste0("MATRIXEQTL/snp"), x = snp , row.names = F, col.names =T,  sep="\t", quote=F )
snp.loc = lRNA.tab.sub[,list(snp=id, chr=V1, pos=floor((V2+V3)/2))]
write.table(file=paste0("MATRIXEQTL/snp.loc"), x = snp.loc, row.names = F, col.names =T,  sep="\t", quote=F )


gene.inx =  colnames(gene.mat) %in% gene.annot$entrezid

gene.mat.sub = t(gene.mat[rownames(cnv.mat), gene.inx])
setkey(gene.annot, entrezid)
gene.annot.sub = (gene.annot[rownames(gene.mat.sub)] )


####################entrezid$$$$$$$$$$$$$
# gene.mat.sub = t(gene.mat)
stopifnot(identical(colnames(gene.mat.sub), rownames(cnv.mat))) 

write.table(file="MATRIXEQTL/exp.loc", x=gene.annot.sub[,list(entrezid, chr, start, end)], quote=F, 
		row.names=F, col.names=F, sep="\t")
exp = cbind( entrezid=entrezid.sub, gene.mat.sub)
write.table(file="MATRIXEQTL/exp", x=exp,
 row.names = F, col.names =T,  sep="\t", quote=F )

source("~/project/tcga/result/24febluad/MATRIXEQTL/sample.cis.r", echo=T) 

#####regressing out copy number variation of gene from gene-expression#####
luad.gene.scna = fread("~/project/tcga/data/luad/LUAD_SCNA.dat")
luad.gene.scna.mat =  as.matrix(luad.gene.scna[,-1,with=F])
rownames(luad.gene.scna.mat) = luad.gene.scna$genes
 

