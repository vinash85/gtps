library(avinash)
data.base.dir = "/cbcb/project2-scratch/vinash85/data/tcga/luad/24FebMatchUmMatch/" // directory to store data
output.base.dir ="~/project/tcga/result/24FebMatchUnmatch/" 
dir.create(output.base.dir)
setwd(output.base.dir)
files.list = fread(paste0(data.base.dir,"file_manifest.txt"))
header = colnames(files.list)
header.new = gsub(header, pattern=" ", replacement="_")
setnames(files.list, header, header.new)

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

