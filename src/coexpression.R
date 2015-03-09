# coexpresion
library(avinash)
setwd("/cbcb/project2-scratch/vinash85/project/1000gene/data/tissues")
load("/cbcb/project2-scratch/vinash85/project/1000gene/data/tissues/gtex.pc.RData")
load("/cbcb/project2-scratch/vinash85/project/1000gene/data/tissues/sample.info.RData")
gtex.pc = gtex.pc[!duplicated(geneSym)]
gtex.mat = as.matrix(gtex.pc[,dt.merge$SAMPID,with=F])
# gtex.mat= gtex.mat[! duplicated(gtex.pc$geneSym), ]
rownames(gtex.mat) = gtex.pc$geneSym
# dt.merge = dt.merge[!duplicated(gtex.pc$geneSym)] 
# na.inx = is.na(gtex.mat[,1])

tissues.type=dt.merge$SMTS
# rownames(gtex.mat) = gtex.pc # Lung kidney prostate 
tissues.sub = c("Lung", "Kidney", "Prostate")
tissues=sort(unique(dt.merge$SMTS))
library(WGCNA)
library(reshape)
tissues.corr = list()
geneNames = rownames(gtex.mat)
for (ii in seq(length(tissues.sub)) ) {
	tissues.curr = tissues.sub[ii]
    gtex.curr = t(gtex.mat[,tissues.type == tissues.curr ])
    tissues.corr = bicor(gtex.curr, nThreads = 8)
    # tissues.corr = tissues.corr[1:100,1:100]
    bb.temp = which( tissues.corr > 0.7 , arr.ind=T)

    bb  = as.data.table(bb.temp)
    bb$bicor <- tissues.corr[bb.temp]
    bb[,gene1 := geneNames[row]]
    bb[,gene2 := geneNames[col]]
    # bb = bb[!duplicated(paste0(gene1, gene2))]
    bb = bb[row > col] 
    # bb = bb[!(gene1 == gene2)]
    bb = bb[,list(gene1, gene2, bicor)]
    write.table(file=paste0(tissues.curr,  "pos.coexp.txt"), x = bb,row.names=F, col.names=T, quote=F)

    bb.temp = which( tissues.corr < -0.7 , arr.ind=T) 
    bb  = as.data.table(bb.temp)
    bb$bicor <- tissues.corr[bb.temp]
    bb[,gene1 := geneNames[row]]
    bb[,gene2 := geneNames[col]]
    # bb = bb[!duplicated(paste0(gene1, gene2))]
    bb = bb[row > col] 
    # bb = bb[!(gene1 == gene2)]
    bb = bb[,list(gene1, gene2, bicor)]
    write.table(file=paste0(tissues.curr,  "neg.coexp.txt"), x = bb,row.names=F, col.names=T, quote=F)

    bb.temp = which( (abs(tissues.corr) < .001)  , arr.ind=T) 
    bb  = as.data.table(bb.temp)
    bb$bicor <- tissues.corr[bb.temp]
    bb[,gene1 := geneNames[row]]
    bb[,gene2 := geneNames[col]]
    # bb = bb[!duplicated(paste0(gene1, gene2))]
    bb = bb[row > col] 
    # bb = bb[!(gene1 == gene2)]
    bb = bb[,list(gene1, gene2, bicor)]
    write.table(file=paste0(tissues.curr,  "zero.coexp.txt"), x = bb,row.names=F, col.names=T, quote=F)
    # diag(tissues.corr[[ii]]) <- 0
     # temp <- lower.tri(tissues.corr[[i]], diag = FALSE)
     # tissues.corr[[ii]] [ is.na(tissues.corr[[ii]])] < - 0 
}
    map = fread("UNIQUE_GENES.txt", header=F)
    map$V1 = as.character(map$V1)
    setkey(map, V1)
for (ii in seq(length(tissues.sub)) ) {
    stopifnot(FALSE)
    aa = fread("IMR90_100KB_Gene_Gene.txt")
    # affy2ensembl1 = data.table(affy2ensembl1)
    # affy2ensembl = affy2ensembl1[!is.na(entrezgene)]
    # affy2ensembl[,entrezgene:=as.character(entrezgene)]
    # setkey(affy2ensembl, entrezgene)
    # aa$V1hgnc=affy2ensembl[as.character(aa$V1)]$hgnc_symbol
    tissues.curr = tissues.sub[ii]
        gtex.curr = t(gtex.mat[,tissues.type == tissues.curr ])
    tissues.corr = bicor(gtex.curr, nThreads = 8)
    aa$V1hg = map[as.character(aa$V1)]$V2
    aa$V2hg = map[as.character(aa$V2)]$V2
    aa = aa[(V1hg %in% colnames(tissues.corr)) & (V2hg %in% colnames(tissues.corr))]
    geneInx = data.table(colnames(tissues.corr), 1:nrow(tissues.corr))
    setkey(geneInx, V1)
    aa$V1Inx = geneInx[as.character(aa$V1hg)]$V2
    aa$V2Inx = geneInx[as.character(aa$V2hg)]$V2
    aa2 = as.matrix(aa[,list(V1Inx, V2Inx)])
    aa$bicor = tissues.corr[aa2]
    aa[, bicor:=ifelse(is.na(bicor), 0, bicor )]
    write.table(file=paste0(tissues.curr,  "pair.txt"), x = aa,row.names=F, col.names=T, quote=F)
}