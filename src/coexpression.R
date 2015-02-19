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