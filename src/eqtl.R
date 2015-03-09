library("avinash")
  load("~/shortcuts/Data/MikeData/MAGnet_combat_eset.Rdata")
  eset.combat
  expression <- exprs(eset.combat)
 load("~/shortcuts/Data/MikeData/magnet_PedTraits.RData")

 h.m.genes = fread("Heart_Muscle.txt", skip=1, header=F)$V2

 library(biomaRt)
 ensembl <- useMart("ensembl", dataset= "hsapiens_gene_ensembl")
 affy2ensembl = getBM(attributes = c("ensembl_gene_id", "affy_hugene_1_0_st_v1"), mart=ensembl)
 coefs <- lasso.eqtl(obj$x, obj$y, obj$mask2[,1], nthreads=2)

affy2ensembl.dt = data.table(affy2ensembl)
affy2ensembl.dt = affy2ensembl.dt[!is.na(affy_hugene_1_0_st_v1)]
setkey(affy2ensembl.dt, ensembl_gene_id)
affy2ensembl.dt = affy2ensembl.dt[h.m.genes]
affy2ensembl.dt = affy2ensembl.dt[!is.na(affy_hugene_1_0_st_v1)]
affy2ensembl.dt=affy2ensembl.dt[affy_hugene_1_0_st_v1 %in% rownames(expression)]
yy = t(expression[as.character(affy2ensembl.dt$affy_hugene_1_0_st_v1),])
SNPcol = grep(colnames(pedtrait), pattern="SNP")
xx = as.matrix(pedtrait[,SNPcol,with=F])
xx[is.na(xx)] <- 0
coefs <- lasso.new.eqtl(xx, yy)

#### control for variance
h.m.var = apply(yy,2 , var)
gene.var = apply(expression,1, var)
controlinx = sample.control(data=gene.var, size = length(h.m.var), target= h.m.var,"" replace = FALSE, prob = seq(0, 1, 0.01))
yy.control = tt(exrpression[controlinx,])
coefs.control <- lasso.new.eqtl(xx, yy.control)