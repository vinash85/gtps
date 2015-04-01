###### Transcriptomic space of donors #######

load("~/shortcuts/Data/MikeData/MAGnet_combat_eset.Rdata") 
eset.combat
expression <- exprs(eset.combat)
source("~/project/gtps/gtps/src/physioProcess.R")
physio.reduced <- physio.subset[,list(sample_name, site.p, age.p, gender.p, race.p, ethnicity.p, patient_weight_kg.p, height_cm.p, heart_weight_grams.p, prior_cabg.p,
				   history_of_afib_aflutter.p, history_of_diabetes.p, history_of_hypertension.p, ace_inhibitor.p, lvedd.p, lvesd.p, lvef.p,
				   mitral_regurgitation.p, tricuspid_regurgitation.p, hemos_rap.p, hemos_systolic.p,hemos_pwcp.p,cardiac_index.p, creatinine_level.p,
				   disease.p, etiology)] 
setkey(physio.reduced, sample_name)



library(missMDA)
physio.reduced.df <- as.data.frame(physio.subset[,list( age.p, gender.p, race.p,  patient_weight_kg.p, height_cm.p, heart_weight_grams.p, prior_cabg.p,
				   history_of_afib_aflutter.p, history_of_diabetes.p, history_of_hypertension.p, ace_inhibitor.p, lvedd.p, lvesd.p, lvef.p,
				   mitral_regurgitation.p, tricuspid_regurgitation.p, cardiac_index.p, creatinine_level.p,disease.p, etiology.Ischemic.p, etiology.Idiopathic.p) ])
imp <-  imputePCA(physio.reduced.df, ncp=10)
physio.imputed  <-  as.data.table(imp$completeObs)
physio.imputed$sample_name = physio.subset$sample_name
physio.imputed$etiology =   physio.subset$etiology

# expression.merge <- merge(x=physio.imputed, y=gene.proj.df, by = "sample_name", sort = F, suffixes= c("", ".gene"))
setkey(physio.imputed, sample_name)
expression = expression[,colnames(expression) %in% physio.subset$sample_name]
expression.m <- sweep(expression, 1, rowMeans(expression), FUN="-")
physio.pedtrait = physio.imputed[colnames(expression)]
# physio.pedtrait$etiology = physio.reduced[physio.pedtrait$sample_name]$etiology


#### GPDI

physio.imputed1 <- physio.imputed[!duplicated(sample_name)]
physio.imputed1 = physio.imputed1
physio.imputed.params <-  as.matrix( physio.imputed1[,list( lvedd.p, lvesd.p, lvef.p,
				   mitral_regurgitation.p, tricuspid_regurgitation.p, cardiac_index.p, creatinine_level.p)])

rownames(physio.imputed.params)  <-  physio.imputed1$sample_name
source("~/project/gtps/gtps/src/fisher.score.R")
physio.weight <- apply(physio.imputed.params,  2,function(tt) fisher.score(tt, physio.imputed1$disease.p) )
gene.weight <- apply(expression ,  1,function(tt) fisher.score(tt, physio.imputed1[colnames(expression)]$disease.p) )
	pca.gene <- PCA(expression.curr.sel, scale.unit=F, 
pca.physio <- PCA(physio.imputed.params, 
	# col.w=physio.weight, 
	ncp=ncol(physio.imputed.params), graph=F)
physio.proj.df <- as.data.table(pca.physio$ind$coord)
eigen.physio <- colnames(pca.physio$ind$coord)
eigen.physio <-  paste0(eigen.physio, ".p")
setnames(physio.proj.df, 1:ncol(physio.proj.df), eigen.physio)
physio.proj.df$sample_name  <-  rownames(pca.physio$ind$coord)
physio.proj.df$etiology  <-  physio.imputed1$etiology
physio.proj.df$disease.p <-  physio.imputed1$disease.p
setkey(physio.proj.df, sample_name)
donor.eigen.physio <- colMeans( physio.proj.df[disease.p==1, eigen.physio, with=F])
physio.proj.df$gpdi <- apply(physio.proj.df[,eigen.physio, with=F], 1 , function(tt) sum((tt - donor.eigen.physio)^2))



########## select genes based on those which are significantly correlated to phenotype
params = c("lvedd.p", "lvesd.p", "lvef.p", "mitral_regurgitation.p", "tricuspid_regurgitation.p", "cardiac_index.p", "creatinine_level.p")
gene.sel = NULL

physio.curr = physio.pedtrait[disease.p==2]
expression.curr = expression[,physio.curr$sample_name]
for (param in params) {
	trait = unlist(physio.curr[,colnames(physio.curr) == param, with =F])
aa = apply(expression.curr,1, function(tt) (cor.test(tt, trait,method="spearman"))$p.value)
gene.sel[[param]] = as.character(names(aa)[ aa<1E-5])
}

gene.sel = unique(unlist(gene.sel))

expression.curr.sel = t(expression.curr[gene.sel,])
pca.gene <- PCA(expression.curr.sel,
		# ncp=ncol(expression.curr.sel),
		ncp=10,
		col.w=gene.weight[colnames(expression.curr.sel)],
		 graph=F)
# pca.gene1 = prinp.comp(expression.curr.join) 
gene.proj.df <- as.data.table(pca.gene$ind$coord)
eigen.gene <- colnames(pca.gene$ind$coord)
gene.proj.df$sample_name  <-  rownames(pca.gene$ind$coord)
setkey(gene.proj.df, sample_name)
expression.curr.merge <- merge(x=physio.proj.df, y=gene.proj.df, by = "sample_name", sort = F, suffixes = c("", ".gene") )
expression.curr.merge <- merge(x=expression.curr.merge, y=physio.pedtrait, by = "sample_name", sort = F, suffixes = c("", ".ped") )
donor.eigen.gene <- colMeans( expression.curr.merge[disease.p==2, eigen.gene, with=F])
eigen.gene.m <- paste0(eigen.gene,".gene")
expression.curr.merge$tdi <- apply(expression.curr.merge[,eigen.gene, with=F], 1 , function(tt) sum((tt - donor.eigen.gene)^2))
expression.curr.merge$gpdi <- physio.proj.df[expression.curr.merge$sample_name]$gpdi
expression.curr.merge$etiology <- physio.imputed[as.character(expression.curr.merge$sample_name)]$etiology

curr.name <- "GPDI"
curr.param <- "gpdi"

aa <-  cor.test(expression.curr.merge$gpdi, expression.curr.merge$tdi, method="spearman")
p <- ggplot(data=expression.curr.merge, aes(x=rank(tdi), y=rank(gpdi)))
p <- p+geom_point(aes(col=factor(etiology)) ) + geom_smooth(method=lm) +
labs(y="Global phenotypic deviation index (Rank) ", x = "Transcriptomic Deviation index (Rank)",   title=paste0(curr.name , "\n correlation =",  format( aa$estimate, digits=3) , " \n p =", format(aa$p.value, digits=3, scientific=T)))   


ggsave(paste0( "tdi_all_imputed_",curr.name, ".value.pdi.pdf"), p)

# create space without weight
# create weighted space