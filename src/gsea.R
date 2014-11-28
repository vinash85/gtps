setkey(pedtrait, sample_name)
gsea.xx <- pedtrait[physio.pedtrait$sample_name, grep(x=colnames(pedtrait), pattern="^X"), with=F ]
gsea.xx <- t(gsea.xx)
colnames(gsea.xx) <-  pedtrait$sample_name
gsea.df <- cbind(NAME=gsub(rownames(gsea.xx), pattern="^X",replacement=""), Description=NA, data.frame(gsea.xx)) 
#http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
# gene expression format
gct.file <- file("magnet.gct", "w")
cat("#1.2\n" , file=gct.file)
cat(nrow(gsea.xx) , "\t", ncol(gsea.xx), "\n",  file=gct.file)
write.table(file=gct.file,x = gsea.df , row.names = F, 
	col.names =T,  sep="\t", quote=F, append=T )
close(gct.file)
# phenotype format
cls.file <- file("magnet.cls", "w")
cat(ncol(gsea.xx), "\t2\t1\n" , file=cls.file)
cat("# 1=donor 2=heartFailure\n" , file=cls.file)
cat(physio.pedtrait$disease.ped, sep="\t",  file=cls.file)
close(cls.file)

# isthemic
physio.ishemic  <- physio.pedtrait[etiology %in% c("donor", "Ischemic")]
gsea.xx <-pedtrait[physio.ishemic$sample_name, grep(x=colnames(pedtrait), pattern="^X"), with=F ]
gsea.xx <- t(gsea.xx)
colnames(gsea.xx) <-  physio.ishemic$sample_name
gsea.df <- cbind(NAME=gsub(rownames(gsea.xx), pattern="^X",replacement=""), Description=NA, data.frame(gsea.xx)) 
#http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
# gene expression format
gct.file <- file("ishemic.gct", "w")
cat("#1.2\n" , file=gct.file)
cat(nrow(gsea.xx) , "\t", ncol(gsea.xx), "\n",  file=gct.file)
write.table(file=gct.file,x = gsea.df , row.names = F, 
	col.names =T,  sep="\t", quote=F, append=T )
close(gct.file)
# phenotype format
cls.file <- file("ishemic.cls", "w")
cat(ncol(gsea.xx), "\t2\t1\n" , file=cls.file)
cat("# 1=donor 2=heartFailure\n" , file=cls.file)
cat(physio.ishemic$disease.ped, sep="\t",  file=cls.file)
close(cls.file)


physio.idiopathic  <- physio.pedtrait[etiology %in% c("donor", "Idiopathic Dilated CMP")]
gsea.xx <-pedtrait[physio.idiopathic$sample_name, grep(x=colnames(pedtrait), pattern="^X"), with=F ]
gsea.xx <- t(gsea.xx)
colnames(gsea.xx) <-  physio.idiopathic$sample_name
gsea.df <- cbind(NAME=gsub(rownames(gsea.xx), pattern="^X",replacement=""), Description=NA, data.frame(gsea.xx)) 
#http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
# gene expression format
gct.file <- file("idiopathic.gct", "w")
cat("#1.2\n" , file=gct.file)
cat(nrow(gsea.xx) , "\t", ncol(gsea.xx), "\n",  file=gct.file)
write.table(file=gct.file,x = gsea.df , row.names = F, 
	col.names =T,  sep="\t", quote=F, append=T )
close(gct.file)
# phenotype format
cls.file <- file("idiopathic.cls", "w")
cat(ncol(gsea.xx), "\t2\t1\n" , file=cls.file)
cat("# 1=donor 2=heartFailureIdiopathic\n" , file=cls.file)
cat(physio.idiopathic$disease.ped, sep="\t",  file=cls.file)
close(cls.file)

#diabetes
physio.diabetes  <- physio.pedtrait[ etiology %in% c("donor")] 
physio.diabetes  <- physio.diabetes[ !is.na(history_of_diabetes.p)] 
gsea.xx <-pedtrait[physio.diabetes$sample_name, grep(x=colnames(pedtrait), pattern="^X"), with=F ]
gsea.xx <- t(gsea.xx)
colnames(gsea.xx) <-  physio.diabetes$sample_name
gsea.df <- cbind(NAME=gsub(rownames(gsea.xx), pattern="^X",replacement=""), Description=NA, data.frame(gsea.xx)) 
#http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
# gene expression format
gct.file <- file("diabetes.gct", "w")
cat("#1.2\n" , file=gct.file)
cat(nrow(gsea.xx) , "\t", ncol(gsea.xx), "\n",  file=gct.file)
write.table(file=gct.file,x = gsea.df , row.names = F, 
	col.names =T,  sep="\t", quote=F, append=T )
close(gct.file)
# phenotype format
cls.file <- file("diabetes.cls", "w")
cat(ncol(gsea.xx), "\t2\t1\n" , file=cls.file)
cat("# 0=nodiabetes 1=diabetes\n" , file=cls.file)
cat(physio.diabetes$history_of_diabetes.p, sep="\t",  file=cls.file)
close(cls.file)


#hypertension
physio.hypertension  <- physio.pedtrait[ etiology %in% c("donor")] 
physio.hypertension  <- physio.hypertension[ !is.na(history_of_hypertension.p)] 
gsea.xx <-pedtrait[physio.hypertension$sample_name, grep(x=colnames(pedtrait), pattern="^X"), with=F ]
gsea.xx <- t(gsea.xx)
colnames(gsea.xx) <-  physio.hypertension$sample_name
gsea.df <- cbind(NAME=gsub(rownames(gsea.xx), pattern="^X",replacement=""), Description=NA, data.frame(gsea.xx)) 
#http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
# gene expression format
gct.file <- file("hypertension.gct", "w")
cat("#1.2\n" , file=gct.file)
cat(nrow(gsea.xx) , "\t", ncol(gsea.xx), "\n",  file=gct.file)
write.table(file=gct.file,x = gsea.df , row.names = F, 
	col.names =T,  sep="\t", quote=F, append=T )
close(gct.file)
# phenotype format
cls.file <- file("hypertension.cls", "w")
cat(ncol(gsea.xx), "\t2\t1\n" , file=cls.file)
cat("# 0=nohypertension 1=hypertension\n" , file=cls.file)
cat(physio.hypertension$history_of_hypertension.p, sep="\t",  file=cls.file)
close(cls.file)

