#Donor 


phyiso.ace_inhibit <-  physio.failure[(( history_of_hypertension.p  ==1) & (ace_inhibitor.p ==1)) ]
phyiso.no_ace_inhibit <-  physio.failure[(( history_of_hypertension.p  ==1) & (ace_inhibitor.p ==0)) ]

which(pedtrait$sample_name %in% phyiso.ace_inhibit$sample_name)
expression.mat  <-  as.matrix(expression)
ace.ped <- expression.mat[which(pedtrait$sample_name %in% phyiso.ace_inhibit$sample_name),]
noace.ped <- expression.mat[which(pedtrait$sample_name %in% phyiso.no_ace_inhibit$sample_name),]
  noace.ped <- pedtrait[which(pedtrait$sample_name %in% phyiso.no_ace_inhibit$sample_name) ]

noace.ped.samp1  <- noace.ped[sample(nrow(noace.ped), .75 * nrow(noace.ped),replace=F),]
noace <-  colMeans(noace.ped.samp1) 
ace.ped.samp1  <- ace.ped[sample(nrow(ace.ped), .75 * nrow(ace.ped),replace=F),]
ace <-  colMeans(ace.ped.samp1)
ace.diff <- sort(noace - ace)
head(ace.diff)

