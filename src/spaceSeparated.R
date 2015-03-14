###### Transcriptomic space of donors #######

load("~/shortcuts/Data/MikeData/MAGnet_combat_eset.Rdata") 
eset.combat
expression <- exprs(eset.combat)
source("physioProcess.R")

# select genes based on those which are significantly correlated to predictors
# create space without weight
# create weighted space