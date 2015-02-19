
#############################################################################
## (1.1) Import all required functions with the following source() command ##
#############################################################################
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/GOHyperGAll.txt")

###########################################################
## (2.1) Using gene-to-GO mappings from geneontology.org ##
###########################################################
readGOorg(myfile = "gene_association.goa_ref_human", colno = c(5,11,9), org = "Human"); gene2GOlist(rootUK=T)
readGOorg(myfile = "gene_association.tair", colno = c(5,11,9), org = "Arabidopsis"); gene2GOlist(rootUK=T)
   # Download the required annotation table from geneontology.org and unzip it. Then point the 'readGOorg()' function to
   # this file name. The two functions generate 4 data frames with the assigned gene-to-GO mappings and 3 lists containing
   # the gene-to-GO-OFFSPRING associations. When the processes are completed, 6 files will be saved in your working directory!
   # They can be reloaded in future R sessions with the 'load' command below. If the argument 'rootUK' is set to TRUE, then
   # the root nodes are treated as terminal nodes to account for the new assignment of unknown genes to the root nodes.

###############################################
## (2.2) Using gene-to-GO mappings from Bioc ##
###############################################
## Note: users should execute either step (2.1) or (2.2), but not both!
sampleDFgene2GO(lib="ath1121501.db"); gene2GOlist(rootUK=T)
   # Similar as above, but the gene-to-GO mappings are obtained from BioC. The generated 4 sample data frame and 3
   # list objects can be reloaded in future R sessions with the 'load' command below.

######################################################################
## (2.3) Obtain AffyID-to-GeneID mappings when working with AffyIDs ##
######################################################################
AffyID2GeneID(map = "ftp://ftp.arabidopsis.org/home/tair/Microarrays/Affymetrix/affy_ATH1_array_elements-2010-12-20.txt")
   # When working with AffyIDs, this function creates a AffyID-to-GeneID mapping data frame using by default the TAIR
   # mappings for the Arabidopsis ATH1 chip. To use the function for th emappings of other chips, one needs to create the
   # corresponding decoding data frame 'affy2locusDF'.

############################################################
## (3.1) Reloading required data objects from local files ##
############################################################
loadData(); load(file="MF_node_affy_list"); load(file="BP_node_affy_list"); load(file="CC_node_affy_list")
   # This step makes future sessions much faster, since it allows to skip the previous data generation steps (2.1-2.3).
   # A sample data set is available here: ArabSampleGOHyperGAll.zip (Jan 2010).

##########################################
## (3.2) Obtain a sample set of GeneIDs ##
##########################################
test_sample <- unique(as.vector(GO_MF_DF[1:40,2])) # When working with GeneIDs.
test_sample <- AffyID2GeneID(affyIDs=affy_sample, probe2gene=1)
   # When working with AffyIDs, one can use the function 'AffyID2GeneID' to obtain for a set of AffyIDs their corresponding
   # GeneIDs from the data frame 'affy2locusDF' (see above). For probe sets that match several loci, only the first locus
   # ID will be used if the argument 'probe2gene' is set to 1. To demo the function, one can use the following sample
affy_sample <- c("266592_at", "266703_at", "266199_at", "246949_at", "267370_at", "267115_s_at", "266489_at", "259845_at", "266295_at", "262632_at") # AffyID sample

##########################################################################
## (4.1) Perform phyper test, goSlim subsetting and plotting of results ##
##########################################################################
GOHyperGAll_result <- GOHyperGAll(gocat="MF", sample=test_sample, Nannot=2); GOHyperGAll_result[1:10,-8]
   # The function 'GOHyperGAll()' performs the phyper test against all nodes in the GO network. It returns raw and
   # adjusted p-values. The Bonferroni correction is used as p-values adjustment method according to Boyle et al, 2004 (online).
   # The argument 'Nannot' defines the minimum number of direct annotations per GO node from the sample set to determine
   # the number of tested hypotheses for the p-value adjustment. The argument 'gocat' can be assigned the values "MF", "BP"
   # and "CC". Omitting the '-8' delimiter will provide the sample keys matching at every GO node.
subset <- GOHyperGAll_Subset(GOHyperGAll_result, sample=test_sample, type="goSlim"); subset[,-8]
   # The function 'GOHyperGAll_Subset()' subsets the GOHyperGAll results by assigned GO nodes or custom goSlim categories.
   # The argument 'type' can be assigned the values "goSlim" or "assigned". The optional argument 'myslimv' can be used to
   # provide a custom goSlim vector. Omitting the '-8' delimiter will show the sample keys matching at every GO node.
pie(subset[subset$SampleMatch>0 ,3], labels=as.vector(subset[subset$SampleMatch>0 ,1]), main=unique(as.vector(subset[subset$SampleMatch>0, 7]))) # Plots pie chart of subsetted results.

##############################################################
## (4.2) Reduce GO Term Redundancy in 'GOHyperGAll_results' ##
##############################################################
simplifyDF <- GOHyperGAll_Simplify(GOHyperGAll_result, gocat="MF", cutoff=0.001, correct=T)
   # The result data frame 'GOHyperGAll_result' often contains several connected GO terms with significant scores which
   # can complicate the interpretation of large sample sets. To reduce the redundancy, the function 'GOHyperGAll_Simplify'
   # subsets the data frame 'GOHyperGAll_result' by a user specified p-value cutoff and removes from it all GO nodes with
   # overlapping children sets (OFFSPRING), while the best scoring nodes remain in the data frame.
data.frame(GOHyperGAll_result[GOHyperGAll_result[,1] %in% simplifyDF[,1], -8], GO_OL_Match=simplifyDF[,2])
   # This command returns the redundancy reduced data set. The column 'GO_OL_Match' provides the number of accessions that
   # match the connected nodes.

################################################
## (4.3) Batch Analysis of Many Gene Clusters ##
################################################
BatchResult <- GOCluster_Report(CL_DF=CL_DF, method="all", id_type="gene", CLSZ=10, cutoff=0.001, gocats=c("MF", "BP", "CC"), recordSpecGO=c("GO:0003674", "GO:0008150", "GO:0005575"))
   # The function 'GOCluster_Report' performs the three GO analyses in batch mode: 'GOHyperGAll', 'GOHyperGAll_Subset'
   # or 'GOHyperGAll_Simplify'. It processes many groups of genes (e.g. gene expression clusters) and returns the results
   # conveniently organized in a single data frame. The gene sets need to be provided in a data frame of the format
   # specified at the end of the GOHyperGAll script. CLSZ: minimum cluster size to consider; method: "all", "slim" or
   # "simplify"; gocat: "MF", "BP" or "CC"; cutoff: adjusted p-value cutoff; recordSpecGO: argument to include one specific
   # GOID in each of the 3 ontologies, e.g: recordSpecGO=c("GO:0003674", "GO:0008150", "GO:0005575").
