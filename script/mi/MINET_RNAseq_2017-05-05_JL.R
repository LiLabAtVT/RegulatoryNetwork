## README: If you want to see how your input data works here,    you run 0->X--.->B->Y--.->D
#-----------------------------
# 0. Setup environment & load dataset
getwd();
workingDir = "/Users/jylee43/Google Drive/coursework_2017SPRING/ProblemSolving/4_ReadDataWork/Shamima";
setwd(workingDir); 
getwd();
fileList=list.files(workingDir, pattern="MeanFPKMGSE")
fileList

#-----------------------------
# 1. load dataset
#InputFileName="MeanFPKMGSE74692.csv"
InputFileName="MeanFPKMGSE80565.csv"
InputFileName=fileList[4]; print(InputFileName)
InputFileName=fileList[5]; print(InputFileName)
ExprData  <- read.table(InputFileName, sep = "," , header = T, na.strings ="", stringsAsFactors= F)
dim(ExprData)
head(ExprData)
exp.data = ExprData[, 2: (dim(ExprData)[2]) ]
rownames(exp.data) = ExprData[,1]

dim(exp.data)
head(exp.data)  

#-----------------------------
# 2. Filter & Get top 25% of high rowVar genes
#--------------------
library(matrixStats)
exp.data_rowvars = cbind(exp.data,rowvars=rowVars(as.matrix(exp.data)))
head(exp.data_rowvars)
head(exp.data_rowvars[order(-exp.data_rowvars[,'rowvars']),])
#topVarGenes <- head(exp.data_rowvars[order(-exp.data_rowvars[,'rowvars']),], (0.1*dim(exp.data_rowvars)[1]) )
topVarGenes <- head(exp.data_rowvars[order(-exp.data_rowvars[,'rowvars']),], 2000 )
head(topVarGenes)
tail(topVarGenes)
dim(topVarGenes)

#head(round(topVarGenes,2))
#topVarGenes =round(topVarGenes,2)

head(topVarGenes[,1: (dim(topVarGenes)[2]-1)])
topVarGenes=topVarGenes[,1: (dim(topVarGenes)[2]-1)]
head(topVarGenes)
tail(topVarGenes)
dim(topVarGenes)
#-----

#OutputFileName=paste0(gsub(".csv","",InputFileName),"_TopVar_", format(Sys.time(), "%Y%m%d_%I%p"),".csv")
OutputFileName=paste0(gsub(".csv","",InputFileName),"_TopVar.csv")
OutputFileName
write.csv(topVarGenes, file = OutputFileName,row.names=T)
#-----------------------------
exp.data=t(topVarGenes)
dim(exp.data)
head(exp.data)  
#-----------------------------
# B. Mutual Information Network. 
# For a given dataset, minet infers the network in two steps. First, the mutual information between all pairs of variables in dataset is computed according to the estimator argument. Then the algorithm given by method considers the estimated mutual informations in order to build the network. 
library(minet)

clr.net <- minet( exp.data, method="clr", estimator="spearman" )
head(clr.net)
arc.net <- minet( exp.data, method="aracne", estimator="spearman" )
head(arc.net)
mrn.net <- minet( exp.data, method="mrnet", estimator="spearman" )  
head(mrn.net)
mnb.net <- minet( exp.data, method="mrnetb", estimator="spearman" ) 
head(mnb.net)

#-----------------------------
# Z. Save nets into RData
# Save module colors and labels for use in subsequent parts
RDataFileName=paste0(gsub(".csv","",InputFileName),"_TopVar_Nets.RData")
RDataFileName
#save(clr.net, arc.net, mrn.net, mnb.net, file = RDataFileName)

#save(clr.net,clr.list, arc.net,arc.list, mrn.net,mrn.list, mnb.net,mnb.list, file = RDataFileName)

# Load network data saved in the second part.
lnames = load(file = RDataFileName);
#The variable lnames contains the names of loaded variables.
lnames

#-----------------------------
# Y. Write network matrix to CSV
write.csv(clr.net, file = paste0(gsub(".csv","",InputFileName),"_minet_clr_net.csv"),row.names=T)
write.csv(arc.net, file = paste0(gsub(".csv","",InputFileName),"_minet_arc_net.csv"),row.names=T)
write.csv(mrn.net, file = paste0(gsub(".csv","",InputFileName),"_minet_mrn_net.csv"),row.names=T)
write.csv(mnb.net, file = paste0(gsub(".csv","",InputFileName),"_minet_mnb_net.csv"),row.names=T)

# Y. Write network List to CSV
# Convert an adj MATRIX to an adj LIST.  
library(PCIT)

clr.list <- getEdgeList(as.matrix(clr.net))
dim(clr.net)	
dim(clr.list)
head(clr.list) 
tail(clr.list) 
head(clr.list[ order(-clr.list[,3]), ])

clr.list_sorted= clr.list[ order(-clr.list[,3]), ]
dim(clr.list_sorted)
#clr.list_sorted_filtered=clr.list_sorted[(clr.list_sorted[,3] >0.1), ]
#dim(clr.list_sorted_filtered)

write.table(format(clr.list_sorted,digit=4), file = paste0(gsub(".csv","",InputFileName),"_minet_clr_list.csv"),row.names=F, quote = FALSE)

#-----
arc.list <- getEdgeList(as.matrix(arc.net))
dim(arc.net)	
dim(arc.list)
head(arc.list) 

arc.list_sorted= arc.list[ order(-arc.list[,3]), ]
head(arc.list_sorted)
write.table(format(arc.list_sorted,digit=4), file = paste0(gsub(".csv","",InputFileName),"_minet_arc_list.csv"),row.names=F, quote = FALSE)
#-----
mrn.list <- getEdgeList(as.matrix(mrn.net))
dim(mrn.net)	
dim(mrn.list)
head(mrn.list) 

mrn.list_sorted= mrn.list[ order(-mrn.list[,3]), ]
head(mrn.list_sorted)
dim(mrn.list_sorted)
write.table(format(mrn.list_sorted,digit=4), file = paste0(gsub(".csv","",InputFileName),"_minet_mrn_list.csv"),row.names=F, quote = FALSE)

#-----
mnb.list <- getEdgeList(as.matrix(mnb.net))
dim(mnb.net)	
dim(mnb.list)
head(mnb.list) 

mnb.list_sorted= mnb.list[ order(-mnb.list[,3]), ]
head(mnb.list_sorted)
dim(mnb.list_sorted)
write.table(format(mnb.list_sorted,digit=4), file = paste0(gsub(".csv","",InputFileName),"_minet_mnb_list.csv"),row.names=F, quote = FALSE)


#-----------------------------
# C. validate: Inference Validation
# validate compares the infered network to the true underlying network for several threshold values and appends the resulting confusion matrices to the returned object.
# syn.net: SynTReN Source Network is the true underlying network used to generate the dataset loaded by data(exp.data) - see exp.data.

data(syn.net)  
clr.tbl <- validate( clr.net, syn.net ) 
head(clr.tbl)
arc.tbl <- validate( arc.net, syn.net ) 
head(arc.tbl)
mrn.tbl <- validate( mrn.net, syn.net ) 
head(mrn.tbl)
mnb.tbl <- validate( mnb.net, syn.net ) 
head(mnb.tbl)

max(fscores(validate( clr.net, syn.net )))
max(fscores(validate( arc.net, syn.net )))
max(fscores(validate( mrn.net, syn.net )))
max(fscores(validate( mnb.net, syn.net )))
#-----------------------------
# D. visualization
library(Rgraphviz)
#graph <- as(clr.net, "graphNEL")
graph <- as(arc.net, "graphNEL")
#graph <- as(mrn.net, "graphNEL")
#graph <- as(mnb.net, "graphNEL")
plot(graph,nodeAttrs=makeNodeAttrs(graph, fontsize=50))

#-----------------------------
# E. Plot PR-Curves max(fscores(mr.tbl)) 
dev <- show.pr(clr.tbl, col="red",type="b") 
dev <- show.pr(arc.tbl, device=dev, col="blue", type="b") 
dev <- show.pr(mrn.tbl, device=dev, col="green", type="b") 
show.pr(mnb.tbl, device=dev, col="black",type="b") 
auc.pr(clr.tbl)
auc.pr(arc.tbl)
auc.pr(mrn.tbl)
auc.pr(mnb.tbl)
