
# sample code for performing the network analysis.
# Partial Correlation method
# required R library: GeneNet
# sample input data: GSE10670_ave_TopVar.csv
# sample output data: pcor_edgelist.csv
#
# How to run the program:
# Assume working direction contain a data folder and a script folder.
# 
# Rscript rf.R [ExpressionMatrix.csv] [InteractionEdges.csv] [number of edges] [OutputFile.csv]
# 
# example
# 
# Rscript ./script/randomforest/rf.R ./data/GSE10670_ave_TopVar.csv ./data/golddata.csv nedges ./result/rf.csv
#
# or simply
# 
# Rscript ./script/randomforest/rf.R ./data/GSE10670_ave_TopVar.csv
# 

 

setwd("/Users/jylee43/Google Drive/coursework_2017SPRING/ProblemSolving/5_2017Summer")
emfn = "./data/GSE10670_ave_TopVar.csv"
iefn = "./data/golddata.csv"
outfn= "./result/rf.csv"

# parse arguements
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "./data/golddata.csv"
  args[3] = 100
  args[4] = "./result/rf.csv"
}


emfn = args[1]
iefn = args[2] 
noutedge = args[3]
outfn = args[4]

# Step 1. load datasets
# load normalized expression data
ExpData<-read.table(emfn, sep=',', as.is=T, header=T, row.names=1)
dim(ExpData)
print(paste0("* Total number of genes from input gene expression matrix   : ", length(ExpData[,1])))

# load interaction edge list data
EdgData  <- read.table(iefn, sep=',', as.is=T, header = T, na.strings ="", stringsAsFactors= F)
print(paste0("* Total number of genes from input edge list                : ", length(unique(c(EdgData[,1],EdgData[,2])))))
print(paste0(" -Total number of transcription factors from input edge list: ", length(unique(EdgData[,1]))))
print(paste0(" -Total number of target genes from input edge list         : ", length(unique(EdgData[,2]))))

# Step 2. filter & Get top 5000 genes with high rowVar genes
if(length(ExpData[,1]) >5000)
{
	library(matrixStats)
	ExpData_rowvars = cbind(ExpData,rowvars=rowVars(as.matrix(ExpData)))
	head(ExpData_rowvars)
	head(ExpData_rowvars[order(-ExpData_rowvars[,'rowvars']),])
	topVarGenes <- head(ExpData_rowvars[order(-ExpData_rowvars[,'rowvars']),], 500 )
	ExpData=topVarGenes[,1: (dim(topVarGenes)[2]-1)]
	print("* Filtered top 5000 genes with high row variance")
}

# Step 3. filter edge list only with gene expression values & write the filtered edge list
# filter out edges without gene expression values for both of genes 
EdgeStatus=as.matrix(EdgData[,1] %in% rownames(ExpData) * EdgData[,2] %in% rownames(ExpData))
EdgData_Status= cbind(EdgData, EdgWgene=EdgeStatus)
EdgFilt= EdgData_Status[which(EdgData_Status[,3]==1), c(1,2) ]

# write the new filtered edge list
write.csv(EdgFilt, "filteredEdge.csv", row.names=FALSE)

# get unique genes from filtered edge list and the number of unique genes
EdgFiltUniqGenes=unique(c(EdgFilt[,1],EdgFilt[,2]))
NumUniqGenes= length(EdgFiltUniqGenes)
print(paste0("* Total number of common genes from two inputs              : ", NumUniqGenes))
print(paste0(" -Total number of transcription factors from common genes.  : ", length(unique(EdgFilt[,1]))))
print(paste0(" -Total number of target genes from input edge list         : ", length(unique(EdgFilt[,2]))))

# Step 4. get an adjacency matrix of interactions between TF and TG
# load igraph package
library(igraph)
# get adjacency matrix of filtered edge list
EdgFiltGraph=graph.data.frame(EdgFilt)
EdgFiltGraphMatrix=get.adjacency(EdgFiltGraph,sparse=FALSE)

# Step 5. extract gene expression only for common genes
# get filtered gene expression matrix
ExpFilt = ExpData[rownames(ExpData) %in% EdgFiltUniqGenes, ]

# match orders of two datasets
indx=match(rownames(EdgFiltGraphMatrix), rownames(ExpFilt))
ExpFiltMatched=ExpFilt[c(indx),]

# write the new filtered gene expression matrix
write.csv(ExpFiltMatched, "filteredExpr.csv")

# Step 6. standardize gene expression values to mean 0 and variance 1 by gene
ExpData.st =(apply(ExpFiltMatched, 1, function(x) { (x - mean(x)) / sd(x) } )) #st for column, genes

# Step 7. Run iRafNet
# load iRafNet package
library(iRafNet)

# variable preparation for iRafNet
data.st=matrix(ExpData.st,dim(ExpData.st)[1], dim(ExpData.st)[2])
W=matrix(EdgFiltGraphMatrix,dim(ExpData.st)[2],dim(ExpData.st)[2])
MTRY=round(sqrt(NumUniqGenes-1))
genes.name=colnames(ExpData.st)
M=100

# run iRafNet and obtain importance score of regulatory relationships
out.iRafNet<-iRafNet(data.st, W, mtry=MTRY, ntree=2000, genes.name)    #We need to decide number of trees
dim(out.iRafNet)
head(out.iRafNet)
	####When I tried 1655 genes, I got errors below, but when I try 271 genes, it works.
	#Error in names(impSD) <- x.col.names : 
	#  'names' attribute [2] must be the same length as the vector [0]
	#In addition: Warning messages:
	#1: In irafnet_onetarget(x = x.sorted, y = as.double(y), importance = TRUE,  :
	#  invalid mtry: reset to within valid range
	#2: In max(ncat) : no non-missing arguments to max; returning -Inf

# run iRafNet for M permuted data sets
out.perm<-Run_permutation(data.st, W, mtry=round(sqrt(p-1)), ntree=1000, genes.name, M)
out.perm.1000<-Run_permutation(data.st, W, mtry=round(sqrt(p-1)), ntree=1000, genes.name, 1000)
dim(out.perm)
head(out.perm)

# derive final networks
final.net<-iRafNet_network(out.iRafNet,out.perm,0.1)
dim(final.net)
head(final.net)

out.iRafNet_above0=out.iRafNet[(out.iRafNet[,3]>0), ]
head(out.iRafNet_above0)

#-----------------------------
# Z. Save nets into RData
# Save module colors and labels for use in subsequent parts
RDataFileName=paste0(gsub(".csv","",InputFileName),"_FilteredSortedFor_", gsub(".tab","", EdgDataFileName), "__iRafNet2.RData")
RDataFileName
save(out.iRafNet, out.perm, data.st, W, MTRY, genes.name,  file = RDataFileName)

#-----------------------------
# Load network data saved in the second part.
lnames = load(file = RDataFileName);
#The variable lnames contains the names of loaded variables.
lnames


# 6. Matrix of true regulations
truth<-out.iRafNet[,seq(1,2)]
truth<-cbind(as.character(truth[,1]),as.character(truth[,2]),as.data.frame(rep(0,,dim(out)[1])));
truth[(truth[,1]=="G2" & truth[,2]=="G1") | (truth[,1]=="G1" & truth[,2]=="G2"),3]<-1

# 6-1. Plot ROC curve and compute AUC
auc<-roc_curve(out,truth)
