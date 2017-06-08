## README: If you want to see how your input data works here, you can run: 

#-----------------------------
# 0. Setup environment & load dataset
getwd();
workingDir = "/Users/jylee43/Google Drive/coursework_2017SPRING/ProblemSolving/4_ReadDataWork/Microarray";
setwd(workingDir); 
getwd();
fileList=list.files(workingDir, pattern="_ave_TopVar.csv")
fileList
#-----------------------------
# 1. load dataset
InputFileName="GSE10670_ave_TopVar.csv"
InputFileName=fileList[1]
InputFileName
ExprData  <- read.table(InputFileName, sep = "," , header = T, na.strings ="", stringsAsFactors= F)
dim(ExprData)
head(ExprData)
#--------
# sort
head(ExprData[,1])
#ExprData[order(ExprData[,1]),] 
ExprDataOrdered= ExprData[order(ExprData[,1]),] 
head(ExprDataOrdered)
#--------
exp.data = ExprDataOrdered[, 2: (dim(ExprDataOrdered)[2]) ]
rownames(exp.data) = ExprDataOrdered[,1]
dim(exp.data)
head(exp.data)

#-----------------------------
# 1. load dataset
DapDataFileName="AllUniq.chr1-5_GEM_events.nS_targets.tab"
DapData  <- read.table(paste0("/Users/jylee43/Google Drive/coursework_2017SPRING/ProblemSolving/4_ReadDataWork/", DapDataFileName), sep = "\t" , header = T, na.strings ="", stringsAsFactors= F)
dim(DapData)
head(DapData)

head(rownames(exp.data) )
head(DapData[1])
head(DapData[2])

DapData[1,1] %in% rownames(exp.data) 
DapData[1,2] %in% rownames(exp.data)
DapData[1,1] %in% rownames(exp.data) && DapData[1,2] %in% rownames(exp.data)

(DapData[2,] %in% rownames(exp.data) )[2] * (DapData[2,] %in% rownames(exp.data) )[1]

head(DapData[,1] %in% rownames(exp.data))
head(DapData[,2] %in% rownames(exp.data))
cbind(head(DapData), (as.matrix(head(DapData[,1] %in% rownames(exp.data)) * head(DapData[,2] %in% rownames(exp.data)))))

DapFiltered=as.matrix(DapData[,1] %in% rownames(exp.data) * DapData[,2] %in% rownames(exp.data))

DapData_Filtered= cbind(DapData, DAP_RNA=DapFiltered)
head(DapData_Filtered)
head(DapData_Filtered[,3])
head( DapData_Filtered[which(DapData_Filtered[,3]==1), ] )
DapFilteredData= DapData_Filtered[which(DapData_Filtered[,3]==1), c(1,2) ]
dim(DapFilteredData)
head(DapFilteredData)
dim(DapData)

#-----------------------------
DapDataFileName
InputFileName
OutputFileName=paste0(gsub(".csv","",InputFileName),"_", gsub(".tab",".csv", DapDataFileName))
OutputFileName
write.csv(DapFilteredData, file = OutputFileName, row.names=FALSE)
#-----------------------------#-----------------------------
format(Sys.time(), "%Y/%m/%d_%H:%M:%S")
DapFilteredList=list()
DapFilteredList
for ( i in c(1:length(DapData[,1])) )
{
	print(i)
	if (DapData[i,1] %in% rownames(exp.data) && DapData[i,2] %in% rownames(exp.data))
	{
	#print(DapData[i,])
	DapFilteredList=rbind(DapFilteredList, DapData[i,])
	}	
}
format(Sys.time(), "%Y/%m/%d_%H:%M:%S")


#-----------------------------#-----------------------------
#http://stackoverflow.com/questions/16584948/how-to-create-weighted-adjacency-list-matrix-from-edge-list
library(igraph)

head(DapFilteredData)
DapFilteredDataList=c(DapFilteredData[,1],DapFilteredData[,2])
NumUniqGenes= length(unique(DapFilteredDataList))
NumUniqGenes

DapFilteredDataGraph=graph.data.frame(DapFilteredData)
DapFilteredDataGraph
DapFilteredDataGraphMatrix=get.adjacency(DapFilteredDataGraph,sparse=FALSE)
dim(DapFilteredDataGraphMatrix)
sum(DapFilteredDataGraphMatrix)

#-----------------------------#-----------------------------

format(Sys.time(), "%Y/%m/%d_%H:%M:%S")
DapFilteredList=list()
DapFilteredList
for ( i in c(1:length(DapData[,1])) )
{
	print(i)
	if (DapData[i,1] %in% rownames(exp.data) && DapData[i,2] %in% rownames(exp.data))
	{
	#print(DapData[i,])
	DapFilteredList=rbind(DapFilteredList, DapData[i,])
	}	
}
#-----------------------------#-----------------------------



#-----------------------------
library(iRafNet)

dim(exp.data)
DapFilteredDataUniqueGenes =unique(DapFilteredDataList)
head(DapFilteredDataUniqueGenes)
length(DapFilteredDataUniqueGenes)

exp.data_Filtered = exp.data[rownames(exp.data) %in% DapFilteredDataUniqueGenes, ]
dim(exp.data_Filtered)
head(exp.data_Filtered)
#-------

head(rownames(DapFilteredDataGraphMatrix))
head(rownames(exp.data_Filtered))

indx=match(rownames(DapFilteredDataGraphMatrix), rownames(exp.data_Filtered))
head(indx)
head(exp.data_Filtered[c(indx),])
exp.data_FilteredMatched=exp.data_Filtered[c(indx),]
head(exp.data_FilteredMatched)
head(rownames(DapFilteredDataGraphMatrix))

OutputFileName=paste0(gsub(".csv","",InputFileName),"_FilteredSortedFor_", gsub("tab","csv", DapDataFileName))
OutputFileName
write.csv(exp.data_FilteredMatched, file = OutputFileName)
#-----------------------------
exp.data=t(exp.data_FilteredMatched)
dim(exp.data)
# 2. Standardize variables to mean 0 and variance 1

exp.data.st1 =(apply(exp.data, 1, function(x) { (x - mean(x)) / sd(x) } )) #st for column, genes
exp.data.st2 =(apply(exp.data.st1, 2, function(x) { (x - mean(x)) / sd(x) } )) #st for column, genes



par(mar=c(8, 4, 4, 4) + 0.1)
par(mfrow = c(1,3));
boxplot(t(exp.data), las = 2)
boxplot((exp.data.st1), las = 2)
boxplot((exp.data.st2), las = 2)


#-----------------------------#-----------------------------

# 2-1. Visualization of Standardized variables
par(mfrow = c(1,2));
boxplot(data, ylim=c(-4, 3))
abline(h=c(-1,1), lty=2,col="blue"); abline(h=c(0), lty=1,col="red"); 
boxplot(data.st, ylim=c(-4, 3))
abline(h=c(-1,1), lty=2,col="blue"); abline(h=c(0), lty=1,col="red"); 
summary(data)
summary(data.st)
#-----------------------------#-----------------------------

# 3. Run iRafNet and obtain importance score of regulatory relationships
library(iRafNet)

data.st=(exp.data.st2)
dim(data.st)

W=DapFilteredDataGraphMatrix
dim(W)

NumUniqGenes
MTRY=round(sqrt(NumUniqGenes-1))
MTRY

genes.name=rownames(data.st)
head(genes.name)

out.iRafNet<-iRafNet(data.st, W, mtry=round(sqrt(NumUniqGenes-1)), ntree=1000, genes.name)


#-----------------------------
str(out.iRafNet)
head(out.iRafNet)
dim(out.iRafNet)
min(out.iRafNet[,3])
max(out.iRafNet[,3])

summary(out.iRafNet[,3])
boxplot(out.iRafNet[,3])
hist(out.iRafNet[,3])
dim(out.iRafNet)
dim(out.iRafNet[(out.iRafNet[,3]>0), ])
out.iRafNet_Filtered=out.iRafNet[(out.iRafNet[,3]>0), ]
head(out.iRafNet_Filtered)

min(out.iRafNet_Filtered[,3])
max(out.iRafNet_Filtered[,3])

summary(out.iRafNet_Filtered[,3])
boxplot(out.iRafNet_Filtered[,3])
hist(out.iRafNet_Filtered[,3])


# 4. Run iRafNet for M permuted data sets
dim(data.st)
head(data.st)

data.st=t(data.st)
dim(W)
out.perm<-Run_permutation(data.st, W, mtry=round(sqrt(NumUniqGenes-1)), ntree=1000, genes.name, 10)

#-----------------------------
# Z. Save nets into RData
# Save module colors and labels for use in subsequent parts
RDataFileName=paste0(gsub(".csv","",InputFileName),"_FilteredSortedFor_", gsub(".tab","", DapDataFileName), "__iRafNet2.RData")
RDataFileName
save(out.perm, out.iRafNet,data.st, W, MTRY,genes.name,  file = RDataFileName)

#-----------------------------
# Load network data saved in the second part.
lnames = load(file = RDataFileName);
#The variable lnames contains the names of loaded variables.
lnames
#-----------------------------
head(out.perm)
dim(out.perm)

head(out.iRafNet)
dim(out.iRafNet)

# 5. Derive final networks     
final.net<-iRafNet_network(out.iRafNet,out.perm, 100)
?iRafNet_network
dim(final.net)
head(final.net)


# 6. Matrix of true regulations
truth<-out.iRafNet[,seq(1,2)]
truth<-cbind(as.character(truth[,1]),as.character(truth[,2]),as.data.frame(rep(0,,dim(out)[1])));
truth[(truth[,1]=="G2" & truth[,2]=="G1") | (truth[,1]=="G1" & truth[,2]=="G2"),3]<-1

# 6-1. Plot ROC curve and compute AUC
auc<-roc_curve(out,truth)

