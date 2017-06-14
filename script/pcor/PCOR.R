# sample code for performing the network analysis.
# Partial Correlation method
# required R library: GeneNet
# sample input data: GSE10670_ave_TopVar.csv
# sample output data: pcor_edgelist.csv
#
# How to run the program:
# Assume working direction contain a data folder and a script folder.
# 
# Rscript PCOR.R [ExpressionMatrix.csv] [GoldStandardEdges.csv] [number of edges] [OutputFile.csv]
# 
# example
# 
# Rscript ./script/pcor/PCOR.R ./data/GSE10670_ave_TopVar.csv ./data/golddata.csv nedges ./result/pcor.csv
#
# or simply
# 
# Rscript ./script/pcor/PCOR.R ./data/GSE10670_ave_TopVar.csv
# 

# parse arguements
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "./data/golddata.csv"
  args[3] = 100
  args[4] = "./results/pcor.csv"
}

inputfn = args[1]
gsfn = args[2] # gold standard file, this is not used in this setting
noutedge = args[3]
outputfn = args[4]

# Step 1. load normalized expression data.
dat<-read.table(inputfn,sep=',',as.is=T,header=T,row.names=1)

# To perform partial correlation analysis
library("GeneNet")

# calculate partial correlation, using "static" method.
pcor=ggm.estimate.pcor(t(dat),method = 'static')

require(reshape)
netmat<-pcor 

# convert adjancy matrix to edgelist 
elist<-melt(netmat[1:30,1:30]) # only use top 30 to reduce number of outputs.
# This needs to be changed to find top ranking edges in a matrix.
elist[,1]<-as.character(elist[,1])
elist[,2]<-as.character(elist[,2])
colnames(elist)<-c('regulator','target','score')


# prepare final output matrix
out_elist<-elist[elist[,'score']<1,] # remove the self connection
out_elist<-data.frame(out_elist,
                      rank=rank(-out_elist[,'score']),# this give the highest score rank 1
                      stringsAsFactors = FALSE)
out_elist<-out_elist[order(out_elist[,'score'],decreasing = TRUE),]

# only write top 100 edges, this can be a new parameter.
write.table(out_elist[1:100,],outputfn,sep=',',row.names=FALSE)

