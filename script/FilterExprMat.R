# 1. load dataset
argv <- commandArgs(trailingOnly = TRUE)
expfn<-argv[1] # expression data
goldsd<-argv[2] # input gold standard file
tfnames<-argv[3] # input list of transcription factor names, downloaded from PlantTFDB.
outprefix<-argv[4] # prefix of output files

# usage:
# Rscript FilterExprMat.R argv1 argv2 argv3 argv4
# Example:
# Rscript script/FilterExprMat.R ./data/TestRunData/GSE10670_ave.csv ./data/golddata.csv ./data/Ath_TF_list GSE10670_ave_topVar
# Rscript script/FilterExprMat.R ./data/TestRunData/GSE46205_ave.csv ./data/golddata.csv ./data/Ath_TF_list GSE46205_ave_topVar
# Rscript script/FilterExprMat.R ./data/TestRunData/GSE5624_ave.csv ./data/golddata.csv ./data/Ath_TF_list GSE5624_ave_topVar



# these are for testing:
#goldsd<-'./data/golddata.csv'
#tfnames<-'./data/Ath_TF_list'
#expfn<-'./data/TestRunData/GSE10670_ave.csv'
#outprefix<-'GSE10670_ave_topVar'

ExprData  <- read.table(expfn, sep = "," , header = T, na.strings ="", 
                        stringsAsFactors= F)

TF2TA<-read.table(goldsd,sep = "," , header = T, na.strings ="", 
                  stringsAsFactors= F)

TFnames<-read.table(tfnames,sep = "\t" , header = T, na.strings ="", 
                    stringsAsFactors= F)

# TFnames has all TF namse from Arabidopsis TF database
# TF2TA has the DAP-seq all TF names as the first column
allTF<-unique(c(TFnames[,2],TF2TA[,1]))
DAPTF<-unique(TF2TA[,1])

# get any data not in DAP TF list.
KeepGeneNames<-rownames(ExprData)[!rownames(ExprData)%in%DAPTF]

EXPnonDAPTF<-ExprData[KeepGeneNames,]
#   > length(DAPTF)
# [1] 387
# > length(allTF)
# [1] 1727
# > dim(EXPnonDAPTF)
# [1] 21331     8

#-----------------------------
# 2. Filter & Get top N% of high rowVar genes
#--------------------
library(matrixStats)
exp.data_rowvars = cbind(EXPnonDAPTF,rowvars=rowVars(as.matrix(EXPnonDAPTF)))
IncludePerc<-0.2 # this controls inclusion percentage.
topVarGenes <- head(exp.data_rowvars[order(-exp.data_rowvars[,'rowvars']),], 
                    (IncludePerc*dim(exp.data_rowvars)[1]) )

# these are genes to keep
tokeep<-rownames(topVarGenes)
dapTF<-unique(TF2TA[TF2TA[,2]%in%tokeep,1]) # get DAP-seq TF list
# > length(dapTF)
# [1] 387

tokeepall<-unique(c(tokeep,dapTF))

# split genes to TF and targets
tokeepTF<-tokeepall[tokeepall%in%allTF]
tokeepOther<-tokeepall[!tokeepall%in%tokeepTF]

outTFmat<-ExprData[tokeepTF,]
outOthermat<-ExprData[tokeepOther,]

# write output:
outTFfn<-paste(outprefix,'_TF.csv',sep='')
outTAfn<-paste(outprefix,'_TA.csv',sep='')
#outGoldFN<-paste(outprefix,'_goldsd.csv',sep='')

write.table(outTFmat,outTFfn,sep=',')
write.table(outOthermat,outTAfn,sep=',')



