# Get all ATH1 array expression data and design matrix
library(GEOquery)
#library(SRAdb)

# get list of GSE id:
GSE<-read.table('data/ExpressionFromArray/GSElist_morethan10sample_033017.csv',
                       as.is=T,sep=',')
dim(GSE)

# Get GSE without expression matrix
# These matrices have design information.
gdsnomat<-list()
for (i in 1:nrow(GSE))
{
  geo<-GSE[i,1]
  gdsnomat[[i]] <- getGEO(geo,GSEMatrix = FALSE)
}

# get tables for GSEs
# column: GSE, GSM, simply summary

outtable<-NULL
for(i in 1:length(gdsnomat))
{
  tmpgds<-gdsnomat[[i]]
  tmpgsm<-GSMList(tmpgds)
  for(j in 1:length(tmpgsm))
  {
    title<-tmpgsm[[j]]@header$title
    gsename<-tmpgsm[[j]]@header$series_id
    gsmname<-tmpgsm[[j]]@header$geo_accession
    print(c(gsename,gsmname,title))
    outtable<-rbind(outtable,c(gsename,gsmname,title))
  }
}
colnames(outtable)<-c('GSE','GSM','Summary')

write.table(outtable,'results/annotation/GSElistAnnotation.csv',sep=',',row.names=FALSE)

