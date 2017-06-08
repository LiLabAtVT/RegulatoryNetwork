library(iRafNet)

# 1. Generate data sets
n<-20  # sample size
p<-5   # number of genes
genes.name<-paste("G",seq(1,p),sep="")  # genes name
#genes.name
data<-matrix(rnorm(p*n),n,p) # generate expression matrix
#data

data[,1]<-data[,2]              # var 1 and 2 interact
#data

W<-abs(matrix(rnorm(p*p),p,p))	# generate weights for regulatory relationships
#W
 
# 2. Standardize variables to mean 0 and variance 1
data.st <- (apply(data, 2, function(x) { (x - mean(x)) / sd(x) } ))
#data.st

# 2-1. Visualization of Standardized variables
par(mfrow = c(1,2));
boxplot(data, ylim=c(-4, 3))
abline(h=c(-1,1), lty=2,col="blue"); abline(h=c(0), lty=1,col="red"); 
boxplot(data.st, ylim=c(-4, 3))
abline(h=c(-1,1), lty=2,col="blue"); abline(h=c(0), lty=1,col="red"); 
summary(data)
summary(data.st)

# 3. Run iRafNet and obtain importance score of regulatory relationships
out.iRafNet<-iRafNet(data.st, W, mtry=round(sqrt(p-1)), ntree=1000, genes.name)
str(out.iRafNet)
#out.iRafNet

# 4. Run iRafNet for M permuted data sets
out.perm<-Run_permutation(data.st, W, mtry=round(sqrt(p-1)), ntree=1000, genes.name, M)
#out.perm

# 5. Derive final networks
final.net<-iRafNet_network(out.iRafNet,out.perm,0.001)
final.net      
final.net<-iRafNet_network(out.iRafNet,out.perm,0.1)
final.net

# 6. Matrix of true regulations
truth<-out.iRafNet[,seq(1,2)]
truth<-cbind(as.character(truth[,1]),as.character(truth[,2]),as.data.frame(rep(0,,dim(out)[1])));
truth[(truth[,1]=="G2" & truth[,2]=="G1") | (truth[,1]=="G1" & truth[,2]=="G2"),3]<-1

# 6-1. Plot ROC curve and compute AUC
auc<-roc_curve(out,truth)

