#This code is for S&P 500 data in the paper
library("igraph")
library("igraphdata")
library(Matrix)
library(LaplacesDemon)
library(rARPACK)
library(nortest)
library(RCurl)
library(emulator)
library(MASS)
set.seed(0)
X=read.csv("spx500.csv",header=T)#S&P 500
X=X[,-1]
X=X[, colSums(is.na(X))==0]
name1=names(X)#the names without 'na'
F=read.csv("factor.csv",header=T)#factors of S&P 500
F=F[,-1]
F=F[2:length(X[,1]),]
F1=F[,-c(4)]
F1=as.matrix(F1)
F2=F[,4]#RF
E=matrix(rep(0,length(F2)*length(name1)),ncol=length(name1))
for (i in 1:length(name1))
{
  y=diff(log(as.numeric(X[,i])))*100-F2
  fit=lm(y~F1)
  E[,i]=fit$residuals
}
corv=cor(E)
Z=corv#correlation matrix
z1=eigen(Z)$vectors[,1]
z2=quantile(abs(z1),0.1)
name1=name1[abs(z1)>z2]
X=X[,abs(z1)>z2]
E=matrix(rep(0,length(F2)*length(name1)),ncol=length(name1))
for (i in 1:length(name1))
{
  y=diff(log(as.numeric(X[,i])))*100-F2
  fit=lm(y~F1)
  E[,i]=fit$residuals
}
corv=cor(E)
X=corv

k=3#spike eigenvalues: K in our paper.
cov1=cov=matrix(rep(0,(k-1)^2),nrow=k-1)
Spi=eigen(X)
d=Spi$values[1:k]
e1=Spi$vectors[,1]
e2=Spi$vectors[,2:k]
e3=cbind(e1,e2)
W=X
for (t in 1:k)
  W=W-d[t]*tcrossprod(e3[,t])
W2=W^2
c=1:k
for (t in 1:k)
{
  c[t]=as.numeric(1/(1/d[t]+sum(W2*e3[,t]^2)/(d[t]^3)))#estimate population eigenvalue
}
W=X
for (t in 1:k)
  W=W-c[t]*tcrossprod(e3[,t])
W2=W^2
sam=function(i1,j1)
{    for (a in 2:k)
  for (b  in 2:k)
  {
    va1=d[1]/d[a]/e3[i1,1]*e3[,a]-e3[i1,a]/e3[i1,1]/e3[i1,1]*e3[,1]
    vb1=d[1]/d[b]/e3[i1,1]*e3[,b]-e3[i1,b]/e3[i1,1]/e3[i1,1]*e3[,1]
    va2=d[1]/d[a]/e3[j1,1]*e3[,a]-e3[j1,a]/e3[j1,1]/e3[j1,1]*e3[,1]
    vb2=d[1]/d[b]/e3[j1,1]*e3[,b]-e3[j1,b]/e3[j1,1]/e3[j1,1]*e3[,1]
    cov1[a-1,b-1]=sum(W2[i1,]*va1*vb1+W2[j1,]*va2*vb2)-W2[i1,j1]*(va1[j1]*vb1[j1]+va2[i1]*vb2[i1])+W2[i1,j1]*(va1[j1]-va2[i1])*(vb1[j1]-vb2[i1])
  }
  cov=cov1/d[1]/d[1]
  
  st=e3[i1,2:k]/e1[i1]-e3[j1,2:k]/e1[j1]
 
  f=t(st)%*%ginv(cov)%*%st#G satistics
  f
}

s=1:length(name1)
l=length(s)
t=0
M=matrix(rep(0, l^2),nrow=l)
for (i in 1:l)
  for (j in 1:l)
  {M[i,j]=sam(s[i],s[j])}#p-value matrix is 1-pchisq(M, df=k-1)
s1=c('TGT','HD','L','AAPL','INTC','MCHP','AEE', 'NEE','EVRG','ADBE','EBAY','AMZN')#the table in the paper
ss=1:length(s1)
for (i in 1:length(ss))
  ss[i]=which(name1==s1[i])
pmat<-1-pchisq(M, df=k-1)#p-values for the selected stocks in s1


