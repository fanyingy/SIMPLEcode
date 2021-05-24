library("igraph")
library("igraphdata")
library(Matrix)
library(LaplacesDemon)
library(rARPACK)
library(nortest)
library(RCurl)
library(emulator)
set.seed(0)
y=read.csv('C://Users//HP//Downloads//political-books-edges.csv')#the edges of political data
y=y+1
X=matrix(rep(0,105^2),nrow=105)
y=as.matrix(y)
X[y[,1:2]]=1
X=X+t(X)
k=2#K
cov1=cov=matrix(rep(0,(k)^2),nrow=k)
Spi=eigs_sym(X,k)
d=Spi$values
e1=Spi$vectors[,1]
e2=Spi$vectors[,2:k]
e3=cbind(e1,e2)
for (t in 1:k)
  e3[,t]=sign(sum(e3[,t]*ev[,t]))*e3[,t]
W=X
for (t in 1:k)
  W=W-d[t]*tcrossprod(e3[,t])
W2=W^2
c=1:k
for (t in 1:k)
{
  c[t]=as.numeric(1/(1/d[t]+sum(W2*e3[,t]^2)/(d[t]^3)))
}
W=X
for (t in 1:k)
  W=W-c[t]*tcrossprod(e3[,t])
W2=W^2
sam=function(i1,j1)
{    
  for (a in 1:k)
  for (b  in 1:k)
  {
    va=e3[,a]
    vb=e3[,b]
    cov1[a,b]=sum(W2[i1,]*va*vb+W2[j1,]*va*vb)-W2[i1,j1]*(va[j1]*vb[j1]+va[i1]*vb[i1])+W2[i1,j1]*(va[j1]-va[i1])*(vb[j1]-vb[i1])
    cov[a,b]=cov1[a,b]/d[a]/d[b]
  }
  
  
  st=e3[i1,]-e3[j1,]
  f=t(st)%*%solve(cov)%*%st
  pchisq(f, df=k)
}

s=1:105
l=length(s)
M=matrix(rep(0, l^2),nrow=l)
for (i in 1:l)
  for (j in 1:l)
  {M[i,j]=1-sam(s[i],s[j])}
s1=c(105,104,59,29,78,77,47,19,50)
P=1-pchisq(M[s1,s1], df=k)#p-values of selected labels s1
