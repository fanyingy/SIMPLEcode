#This code is for political data in our paper with statistics G_ij
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
cov1=cov=matrix(rep(0,(k-1)^2),nrow=k-1)
Spi=eigs_sym(X,k)
d=Spi$values
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
  c[t]=as.numeric(1/(1/d[t]+sum(W2*e3[,t]^2)/(d[t]^3)))
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

  f=t(st)%*%solve(cov)%*%st
  f
}
s=1:105
l=length(s)
t=0
M=matrix(rep(0, l^2),nrow=l)
for (i in 1:l)
  for (j in 1:l)
  {M[i,j]=sam(s[i],s[j])}
s1=c(105,104,59,29,78,77,47,19,50)
P=1-pchisq(M[s1,s1], df=k-1)#p-values of selected labels s1
