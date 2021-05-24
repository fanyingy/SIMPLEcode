#Mixed membership model with estimated K or true K for T_{ij}
set.seed(3000)
library(Matrix)
library(LaplacesDemon)
library(rARPACK)
library(nortest)
library(emulator)
po<-function(W,n){  #x is the matrix£¬n is the power
  t=W%*%W
  i=2
  while (i<n)
  {i=i+1
  t=t%*%W
  }
  t}
R=function(b,x,y,t)
{
  a=0
  for (j in 2:l)
    a=a-b[,((j-2)*n+1):((j-1)*n)]/(t^{j+1})
  c=-sum(x*y)/t+quad.3form(a,x,y)
  as.numeric(c)
}
P=function(b,x,y,t)
{
  a=0
  for (j in 2:l)
    a=a-b[,((j-2)*n+1):((j-1)*n)]/(t^j)
  a=-sum(x*y)+quad.3form(a,x,y)
  as.numeric(a)
}
TP=function(me,x,y,t)
{
  a=0
  for (j in 2:l)
    a=a+(j+1)*me[,((j-2)*n+1):((j-1)*n)]/(t^j)
  a=as.numeric(sum(x*y)+quad.3form(a,x,y))
  1/a
}
A=function(me,D,x,y,z,t)
{
  
  a=P(me,x,y,t)
  b=P(me,x,z,t)
  c=P(me,z,y,t)
  d=P(me,z,z,t)
  e=as.numeric(a-b%*%solve(t*solve(D)+d)%*%c)
  return(list(p1=a,p2=b,p3=c,p4=d,p5=e))
}

be=function(me,D,u,y,t)
{
  u-y%*%solve(solve(D)+R(me,y,y,t))%*%R(me,y,u,t)
}

mew=function(l,n,k)
{
  b=matrix(rep(0,(l-1)*n^2),nrow=n)
  for (i in 1:1000)
  {
    X <- sample_sbm(pm,n,k)
    W=X-E1
    for (j in 1:(l-1))
      b[,((j-1)*n+1):(j*n)]=po(W,j+1)+b[,((j-1)*n+1):(j*n)]
  }
  b/1000
}


sample_sbm=function(pm,n,k)
{
  W=matrix(rbern(n^2,pm),nrow=n)
  W[lower.tri(W)]=0
  W=W+t(W)
  diag(W)=diag(W)/2
  W
}




sa=function(n,X)#Estimate K
{
  Spi=eigs_sym(X,5)
  d=Spi$values
  a=colSums(X)
  cri=sqrt(2)*sqrt(max(a)*log(n))
  sum(abs(d)>cri)
  
}
samp=function(n,k0,PM,i1,j1,j2)#Calculate T_{ij}, K0 is the true K, n is dimension of the matrix, PM is the expectation of the adjacency matrix, we test the equality of i1, j1 and i1,j2
{ 
  el=eigen(PM)
  e3=el$vectors[,1]
  e4=el$vectors[,2:5]
  W12=PM
  ev=el$vectors[,1:5]
  d=el$values
  f=matrix(rep(0,2*o),nrow=2)
  for (j in 1:o)
  {
    X <- sample_sbm(PM,n,k0)
      k=max(sa(n,X),1)#Estimated K
      #   k=k0# True K
    thr=qchisq(0.95, df=k)
    cov1=cov=matrix(rep(0,k^2),nrow=k)
    Spi=eigs_sym(X,k)
    d=Spi$values
    e1=Spi$vectors[,1]
    e1=sign(sum(e1*ev[,1]))*e1
    if (k>1)
    {
    e2=Spi$vectors[,2:k]
    e3=cbind(e1,e2)
    } else {e3=cbind(e1,e1)}
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
    for (a in 1:k)
      for (b  in 1:k)
      {
        va=e3[,a]
        vb=e3[,b]
        cov1[a,b]=sum(W2[i1,]*va*vb+W2[j1,]*va*vb)-W2[i1,j1]*(va[j1]*vb[j1]+va[i1]*vb[i1])+W2[i1,j1]*(va[j1]-va[i1])*(vb[j1]-vb[i1])
        cov[a,b]=cov1[a,b]/d[a]/d[b]
      }
    
    if (k>1)
    {
    st=e3[i1,]-e3[j1,]
    } else {st=e3[i1,1]-e3[j1,1]}
   f[1,j]=t(st)%*%solve(cov)%*%st>thr 
    
    for (a in 1:k)
      for (b  in 1:k)
      {
        va=e3[,a]
        vb=e3[,b]
        cov1[a,b]=sum(W2[i1,]*va*vb+W2[j2,]*va*vb)-W2[i1,j2]*(va[j2]*vb[j2]+va[i1]*vb[i1])+W2[i1,j2]*(va[j2]-va[i1])*(vb[j2]-vb[i1])
        cov[a,b]=cov1[a,b]/d[a]/d[b]
      }
    
    if (k>1)
    {st=e3[i1,]-e3[j2,]} else {st=e3[i1,1]-e3[j2,1]}
   f[2,j]=t(st)%*%solve(cov)%*%st>thr  }
  f
}
samk=function(n,k0,PM)
{
  f=rep(0,o)
  
  for (j in 1:o)
  {X <- sample_sbm(PM,n,k0)
  f[j]=sa(n,X)
  }
  f
}
n=3000; k0=3; o=500;#o is the number of repetitions, k0 is the true K 
x=0.2;rho=0.2;n0=500;n1=(n-3*n0)/4#x=0.2 is equal to the first entry of membership vector a_1
r=(2:9)/10#r represents \theta in the paper
le=length(r)
mea=matrix(rep(0,2*le*o),ncol=o)
qq=c(2,n)
pii=matrix(rep(0,7*k0),ncol=k0)
pii[1,]=c(x,1-2*x,x)
pii[2,]=c(1-2*x,x,x)
pii[3,]=c(x,x,1-2*x)
pii[4,]=c(1/3,1/3,1/3)
pii[5,]=c(1,0,0)
pii[6,]=c(0,0,1)
pii[7,]=c(0,1,0)
pm=matrix(rep(rho,k0^2),nrow=k0)
pm[1,k0]=pm[k0,1]=rho/2
diag(pm)=1
PI=matrix(rep(0,n*k0),ncol=k0)
for (i in 1:4)
{PI[(1+(i-1)*n1):(i*n1),]=rep(1,n1)%*%t(pii[i,])}
for (i in 1:3)
{PI[(n-3*n0+1+(i-1)*n0):(n-3*n0+i*n0),]=rep(1,n0)%*%t(pii[i+4,])}
PM1=PI%*%pm%*%t(PI)
for (j1 in 1:le)
{ 
  PM=r[j1]*PM1
  s=samp(n,k0,PM,1,qq[1],qq[2])
  j2=j1
  mea[j2,]=s[1,]
  j2=j1+le
  mea[j2,]=s[2,]
  print(j2)
}
mee=matrix(rep(0,4*le),nrow=4)

for (i in 1:le)
{
  mee[1,i]=sum(mea[i,]==0)/o
  mee[2,i]=sum(mea[i+le,]==0)/o
  mee[3,i]=mean(mea[i,])
  mee[4,i]=var(mea[i,])
}
#the first two rows of mee are the p-values of T_{i1j1} and T_{i2j2} 