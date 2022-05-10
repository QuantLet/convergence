elimna<-function(m){
#
# remove any rows of data having missing values
#
DONE=FALSE
if(is.list(m) && is.matrix(m)){
z=pool.a.list(m)
m=matrix(z,ncol=ncol(m))
DONE=TRUE
}
if(!DONE){
if(is.list(m) && is.matrix(m[[1]])){
for(j in 1:length(m))m[[j]]=na.omit(m[[j]])
e=m
DONE=TRUE
}}
if(!DONE){
if(is.list(m) && is.null(dim(m))){ #!is.matrix(m))
for(j in 1:length(m))m[[j]]=as.vector(na.omit(m[[j]]))
e=m
DONE=TRUE
}}
if(!DONE){
#if(!is.list(m)){
#if(is.null(dim(m)))
m<-as.matrix(m)
ikeep<-c(1:nrow(m))
for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
e<-m[ikeep[ikeep>=1],]
#}
}
e
}

disc2comSK<-function(x, y, alpha=.05, nboot = 500, SEED = TRUE){
#
#  Comparing two independent variables in terms of their probability function.
#  A global test of P(X=x)=P(Y=x) for all x.
#  Appears  to have no advantage over a chi-square test done by the R function disc2com
#
#  The R function binband tests this hypothesis for each x.
#
library(mc2d)
x=elimna(x)
y=elimna(y)
if(SEED)set.seed(2)
x=elimna(x)
y=elimna(y)
vals=sort(unique(c(x,y)))
n1=length(x)
n2=length(y)
K=length(vals)
C1=NULL
C2=NULL
HT=NULL
for(i in 1:K){
C1[i]=sum(x==vals[i])
C2[i]=sum(y==vals[i])
HT[i]=(C1[i]+C2[i])/(n1+n2)
}
p1hat=C1/n1
p2hat=C2/n2
test=sum((p1hat-p2hat)^2)
tv=NULL
TB=NA
VP=NA
for(ib in 1:nboot){
xx=rmultinomial(n1,1,HT)
yy=rmultinomial(n2,1,HT)
B1=NA
B2=NA
BP=NA
for(i in 1:K){
B1[i]=sum(xx[,i])
B2[i]=sum(yy[,i])
}
B1hat=B1/n1
B2hat=B2/n2
TB[ib]=sum((B1hat-B2hat)^2)
}
pv=1-mean(test>TB)-.5*mean(test==TB)
list(test=test,p.value=pv)
}


twobinom<-function(r1=sum(elimna(x)),n1=length(elimna(x)),r2=sum(elimna(y)),n2=length(elimna(y)),x=NA,y=NA,alpha=.05){
#
# Test the hypothesis that two independent binomials have equal
# probability of success using the Storer--Kim method.
#
# r1=number of successes in group 1
# n1=number of observations in group 1
#
n1p<-n1+1
n2p<-n2+1
n1m<-n1-1
n2m<-n2-1
chk<-abs(r1/n1-r2/n2)
x<-c(0:n1)/n1
y<-c(0:n2)/n2
phat<-(r1+r2)/(n1+n2)
m1<-outer(x,y,"-")
m2<-matrix(1,n1p,n2p)
flag<-(abs(m1)>=chk)
m3<-m2*flag
b1<-1
b2<-1
xv<-c(1:n1)
yv<-c(1:n2)
xv1<-n1-xv+1
yv1<-n2-yv+1
dis1<-c(1,pbeta(phat,xv,xv1))
dis2<-c(1,pbeta(phat,yv,yv1))
pd1<-NA
pd2<-NA
for(i in 1:n1)pd1[i]<-dis1[i]-dis1[i+1]
for(i in 1:n2)pd2[i]<-dis2[i]-dis2[i+1]
pd1[n1p]<-phat^n1
pd2[n2p]<-phat^n2
m4<-outer(pd1,pd2,"*")
test<-sum(m3*m4)
list(p.value=test,p1=r1/n1,p2=r2/n2,est.dif=r1/n1-r2/n2)
}


splotg2<-function(x,y,op=TRUE,xlab="X",ylab="Rel. Freq."){
#
# Frequency plot
#  op=TRUE, lines are drawn.
#
x<-x[!is.na(x)]
temp<-sort(unique(x))
freqx<-NA
for(i in 1:length(temp)){
freqx[i]<-sum(x==temp[i])
}
freqx<-freqx/length(x)
y<-y[!is.na(y)]
tempy<-sort(unique(y))
freqy<-NA
for(i in 1:length(tempy)){
freqy[i]<-sum(y==tempy[i])
}
freqy<-freqy/length(y)
plot(c(temp,tempy),c(freqx,freqy),type="n",xlab=xlab,ylab=ylab)
points(temp,freqx)
points(tempy,freqy,pch="o")
if(op){
lines(temp,freqx)
lines(tempy,freqy,lty=2)
}
}



binband<-function(x, y, KMS = FALSE, alpha = .05, ADJ.P = FALSE, 
                  plotit = TRUE, op = TRUE, xlab = "X", ylab = "Rel. Freq.", 
				  ADJ.CI = TRUE){
#
#  Comparing two independent variables in terms of their probability function.
#  For each value that occurs, say x, test P(X=x)=P(Y=x)
#  So this method is useful when dealing with highly discrete data.
#
#  If KMS=T, use Kulinskaya, Morgenthaler and Staudte (2010)
#   method for comparing binomials
# Kulinskaya, E., Morgenthaler, S. and Staudte, R. (2010). 
# Variance Stabilizing the Difference of two Binomial
#  Proportions. {\em American Statistician, 64}, 
#  350--356 DOI:10.1198/tast.2010.09096

#  Otherwise use Storer and Kim.
#
#   ADJ.P=T means that critical p-value is adjusted to control FWE when the sample
#   size is small (<50).
#
#
#  Hochberg's method is used to determine critical p-values so that FWE=alpha
#
x=elimna(x)
y=elimna(y)
vals=sort(unique(c(x,y)))
ncon=length(vals)
n1=length(x)
n2=length(y)
p.values=NA
adj=1
cv=1
if(!KMS){
output=matrix(NA,ncol=6,nrow=length(vals))
dimnames(output)=list(NULL,c("Value","p1.est","p2.est","p1-p2","p.value","p.crit"))
}
if(KMS){
output=matrix(NA,ncol=8,nrow=length(vals))
dimnames(output)=list(NULL,c("Value","p1.est","p2.est","p1-p2","ci.low","ci.up","p.value",
"p.crit"))
}
for(i in 1:length(vals)){
x1=sum(x==vals[i])
y1=sum(y==vals[i])
if(!KMS){
output[i,5]=twobinom(x1,n1,y1,n2)$p.value
output[i,2]=x1/n1
output[i,3]=y1/n2
output[i,1]=vals[i]
output[i,4]=output[i,2]-output[i,3]
}
if(KMS){
temp=bi2KMSv2(x1,n1,y1,n2)
output[i,1]=vals[i]
output[i,5]=temp$ci[1]
output[i,6]=temp$ci[2]
output[i,2]=x1/n1
output[i,3]=y1/n2
output[i,4]=output[i,2]-output[i,3]
output[i,7]=temp$p.value
}}
# Determine adjusted  critical p-value using Hochberg method
ncon=length(vals)
dvec=alpha/c(1:ncon)
if(ADJ.P){
mn=max(c(n1,n2))
cv=1
if(ncon!=2){
if(mn>50){
cv=2-(mn-50)/50
if(cv<1)cv=1
}
if(mn<=50)cv=2
}
if(KMS){
flag=(output[,7]<=2*alpha)
output[flag,8]=output[flag,8]/cv
}
if(!KMS){
cv=1
flag=(output[,5]<=2*alpha)
if(min(c(n1,n2))<20 && n1!=n2 && ncon>=5)cv=2
output[flag,5]=output[flag,5]/cv
}}
if(KMS){
temp2=order(0-output[,7])
output[temp2,8]=dvec
}
if(!KMS){
temp2=order(0-output[,5])
output[temp2,6]=dvec
}
if(plotit)splotg2(x,y,op=op, xlab=xlab, ylab=ylab)
output
}
