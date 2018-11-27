a<-c(18,12,27,7,14)
b<-c(162,26,121,21,353)
c<-c(25,123,104,3,7)
d<-c(252,431,475,25,359)

OR<-numeric(5)
for(i in 1:5)
  OR[i]=a[i]*d[i]/(c[i]*b[i])

asum=sum(a)
bsum=sum(b)
csum=sum(c)
dsum=sum(d)

OR_unadjusted=asum*dsum/(bsum*csum)
OR_unadjusted

logoR_unadjusted=log(OR_unadjusted)
varlog_unadjusted=1/asum+1/bsum+1/csum+1/dsum
sdlog_unadjusted=sqrt(varlog_unadjusted)

logupper_unadjusted=logoR_unadjusted+1.96*sdlog_unadjusted
loglower_unadjusted=logoR_unadjusted-1.96*sdlog_unadjusted

exp(logupper_unadjusted)
exp(loglower_unadjusted)

vsub=numeric(5)
for (i in 1:5) 
  vsub[i]=b[i]*c[i]/(a[i]+b[i]+c[i]+d[i])

OR_adjusted=0
for(i in 1:5)
  OR_adjusted=OR_adjusted+vsub[i]/sum(vsub)*OR[i]

OR_adjusted

logoR_adjusted=log(OR_adjusted)

S1=0
S2=0
S3=0
S4=0
S5=0

for(i in 1:5)
{
  S1=S1+a[i]*d[i]/(a[i]+b[i]+c[i]+d[i])
  S2=S2+b[i]*c[i]/(a[i]+b[i]+c[i]+d[i])
  S3=S3+(a[i]+d[i])*a[i]*d[i]/(a[i]+b[i]+c[i]+d[i])^2
  S4=S4+(b[i]+c[i])*b[i]*c[i]/(a[i]+b[i]+c[i]+d[i])^2
  S5=S5+((a[i]+d[i])*b[i]*c[i]+(b[i]+c[i])*a[i]*d[i])/(a[i]+b[i]+c[i]+d[i])^2
}

varlog_adjusted=S3/(2*S1^2)+S5/(2*S1*S2)+S4/(2*S2^2)

sdlog_adjusted=sqrt(varlog_adjusted)

logupper_adjusted=logoR_adjusted+1.96*sdlog_adjusted
loglower_adjusted=logoR_adjusted-1.96*sdlog_adjusted

exp(logupper_adjusted)
exp(loglower_adjusted)
