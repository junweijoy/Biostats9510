library(binom)
library(ggplot2)

Size=numeric(180)
Prop=numeric(180)
Level=numeric(180)
Err.rate=numeric(180)
Method=numeric(180)
CP=numeric(180)
Len=numeric(180)

n<-c(10,25,50,100,250)
p<-c(0.05,0.1,0.2,0.3,0.4,0.5)
alpha<-c(0.05,0.01)

set.seed(55)
## Problem 2

#LRT test
for(i in 1:length(n))
  for(j in 1:length(p))
    for(k in 1:length(alpha))
    {
      true=p[j]
      critical=qchisq(1-alpha[k],df=1)
      m=0
      for(s in 1:1000)
      {
        x=rbinom(1,n[i],true)
        estimate=x/n[i]
        gam=((true/estimate)^x)*(((1-true)/(1-estimate))^(n[i]-x))
        loggam=(-2)*log(gam)
        if (loggam > critical)
          m=m+1
      }
      Err.rate[(i-1)*12+(j-1)*2+k]=m/1000
      Method[(i-1)*12+(j-1)*2+k]="LRT"
      Size[(i-1)*12+(j-1)*2+k]=n[i]
      Prop[(i-1)*12+(j-1)*2+k]=p[j]
      Level[(i-1)*12+(j-1)*2+k]=alpha[k]
    }

#Score test

for(i in 1:length(n))
  for(j in 1:length(p))
    for(k in 1:length(alpha))
    {
      true=p[j]
      critical=qchisq(1-alpha[k],df=1)
      m=0
      fisher=n[i]/true+n[i]/(1-true)
      for(s in 1:1000)
      {
        x=rbinom(1,n[i],true)
        score=x/true-(n[i]-x)/(1-true)
        if ((score^2)/fisher > critical)
          m=m+1
      }
      Err.rate[(i-1)*12+(j-1)*2+k+60]=m/1000
      Method[(i-1)*12+(j-1)*2+k+60]="Score"
      Size[(i-1)*12+(j-1)*2+k+60]=n[i]
      Prop[(i-1)*12+(j-1)*2+k+60]=p[j]
      Level[(i-1)*12+(j-1)*2+k+60]=alpha[k]
    }

#Wald test
for(i in 1:length(n))
  for(j in 1:length(p))
    for(k in 1:length(alpha))
    {
      true=p[j]
      critical=qchisq(1-alpha[k],df=1)
      m=0
      for(s in 1:1000)
      {
        x=rbinom(1,n[i],true)
        estimate=x/n[i]
        fisher=n[i]/estimate+n[i]/(1-estimate)
        if ((estimate-true)^2*fisher > critical)
          m=m+1
      }
      Err.rate[(i-1)*12+(j-1)*2+k+120]=m/1000
      Method[(i-1)*12+(j-1)*2+k+120]="Wald"
      Size[(i-1)*12+(j-1)*2+k+120]=n[i]
      Prop[(i-1)*12+(j-1)*2+k+120]=p[j]
      Level[(i-1)*12+(j-1)*2+k+120]=alpha[k]
    }


result=data.frame(Size,Prop,Level,Method,Err.rate)
ggplot(result,aes(factor(Size),Err.rate,color=Method))+geom_jitter(width=0.3,height = 0,size=2.5)+ylab("Type 1 Error Rate")+xlab("Sample Size")+facet_wrap(~Prop+Level)

ggplot(result,aes(factor(Prop),Err.rate,color=Method))+geom_jitter(width=0.3,height = 0,size=2.5)+ylab("Type 1 Error Rate")+xlab("True P")+facet_wrap(~Size+Level)



##Problem 3
set.seed(55)

#LRT test


for(i in 1:length(n))
  for(j in 1:length(p))
    for(k in 1:length(alpha))
    {
      true=p[j]
      critical=qchisq(1-alpha[k],df=1)
      l=0
      c=0
      for(s in 1:1000)
      {
        x=rbinom(1,n[i],true)
        estimate=x/n[i]
        low=binom.lrt(x,n[i])$lower
        up=binom.lrt(x,n[i])$upper
        if(true>=low & true<=up) c=c+1
        l=l+up-low
      }
       CP[(i-1)*12+(j-1)*2+k]=c/1000
       Len[(i-1)*12+(j-1)*2+k]=l/1000
    }

#Score test

for(i in 1:length(n))
  for(j in 1:length(p))
    for(k in 1:length(alpha))
    {
      true=p[j]
      c=0
      l=0
      for(s in 1:1000)
      {
        x=rbinom(1,n[i],true)
        low=binom.confint(x,n[i],methods = "wilson")$lower
        up=binom.confint(x,n[i],methods = "wilson")$upper
        if(true>=low & true<=up)
          c=c+1
        l=l+up-low
      }
      CP[(i-1)*12+(j-1)*2+k+60]=c/1000
      Len[(i-1)*12+(j-1)*2+k+60]=l/1000
    }


#Wald test
for(i in 1:length(n))
  for(j in 1:length(p))
    for(k in 1:length(alpha))
    {
      true=p[j]
      critical=qnorm(1-alpha[k]/2)
      c=0
      l=0
      for(s in 1:1000)
      {
        x=rbinom(1,n[i],true)
        estimate=x/n[i]
        se=sqrt(estimate*(1-estimate)/n[i])
        low=estimate-critical*se
        up=estimate+critical*se
        if(true>=low & true<=up)
          c=c+1
        l=l+up-low
      }
      CP[(i-1)*12+(j-1)*2+k+120]=c/1000
      Len[(i-1)*12+(j-1)*2+k+120]=l/1000
    }

result2=data.frame(Size,Prop,Level,Method,CP,Len)

ggplot(data=result2,aes(factor(Size),CP))+geom_jitter(size=2.5,height = 0,width = 0.3)+geom_hline(yintercept = 0.95,color="red")+geom_hline(yintercept=0.99,color="blue")+facet_grid(Level~Method)+xlab("Sample Size")+ylab("Coverage Probility")

ggplot(data=result2,aes(x=factor(Size),y=Len))+geom_jitter(size=2.5,width = 0.3,height = 0)+geom_smooth(se=FALSE)+facet_grid(Level~Method)+xlab("Sample Size")+ylab("Length")


      
      
