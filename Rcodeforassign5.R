##################### Problem 7.4###################
y=c(rep(1,44),rep(0,1543),rep(1,21),rep(0,460),rep(1,60),rep(0,1061),rep(1,19),rep(0,1123),rep(1,5),rep(0,226),rep(1,16),rep(0,841))
smoke=c(rep(1,3189),rep(0,2230))
stratum=c(rep(1,1587),rep(2,481),rep(3,1121),rep(1,1142),rep(2,231),rep(3,857))

Data<-data.frame(smoke,stratum,y)

loglink<-glm(y~.,family = binomial(link="log"),data=Data)
summary(loglink)

logloglink<-glm(y~.,family = binomial(link = cloglog),data=Data)
summary(logloglink)

probitlink<-glm(y~.,family=binomial(link=probit),data=Data)
summary(probitlink)


#####################Problem 7.10#############
DCCT<-read.csv("C:/Users/37099/OneDrive - The University of Western Ontario/course/Biostats 9510/assignment5/renal1.csv",sep=",",header = FALSE)
yearsdm=DCCT[,4]/12
DCCT=DCCT[,-4]
DCCT=cbind(micro24=DCCT[,1],int=DCCT[,2],hbael=DCCT[,3],yearsdm,sbp=DCCT[,4],female=DCCT[,5])
DCCT=data.frame(DCCT)
model <-glm(micro24~int+hbael+yearsdm+sbp+female,family = binomial,data=DCCT) 

score=estfun(model) #score vector 
u=t(score)%*%score 
va=vcov(model)#covariance matrix, inverse matrix of observed information matrix
#In logistic regression, observed I is equal to the expected I

rva=va%*%u%*%va #Sigma

var1=diag(rva) #variance of parameter
beta=model$coefficients
se=sapply(var1, sqrt)#standard error

upper=beta+1.96*se
lower=beta-1.96*se#CI
wald=beta^2/var1 #wald statistic

m1=sum(DCCT[,1])
piest=m1/length(DCCT[,1]) 
sub=DCCT[DCCT$micro24==1,]
l=rep(0,6)
mat=matrix(0,nrow=6,ncol=6)
for(i in 1:172)
{
  u0=rep(0,6)
  for(j in 1:6)
  {
    if(j==1)
      u0[j]=DCCT[i,1]-piest
     else
       u0[j]=(DCCT[i,1]-piest)*DCCT[i,j]
  }
  mat=mat+u0%*%t(u0)
}  

mat #J(theta0)

i0=matrix(0,nrow = 6,ncol=6)
for(i in 1:6)
  for(j in 1:6)
  {
    if(i==1 && j==1)
      i0[i,j]=172*piest*(1-piest)
    if(i==1 && j!=1)
      i0[i,j]=172*mean(DCCT[,j])*piest*(1-piest)
    if(i!=1 && j==1)
      i0[i,j]=172*mean(DCCT[,i])*piest*(1-piest)
    if(i!=j && (i!=1 && j!=1))
      i0[i,j]=piest*(1-piest)*sum(DCCT[,i]*DCCT[,j])
    if(i==j && j!=1)
      i0[i,j]=piest*(1-piest)*sum((DCCT[,i])^2) 
  }#information matrix under null hypothesis

unull=rep(0,6)
l1=mean(sub[,2]-mean(DCCT[,2]))
l2=mean(sub[,3]-mean(DCCT[,3]))
l3=mean(sub[,4]-mean(DCCT[,4]))
l4=mean(sub[,5]-mean(DCCT[,5]))
l5=mean(sub[,6]-mean(DCCT[,6]))
unull=c(0,l1,l2,l3,l4,l5)*m1 #score vector under null hypothesis

scoresquare=t(unull)%*%ginv(mat)%*%unull #score test statistic

Sigma=ginv(i0)%*%mat%*%ginv(i0) #Sigma_0


##################### Problem 3####################
bresmack<-read.csv("C:/Users/37099/OneDrive - The University of Western Ontario/course/Biostats 9510/assignment5/renal1.csv",sep=",")
conjest=ifelse(bresmack$dose==0,0,1)
bresmack$obese[bresmack$obese==9]=0

#(a)
con1<-clogit(case~estrogen+strata(caseset),data=bresmack)
summary(con1)

con2<-clogit(case~estrogen+gbdx+hyper+obese+nonest+strata(caseset),data=bresmack)
summary(con2)

#(b)
new=bresmack$dose
new[new==0]="none"
new[new==1]="low"
new[new==2]="middle"
new[new>2 & new<9]="high"
new[new==9]=NA
new<-C(as.factor(new),base=4)
bresmack$dose=new


con3<-clogit(case~dose+strata(caseset),data=bresmack)
summary(con3)

con4<-clogit(case~dose+ gbdx+ hyper+ obese + nonest + strata(caseset),data=bresmack)
summary(con4)

#(c)
new1<-C(as.factor(new),base=2)
bresmack$dose=new1
con5<-clogit(case~dose+strata(caseset),data=bresmack)
summary(con5)

con6<-clogit(case~dose+ gbdx+ hyper+ obese + nonest + strata(caseset),data=bresmack)
summary(con6)

#(d)

bresmack=data.frame(bresmack,conjest)
bresmack$dur[bresmack$dur==99]=NA
con7<-clogit(case~conjest:dur+strata(caseset),data=bresmack)
summary(con7)

con8<-clogit(case~conjest:dur+gbdx+ hyper+ obese + nonest + strata(caseset),data=bresmack)
summary(con8)
