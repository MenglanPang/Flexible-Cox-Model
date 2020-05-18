##Example for Flexible Cox model using septic shock dataset
library(survival)
library(splines)
sepsis<-read.csv("sepsis.csv")
source("CoxFlex.R")
source("CoxFlex_BH.R")
Vars<-c("immunodep","knauss","urinaire","cirrhose","infection","age","sofa")

sepsis$ID<-1:nrow(sepsis)
sepsis<-data.frame(sepsis[,c("ID","time.obs","event",Vars)])

#####Fit the flexible Cox model###
cox_fit<-last_prog(data=sepsis, Type=c("time.obs","event"),
                         variables=Vars, 
                         NL=c(rep(0,5),rep(1,2)),TD=c(0,0,0,rep(1,4)),
                         m=1, p=2, knots=-999)

#fit: object after fitting the flexible cox model using the last_prog function
#Data: Dataset that is used to fit the flexible cox model
#nknot.bh: number of the interior knots of the splines for modeling the baseline hazard
#degree.bh: degree of the splines for modeling the baseline hazard
#knot_time: "eventtime" or "alltime"; whether the knots are placed based on the quantiles of the uncensored event time or the overal observed (event or censoring) time 
#ndivision: the number of intervals defining the granularity of the numerical computation of integrals; 
##default is "maxtime" such that the number of interval=100*floor(Observed event/censoring time) 

#####Fit the extended flexible Cox model (estimation of the baseline hazard)###
FlexCoxbh<-CoxFlex_BH(fit=cox_fit,Data=sepsis,nknot.bh=2,degree.bh=3,knot_time="eventtime",ndivision=600)
save(FlexCoxbh,file="sepsis_flexcox.rda")

cov.val<-rep(NA,7)
for (i in 1:5){
  cov.val[i]<-0
}
cov.val[6]<-min(sepsis$age)
cov.val[7]<-min(sepsis$sofa)

par(mfrow=c(2,3),mgp=c(2,1,0),mar=c(3.5,3.5,2,1),
    oma=c(0.25,0.25,0.25,0.25))

####Estimation of hazard function#####

#fit: object after fitting the extended flexible cox model using the "CoxFlex_BH" function
#time: a vector of time values at which the hazard function is estimated
#cov: the specific covariate pattern for which hazard function is estimated
#Data: Dataset that was used to fit the flexible cox model

Hazard.est<-HazardEst(fit=FlexCoxbh,time=seq(0.3,60,0.1),cov=cov.val,Data=sepsis)
plot(0,0,type="n",xlim=c(0,60),ylim=c(0,0.006),xlab="Time (day)",ylab="Hazard",cex.main=0.75,main="Hazard")
lines(Hazard.est$time,Hazard.est$hazard)
Data<-sepsis
Time.Obs<-"time.obs"
Delta<-"event"

Time.obs<-Data[Data[,Delta]==1,Time.Obs]
Time.obs.freq<-table(Time.obs)
Time.obs.uniq<-as.numeric(names(Time.obs.freq))
for (j in 1:length(Time.obs.uniq)){
  rug(Time.obs.uniq[j],ticksize=0.005*Time.obs.freq[j])
}
mtext('(a)', side=1, line=2, at=1,cex=0.75)

###Estimation of TD effects###

plot(0,0,type="n",xlim=c(0,60),ylim=c(0,6),ylab=expression(e^beta(t)),xlab="Time(Day)",cex.axis=0.8,cex.lab=0.8,cex.main=0.75,main="TD effect")

#model.FlexSurv: object after fitting the extended flexible cox model using the "CoxFlex_BH" function
#variable: variable name whose TD/NL effect is to be estimated
#TD: 1 if a TD is to be estimated; otherwise 0
#NL: 1 if a NL is to be estimated; otherwise 0
#est.range: the range of time where the TD effect is to be estimated;
#           the range of the variable where the NL is to be estimated

TD_effect<-est.FlexSurv(model.FlexSurv=FlexCoxbh, variable="cirrhose", TD=1, NL=0,est.range=c(0,60))
lines(TD_effect$axis,exp(TD_effect$estimate),lty=1,lwd=1)

TD_effect<-est.FlexSurv(model.FlexSurv=FlexCoxbh, variable="infection", TD=1, NL=0,est.range=c(0,60))
lines(TD_effect$axis,exp(TD_effect$estimate),lty=2,lwd=1)

legend("topright",c("Cirrhosis","Infection type(nosocomial)"),
       lty=1:2,
       #fill=c(rgb(0.192,0.192,0.192,0.1),rgb(0.128,0.128,0.128,0.3)),
       cex=0.65)
for (j in 1:length(Time.obs.uniq)){
  rug(Time.obs.uniq[j],ticksize=0.005*Time.obs.freq[j])
}
#abline(h=0,lty=2,col="grey")
mtext('(b)', side=1, line=2, at=1,cex=0.75)

plot(0,0,type="n",xlim=c(20,80),ylim=c(-0.8,0.6),xlab="Age",ylab=expression(g(X)),cex.axis=0.8,cex.lab=0.8,cex.main=0.75,main="NL effect")

###Estimation of NL effects###

#model.FlexSurv: object after fitting the extended flexible cox model using the "CoxFlex_BH" function
#variable: variable name whose TD/NL effect is to be estimated
#TD: 1 if a TD is to be estimated; otherwise 0
#NL: 1 if a NL is to be estimated; otherwise 0
#ref.value.NL: the reference value to which the NL effect is relatively estimated 
#est.range: the range of time values where the TD effect is to be estimated;
#           the range of the variable vallues where the NL is to be estimated
NL_effect<-est.FlexSurv(model.FlexSurv=FlexCoxbh, variable="age", TD=0, NL=1,ref.value.NL=66,est.range=c(20,80))

lines(NL_effect$axis,NL_effect$estimate,lty=1,lwd=1)
X.freq<-table(sepsis[,"age"])
X.uniq<-as.numeric(names(X.freq))
for (j in 1:length(X.uniq)){
  rug(X.uniq[j],ticksize=X.freq[j]/50)
  
}
mtext('(c)', side=1, line=2, at=20,cex=0.75)


plot(0,0,type="n",xlim=c(0,60),ylim=c(0.45,5),ylab=expression(beta(t)),xlab="Time (Day)",cex.axis=0.8,cex.lab=0.8,cex.main=0.75,main="TD effect")

TD_effect<-est.FlexSurv(model.FlexSurv=FlexCoxbh, variable="age", TD=1, NL=0,est.range=c(0,60))
lines(TD_effect$axis,TD_effect$estimate,lty=1,lwd=1)

for (j in 1:length(Time.obs.uniq)){
  rug(Time.obs.uniq[j],ticksize=0.005*Time.obs.freq[j])
}
legend("bottomright","Age",lty=1,cex=0.65)
mtext('(d)', side=1, line=2, at=1,cex=0.75)

plot(0,0,type="n",xlim=c(3,20),ylim=c(-2,2),xlab="SOFA",ylab=expression(g(X)),cex.axis=0.8,cex.lab=0.8,cex.main=0.75,main="NL effect")
NL_effect<-est.FlexSurv(model.FlexSurv=FlexCoxbh, variable="sofa", TD=0, NL=1,ref.value.NL=11,est.range=c(3,20))

lines(NL_effect$axis,NL_effect$estimate,lty=1,lwd=1)

X.freq<-table(sepsis[,"sofa"])
X.uniq<-as.numeric(names(X.freq))
for (j in 1:length(X.uniq)){
  rug(X.uniq[j],ticksize=X.freq[j]/1000)
  
}
mtext('(e)', side=1, line=2, at=4,cex=0.75)

plot(0,0,type="n",xlim=c(0,60),ylim=c(0,2),ylab=expression(beta(t)),xlab="Time (Day)",cex.axis=0.8,cex.lab=0.8,cex.main=0.75,main="TD effect")

TD_effect<-est.FlexSurv(model.FlexSurv=FlexCoxbh, variable="sofa", TD=1, NL=0,est.range=c(0,60))

lines(TD_effect$axis,TD_effect$estimate,lty=1,lwd=1)
for (j in 1:length(Time.obs.uniq)){
  rug(Time.obs.uniq[j],ticksize=0.005*Time.obs.freq[j])
}
legend("bottomright","SOFA",lty=1,cex=0.65)
mtext('(f)', side=1, line=2, at=1,cex=0.75)


####Estimation of survival curve#####
par(mfrow=c(1,1))
#fit: object after fitting the extended flexible cox model using the "CoxFlex_BH" function
#time: a vector of time values at which the survivl curve is estimated
#cov: the specific covariate pattern for which survival curve is estimated
#Data: Dataset that was used to fit the flexible cox model
Surv.est<-SurvEst(fit=FlexCoxbh,time=seq(0.3,60,0.1),cov=cov.val,Data=sepsis)
plot(0,0,type="n",xlim=c(0,60),ylim=c(0,1),xlab="Time (day)",ylab="Survival",cex.main=0.75,main="Survival")
lines(Surv.est$time,Surv.est$survival)


###linear predictor
#fit: object after fitting the extended flexible cox model using the "CoxFlex_BH" function
#Data: Dataset that was used to fit the flexible cox model
#time: a specific time point at which the linear preditor is computed
#newdata: a new dataset that contains the variables for which the linear predictor is computed
predlp<-pred.lp(fit=FlexCoxbh,Data=sepsis,time=10,newdata=sepsis)

#Loglikelihood
#fit: object after fitting the extended flexible cox model using the "CoxFlex_BH" function
#Data: Dataset that was used to fit the flexible cox model
#newdata: a new dataset that contains the variables for which the loglikelihood is computed
loglik<-CoxLoglike(fit=FlexCoxbh,Data=sepsis,newdata=sepsis)


