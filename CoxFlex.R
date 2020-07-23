###
### Flexible modeling program
###

# Original code by Willy Wynant
# New version by Yishu Wang
# Date of current version: February 26, 2018 

# Small modifications done by Marie-Eve Beauchamp (MEB)
# Last update: September 17, 2018


# Reference:
# Wynant W, Abrahamowicz M. Impact of the model building strategy on the inference 
# about time-dependent and # non-linear covariate effects in survival analysis. 
# Stat Med. 2014 Aug;33(19):3318-37.

# Functions in this program:
  # backward_selection2 
  # CoxFlex
  # DvlpMatrix
  # last_prog
  # lines.FlexSurv
  # plot.FlexSurv
  # spli


library(survival)
library(splines)


backward_selection2 <- function(data,Type,variables,continuous, TD, NL,m=1,p=2,alpha_back=0.05,knots=-999){

  V<-length(variables) # number of total variables
  
  covB<-variables[order(1-continuous)] 
  Ftd<-TD[order(1-continuous)]
  Fnl<-NL[order(1-continuous)]
  
  TDf<-rep(1,V) 
  NLf<-c(rep(1,sum(continuous)),rep(0,sum(1-continuous)))
  nNLf<-sum(continuous) 
  nNLfNOT=V-nNLf 

  m1nl<-matrix(nrow=nNLf,ncol=nNLf,1)
  diag(m1nl)<-0
  m2nl<-matrix(nrow=nNLf,ncol=nNLfNOT,0)
  m3nl<-cbind(m1nl,m2nl)
  m4nl<-matrix(nrow=V,ncol=V,c(rep(1,nNLf),rep(0,nNLfNOT)),byrow=T)
  mNL<-rbind(NLf,m3nl,m4nl)
  
  ## TD matrix
  m1td<-matrix(nrow=nNLf,ncol=V,1)
  m2td<-matrix(nrow=nNLfNOT,ncol=nNLf,1)
  m3td<-matrix(nrow=nNLfNOT,ncol=nNLfNOT,1)
  diag(m3td)<-0
  m4td<-cbind(m2td,m3td)
  m5td<-matrix(nrow=nNLf,ncol=nNLf,1)
  diag(m5td)<-0
  m6td<-matrix(nrow=nNLf,ncol=nNLfNOT,1)
  m7td<-cbind(m5td,m6td)
  mTD<-rbind(TDf,m1td,m4td,m7td)
  
 
  for (i in 1:V){
    if(Ftd[i]==1) mTD[,i]<-1
    else if(Ftd[i]==-1) mTD[,i]<-0
    
    if(Fnl[i]==1) mNL[,i]<-1
    else if(Fnl[i]==-1) mNL[,i]<-0
  }
  
  PLback<-rep(0,V+nNLf+1) 
  DFback<-rep(-999,V+nNLf+1) 
  
  for (i in 1:(V+nNLf+1)){
    res<-CoxFlex(data,Type,variables=covB,TD=mTD[i,],NL=mNL[i,],m,p,knots=-999)
    PLback[i]<-res$Partial_Log_Likelihood 
    DFback[i]<-res$Number_of_parameters
  }
  
  
  Mind<-matrix(nrow=4,ncol=V,c(rep(1,2*V),rep(1,nNLf),rep(0,nNLfNOT),rep(2,nNLf),rep(1,nNLfNOT)),byrow=TRUE)
  for (i in 1:V){
    if(Ftd[i]==-1) Mind[2,i]<-0 # 2nd index TD 
    if(Fnl[i]==-1) Mind[3,i]<-0 # 3rd index NL 
  }
  Mind[4,]<-Mind[2,]+Mind[3,]
  
  
  pB<-rep(0,(V+nNLf)) # stores the test results of likelihood ratio tests for each model 
  
  for (j in 1:(V+nNLf)){
    if(DFback[1]-DFback[j+1]!=0)
      pB[j]<-1-pchisq(-2*(PLback[j+1]-PLback[1]),DFback[1]-DFback[j+1])
  }
  
  
  MpB<-matrix(nrow=5,ncol=V,0)
  
  MpB[4,1:nNLf]<-pB[1:nNLf] 
  if(sum(1-continuous)!=0){
  MpB[2,(nNLf+1):V]<-pB[(nNLf+1):V] 
  }
  MpB[5,1:nNLf]<-pB[(V+1):(V+nNLf)] 
  
  
  a<-which.max(MpB) 
  nc<-floor(a/5)+((a-5*floor(a/5))!=0) 
  nr<-a-5*floor(a/5)+5*((a-5*floor(a/5))==0)

  indexrow<-1                                       
  
  
  MatResPvalue<-matrix(nrow=2*length(variables)+sum(continuous),ncol=1)   
  MatResVariables<-matrix(nrow=2*length(variables)+sum(continuous),ncol=1) 
  
  if(MpB[nr,nc]<alpha_back){ 
    mbase<-CoxFlex(data,Type,variables=covB,TD=mTD[1,],NL=mNL[1,],m=1,p=2,knots=-999)
  }
  
  while(MpB[nr,nc]>=alpha_back) {     
    
    MatResPvalue[indexrow,1]<-MpB[nr,nc]       
    MatResVariables[indexrow,1]<-paste(nr,nc,sep="") 
    indexrow<-indexrow+1                              
    
    # print("indexrow",\n)
    
    if (nr==1) Mind[1,nc]<-0 
    if (nr==2 | nr==5) Mind[2,nc]<-0 
    if (nr==3 | nr==4) Mind[3,nc]<-0 
    Mind[4,]<-Mind[2,]+Mind[3,] 
    
    mbase<-CoxFlex(data,Type,variables=covB[Mind[1,]==1],TD=Mind[2,][Mind[1,]==1],NL=Mind[3,][Mind[1,]==1],m=1,p=2,knots=-999)
    
    MpB<-matrix(nrow=5,ncol=V,rep(0,V*5))
    
    
    for (k in 1:V){ 
      
      if (Mind[2,k]==1 & Mind[3,k]==0) { 
        MindNew<-Mind
        if(Ftd[k]!=1) {
          MindNew[2,k]<-0 
          mnew<-CoxFlex(data,Type,variables=covB[MindNew[1,]==1],TD=MindNew[2,][MindNew[1,]==1],NL=MindNew[3,][MindNew[1,]==1],m=1,p=2,knots=-999)
          MpB[2,k]<-1-pchisq(-2*(mnew$Partial_Log_Likelihood-mbase$Partial_Log_Likelihood),mbase$Number_of_parameters-mnew$Number_of_parameters)
         
        }
        else MpB[2,k]<-0
      }
      
      if (Mind[2,k]==1 & Mind[3,k]==1) {
        MindNew<-Mind
        if(Ftd[k]!=1){
          MindNew[2,k]<-0 
          mnew<-CoxFlex(data,Type,variables=covB[MindNew[1,]==1],TD=MindNew[2,][MindNew[1,]==1],NL=MindNew[3,][MindNew[1,]==1],m=1,p=2,knots=-999)
          MpB[5,k]<-1-pchisq(-2*(mnew$Partial_Log_Likelihood-mbase$Partial_Log_Likelihood),mbase$Number_of_parameters-mnew$Number_of_parameters)
          
        }
        else MpB[5,k]<-0
      }
      
      if (Mind[3,k]==1 & Mind[2,k]==0) { #if NL + non-TD
        MindNew<-Mind
        if(Fnl[k]!=1){
          MindNew[3,k]<-0 
          mnew<-CoxFlex(data,Type,variables=covB[MindNew[1,]==1],TD=MindNew[2,][MindNew[1,]==1],NL=MindNew[3,][MindNew[1,]==1],m=1,p=2,knots=-999)
          MpB[3,k]<-1-pchisq(-2*(mnew$Partial_Log_Likelihood-mbase$Partial_Log_Likelihood),mbase$Number_of_parameters-mnew$Number_of_parameters)
         
        }
        else MpB[3,k]<-0
      }
      
      if (Mind[3,k]==1 & Mind[2,k]==1) { 
        MindNew<-Mind
        if(Fnl[k]!=1){
          MindNew[3,k]<-0 
          mnew<-CoxFlex(data,Type,variables=covB[MindNew[1,]==1],TD=MindNew[2,][MindNew[1,]==1],NL=MindNew[3,][MindNew[1,]==1],m=1,p=2,knots=-999)
          MpB[4,k]<-1-pchisq(-2*(mnew$Partial_Log_Likelihood-mbase$Partial_Log_Likelihood),mbase$Number_of_parameters-mnew$Number_of_parameters)
          
        }
        else MpB[4,k]<-0
      }
      
      if (Mind[3,k]==0 & Mind[2,k]==0 & Mind[1,k]==1) { 
        MindNew<-Mind
        if(Fnl[k]!=1 & Fnl[k]!=-1 & Ftd[k]!=1 & Ftd[k]!=-1){
          MindNew[1,k]<-0 # remove this variable 
          mnew<-CoxFlex(data,Type,variables=covB[MindNew[1,]==1],TD=MindNew[2,][MindNew[1,]==1],NL=MindNew[3,][MindNew[1,]==1],m=1,p=2,knots=-999)
          MpB[1,k]<-1-pchisq(-2*(mnew$Partial_Log_Likelihood-mbase$Partial_Log_Likelihood),mbase$Number_of_parameters-mnew$Number_of_parameters)
         
        }
        else MpB[1,k]<-0
      }
    } 
    
    a<-which.max(MpB)
    nc<-floor(a/5)+((a-5*floor(a/5))!=0)
    nr<-a-5*floor(a/5)+5*((a-5*floor(a/5))==0)
  }
  
  list(final_model=mbase)
}





CoxFlex <- function(data, Type, variables, TD, NL, m, p, knots){
  
  res_principal<-last_prog(data, Type, variables, TD, NL, m, p, knots)
  
  TDtestP<-rep(-999,length(variables))
  TDtestN<-rep(-999,length(variables))
  TDtest<-rep(-999,length(variables))
  NLtestP<-rep(-999,length(variables))
  NLtestN<-rep(-999,length(variables))
  NLtest<-rep(-999,length(variables))
  
  for (kt in 1:length(variables)){
    if (TD[kt]==1){
      TDnew<-TD
      TDnew[kt]<-0
      resSEC<-last_prog(data, Type, variables, TDnew, NL, m, p, knots)
      TDtestP[kt]<-resSEC$Partial_Log_Likelihood
      TDtestN[kt]<-resSEC$Number_of_parameters
      TDtest[kt]<-1-pchisq(-2*(TDtestP[kt]-res_principal$Partial_Log_Likelihood),res_principal$Number_of_parameters-TDtestN[kt])
    }
    if (NL[kt]==1){
      NLnew<-NL
      NLnew[kt]<-0
      resSEC<-last_prog(data, Type, variables, TD, NLnew, m, p,knots)
      NLtestP[kt]<-resSEC$Partial_Log_Likelihood
      NLtestN[kt]<-resSEC$Number_of_parameters
      NLtest[kt]<-1-pchisq(-2*(NLtestP[kt]-res_principal$Partial_Log_Likelihood),res_principal$Number_of_parameters-NLtestN[kt])
    }
  }
 
  
  #   cat("\n")
  #   cat("Call: \n")
  #   cat("Cox(formula=Surv(")
  #   
  #   if (length(Type)==3){
  #     cat(Type[1])
  #     cat(",")
  #     cat(Type[2])
  #     cat(",")
  #     cat(Type[3])}
  #   
  #   if (length(Type)==2){
  #     cat(Type[1])
  #     cat(",")
  #     cat(Type[2])
  #   }
  #   
  #   
  #   cat(")~")
  #   for (yu in 1:(length(variables)-1)){
  #     if (TD[yu]==0 & NL[yu]==0){
  #       cat(variables[yu]) 
  #       cat("+")
  #     }
  #     if (TD[yu]==0 & NL[yu]==1){
  #       cat("NL(") 
  #       cat(variables[yu]) 
  #       cat(")+")
  #     }
  #     if (TD[yu]==1 & NL[yu]==0){
  #       cat("TD(") 
  #       cat(variables[yu]) 
  #       cat(")+")
  #     }
  #     if (TD[yu]==1 & NL[yu]==1){
  #       cat("NL(") 
  #       cat(variables[yu]) 
  #       cat(")+")
  #       cat("TD(") 
  #       cat(variables[yu]) 
  #       cat(")+")}
  #   }
  #   if (TD[length(variables)]==0 & NL[length(variables)]==0){
  #     cat(variables[length(variables)]) 
  #     cat(")")
  #   }
  #   if (TD[length(variables)]==0 & NL[length(variables)]==1){
  #     cat("NL(") 
  #     cat(variables[length(variables)]) 
  #     cat("))")
  #   }
  #   if (TD[length(variables)]==1 & NL[length(variables)]==0){
  #     cat("TD(") 
  #     cat(variables[length(variables)]) 
  #     cat("))")
  #   }
  #   if (TD[length(variables)]==1 & NL[length(variables)]==1){
  #     cat("NL(") 
  #     cat(variables[length(variables)]) 
  #     cat(")+")
  #     cat("TD(") 
  #     cat(variables[length(variables)]) 
  #     cat("))")
  #   }
  #   cat("\n")
  #   cat("Using splines of degree ")
  #   cat(p)
  # #   cat(" and ")
  # #   cat(m)
  # #   if (m==1) {cat(" knot")} else {cat(" knots")}
  # #   cat("\n")
  # #   cat("\n")
  # #   
  
  Rescox<-matrix(ncol=6,nrow=sum(2*(NL+TD==2)+1*(NL+TD!=2)) )
  colnames(Rescox)<-c("","coef","exp(coef)","se(coef)","z","p")
  rownames(Rescox)<-rep(c(""),sum(2*(NL+TD==2)+1*(NL+TD!=2)))
  
  
  indexrescox<-0
  for (yu in 1:length(variables)){
    indexrescox<-indexrescox+1
    if (TD[yu]==0 & NL[yu]==0){
      Rescox[indexrescox,1]<-variables[yu]
      Rescox[indexrescox,2]<-round(as.numeric(res_principal$coefficients[yu]),3)
      Rescox[indexrescox,3]<-round(as.numeric(exp(res_principal$coefficients[yu])),3)
      Rescox[indexrescox,4]<-round(as.numeric(res_principal$Standard_Error[yu]),3)
      Rescox[indexrescox,5]<-round(as.numeric(res_principal$coefficients[yu]/res_principal$Standard_Error[yu]),3)
      Rescox[indexrescox,6]<-round(as.numeric(1-pchisq((res_principal$coefficients[yu]/res_principal$Standard_Error[yu])^2,df=1)) ,3)
    }
    if (TD[yu]==0 & NL[yu]==1){
      Rescox[indexrescox,1]<-paste("NL(",variables[yu],")",sep="")
      Rescox[indexrescox,2]<-"---"
      Rescox[indexrescox,3]<-"splines"
      Rescox[indexrescox,4]<-"---"
      Rescox[indexrescox,5]<-"---"
      Rescox[indexrescox,6]<-round(NLtest[yu],3)
    }
    if (TD[yu]==1 & NL[yu]==0){
      Rescox[indexrescox,1]<-paste("TD(",variables[yu],")",sep="")
      Rescox[indexrescox,2]<-"---"
      Rescox[indexrescox,3]<-"splines"
      Rescox[indexrescox,4]<-"---"
      Rescox[indexrescox,5]<-"---"
      Rescox[indexrescox,6]<-round(TDtest[yu],3)
    }
    if (TD[yu]==1 & NL[yu]==1){
      Rescox[indexrescox,1]<-paste("NL(",variables[yu],")",sep="")
      Rescox[indexrescox,2]<-"---"
      Rescox[indexrescox,3]<-"splines"
      Rescox[indexrescox,4]<-"---"
      Rescox[indexrescox,5]<-"---"
      Rescox[indexrescox,6]<-round(NLtest[yu],3)
      indexrescox<-indexrescox+1
      Rescox[indexrescox,1]<-paste("TD(",variables[yu],")",sep="")
      Rescox[indexrescox,2]<-"---"
      Rescox[indexrescox,3]<-"splines"
      Rescox[indexrescox,4]<-"---"
      Rescox[indexrescox,5]<-"---"
      Rescox[indexrescox,6]<-round(TDtest[yu],3)
    }
  }
  #   print(Rescox,quote=F)
  #   cat("\n")
  #   
  #   cat("Partial log-likelihood: ")
  #   cat(res_principal$Partial_Log_Likelihood)
  #   cat("\n")
  #   cat("\n")
  #   cat("Number of events: ")
  #   cat(res_principal$Number_events)
  #   cat("\n")
  #   cat("\n")
  #   cat("Number of parameters to estimate in the model: ")
  #   cat(res_principal$Number_of_parameters)
  #   cat("\n")
  #   cat("\n")
  
  # list(model = res_principal, coef=suppressWarnings(as.numeric(Rescox[,2])),var=suppressWarnings((as.numeric(Rescox[,2]))^2),pvalue=suppressWarnings(as.numeric(Rescox[,6])))
  
  list(Partial_Log_Likelihood = res_principal$Partial_Log_Likelihood, 
       Number_of_parameters = res_principal$Number_of_parameters, 
       Number_events=res_principal$Number_events, 
       Number_knots = res_principal$Number_knots, 
       Degree_of_splines = res_principal$Degree_of_splines, 
       knots_covariates = res_principal$knots_covariates, 
       knots_time = res_principal$knots_time, 
       coefficients = res_principal$coefficients, 
       Standard_Error=res_principal$Standard_Error,
       coefficients_splines_NL = res_principal$coefficients_splines_NL,
       coefficients_splines_TD = res_principal$coefficients_splines_TD,
       variables=res_principal$variables,
       # Modif by MEB: correction of the wrong variance in the next line
       #coef=suppressWarnings(as.numeric(Rescox[,2])),var=suppressWarnings((as.numeric(Rescox[,2]))^2),pvalue=suppressWarnings(as.numeric(Rescox[,6])))
       coef=suppressWarnings(as.numeric(Rescox[,2])),var=suppressWarnings((as.numeric(Rescox[,4]))^2),pvalue=suppressWarnings(as.numeric(Rescox[,6])))

}


DvlpMatrix <- function(data, listeT, ncol, TypeNEW) {
  data <- matrix(data, ncol = ncol) # data is a list of full data on each individual, first transform this list to the matrix 
  if (max(data[, TypeNEW[2]]) < min(listeT[listeT != 0])) {  
    XX <- data  ## if max(stop time) is less then min(event time), then data remain the same
  }
  else { ## if max(stop time of certain individual) > min(event time)
    aindex <- rep(0, (sum(listeT <= max(data[, TypeNEW[2]])) - 
                        1))
    # apb <- rep(0, (sum(listeT <= max(data[, TypeNEW[2]])) - 
    #                  1))
    for (i in 1:(sum(listeT <= max(data[, TypeNEW[2]])) -  ## for the length of unique event time 
                 1)) {
      for (j in 1:(dim(data)[1])) { #for each row of certain individual, if start time < ith order of event time<=stop time of row j, stores
        # the largest row satisfy this condition in aindex[i]
        if (as.numeric(data[j, TypeNEW[1]]) < as.numeric(listeT[1 + i]) & as.numeric(data[j, TypeNEW[2]]) >= as.numeric(listeT[1 + i]))
          aindex[i] <- j
      }
    }   ## xx first coloum repeat ID, coloum for start day assigns c(0,unique event time-last one), colum for stop day assigns unique event time 
    XX <- matrix(nrow = sum(listeT <= max(data[, TypeNEW[2]])) - 
                   1, ncol = ncol)
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), 1] <- rep(data[1, 
                                                                      1], sum(listeT <= max(data[, TypeNEW[2]])) - 1)
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), TypeNEW[1]] <- listeT[1:(sum(listeT <= 
                                                                                      max(data[, TypeNEW[2]])) - 1)]
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), TypeNEW[2]] <- listeT[2:(sum(listeT <= 
                                                                                      max(data[, TypeNEW[2]])))]
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), TypeNEW[3]] <- c(rep(0,(sum(listeT <= max(data[, TypeNEW[2]])) - 2)), data[dim(data)[1],TypeNEW[3]])
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), -c(1, 
                                                          TypeNEW[1], TypeNEW[2], TypeNEW[3])] <- as.matrix(data[aindex, 
                                                                                                                 -c(1, TypeNEW[1], TypeNEW[2], TypeNEW[3])])
    
  }
  X <- XX

  list(X)
}





last_prog <- function (data, Type, variables, TD, NL, m, p, knots=-999) {
  
  if (length(Type)==2){Type2<-Type
  data$StartV0<-rep(0,dim(data)[1])
  Type<-c("StartV0",Type2[1],Type2[2])
  } 
  
  i1 <- sum((NL + TD) == 0)   
  i2 <- sum(((NL == 1) & (TD == 0))) 
  i3 <- sum(((NL == 0) & (TD == 1)))  
  i4 <- sum((NL + TD) == 2) 
  nonpara <- TD + NL 
  
  
  V <- length(variables) # number of total covariates
  
  variablesNEW<-match(variables,names(data)) 
  TypeNEW<-match(Type,names(data)) 
  
  
  listeprobaquantile <- seq(1, m)/(m + 1) 
  knotsNEW <- matrix(nrow = V + 1, ncol = p + 1 + m + p + 1) 
  

  if (is.matrix(knots)==FALSE){ ## if not user defined knots 
    if (is.numeric(knots)==TRUE & knots==-999) { ## if the knots is set as default 
      for (i in 1:V) {
        knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p + 1), ## place the first p+1 exterior knots at the min(variables[i])
                           quantile(data[, variables[i]], probs = listeprobaquantile), ## place the interior knots at quantile 
                           seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1)) ### place the last p+1 exterior knots equally spaced between (max(variables[i]), ...+p)
      }
      knotsNEW[V + 1, ] <- c(rep(0, p + 1), ## place the exterior knots for time at 0
                             quantile(data[data[, TypeNEW[3]] == 1, TypeNEW[2]], probs = listeprobaquantile), ### place the interior knots at quantile
                             seq(max(data[data[,TypeNEW[3]] == 1, TypeNEW[2]]), max(data[data[, TypeNEW[3]] ==1, TypeNEW[2]]) + p, 1)) 
    }                                                                          
  }
 
  if (is.matrix(knots)==TRUE){ ## if user defined interior knots
    
    if (dim(knots)[1]!=(length(variables)+1) | dim(knots)[2]!=m) stop("Error Message: variable knots should be a matrix of dimension (length(variables)+1)*m")
    
    for (i in 1:V) {
      if (is.na(knots[i,1])==TRUE){  knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p + 1),   
                                                        quantile(data[, variables[i]], probs = listeprobaquantile), 
                                                        seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1))
      ## if user defined knots for certain variable is NA, then use default knots setup
      } else {
        
        knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p + 1),   
                           knots[i,], 
                           seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1))
        ## else use user defined interior knots 
      }
    }
    
    if (is.na(knots[(V+1),1])==TRUE){  knotsNEW[(V+1), ] <- c(rep(0, p + 1), 
                                                              quantile(data[data[, TypeNEW[3]] == 1, TypeNEW[2]], probs = listeprobaquantile), 
                                                              seq(max(data[data[,TypeNEW[3]] == 1, TypeNEW[2]]), max(data[data[, TypeNEW[3]] ==1, TypeNEW[2]]) + p, 1))
     } else {
      ## else use default knots for time
      knotsNEW[(V+1), ] <- c(rep(0, p + 1), 
                             knots[(V+1),], 
                             seq(max(data[data[,TypeNEW[3]] == 1, TypeNEW[2]]), max(data[data[, TypeNEW[3]] ==1, TypeNEW[2]]) + p, 1))
      
    }
  }
  
 
  data <- as.matrix(data)
  listeT <- c(0, sort(unique(data[data[, TypeNEW[3]] == 1, TypeNEW[2]]))) 
  ncol <- dim(data)[2] # number of coloum 

  
  X <- split(data, data[, 1]) 
  matX <- sapply(X, DvlpMatrix, listeT = listeT, ncol = ncol, TypeNEW=TypeNEW) 
  QWR <- do.call(rbind, matX) 
  
  nbNL <- sum(NL)
  if (nbNL != 0) {  
    for (i in 1:nbNL) { 
      QWR <- cbind(QWR,  splineDesign(knotsNEW[seq(1,V, 1)[NL == 1][i],], x=QWR[, variablesNEW[NL == 1][i]],ord=p+1)[,-1])
    }
  }
  
  QWR <- cbind(QWR, splineDesign(knotsNEW[V+1,], x=QWR[, TypeNEW[2]],ord=p+1))
  
  tt <- paste("modX<-coxph(Surv(QWR[,", TypeNEW[1], "],QWR[,", 
              TypeNEW[2], "],QWR[,", TypeNEW[3], "])~", sep = "")
 
  
  Nn <- 0
  nbpara <- sum((NL == 0 & TD == 0)) 
  
  if (nbpara != 0) { 
    for (k in 1:nbpara) {
      tt <- paste(tt, "QWR[,", variablesNEW[NL == 0 & TD == 
                                              0][k], "]+", sep = "")
      Nn <- Nn + 1
    }
  }
  
  
  nbonlyNL <- sum((NL == 1 & TD == 0)) 
  if (nbonlyNL != 0) {
    onlyNL <- match(variablesNEW[NL == 1 & TD == 0], variablesNEW[NL == 
                                                                    1]) 
    for (k in 1:nbonlyNL) {
      covp<-paste( "QWR[,", (dim(data)[2] + (onlyNL[k]-1)*(m+p)+1):(dim(data)[2] + onlyNL[k]*(m+p))  , "]", sep = "")
      tt<-paste(tt, paste(c(covp,""), collapse= "+") )
      Nn <- Nn +m+p
    }
  }
  
  
  nbonlyTD <- sum((NL == 0 & TD == 1)) # no NL only TD
  if (nbonlyTD != 0) {
    for (k in 1:nbonlyTD) {
      flag<-dim(QWR)[2]+1
      QWR <- cbind(QWR, QWR[, variablesNEW[NL == 0 & TD ==1][k]] * QWR[, (dim(data)[2] +(m+p)*nbNL 
                                                                          +1): (dim(data)[2] +(m+p)*(nbNL+1)+1)]) 
      covp<-paste("QWR[,", flag: dim(QWR)[2], "]", sep = "" )
      tt<-paste(tt, paste(c(covp,""), collapse= "+"))  
      Nn <- Nn + m+p+1
    }#
  }
  
  tt2 <- tt
  
  
  nbNLTD <- sum((NL == 1 & TD == 1)) 
  NOTonlyNL <- match(variablesNEW[NL == 1 & TD == 1], variablesNEW[NL ==  
                                                                     1])
  if (nbNLTD != 0) { 
    for (k in 1:nbNLTD) {
      covp<-paste("QWR[,",(dim(data)[2] + (NOTonlyNL[k]-1)*(m+p)+1):(dim(data)[2] + NOTonlyNL[k]*(m+p)) , "]", sep = "" )
      tt2<-paste(tt2, paste(c(covp,""), collapse= "+")) 
    }
    vrais <- c()
    modX <- coxph(eval(parse(text = substr(tt2, 13, nchar(tt2) - 
                                             1))), method = "efron")
    vrais <- c(vrais, modX$loglik[2])
    
    
    mod<-modX ## estimates of NL effects
    tt2<-tt
    VV<-matrix(ncol=nbNLTD*(m+p+1),nrow=dim(QWR)[1])
    
    for (k in 1:nbNLTD){
      
      VV[,((1+m+p)*(k-1)+1):((1+m+p)*k)]<-QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                              +(NOTonlyNL[k]-1)*(m+p)+1):(dim(data)[2]+NOTonlyNL[k]*(m+p))]%*%mod$coef[(Nn+(k-1)*(m+p)+1):(Nn+k*(m+p))])
      covp<-paste( "VV[,", ((m+p+1)*(k-1)+1):((k)*(m+p+1)), "]", sep = "")
      tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
      
    }
    modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
    vrais<-c(vrais,modX$loglik[2])
   
    
    diff<-1
    #w<-1
    #w<-0
    while(diff>0.00001){ ## ACE algorithm until converge 
      #cat(w,"\n")
      
      #if (w%%2==0){
      mod<-modX
      tt2<-tt
      VV<-matrix(ncol=nbNLTD*(m+p),nrow=dim(QWR)[1])
      for (k in 1:nbNLTD){
        VV[,((m+p)*(k-1)+1):((m+p)*(k))]<-QWR[,(dim(data)[2]+(NOTonlyNL[k]-1)*(m+p)+1):(dim(data)[2]+NOTonlyNL[k]*(m+p))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                          +(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))]%*%mod$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+(k)*(m+p+1))])
        covp<-paste( "VV[,", ((m+p)*(k-1)+1):((k)*(m+p)), "]", sep = "")
        tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
      }
      
      modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
      vrais<-c(vrais,modX$loglik[2])
     

      mod<-modX
      tt2<-tt
      VV<-matrix(ncol=nbNLTD*(m+p+1),nrow=dim(QWR)[1])
      
 
      for (k in 1:nbNLTD){
        VV[,((1+m+p)*(k-1)+1):((1+m+p)*(k))]<-QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                  +(NOTonlyNL[k]-1)*(m+p)+1):(dim(data)[2]+NOTonlyNL[k]*(m+p))]%*%mod$coef[(Nn+(k-1)*(m+p)+1):(Nn+(k)*(m+p))])
        covp<-paste( "VV[,", ((m+p+1)*(k-1)+1):((k)*(m+p+1)), "]", sep = "")
        tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
      }
     
      
      modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
      vrais<-c(vrais,modX$loglik[2])
      diff<-abs(vrais[length(vrais)]-vrais[length(vrais)-2])
      #w<-w+1
      #}
      
    }
    
 
  }    else {
    modX <- coxph(eval(parse(text = substr(tt2, 13, nchar(tt2) - 
                                             1))), method = "efron")
    vrais <- c(modX$loglik[2])
  }
  rm(QWR, X, matX)
  gc() # Added by MEB
  
  MAT <- matrix(ncol = V, nrow = 1 + m + p + m + p + 1) 
  if (i1 != 0) { #if no NL and TD
    for (j in 1:i1) {
      MAT[1, j] <- modX$coef[j]
    }
  }
  if (i2 != 0) {
    for (j in 1:i2) {
      MAT[2:(m + p + 1), i1 + j] <- modX$coef[(i1 + (j - 
                                                       1) * (m + p) + 1):(i1 + (j - 1) * (m + p) + m + 
                                                                            p)]
    }
  }
  if (i3 != 0) {
    for (j in 1:i3) {
      MAT[(m + p + 2):(2 * m + 2 * p + 2), i1 + i2 + j] <- modX$coef[(i1 + 
                                                                        i2 * (m + p) + (j - 1) * (m + p + 1) + 1):(i1 + 
                                                                                                                     i2 * (m + p) + (j - 1) * (m + p + 1) + m + p + 
                                                                                                                     1)]
    }
  }
  if (i4 != 0) {
    for (j in 1:i4) {
      MAT[2:(m + p + 1), i1 + i2 + i3 + j] <- mod$coef[(i1 + 
                                                          i2 * (m + p) + i3 * (m + p + 1) + (j - 1) * (m + 
                                                                                                         p) + 1):(i1 + i2 * (m + p) + i3 * (m + p + 1) + 
                                                                                                                    (j - 1) * (m + p) + m + p)]
      MAT[(m + p + 2):(2 * (m + p + 1)), i1 + i2 + i3 + 
            j] <- modX$coef[(i1 + i2 * (m + p) + i3 * (m + 
                                                         p + 1) + (j - 1) * (m + p + 1) + 1):(i1 + i2 * 
                                                                                                (m + p) + i3 * (m + p + 1) + (j - 1) * (m + p + 
                                                                                                                                          1) + m + p + 1)]
    }
  }
  
  
  MATse <- matrix(ncol = V, nrow = 1 + m + p + m + p + 1) 
  if (i1 != 0) {
    for (j in 1:i1) {
      MATse[1, j] <- sqrt(diag(modX$var)[j])
    }
  }
  if (i2 != 0) {
    for (j in 1:i2) {
      MATse[2:(m + p + 1), i1 + j] <- sqrt(diag(modX$var)[(i1 + (j - 
                                                                   1) * (m + p) + 1):(i1 + (j - 1) * (m + p) + m + 
                                                                                        p)])
    }
  }
  if (i3 != 0) {
    for (j in 1:i3) {
      MATse[(m + p + 2):(2 * m + 2 * p + 2), i1 + i2 + j] <- sqrt(diag(modX$var)[(i1 + 
                                                                                    i2 * (m + p) + (j - 1) * (m + p + 1) + 1):(i1 + 
                                                                                                                                 i2 * (m + p) + (j - 1) * (m + p + 1) + m + p + 
                                                                                                                                 1)])
    }
  }
  if (i4 != 0) {
    for (j in 1:i4) {
      MATse[2:(m + p + 1), i1 + i2 + i3 + j] <- sqrt(diag(mod$var)[(i1 + 
                                                                      i2 * (m + p) + i3 * (m + p + 1) + (j - 1) * (m + 
                                                                                                                     p) + 1):(i1 + i2 * (m + p) + i3 * (m + p + 1) + 
                                                                                                                                (j - 1) * (m + p) + m + p)])
      MATse[(m + p + 2):(2 * (m + p + 1)), i1 + i2 + i3 + 
              j] <- sqrt(diag(modX$var)[(i1 + i2 * (m + p) + i3 * (m + 
                                                                     p + 1) + (j - 1) * (m + p + 1) + 1):(i1 + i2 * 
                                                                                                            (m + p) + i3 * (m + p + 1) + (j - 1) * (m + p + 
                                                                                                                                                      1) + m + p + 1)])
    }
  }
  
  
  var_order <- c(variablesNEW[NL == 0 & TD == 0], variablesNEW[NL == 
                                                                 1 & TD == 0], variablesNEW[NL == 0 & TD == 1], variablesNEW[NL == 
                                                                                                                               1 & TD == 1])
  
  coefficients<-MAT[1,match(variablesNEW,var_order)] 
  names(coefficients)<-variables
  
  se_coef<-MATse[1,match(variablesNEW,var_order)] 
  
  coefficients_splines_NL<-as.matrix(MAT[2:(m+p+1),match(variablesNEW,var_order)])
  coefficients_splines_NL<-rbind(rep(0,V),coefficients_splines_NL)
  coefficients_splines_NL[1,(NL==0)]<-NA
  colnames(coefficients_splines_NL)<-variables
  
  coefficients_splines_TD<-as.matrix(MAT[(m+p+2):(2*(m+p+1)),match(variablesNEW,var_order)])
  colnames(coefficients_splines_TD)<-variables

  knots_covariates<-knotsNEW[1:V,]
  if (V>1) {rownames(knots_covariates)<-variables}
  if (V>1) {knots_covariates[(NL==0),]<-rep(NA,p + 1 + m + p + 1)} else { if (NL==0) {knots_covariates<-rep(NA,p + 1 + m + p + 1)} }
  knots_time<-knotsNEW[V+1,]
  
  nEvents<-sum(data[,TypeNEW[3]]==1) 
  
  rm(data, modX)
  gc()
  
  nknot.NL<-NL
  degree.NL<-NL
  nknot.NL[NL==1]<-m
  degree.NL[NL==1]<-p
  
  nknot.TD<-TD
  degree.TD<-TD
  nknot.TD[TD==1]<-m
  degree.TD[TD==1]<-p
  
  list(Partial_Log_Likelihood = vrais[length(vrais)], Number_of_parameters = nbpara + 
         nbonlyNL * (m + p) + nbonlyTD * (1 + m + p) + nbNLTD * 
         (m + p + m + p + 1), Number_events=nEvents, Number_knots = m, Degree_of_splines = p, 
       knots_covariates = knots_covariates, 
       knots_time = knots_time, 
       coefficients = coefficients, Standard_Error=se_coef,coefficients_splines_NL = coefficients_splines_NL,coefficients_splines_TD = coefficients_splines_TD,variables=variables,
       NL=NL,TD=TD,nknot.NL=nknot.NL,nknot.TD=nknot.TD,degree.NL=degree.NL,degree.TD=degree.TD,Time.Obs=Type2[1],Delta=Type2[2])
}




lines.FlexSurv <- function(model.FlexSurv,variable,TD,NL,TimePoint=-999,ref.value.NL=0,...){
  
  variableNEW<-match(variable,model.FlexSurv$variables)

  # NOTE: line added by MEB because was causing error messages when only one variable was modeled as NL 
  #       (model.FlexSurv$knots_covariates treated as a matrix below at 2 places)
  if (is.vector(model.FlexSurv$knots_covariates)){model.FlexSurv$knots_covariates <- matrix(model.FlexSurv$knots_covariates, nrow=1)}
  
  tdestim1<-function(x){
    tdfct<-0
    for (k in 1:(model.FlexSurv$Number_knots+model.FlexSurv$Degree_of_splines+1)){
      tdfct<-tdfct+model.FlexSurv$coefficients_splines_TD[k,variableNEW]*spli(x,k,model.FlexSurv$Degree_of_splines,model.FlexSurv$knots_time)
    }
    return(tdfct)
  }
  
  nlestim1<-function(y){
    nlfct<-0
    for (k in 1:(model.FlexSurv$Number_knots+model.FlexSurv$Degree_of_splines+1)){
      nlfct<-nlfct+model.FlexSurv$coefficients_splines_NL[k,variableNEW]*spli(y,k,model.FlexSurv$Degree_of_splines,model.FlexSurv$knots_covariates[variableNEW,])
    }
    return(nlfct)
  }
  
  if(TD==1){
    # Modif by MEB: seq by 1 instead of 0.01 because this is precise enough to have a smooth line in most analyses,
      # e.g. 100 units of follow-up imply 100 points. Seq by 0.01 was time consuming and create very large graphs
      # with normal follow-up, e.g. 365 days implied 36,500 points and 10 years or daily units implied 365,000 points.
  #axist<-seq(0,model.FlexSurv$knots_time[model.FlexSurv$Degree_of_splines+2+model.FlexSurv$Number_knots],0.01)
  axist<-seq(0,model.FlexSurv$knots_time[model.FlexSurv$Degree_of_splines+2+model.FlexSurv$Number_knots])
  }
  if(NL==1){
  axisx<-seq(model.FlexSurv$knots_covariates[variableNEW,1],model.FlexSurv$knots_covariates[variableNEW,model.FlexSurv$Degree_of_splines+2+model.FlexSurv$Number_knots],0.01)
  }
  
  if (TD==1 & NL==0){
    lines(axist,tdestim1(axist),...)
  }
  
  if (TD==0 & NL==1){
    lines(axisx,nlestim1(axisx)-nlestim1(ref.value.NL),...)
  }
  
  if (TD==1 & NL==1 & TimePoint!=-999){
    lines(axisx,(nlestim1(axisx)-nlestim1(ref.value.NL))*tdestim1(TimePoint),...)
  }
  
}


####Add by Menglan
##Start here
est.FlexSurv <- function(model.FlexSurv,variable,TD,NL,TimePoint=-999,ref.value.NL=0,est.range){
  
  variableNEW<-match(variable,model.FlexSurv$variables)
  
  # NOTE: line added by MEB because was causing error messages when only one variable was modeled as NL 
  #       (model.FlexSurv$knots_covariates treated as a matrix below at 2 places)
  if (is.vector(model.FlexSurv$knots_covariates)){model.FlexSurv$knots_covariates <- matrix(model.FlexSurv$knots_covariates, nrow=1)}
  
  tdestim1<-function(x){
    tdfct<-0
    for (k in 1:(model.FlexSurv$Number_knots+model.FlexSurv$Degree_of_splines+1)){
      tdfct<-tdfct+model.FlexSurv$coefficients_splines_TD[k,variableNEW]*spli(x,k,model.FlexSurv$Degree_of_splines,model.FlexSurv$knots_time)
    }
    return(tdfct)
  }
  
  nlestim1<-function(y){
    nlfct<-0
    for (k in 1:(model.FlexSurv$Number_knots+model.FlexSurv$Degree_of_splines+1)){
      nlfct<-nlfct+model.FlexSurv$coefficients_splines_NL[k,variableNEW]*spli(y,k,model.FlexSurv$Degree_of_splines,model.FlexSurv$knots_covariates[variableNEW,])
    }
    return(nlfct)
  }
  
  if(TD==1){
    # Modif by MEB: seq by 1 instead of 0.01 because this is precise enough to have a smooth line in most analyses,
    # e.g. 100 units of follow-up imply 100 points. Seq by 0.01 was time consuming and create very large graphs
    # with normal follow-up, e.g. 365 days implied 36,500 points and 10 years or daily units implied 365,000 points.
    #axist<-seq(0,model.FlexSurv$knots_time[model.FlexSurv$Degree_of_splines+2+model.FlexSurv$Number_knots],0.01)
    axist<-seq(est.range[1],est.range[2],0.1)
  }
  if(NL==1){
    axisx<-seq(est.range[1],est.range[2],0.01)
  }
  
  if (TD==1 & NL==0){
    axis<-axist
    estimate<-tdestim1(axist)
  }
  
  if (TD==0 & NL==1){
    axis<-axisx
    estimate<-nlestim1(axisx)-nlestim1(ref.value.NL)
  }
  
  if (TD==1 & NL==1 & TimePoint!=-999){
    axis<-axisx
    estimate<-(nlestim1(axisx)-nlestim1(ref.value.NL))*tdestim1(TimePoint)
  }
  return(list("axis"=axis,"estimate"=estimate))
}
###Finish Here###


plot.FlexSurv <- function(model.FlexSurv,variable,TD,NL,TimePoint=-999,ref.value.NL=0,...){
  
  variableNEW<-match(variable,model.FlexSurv$variables)
  
  # NOTE: line added by MEB because was causing error messages when only one variable was modeled as NL 
  #       (model.FlexSurv$knots_covariates treated as a matrix below at 2 places)
  if (is.vector(model.FlexSurv$knots_covariates)){model.FlexSurv$knots_covariates <- matrix(model.FlexSurv$knots_covariates, nrow=1)}

  tdestim1<-function(x){
    tdfct<-0
    for (k in 1:(model.FlexSurv$Number_knots+model.FlexSurv$Degree_of_splines+1)){
      tdfct<-tdfct+model.FlexSurv$coefficients_splines_TD[k,variableNEW]*spli(x,k,model.FlexSurv$Degree_of_splines,model.FlexSurv$knots_time)
    }
    return(tdfct)
  }
  
  nlestim1<-function(y){
    nlfct<-0
    for (k in 1:(model.FlexSurv$Number_knots+model.FlexSurv$Degree_of_splines+1)){
      nlfct<-nlfct+model.FlexSurv$coefficients_splines_NL[k,variableNEW]*spli(y,k,model.FlexSurv$Degree_of_splines,model.FlexSurv$knots_covariates[variableNEW,])
    }
    return(nlfct)
  }
  
  if(TD==1){
    # Modif by MEB: seq by 1 instead of 0.01 because this is precise enough to have a smooth line in most analyses,
      # e.g. 100 units of follow-up imply 100 points. Seq by 0.01 was time consuming and create very large graphs
      # with normal follow-up, e.g. 365 days implied 36,500 points and 10 years or daily units implied 365,000 points.
    #axist<-seq(0,model.FlexSurv$knots_time[model.FlexSurv$Degree_of_splines+2+model.FlexSurv$Number_knots],0.01)
    axist<-seq(0,model.FlexSurv$knots_time[model.FlexSurv$Degree_of_splines+2+model.FlexSurv$Number_knots])
  }

  if(NL==1){
    axisx<-seq(model.FlexSurv$knots_covariates[variableNEW,1],model.FlexSurv$knots_covariates[variableNEW,model.FlexSurv$Degree_of_splines+2+model.FlexSurv$Number_knots],0.01)
  }
  
  if (TD==1 & NL==0){
    plot(axist,tdestim1(axist),...)
  }
    
  if (TD==0 & NL==1){
    plot(axisx,nlestim1(axisx)-nlestim1(ref.value.NL),...)
  }  
  
  if (TD==1 & NL==1 & TimePoint==-999){
    op<-par(mfrow=c(2,1))
    plot(axist,tdestim1(axist),...)
    plot(axisx,nlestim1(axisx)-nlestim1(ref.value.NL),...)
    par(op)
  }
  
  if (TD==1 & NL==1 & TimePoint!=-999){
    plot(axisx,(nlestim1(axisx)-nlestim1(ref.value.NL))*tdestim1(TimePoint),...)
  }
}





spli <- function(x, j, p, knots) {
        if (p == 0) {
            b <- ifelse(x >= knots[j] & x < knots[j + 1], 1, 
                0)
            return(b)
        }
        else {
            a1 <- ifelse(rep(knots[j] != knots[j + p], length(x)), 
                (x - knots[j])/(knots[j + p] - knots[j]), 0)
            a2 <- ifelse(rep(knots[j + p + 1] != knots[j + 1], 
                length(x)), (knots[j + p + 1] - x)/(knots[j + 
                p + 1] - knots[j + 1]), 0)
            return(a1 * spli(x, j, p - 1, knots) + a2 * spli(x, 
                j + 1, p - 1, knots))
        }
    }

