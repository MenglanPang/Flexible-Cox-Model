### Program for flexible Cox model in time-to-event analysis
### 
### Created by Menglan 
#### Estimating the baseline hazard using regression B-splines
#### for the extended flexible Cox TD/NL model
#### Call "CoxFlex.R"

get_args_for <- function(fun, env = parent.frame(), inherits = FALSE, ..., dots) {
  potential <- names(formals(fun))
  
  if ("..." %in% potential) {
    if (missing(dots)) {
      # return everything from parent frame
      return(as.list(env))
    }
    else if (!is.list(dots)) {
      stop("If provided, 'dots' should be a list.")
    }
    
    potential <- setdiff(potential, "...")
  }
  
  # get all formal arguments that can be found in parent frame
  args <- mget(potential, env, ..., ifnotfound = list(NULL), inherits = inherits)
  # remove not found
  args <- args[sapply(args, Negate(is.null))]
  # return found args and dots
  c(args, dots)
}

CoxFlex_BH<-function(fit,Data,nknot.bh,degree.bh,knot_time="eventtime",ndivision="maxtime",init="exp"){
  Var<-fit$variables
  NL<-fit$NL
  TD<-fit$TD
  nknot.NL<-fit$nknot.NL
  nknot.TD<-fit$nknot.TD
  degree.NL<-fit$degree.NL
  degree.TD<-fit$degree.TD
  Time.Obs<-fit$Time.Obs
  Delta<-fit$Delta
  
  num.event<-sum(Data[,Delta])
  
  if (init=="exp"){
  ##get the initial value of gamma from an exponential AFT model (constant baseline hazard)
  exp_formula<-paste0("Surv(",Time.Obs,",",Delta,")~",paste0(Var,collapse = "+"))

  exponential<-survreg(formula(exp_formula),dist="exponential",data=Data)
  log_rate0<--exponential$coef[1]
  } else {log_rate0<-init}

  #Observed time
  Tt<-as.matrix(Data[,Time.Obs])
  
  #censoring indicator
  delta<-Data[,Delta]
  
  #number of covariates included in the model
  nvar<-length(Var)
  
  #number of paramter for modeling each covariate X 
  df.X<-nknot.NL+degree.NL+1
  #number of paramter for moding beta(t) for each covariate
  df.bt<-nknot.TD+degree.TD+1
  df.bt[NL==0 & TD==0]<-1
  
  
  alpha.knots<-vector("list",length=nvar)
  beta.knots<-vector("list",length=nvar)
  Ax<-vector("list",length=nvar)
  Bst<-vector("list",length=nvar)
  
  ##Create the spline basis functions for X and for t;
  #depending on the specification of NL and TD effect
  
  ###Bst are the spline basis functions for time; bs(t)
  ###Ax are the spline basis functions for X; bs(X)
  ###TD=1 and NL=1, Bs(t)*g(x)=Bst*Ax; both alpha and beta are being estimted; TD and NL effect for X
  ###TD=1 and NL=0, Bs(t)=Bst<-Bst*X; g(x)=Ax<-1; only beta is being estimated; TD effect for X
  ###TD=0 and NL=1, Bs(t)=Bst<-1; g(x)=Ax; only alpha is being estimated; NL effect for X
  ###TD=0 and NL=0; Bs(t)=Bst<-X; g(x)=Ax<-1; only beta is being estimated; constant effect for X
  
  
  for (i in 1:nvar){
    if (NL[i]==1){
      alpha.knots[[i]]<-fit$knots_covariates[i,]
      Ax[[i]]<-as.matrix(cbind(spli(Data[,Var[i]],1,2,alpha.knots[[i]]),
                               spli(Data[,Var[i]],2,2,alpha.knots[[i]]),
                               spli(Data[,Var[i]],3,2,alpha.knots[[i]]),
                               spli(Data[,Var[i]],4,2,alpha.knots[[i]])))
    } else {
      alpha.knots[[i]]<-NA;Ax[[i]]<-NA}
    
    if (TD[i]==1 & NL[i]==1){
      beta.knots[[i]]<-fit$knots_time
      Bst[[i]]<-as.matrix(cbind(spli(Tt,1,2,beta.knots[[i]]),
                                spli(Tt,2,2,beta.knots[[i]]),
                                spli(Tt,3,2,beta.knots[[i]]),
                                spli(Tt,4,2,beta.knots[[i]])))
      } else if (TD[i]==1 & NL[i]==0) {
      
        beta.knots[[i]]<-fit$knots_time
        Bst[[i]]<-as.matrix(cbind(spli(Tt,1,2,beta.knots[[i]]),
                                  spli(Tt,2,2,beta.knots[[i]]),
                                  spli(Tt,3,2,beta.knots[[i]]),
                                  spli(Tt,4,2,beta.knots[[i]])))*Data[,Var[i]]
    }else if (NL[i]==1 & TD[i]==0) 
    {
      beta.knots[[i]]<-NA;Bst[[i]]<-NA} else if (NL[i]==0 & TD[i]==0) {
      beta.knots[[i]]<-NA;Bst[[i]]<-matrix(Data[,Var[i]],ncol=1)}
  }
  
  if (knot_time=="alltime"){
    bh.knots<-c(rep(0,3+1),quantile(Tt,probs=c(1/3,2/3)),seq(max(Tt),max(Tt)+3,1))
  } else if (knot_time=="eventtime"){
    bh.knots<-c(rep(0,3+1),quantile(Tt[delta==1],probs=c(1/3,2/3)),seq(max(Tt[delta==1]),max(Tt[delta==1])+3,1))
  } 
  #calculate df for gamma
  df.gamma<-nknot.bh+degree.bh+1
  
  #assign initial value to gamma
  gammaVec<-rep(log_rate0,df.gamma)
  names(gammaVec)<-NULL
  
  df.X.pos<-df.X
  df.X.pos[is.na(df.X)]<-0
  df.X.pos<-cumsum(df.X.pos)
  
  df.bt.pos<-df.bt
  df.bt.pos[is.na(df.bt)]<-0
  df.bt.pos<-cumsum(df.bt.pos)
  
  #calculate df for alpha corresponding to each covariate X 
  pos<-0
  for (i in 1:nvar){
    if (NL[i]==1){
      pos<-pos+1
      pos<-df.X.pos[i]
    } 
  }  
  #total df for NL effects
  df.alpha<-pos

  
  #calculate df for beta corresponding to each covariate X 
  pos<-0
  for (i in 1:nvar){
    if (TD[i]==1 | (TD[i]==0 & NL[i]==0)){
      pos<-pos+1
      pos<-df.bt.pos[i]
    } 
  } 
  #total df for TD effects
   df.beta<-pos
  
   starttime<-proc.time()[3]
   
   environment(AFT.logLik.gamma) <- environment()
   environment(AFT.logLik.gamma.der) <- environment()
   fn<-AFT.logLik.gamma;gr <- AFT.logLik.gamma.der
   
   par<-gammaVec
   
   method <- "BFGS";
   lower <- -Inf; upper <- Inf;
   control <- list(maxit = 5000,fnscale=-1); hessian <- FALSE;
   
   arg_list <- get_args_for(optim, dots = list())
   update_gamma<- do.call(optim, arg_list)
   
  #update_gamma<-optim(gammaVec,fn=AFT.logLik.gamma,gr = AFT.logLik.gamma.der, 
  #                      fit,beta.knots,
  #                      nvar,NL,TD,df.gamma,Ax,Bst,df.bt,df.X,Var,
  #                      Tt,delta,bh.knots,upper,Data,
  #                      method =method,control=list(maxit = 1000,fnscale=-1))
    
    gammaVec<-update_gamma$par
    cat("gamma=",gammaVec,"\n")
    cat("gamma_con=",update_gamma$convergence,"\n")
    cat("gamma_likelihood=",update_gamma$value,"\n")
    cat("time_gamma",proc.time()[3]-starttime,"\n")
    loglikelihood<-update_gamma$value
  #calculate the total df in the estimation
  df.all<-fit$Number_of_parameters+df.gamma
  #calculate the time costs
  elapsedtime<-proc.time()[3]-starttime
  #output
  fit$spline_coef_bh<-gammaVec
  fit$Loglik<-loglikelihood
  fit$bh.knots<-bh.knots
  fit$degree.bh<-degree.bh
  fit$nknot.bh<-nknot.bh
  fit$df<-df.all
  fit$NL<-NL
  fit$TD<-TD
  fit$nknot.NL<-nknot.NL
  fit$nknot.TD<-nknot.TD
  fit$degree.NL<-degree.NL
  fit$degree.TD<-degree.TD
  fit$runtime<-unname(elapsedtime)
  fit$ndivision<-ndivision
  fit$num.events<-num.event
  return(fit)
}

############################# Loglikelihood Function ##############################################
#loglikelihood function for estimating gamma, required by the optim function
#gamma, has to be a vector and has to be the first input in the function

AFT.logLik.gamma<-function(gammaVec
                           #,fit,beta.knots,nvar,NL,TD,
                           #df.gamma,Ax,Bst,df.bt,df.X,Var,
                           #Tt,delta,bh.knots,upper,Data
                           ){
  gx<-vector("list",length=nvar)
  Bt<-vector("list",length=nvar)
  Btgx<-0
  
  alpha<-vector("list",length=nvar)
  beta<-vector("list",length=nvar)
  
  for (i in 1:nvar){
    if (NL[i]==1 & TD[i]==1){
    alpha[[i]]<-fit$coefficients_splines_NL[,i]
    beta[[i]]<-fit$coefficients_splines_TD[,i]
    } else if (NL[i]==1 & TD[i]==0){
    alpha[[i]]<-fit$coefficients_splines_NL[,i]
    beta[[i]]<-NA  
    } else if (NL[i]==0 & TD[i]==1){
      alpha[[i]]<-NA
      beta[[i]]<-fit$coefficients_splines_TD[,i] 
    }else if (NL[i]==0 & TD[i]==0){
      alpha[[i]]<-NA
      beta[[i]]<-fit$coefficients[i] 
    }
  }
  
  for (i in 1:nvar){
    if (is.na(Ax[[i]][1])) {Ax[[i]]<-matrix(1)}
    if (is.na(alpha[[i]][1])) {alpha[[i]]<-1}
    if (is.na(Bst[[i]][1])) {Bst[[i]]<-matrix(1)}
    if (is.na(beta[[i]][1])) {beta[[i]]<-1}
    if (is.na(df.bt[i])) {df.bt[i]<-1}
    if (is.na(df.X[i])) {df.X[i]<-1}
    
    Bt[[i]]<-matrix(Bst[[i]],ncol=df.bt[i],byrow=FALSE)%*%matrix(beta[[i]],ncol=1)
    gx[[i]]<-matrix(Ax[[i]],ncol=df.X[i],byrow=FALSE)%*%matrix(alpha[[i]],ncol=1)
    Btgx<-Btgx+drop(Bt[[i]])*drop(gx[[i]])
  }
  
  gamma<-gammaVec
  
  ### Compute the spline basis for the full dataset, which is need for log likelihood calculation
  BsW<-as.matrix(cbind(spli(Tt,1,3,bh.knots),
                       spli(Tt,2,3,bh.knots),
                       spli(Tt,3,3,bh.knots),
                       spli(Tt,4,3,bh.knots),
                       spli(Tt,5,3,bh.knots),
                       spli(Tt,6,3,bh.knots)))
  
  logL1<-delta%*%(Btgx+BsW%*%gamma)
  
  environment(integrat) <- environment()
  logL2<-integrat("none"
                  #,df.bt,gx,beta,gamma,df.gamma,
                  #Var,Tt,nvar,NL,TD,beta.knots,
                  #bh.knots,upper,Data
                  )
  LogLik<-logL1-sum(unlist(logL2))
  
  return(list(LogLik=LogLik))
}

#first derivative of loglikelihood rwt gamma, optional the optim function
#gamma, has to be a vector and has to be the first input in the function

AFT.logLik.gamma.der<-function(gammaVec
                               #,fit,beta.knots,nvar,NL,TD,
                               #df.gamma,Ax,Bst,df.bt,df.X,Var,
                               #Tt,delta,bh.knots,upper,Data
                               ){
  gx<-vector("list",length=nvar)
  Bt<-vector("list",length=nvar)
  Btgx<-0
  
  alpha<-vector("list",length=nvar)
  beta<-vector("list",length=nvar)
  
  for (i in 1:nvar){
    if (NL[i]==1 & TD[i]==1){
      alpha[[i]]<-fit$coefficients_splines_NL[,i]
      beta[[i]]<-fit$coefficients_splines_TD[,i]
    } else if (NL[i]==1 & TD[i]==0){
      alpha[[i]]<-fit$coefficients_splines_NL[,i]
      beta[[i]]<-NA  
    } else if (NL[i]==0 & TD[i]==1){
      alpha[[i]]<-NA
      beta[[i]]<-fit$coefficients_splines_TD[,i] 
    }else if (NL[i]==0 & TD[i]==0){
      alpha[[i]]<-NA
      beta[[i]]<-fit$coefficients[i] 
    }
  }
  
  
  for (i in 1:nvar){
    if (is.na(Ax[[i]][1])) {Ax[[i]]<-matrix(1)}
    if (is.na(alpha[[i]][1])) {alpha[[i]]<-1}
    if (is.na(Bst[[i]][1])) {Bst[[i]]<-matrix(1)}
    if (is.na(beta[[i]][1])) {beta[[i]]<-1}
    if (is.na(df.bt[i])) {df.bt[i]<-1}
    if (is.na(df.X[i])) {df.X[i]<-1}
    
    Bt[[i]]<-matrix(Bst[[i]],ncol=df.bt[i],byrow=FALSE)%*%matrix(beta[[i]],ncol=1)
    gx[[i]]<-matrix(Ax[[i]],ncol=df.X[i],byrow=FALSE)%*%matrix(alpha[[i]],ncol=1)
    Btgx<-Btgx+drop(Bt[[i]])*drop(gx[[i]])
  }
  
  gamma<-gammaVec
  
  ### Compute the spline basis for the full dataset, which is need for log likelihood calculation
  BsW<-as.matrix(cbind(spli(Tt,1,3,bh.knots),
                       spli(Tt,2,3,bh.knots),
                       spli(Tt,3,3,bh.knots),
                       spli(Tt,4,3,bh.knots),
                       spli(Tt,5,3,bh.knots),
                       spli(Tt,6,3,bh.knots)))
  
  environment(integrat) <- environment()
  Interg<-integrat("gamma"
                   #,df.bt,gx,beta,gamma,df.gamma,
                   #Var,Tt,nvar,NL,TD,beta.knots,
                   #bh.knots,upper,Data
                   )
  
  Interg.matrix<-matrix(unlist(Interg),nrow=length(Tt))
  loglik.der<-drop(delta%*%BsW-apply(Interg.matrix,2,sum))
  
  return(LogLik.der=loglik.der)
}

###hazard function at u: a matrix of time, each row corresponds to an individual, and each column corresponds to the small interval dividing the entire range 
#the function to be integrated
###a different function concerning the estimation of gamma and the calculation of the loglikelihood function
integrand<-function(u,wrt
                    #,df.bt,gx,beta,gamma,df.gamma,
                    #Var,Tt,nvar,NL,TD,beta.knots,
                    #bh.knots,Data
                    ){
  n.sam<-length(Tt)
  exp2<-matrix(NA,nrow=n.sam,ncol=ncol(u))
  
  BSt<-vector("list",length=ncol(u))   #spline basis for beta(t)
  
  BT<-vector("list",length=nvar)  #beta(t)=BSt%beta for each covariate
  for (i in 1:nvar){
    BT[[i]]<-matrix(NA,nrow=n.sam,ncol=ncol(u))
  }
  
  for (i in 1:ncol(u)){
    BSt[[i]]<-vector("list",nvar)  #BSt is a list has length ncol(u); each element BSt[[i]] is also a list that contains nvar elements for the spline basis corresponding to each TD effect
  }  
  ###new April 24  TD=1 & NL=0
  BTGX<-matrix(0,nrow=n.sam,ncol=ncol(u))  #Btgx for each time u
  for (i in 1:ncol(u)){
    for (j in 1:nvar){
      if (TD[j]==1 & NL[j]==1){
        BSt[[i]][[j]]<-as.matrix(cbind(spli(u[,i],1,2,beta.knots[[j]]),
                                       spli(u[,i],2,2,beta.knots[[j]]),
                                       spli(u[,i],3,2,beta.knots[[j]]),
                                       spli(u[,i],4,2,beta.knots[[j]])))
                } else if (TD[j]==1 & NL[j]==0) {
        BSt[[i]][[j]]<-as.matrix(cbind(spli(u[,i],1,2,beta.knots[[j]]),
                                       spli(u[,i],2,2,beta.knots[[j]]),
                                       spli(u[,i],3,2,beta.knots[[j]]),
                                       spli(u[,i],4,2,beta.knots[[j]])))*Data[,Var[j]]
      }else if (NL[j]==1 & TD[j]==0){
        BSt[[i]][[j]]<-matrix(1)} else if (NL[j]==0 & TD[j]==0) {
          BSt[[i]][[j]]<-matrix(Data[,Var[j]],ncol=1)
        }
      BT[[j]][,i]<-drop(matrix(BSt[[i]][[j]],ncol=df.bt[j],byrow=FALSE)%*%matrix(beta[[j]],ncol=1))
      BTGX[,i]<-BTGX[,i]+drop(BT[[j]][,i])*drop(gx[[j]])
    }
  }

  if (wrt=="gamma") {
    Return.ls<-vector("list",df.gamma)
    for (i in 1:df.gamma) {Return.ls[[i]]<-matrix(NA,nrow=n.sam,ncol=ncol(u))}
  }
  
  if (wrt=="none"){
    Return.ls<-vector("list",1)
  }
  
  for (i in 1:ncol(u)){
    #spline basis functions for the baseline hazard
    Bs.bh<-as.matrix(cbind(spli(u[,i],1,3,bh.knots),
                           spli(u[,i],2,3,bh.knots),
                           spli(u[,i],3,3,bh.knots),
                           spli(u[,i],4,3,bh.knots),
                           spli(u[,i],5,3,bh.knots),
                           spli(u[,i],6,3,bh.knots)))
    exp2[,i]<-exp(Bs.bh%*%gamma)
    #return values that corresponds to the derivative of the log-likelihood wrt to gamma, or 
    #or calculate the log-likelihood
  if (wrt=="gamma"){
      for (k in 1:df.gamma){
        Return.ls[[k]][,i]<-exp(BTGX[,i])*exp2[,i]*Bs.bh[,k]}
    } 
  }
  if (wrt=="none"){
    Return.ls[[1]]<-exp(BTGX)*exp2
  }
  
  return(Return.ls)
}


##use numeric integration to compute the cumulative hazard
#take advantage of matrix multiplication
integrat<-function(wrt
                   #,df.bt,gx,beta,gamma,df.gamma,
                   #Var,Tt,nvar,NL,TD,beta.knots,
                   #bh.knots,upper,Data
                   ){
  #integrate from 0 to the observed time 
  bound<-cbind(0,Tt)
  #divide the entire range to many small intervals;
  #the number of intervals is to be specified by user (ndivision argument)
  #it defines the granularity of the calculation
  #the larger the number is, the longer the estimation takes
  #we may take into account the the maximum of the event time and data granularity,i.e., how precise the event is being recorded
  #by default, the number of interval=100*floor(maximum of the observed time)
  if (ndivision=="maxtime"){
    max_obsT<-floor(max(Tt))  
    num_divide<-max_obsT*100
  } else{
    num_divide<-ndivision
  }
  
  xmatrix<-t(apply(bound,1,function(x) {seq(x[1],x[2],length=num_divide)}))
  step<-apply(xmatrix,1,function(x) (x[2]-x[1]))
  
  xmatrix<-(xmatrix+step/2)[,-ncol(xmatrix)]
  #compute the value at the middle point of each interval
  environment(integrand) <- environment()
  yvalue<-integrand(xmatrix,wrt
                    #,df.bt,gx,beta,gamma,df.gamma,
                    #Var,Tt,nvar,NL,TD,beta.knots,
                    #bh.knots,Data
                    )
  #numerical integration by matrix multiplication
  value<-lapply(yvalue,function(x) apply(x*step,1,sum))
  
  return(value)
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



############################# Estimation of Hazard Function #####################################
HazardEst<-function(fit,time,cov,Data,Var){
  NL<-fit$NL
  TD<-fit$TD
  
  Var<-fit$variables
  Time.Obs<-fit$Time.Obs
  Delta<-fit$Delta
  
  
  Tt<-Data[,Time.Obs]
  delta<-Data[,Delta]

  nvar<-length(Var)
  n.row<-nrow(Data)
  res.gamma<-fit$spline_coef_bh
  
  bh.knots<-fit$bh.knots
  
  alpha.knots<-vector("list",length=nvar)
  beta.knots<-vector("list",length=nvar)
  
  alpha<-vector("list",length=nvar)
  beta<-vector("list",length=nvar)
  
  Cov.Ax<-vector("list",length=nvar)
  Time.Bst<-vector("list",length=nvar)
  
  
  for (i in 1:nvar){
    if (NL[i]==1){
      alpha.knots[[i]]<-fit$knots_covariates[i,]
      Cov.Ax[[i]]<-as.matrix(cbind(spli(unlist(cov[i]),1,2,alpha.knots[[i]]),
                                   spli(unlist(cov[i]),2,2,alpha.knots[[i]]),
                                   spli(unlist(cov[i]),3,2,alpha.knots[[i]]),
                                   spli(unlist(cov[i]),4,2,alpha.knots[[i]])))
    } else {
      alpha.knots[[i]]<-NA;Cov.Ax[[i]]<-NA}
    
    if (TD[i]==1 & NL[i]==1){
      beta.knots[[i]]<-fit$knots_time
      alpha[[i]]<-fit$coefficients_splines_NL[,i]
      beta[[i]]<-fit$coefficients_splines_TD[,i]
      Time.Bst[[i]]<-as.matrix(cbind(spli(time,1,2,beta.knots[[i]]),
                                     spli(time,2,2,beta.knots[[i]]),
                                     spli(time,3,2,beta.knots[[i]]),
                                     spli(time,4,2,beta.knots[[i]])))
    } else if (TD[i]==1 & NL[i]==0) {
      alpha[[i]]<-NA
      beta[[i]]<-fit$coefficients_splines_TD[,i] 
      beta.knots[[i]]<-fit$knots_time
      Time.Bst[[i]]<-as.matrix(cbind(spli(time,1,2,beta.knots[[i]]),
                                     spli(time,2,2,beta.knots[[i]]),
                                     spli(time,3,2,beta.knots[[i]]),
                                     spli(time,4,2,beta.knots[[i]])))*cov[i]
    }else if (NL[i]==1 & TD[i]==0) 
    {
      alpha[[i]]<-fit$coefficients_splines_NL[,i]
      beta[[i]]<-NA  
      beta.knots[[i]]<-NA;Time.Bst[[i]]<-NA} else if (NL[i]==0 & TD[i]==0) {
        alpha[[i]]<-NA
        beta[[i]]<-fit$coefficients[i] 
        beta.knots[[i]]<-NA; Time.Bst[[i]]<-matrix(cov[i],ncol=1)}
  }
  
  Cov.Btgx<-0
  Cov.gx<-vector("list",length=nvar)
  Cov.Bt<-vector("list",length=nvar)
  
  for (i in 1:nvar){
    if (is.na(Cov.Ax[[i]][1])) {Cov.Ax[[i]]<-matrix(1)}
    if (is.na(alpha[[i]][1])) {alpha[[i]]<-1}
    if (is.na(Time.Bst[[i]][1])) {Time.Bst[[i]]<-matrix(1)}
    if (is.na(beta[[i]][1])) {beta[[i]]<-1}
    
    Cov.Bt[[i]]<-matrix(as.numeric(Time.Bst[[i]]),nrow=length(time),byrow=FALSE)%*%matrix(beta[[i]],ncol=1)
    Cov.gx[[i]]<-matrix(as.numeric(Cov.Ax[[i]]),nrow=1,byrow=FALSE)%*%matrix(alpha[[i]],ncol=1)
    Cov.Btgx<-Cov.Btgx+drop(Cov.Bt[[i]])*drop(Cov.gx[[i]])
  }
  
  
  
  BsW<-as.matrix(cbind(spli(time,1,3,bh.knots),
                       spli(time,2,3,bh.knots),
                       spli(time,3,3,bh.knots),
                       spli(time,4,3,bh.knots),
                       spli(time,5,3,bh.knots),
                       spli(time,6,3,bh.knots)))
  
  flex.hazard<-exp(c(Cov.Btgx))*exp(tcrossprod(BsW,matrix(res.gamma,nrow=1)))
  
  
  return(list("time"=time,"hazard"=flex.hazard))
}


############################# Estimation of Survival Function #####################################
S.integrand<-function(u,fit,alpha.knots,beta.knots,alpha,beta,bh.knots,cov,nvar,NL,TD,Var,Data)
{
  res.gamma<-fit$spline_coef_bh
  Cov.Ax<-vector("list",length=nvar)
  Time.Bst<-vector("list",length=nvar)
  
  
  for (i in 1:nvar){
    
    if (NL[i]==1){
      Cov.Ax[[i]]<-as.matrix(cbind(spli(unlist(cov[i]),1,2,alpha.knots[[i]]),
                                   spli(unlist(cov[i]),2,2,alpha.knots[[i]]),
                                   spli(unlist(cov[i]),3,2,alpha.knots[[i]]),
                                   spli(unlist(cov[i]),4,2,alpha.knots[[i]])))
    } else {
      Cov.Ax[[i]]<-NA}
    
    if (TD[i]==1 & NL[i]==1){
      Time.Bst[[i]]<-as.matrix(cbind(spli(u,1,2,beta.knots[[i]]),
                                     spli(u,2,2,beta.knots[[i]]),
                                     spli(u,3,2,beta.knots[[i]]),
                                     spli(u,4,2,beta.knots[[i]])))
    } else if (TD[i]==1 & NL[i]==0) {
      Time.Bst[[i]]<-as.matrix(cbind(spli(u,1,2,beta.knots[[i]]),
                                     spli(u,2,2,beta.knots[[i]]),
                                     spli(u,3,2,beta.knots[[i]]),
                                     spli(u,4,2,beta.knots[[i]])))*unlist(cov[i])
    }else if (NL[i]==1 & TD[i]==0) 
    {
      Time.Bst[[i]]<-NA} else if (NL[i]==0 & TD[i]==0) {
        Time.Bst[[i]]<-matrix(cov[i],ncol=1)}
  }
  
  Cov.Btgx<-0
  Cov.gx<-vector("list",length=nvar)
  Cov.Bt<-vector("list",length=nvar)
  
  for (i in 1:nvar){
    if (is.na(Cov.Ax[[i]][1])) {Cov.Ax[[i]]<-matrix(1)}
    if (is.na(alpha[[i]][1])) {alpha[[i]]<-1}
    if (is.na(Time.Bst[[i]][1])) {Time.Bst[[i]]<-matrix(1)}
    if (is.na(beta[[i]][1])) {beta[[i]]<-1}
    
    Cov.Bt[[i]]<-matrix(as.numeric(Time.Bst[[i]]),nrow=1,byrow=FALSE)%*%matrix(beta[[i]],ncol=1)
    Cov.gx[[i]]<-matrix(as.numeric(Cov.Ax[[i]]),nrow=1,byrow=FALSE)%*%matrix(alpha[[i]],ncol=1)
    Cov.Btgx<-Cov.Btgx+drop(Cov.Bt[[i]])*drop(Cov.gx[[i]])
  }
  
  
  BsW<-as.matrix(cbind(spli(u,1,3,bh.knots),
                       spli(u,2,3,bh.knots),
                       spli(u,3,3,bh.knots),
                       spli(u,4,3,bh.knots),
                       spli(u,5,3,bh.knots),
                       spli(u,6,3,bh.knots)))
  
  flex.hazard<-exp(c(Cov.Btgx))*exp(tcrossprod(BsW,matrix(res.gamma,nrow=1)))
  return(flex.hazard)
}


SurvEst<-function(fit,time,cov,Data){
  NL<-fit$NL
  TD<-fit$TD
  
  Var<-fit$variables
  Time.Obs<-fit$Time.Obs
  Delta<-fit$Delta
  
  
  Tt<-Data[,Time.Obs]
  delta<-Data[,Delta]
  
  bh.knots<-fit$bh.knots
  Tt<-Data[,Time.Obs]
  delta<-Data[,Delta]
  
  nvar<-length(Var)
  n.row<-nrow(Data)
  
  alpha.knots<-vector("list",length=nvar)
  beta.knots<-vector("list",length=nvar) 
  
  alpha<-vector("list",length=nvar)
  beta<-vector("list",length=nvar)
  
  
  for (i in 1:nvar){
    if (NL[i]==1){
      alpha.knots[[i]]<-fit$knots_covariates[i,]
    } else {
      alpha.knots[[i]]<-NA}
    
    if (TD[i]==1 & NL[i]==1){
      alpha[[i]]<-fit$coefficients_splines_NL[,i]
      beta[[i]]<-fit$coefficients_splines_TD[,i]
      beta.knots[[i]]<-fit$knots_time
      
    } else if (TD[i]==1 & NL[i]==0) {
      alpha[[i]]<-NA
      beta[[i]]<-fit$coefficients_splines_TD[,i] 
      beta.knots[[i]]<-fit$knots_time
      
    }else if (NL[i]==1 & TD[i]==0) 
    {
      alpha[[i]]<-fit$coefficients_splines_NL[,i]
      beta[[i]]<-NA
      beta.knots[[i]]<-NA} else if (NL[i]==0 & TD[i]==0) {
        alpha[[i]]<-NA
        beta[[i]]<-fit$coefficients[i] 
        beta.knots[[i]]<-NA}
  }
  
  
  for (i in 1:nvar){
    
    if (is.na(alpha[[i]][1])) {alpha[[i]]<-1}
    
    if (is.na(beta[[i]][1])) {beta[[i]]<-1}
    
  }
  
  flex.S<-exp(-1*sapply(time,function(x) integrate(Vectorize(S.integrand,"u"),
                                                 lower=0,upper=x,fit=fit,alpha.knots=alpha.knots,beta.knots=beta.knots,alpha=alpha,beta=beta,bh.knots=bh.knots,cov=cov,nvar=nvar,
                                                 NL=NL,TD=TD,Var=Var,Data=Data)$value))
  
  return(list("time"=time,"survival"=flex.S))
}


###linear predictor
pred.lp<-function(fit,Data,time,newdata){
  NL<-fit$NL
  TD<-fit$TD
  
  Var<-fit$variables
  Time.Obs<-fit$Time.Obs
  Delta<-fit$Delta
  
  n.row<-nrow(newdata)
  Tt<-Data[,Time.Obs]
  delta<-Data[,Delta]
  tt<-rep(time,n.row)
  
  nvar<-length(Var)
  
  
  
  alpha.knots<-vector("list",length=nvar)
  beta.knots<-vector("list",length=nvar)
  
  alpha<-vector("list",length=nvar)
  beta<-vector("list",length=nvar)
  
  Cov.Ax<-vector("list",length=nvar)
  Time.Bst<-vector("list",length=nvar)
  
  
  for (i in 1:nvar){
    if (NL[i]==1){
      alpha.knots[[i]]<-fit$knots_covariates[i,]
      Cov.Ax[[i]]<-as.matrix(cbind(spli(newdata[,Var[i]],1,2,alpha.knots[[i]]),
                                   spli(newdata[,Var[i]],2,2,alpha.knots[[i]]),
                                   spli(newdata[,Var[i]],3,2,alpha.knots[[i]]),
                                   spli(newdata[,Var[i]],4,2,alpha.knots[[i]])))
    } else {
      alpha.knots[[i]]<-NA;Cov.Ax[[i]]<-NA}
    
    if (TD[i]==1 & NL[i]==1){
      beta.knots[[i]]<-fit$knots_time
      alpha[[i]]<-fit$coefficients_splines_NL[,i]
      beta[[i]]<-fit$coefficients_splines_TD[,i]
      Time.Bst[[i]]<-as.matrix(cbind(spli(tt,1,2,beta.knots[[i]]),
                                     spli(tt,2,2,beta.knots[[i]]),
                                     spli(tt,3,2,beta.knots[[i]]),
                                     spli(tt,4,2,beta.knots[[i]])))
    } else if (TD[i]==1 & NL[i]==0) {
      alpha[[i]]<-NA
      beta[[i]]<-fit$coefficients_splines_TD[,i] 
      beta.knots[[i]]<-fit$knots_time
      Time.Bst[[i]]<-as.matrix(cbind(spli(tt,1,2,beta.knots[[i]]),
                                     spli(tt,2,2,beta.knots[[i]]),
                                     spli(tt,3,2,beta.knots[[i]]),
                                     spli(tt,4,2,beta.knots[[i]])))*newdata[,Var[i]]
    }else if (NL[i]==1 & TD[i]==0) 
    {
      alpha[[i]]<-fit$coefficients_splines_NL[,i]
      beta[[i]]<-NA  
      beta.knots[[i]]<-NA;Time.Bst[[i]]<-NA} else if (NL[i]==0 & TD[i]==0) {
        alpha[[i]]<-NA
        beta[[i]]<-fit$coefficients[i] 
        beta.knots[[i]]<-NA; Time.Bst[[i]]<-matrix(newdata[,Var[i]],ncol=1)}
  }
  
  Cov.Btgx<-0
  Cov.gx<-vector("list",length=nvar)
  Cov.Bt<-vector("list",length=nvar)
  
  for (i in 1:nvar){
    if (is.na(Cov.Ax[[i]][1])) {Cov.Ax[[i]]<-matrix(1)}
    if (is.na(alpha[[i]][1])) {alpha[[i]]<-1}
    if (is.na(Time.Bst[[i]][1])) {Time.Bst[[i]]<-matrix(1)}
    if (is.na(beta[[i]][1])) {beta[[i]]<-1}
    
    Cov.Bt[[i]]<-matrix(as.numeric(Time.Bst[[i]]),nrow=n.row,byrow=FALSE)%*%matrix(beta[[i]],ncol=1)
    Cov.gx[[i]]<-matrix(as.numeric(Cov.Ax[[i]]),nrow=n.row,byrow=FALSE)%*%matrix(alpha[[i]],ncol=1)
    Cov.Btgx<-Cov.Btgx+drop(Cov.Bt[[i]])*drop(Cov.gx[[i]])
  }
  
  return(list("time"=time,"lp"=Cov.Btgx))
}

CoxLoglike<-function(fit,Data,newdata){
  Var<-fit$variables
  Time.Obs<-fit$Time.Obs
  Delta<-fit$Delta
  ndivision<-fit$ndivision
  
  NL<-fit$NL
  nknot.NL<-fit$nknot.NL
  degree.NL<-fit$degree.NL
  TD<-fit$TD
  nknot.TD<-fit$nknot.TD
  degree.TD<-fit$degree.TD
  bh.knots<-fit$bh.knots
  df.X<-nknot.NL+degree.NL+1
  n.row<-nrow(Data)
  n.row.new<-nrow(newdata)
  #Observed time
  Tt<-as.matrix(Data[,Time.Obs])
  
  newTt<-as.matrix(newdata[,Time.Obs])
  #censoring indicator
  delta.new<-newdata[,Delta]
  
  #number of covariates included in the model
  nvar<-length(Var)
  
  #number of paramter for moding beta(t) for each covariate
  df.bt<-nknot.TD+degree.TD+1
  df.bt[NL==0 & TD==0]<-1
  
  df.bt.pos<-df.bt
  df.bt.pos[is.na(df.bt)]<-0
  df.bt.pos<-cumsum(df.bt.pos)
  
  #total df for TD effects
  df.beta<-max(df.bt.pos)
  
  df.gamma<-fit$nknot.bh+fit$degree.bh+1
  
  Ax.new<-vector("list",length=nvar)
  Bst.new<-vector("list",length=nvar)
  
  alpha.knots<-vector("list",length=nvar)
  beta.knots<-vector("list",length=nvar)
  
  alpha<-vector("list",length=nvar)
  beta<-vector("list",length=nvar)
  
  for (i in 1:nvar){
    if (NL[i]==1){
      alpha.knots[[i]]<-fit$knots_covariates[i,]
      Ax.new[[i]]<-as.matrix(cbind(spli(newdata[,Var[i]],1,2,alpha.knots[[i]]),
                                   spli(newdata[,Var[i]],2,2,alpha.knots[[i]]),
                                   spli(newdata[,Var[i]],3,2,alpha.knots[[i]]),
                                   spli(newdata[,Var[i]],4,2,alpha.knots[[i]])))
    } else {
      alpha.knots[[i]]<-NA;Ax.new[[i]]<-NA}
    
    if (TD[i]==1 & NL[i]==1){
      alpha[[i]]<-fit$coefficients_splines_NL[,i]
      beta[[i]]<-fit$coefficients_splines_TD[,i]
      
      beta.knots[[i]]<-fit$knots_time
      Bst.new[[i]]<-as.matrix(cbind(spli(newTt,1,2,beta.knots[[i]]),
                                    spli(newTt,2,2,beta.knots[[i]]),
                                    spli(newTt,3,2,beta.knots[[i]]),
                                    spli(newTt,4,2,beta.knots[[i]])))
    } else if (TD[i]==1 & NL[i]==0) {
      alpha[[i]]<-NA
      beta[[i]]<-fit$coefficients_splines_TD[,i] 
      
      beta.knots[[i]]<-fit$knots_time
      Bst.new[[i]]<-as.matrix(cbind(spli(newTt,1,2,beta.knots[[i]]),
                                    spli(newTt,2,2,beta.knots[[i]]),
                                    spli(newTt,3,2,beta.knots[[i]]),
                                    spli(newTt,4,2,beta.knots[[i]])))*newdata[,Var[i]]
    }else if (NL[i]==1 & TD[i]==0) 
    {
      alpha[[i]]<-fit$coefficients_splines_NL[,i]
      beta[[i]]<-NA  
      
      beta.knots[[i]]<-NA;Bst.new[[i]]<-NA} else if (NL[i]==0 & TD[i]==0) {
        alpha[[i]]<-NA
        beta[[i]]<-fit$coefficients[i] 
        beta.knots[[i]]<-NA;Bst.new[[i]]<-matrix(newdata[,Var[i]],ncol=1)}
  }
  
  gx.new<-vector("list",length=nvar)
  Bt.new<-vector("list",length=nvar)
  Btgx.new<-0
  
  for (i in 1:nvar){
    if (is.na(Ax.new[[i]][1])) {Ax.new[[i]]<-matrix(1)}
    if (is.na(alpha[[i]][1])) {alpha[[i]]<-1}
    if (is.na(Bst.new[[i]][1])) {Bst.new[[i]]<-matrix(1)}
    if (is.na(beta[[i]][1])) {beta[[i]]<-1}
    if (is.na(df.bt[i])) {df.bt[i]<-1}
    if (is.na(df.X[i])) {df.X[i]<-1}
    
    Bt.new[[i]]<-matrix(Bst.new[[i]],ncol=df.bt[i],byrow=FALSE)%*%matrix(beta[[i]],ncol=1)
    gx.new[[i]]<-matrix(Ax.new[[i]],ncol=df.X[i],byrow=FALSE)%*%matrix(alpha[[i]],ncol=1)
    Btgx.new<-Btgx.new+drop(Bt.new[[i]])*drop(gx.new[[i]])
  }
  
  
  gamma<-fit$spline_coef_bh
  
  
  BsW.new<-as.matrix(cbind(spli(newTt,1,3,bh.knots),
                           spli(newTt,2,3,bh.knots),
                           spli(newTt,3,3,bh.knots),
                           spli(newTt,4,3,bh.knots),
                           spli(newTt,5,3,bh.knots),
                           spli(newTt,6,3,bh.knots)))
  
  logL1<-delta.new%*%(Btgx.new+BsW.new%*%gamma)
  
  gx<-gx.new
  Tt<-newTt
  Data<-newdata
  
  environment(integrat) <- environment()
  
  logL2<-integrat("none"
                  #,df.bt,gx.new,beta,gamma,df.gamma,
                  #Var,newTt,nvar,NL,TD,beta.knots,
                  #bh.knots,upper,newdata
                  )
  LogLik<-logL1-sum(unlist(logL2))
  
  return(list(LogLik=LogLik))
}
