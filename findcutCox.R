findcutCox<-function(target,event,time,confound=NULL,numcut,initial_rr=NULL,initial_cut=NULL,initial_domain=NULL,numgen,gap=NULL){
  libraries1=c("rgenoud","survival","foreach","doParallel","doRNG","xtable")
  lapply(libraries1, library, quietly = TRUE, character.only = TRUE)
  confound <<- confound
  if(is.null(confound)){
    userdata<-na.omit(data.frame(target,event,time))
  }else{
    userdata<-na.omit(data.frame(target,event,time,confound))
    confound<-data.frame(userdata[,4:dim(userdata)[2]])
  }
  target<-userdata$target
  censor<-userdata$event
  Ti<-userdata$time
  N<-dim(userdata)[1]
  numcut<<-numcut
  
  gen<-maxloglik(target,numcut,Ti,censor,confound,numgen,domain_range(target,initial_domain,numcut,confound),initial(target,initial_rr,initial_cut,numcut),gap)
  group<-cut(target,breaks=c(0,gen$par_corr[1:numcut],max(target)),labels = 0:numcut,quantile=FALSE)
  if(is.null(confound)){
    fit<-coxph(Surv(Ti,censor)~factor(group))
    pvalue<-as.numeric(summary(fit)$coefficients[,5])
    # Wald Test
    varcov = fit$var
    overallstatisice.cross <- as.numeric(summary(fit)[["logtest"]]["test"])
    cutpvalue<-as.numeric(summary(fit)$coefficients[,5])
    beta<-as.numeric(fit$coefficients)
    cut<-gen$par_corr
    return(list(like=fit$loglik[2],cutpvalue=cutpvalue,beta=beta,
                overallstatisice.cross=overallstatisice.cross,
                chiWald=varcov,cut=cut))
  }else{
    if(length(levels(group))==1){
      group<-group
    }else{
      group<-factor(group)
    }
    fit<-coxph(Surv(Ti,censor)~group+as.matrix(confound))
    fit.re <- coxph(Surv(Ti,censor)~as.matrix(confound))
    #Wald Test
    varcov = fit$var
    overallstatisic <- fit$loglik[2]-fit.re$loglik[2]
    cutpvalue<-as.numeric(summary(fit)$coefficients[,5])
    beta<-as.numeric(fit$coefficients)
    cut<-gen$par_corr
    return(list(loglik=fit$loglik[2],cutpvalue=cutpvalue,beta=beta,
                overallstatisic=overallstatisic,
                chiWald=varcov,cut=cut))
  }
}

aictest<-function(target,event,time,confound=NULL,totalcut=totalcut,initial_rr=NULL,initial_cut=NULL,initial_domain=NULL,numgen=numgen,gap=NULL,numboot){
  confound <<- confound
  if(is.null(confound)){
    userdata<-na.omit(data.frame(target,event,time))
  }else{
    userdata<-na.omit(data.frame(target,event,time,confound))
    confound<-userdata[,4:dim(userdata)[2]]
  }
  target<-userdata$target
  censor<-userdata$event
  Ti<-userdata$time
  N<-dim(userdata)[1]
  aicfunction<-function(numcut,numgen,numboot){
    set1 <- numboot[1:round(N/2)]
    set2 <- numboot[(round(N/2)+1):N] #cross validation
    numcut <<- numcut
    
    generate1 <- maxloglik(target[set1],numcut,Ti[set1],censor[set1],confound[set1],numgen,domain_range(target[set1],initial_domain,numcut,confound[set1]),initial(target[set1],initial_rr,initial_cut,numcut),gap)
    generate2 <- maxloglik(target[set2],numcut,Ti[set2],censor[set2],confound[set2],numgen,domain_range(target[set2],initial_domain,numcut,confound[set2]),initial(target[set2],initial_rr,initial_cut,numcut),gap)
    cut1 <- generate1$par_corr
    cut2 <- generate2$par_corr
    group1 <- cut(target[set1],breaks = c(0,cut2[1:numcut],max(target[set1])),labels = c(0:numcut),quantile=FALSE)
    group2 <- cut(target[set2],breaks = c(0,cut1[1:numcut],max(target[set2])),labels = c(0:numcut),quantile=FALSE)
    if(is.null(confound)){
      fit <- coxph(Surv(time[c(set1,set2)],censor[c(set1,set2)])~factor(c(group1,group2))+strata(rep(c(1,2),c(round(N/2),N-round(N/2)))))
      overallstatistic.cross <- as.numeric(summary(fit)[["logtest"]]["test"])
    }else{
      fit <- coxph(Surv(time[c(set1,set2)],censor[c(set1,set2)])~factor(c(group1,group2))+strata(rep(c(1,2),c(round(N/2),N-round(N/2))))+as.matrix(confound[c(set1,set2)]))
      fit.re <- coxph(Surv(time[c(set1,set2)],censor[c(set1,set2)])~strata(rep(c(1,2),c(round(N/2),N-round(N/2))))+as.matrix(confound[c(set1,set2)]))
      overallstatistic.cross <- fit$loglik[2]-fit.re$loglik[2]
    }
    #Wald Test
    varcov = fit$var
    cutpvalue<-as.numeric(summary(fit)$coefficients[,5])
    beta<-as.numeric(fit$coefficients)
    return(list(aic=(-2*fit$loglik[2]+2*numcut) ,loglik=fit$loglik[2],cutpvalue=cutpvalue,beta=beta,
                cut1=cut1,cut2=cut2,overallstatistic.cross=overallstatistic.cross,
                chiWald=varcov))
  }
  catchwrongaic <- function (numcut,baby,rand) {
    out <- tryCatch(aicfunction(numcut,baby,rand), 
                    error = function(e) mget(c("aic","loglik","pvalue","beta","cut1","cut2","overall","wald"),new.env(),
                                             ifnotfound=as.list(list(c("wrong"),c("wrong"),rep(c("wrong"),numcut),
                                                                     rep(c("wrong"),numcut),
                                                                     rep(c("wrong"),2*numcut),
                                                                     c("wrong"),c("wrong"))))
    )
    return(out)
  }
  return(sapply(c(1:totalcut),aicfunction,numgen=numgen,numboot=numboot))
}

maxloglik<-function(target,numcut,time,censor,confound=NULL,baby,domain,initial,gap=NULL){
  time_global<<-time
  censor_global<<-censor
  target_global<<-target
  nvars<<-2*numcut+dim(data.frame(confound))[2]
  if(is.null(confound)){
    confound_global <<- NULL
  }else{
    confound_global<<-as.matrix(confound)
  }
  numcut_global<<-numcut
  target_max<<-max(target)
  if(is.null(gap)){
    gap_global <<- quantile(sort(diff(sort(na.omit(target)))),probs = 0.5)
  }else{
    gap_global<<-gap
  }
  ccc<-genoud(obj, nvars, max=TRUE, pop.size=100, max.generations=baby, wait.generations=10,
              hard.generation.limit=TRUE, starting.values=initial, MemoryMatrix=TRUE,
              Domains=domain, solution.tolerance=0.001,print.level = 0,
              gr=NULL, boundary.enforcement=2, lexical=FALSE, gradient.check=TRUE)
  ccc$par_corr<-ccc$par #the coefficients of genoud
  ccc$par_corr[1:numcut]<-sort(ccc$par[1:numcut])+seq(0,gap_global*(numcut-1),by=gap_global) #sort cutpoint
  return(ccc)
}

initial<-function(target,initial_rr=NULL,initial_cut=NULL,numcut){
  if(is.null(initial_cut)||is.null(initial_rr)){
    incut <- quantile(target,probs = seq(0,1,1/(numcut+1)))
    initial <- c(incut[2:(numcut+1)],initial_rr)
  }else{
    incut<-array(matrix(c(initial_cut),ncol=3,byrow=TRUE),dim=c(1,3))
    inrr<-array(matrix(c(initial_rr),ncol=3,byrow=TRUE),dim=c(1,3))
    initial <- c(as.numeric(unlist(incut[1,numcut])),as.numeric(unlist(inrr[1,numcut])))
  }
  return(initial)
}

domain_range<-function(target,initial_domain=NULL,numcut,confound=NULL){
  if(is.null(confound)){
    vars<-numcut
  }else{
    vars<-numcut+dim(data.frame(confound))[2]
  }
  if(is.null(initial_domain)){
    cutdomain <- c(quantile(target,probs = 0.1),quantile(target,probs = 0.99))
    rrdomain <- rep(c(-10,10),vars)
    ran<-if(numcut==1){
      c(1,2)
    }else if(numcut==2){
      c(1,2,1,2)
    }else{rep(1:2,times=numcut)}
    indomain <- matrix(c(cutdomain[ran],rrdomain),ncol = 2,byrow = TRUE)
  }else{
    indomain <- initial_domain[[numcut]]
  }
  return(indomain)
}

obj <- function(xx){  
  cutoff <- xx[1:numcut_global] #cutpoint
  cut_design <- cut(target_global,breaks=c(0,sort(cutoff)+seq(0,gap_global*(length(cutoff)-1),by=gap_global),target_max),quantile=FALSE,labels=c(0:numcut_global))
  beta <- xx[(numcut+1):nvars]  #coefficients of parameters
  if(is.null(confound)){
    beta1 <- coxph(Surv(time_global,censor_global)~cut_design, init = beta,iter.max=0)#iteration is zero
  }else{
    beta1 <- coxph(Surv(time_global,censor_global)~cut_design+confound_global, init = beta,iter.max=0)#iteration is zero
  }
  return(beta1$loglik[2])
}