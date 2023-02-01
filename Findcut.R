findcut= function(factor=NULL,outcome=NULL,cutnum=NA,datatype=c("survival","logistic"),nmin=20,segment=100)
{
  if (missing(factor))
    stop("The argument factor is missing")
  if (missing(outcome))
    stop("The argument outcome is missing")
  if (missing(datatype))
    stop("The argument datatype is missing")
  if (datatype=="survival"){
    if(!is.matrix(outcome)){
      stop("The outcome must be matrix")
    }
    if(dim(outcome)[2]!=2){
      stop("The outcome's column dimensions must be two ")
    }
  }
  if (datatype=="logistic"){
    if (!is.vector(outcome))
      stop("The argument outcome must be vector")
  }
  #if (cutnum>4)
   # stop("The argument cutnum must be less than or equal to 4")
  z=qnorm(1-(1-0.95)/2)#caculate 95% confidence interval(cut1 need)
  if (datatype == "survival"){
    delta<-data.frame(outcome)
    colnames(delta)=c("event","time")
    userdata=na.omit(data.frame(delta$event,delta$time,factor))  #Collating of survival data
    colnames(userdata)=c("event","time","factor")
    index=order(userdata$factor)
    userdata=userdata[index,]
    n=dim(userdata)[1] #number of subject
    range=range(userdata$factor) #range of continous predictor
    cutunit=diff(range)/segment
    k=1 #count for result
    start <- c()
    end <- c()
    group <- rep(0,cutnum)#dummy variable
  ###############################
     if(cutnum==1){
       start[1]=userdata$factor[nmin]
       end[1]=userdata$factor[n-nmin*cutnum]
       result=matrix(NA,ncol=11,nrow=100000)
       colnames(result)=c("Cut1","Log-rank test","Gehan-Wilcoxon test",
                             "PH_assumption1","HR1","HRlow","HRup","P1","Likelihood ratio test","Waldtest","Scoretest")

       }else{
         for(i in 2:cutnum){
         start[1]=userdata$factor[nmin]
         end[1]=userdata$factor[n-nmin*cutnum]
         start[i]=userdata$factor[which(userdata$factor>start[i-1])[nmin]] #confirmed group2 has 20 person and where cut2 to start cutting
         end[i]=userdata$factor[n-nmin*(cutnum-i+1)]
         result=matrix(NA,ncol=4*cutnum+5,nrow=100000)
         colnames(result)=c(paste("Cut",1:cutnum,sep=""),"Log-rank test","Gehan-Wilcoxon test",
                                 paste("PH_assumption",1:cutnum,sep=""),paste("HR",1:cutnum,sep=""),paste("P",1:cutnum,sep=""),"Likelihood ratio test","Waldtest","Scoretest")}
         }
      process <- function(userdata,start,end,i,cutnum,group,n,cutunit,nmin,result,k){
        if(i==cutnum)
          {
          while(start[i]<end[i])
            {
            if((i-1)==0){group[i]=sum(userdata$factor<=start[i])} else {group[i]=sum(userdata$factor<=start[i])-sum(group[1:(i-1)])}
            j <- 0:cutnum
            factor_status=c(rep(j[1:(tail(j+1,1)-1)],group[1:(tail(j+1,1)-1)]),rep(tail(j,1),n-sum(group[1:tail(j,1)])))
            coxfit=coxph(Surv(userdata$time,userdata$event) ~ factor(factor_status))   #?Xt1?F?n?LO-?PAI?O??
            model= summary(coxfit)
            coef = model$coefficients
            z=qnorm(1-(1-0.95)/2)
            lrtest_fit= survdiff(Surv(userdata$time,userdata$event)~ factor(factor_status)) #log-rank test
            wilcoxon_fit=survdiff(Surv(userdata$time,userdata$event)~factor(factor_status),rho=1)  #wilcoxon test
            if(i==1){
              result[k,]=c(start[1],1-pchisq(lrtest_fit$chisq,1), 1-pchisq(wilcoxon_fit$chisq,1)
                           ,cox.zph(coxfit)$table[,3],coef[,2],exp(coef[1] - z * coef[3]),exp(coef[1] + z * coef[3]),coef[5]
                           ,model[["logtest"]]["pvalue"],model[["waldtest"]]["pvalue"],model[["sctest"]]["pvalue"])
            }else{
              result[k,]=c(start[1:cutnum],1-pchisq(lrtest_fit$chisq,cutnum), 1-pchisq(wilcoxon_fit$chisq,cutnum)
                           ,cox.zph(coxfit)$table[,3][1:cutnum],coef[,2],coef[,5],model[["logtest"]]["pvalue"],model[["waldtest"]]["pvalue"],model[["sctest"]]["pvalue"])
            }
            start[i]=start[i]+cutunit
            k=k+1
          }
          if(cutnum==1){
            return(result)
          }else{
            return(list(start=start,group=group,result=result,k=k))
          }
        }
        else
        {
          while(start[i]<end[i])
          {
            if(i==1){group[i]=sum(userdata$factor<=start[i])} else {group[i]=sum(userdata$factor<=start[i])-sum(group[1:(i-1)])}
            processin=process(userdata,start,end,i+1,cutnum,group,n,cutunit,nmin,result,k)
            group=processin$group
            start=processin$start
            result=processin$result
            k=processin$k;
            start[i] = start[i]+cutunit
            start[i+1] = userdata$factor[sum(userdata$factor<=start[i])+nmin]
          }
          if(i!=1)
          {
            return(list(start=start,group=group,result=result,k=k))
          }
          else
          {
            return(result)
          }
        }
      }
      allcut=na.omit(data.frame(process(userdata,start,end,i=1,cutnum,group,n,cutunit,nmin,result,k)))
      if(cutnum==1){
        par(mfrow=c(1,2))
        plotdata=data.frame(cut=allcut[,"Cut1"],HR=allcut[,"HR1"],HRlow=allcut[,"HRlow"],HRup=allcut[,"HRup"],p=allcut[,"Log.rank.test"])
        plot(plotdata$cut,plotdata$HR,type="l",lwd=2,xlab="Cutoff point",ylab="Hazard ratio",main="Cutoff point and Hazard Ratio",ylim=c(min(allcut[,"HRlow"]),max(allcut[,"HR1"])))
        lines(plotdata$cut,plotdata$HRlow,lty=2)
        lines(plotdata$cut,plotdata$HRup,lty=2)
        plot(plotdata$cut,plotdata$p,type="l",lwd=2,xlab="Cutoff point",ylab="p-value",main="Cutoff point and p-value")
        phokin=which(allcut[,"PH_assumption1"]>0.05)
     }else{
       ph <- paste("PH_assumption",1:cutnum,sep="")
       phok<-matrix(data=NA,ncol=10000,nrow=4)
      for(i in 1:cutnum)
      {
        maxph <- length(which(allcut[,ph[i]]>0.05))
        for(x in 1:maxph){
          phok[i,x]=which(allcut[,ph[i]]>0.05)[x]
        }
      }
      phokin <- as.vector(phok[1:cutnum,])
      for(i in 1:(cutnum-1)){
        phokin=na.exclude(intersect(intersect(phok[i,],phokin),phok[i+1,]))
      }
      }
################################################################################
      bestcut <- matrix(NA,ncol=cutnum,nrow=1)
      colnames(bestcut)=c(paste("Cut",1:cutnum,sep=""))
      cutpoint <- c(allcut[which.min(allcut[,"Log.rank.test"]),][1:cutnum],allcut[which.min(allcut[,"Gehan.Wilcoxon.test"]),][1:cutnum],
                      allcut[phokin[which.min(allcut[phokin,"Likelihood.ratio.test"])],][1:cutnum],allcut[phokin[which.min(allcut[phokin,"Waldtest"])],][1:cutnum],
                      allcut[phokin[which.min(allcut[phokin,"Scoretest"])],][1:cutnum])
	cutmatrix <- matrix(as.numeric(cutpoint),ncol=cutnum,byrow=T)
	count=c()
	for(i in 1:5){
	  count[i]=paste(cutmatrix[i,],collapse = "")
	  }
	tablenumber <- data.frame(table(count))
	bestcut[1,] <- cutmatrix[which(count==tablenumber[which.max(tablenumber$Freq),1])[1],]
    if(cutnum==1){
      output=list(allcut=allcut,
                  logranktest=allcut[which.min(allcut[,"Log.rank.test"]),],
                  wilcoxon=allcut[which.min(allcut[,"Gehan.Wilcoxon.test"]),],
                  logtest=allcut[phokin[which.min(allcut[phokin,"Likelihood.ratio.test"])],],
                  waldtest=allcut[phokin[which.min(allcut[phokin,"Waldtest"])],],
                  scoretest=allcut[phokin[which.min(allcut[phokin,"Scoretest"])],],
                  HRabsmax=allcut[phokin[which.max(abs(allcut[phokin,"HR1"]-1))],])
                  #,bestcut=bestcut)
      }
    else {
      output=list(allcut=allcut,
                  logranktest=allcut[which.min(allcut[,"Log.rank.test"]),],
                  wilcoxon=allcut[which.min(allcut[,"Gehan.Wilcoxon.test"]),],
                  logtest_likelihood.ratio.test=allcut[phokin[which.min(allcut[phokin,"Likelihood.ratio.test"])],],
                  waldtest=allcut[phokin[which.min(allcut[phokin,"Waldtest"])],],
                  scoretest=allcut[phokin[which.min(allcut[phokin,"Scoretest"])],])
                  #,bestcut=bestcut)
       	    }
    return(output)
    }#datatype=survival end
  #----------------------------------------logit regression---------------------------------------------------------#
  if (datatype == "logistic"){
    userdata=na.omit(data.frame(outcome,factor))
    index=order(userdata$factor)
    userdata=userdata[index,]
    n=dim(userdata)[1] #number of subject
    range=range(userdata$factor) #range of continous predictor
    cutunit=diff(range)/segment
    logist_red <- glm(userdata$outcome ~1, family=binomial())# no other variable so the natural prediction is the mean of userdata$outcome ,reduce model
    k=1 #count for result
    start <- c()
    end <- c()
    group <- rep(0,cutnum)  #dummy variable
    #----------------------------------------cut---------------------------------------------------------#
    if(cutnum==1){
      start[1]=userdata$factor[nmin]
      end[1]=userdata$factor[n-nmin*cutnum]
      result=matrix(NA,ncol=13,nrow=100000)
      colnames(result)=c("Cut1","OR1","OR_up","OR_low","P1","Likelihood ratio test","Waldtest","Scoretest","Specificity","Sensitivity","Youden",
                         "Fisher's Exact Test","AUC")
        }else{
      for(i in 2:cutnum){
        start[1]=userdata$factor[nmin]
        end[1]=userdata$factor[n-nmin*cutnum]
        start[i]=userdata$factor[which(userdata$factor>start[i-1])[nmin]] #confirmed group2 has 20 person and where cut2 to start cutting
        end[i]=userdata$factor[n-nmin*(cutnum-i+1)]
        result=matrix(NA,ncol=3*cutnum+11,nrow=100000)
        colnames(result)=c(paste("Cut",1:cutnum,sep=""),paste("OR",1:cutnum,sep=""),paste("P",1:cutnum,sep=""),"Likelihood ratio test","Waldtest","Scoretest","Specificity","Sensitivity","Lower","Higher","Youden","euclidean","manhattan","AUC")
    		 }
       }
           process <- function(userdata,start,end,i,cutnum,group,n,cutunit,nmin,result,k){
        if(i==cutnum)
        {
          while(start[i]<end[i])
          {
            if(i==1){group[i]=sum(userdata$factor<=start[i])} else {group[i]=sum(userdata$factor<=start[i])-sum(group[1:(i-1)])}
            j <- 0:cutnum
            factor_status=c(rep(j[1:(tail(j+1,1)-1)],group[1:cutnum]),rep(tail(j,1),n-sum(group[1:tail(j,1)])))
            logist_fit =glm(userdata$outcome ~factor(factor_status), family=binomial())
            model= summary(logist_fit)
            coef = model$coefficients
            table<- table(userdata$outcome, factor_status)
            Fisher.test=fisher.test(table)$p.value
            roc <- roc(userdata$outcome,factor_status)
            spe=roc$specificities*100
            sen=roc$sensitivities*100
            euc <- sqrt((rep(100,cutnum+2)-spe)^2+(rep(100,cutnum+2)-sen)^2)
            manhan <- abs(rep(100,cutnum+2)-spe)+abs(sen-rep(100,cutnum+2))
            LRT=anova(logist_fit,logist_red,test="Chisq")
            Scoretest=anova(logist_fit,logist_red,test="Rao")
            z=qnorm(1-(1-0.95)/2)
            if(i==1){
              Waldtest=wald.test(b = coef(logist_fit), Sigma = vcov(logist_fit), Terms = 2)
              result[k,]=c(start[1],exp(coef[2]),exp((coef[2]+z*coef[4])),exp((coef[2]-z*coef[4])),coef[8],LRT$P[2],Waldtest$result$chi2[3],Scoretest$P[2],
                spe[2],sen[2],sen[2]+spe[2]-1,Fisher.test,as.numeric(roc$auc))
          }else{
            Waldtest=wald.test(b = coef(logist_fit), Sigma = vcov(logist_fit), Terms = 2:(cutnum+1))
            maxroc=which.max(sen+spe)
            mineuc=which.min(euc)
            minman=which.min(manhan)
            result[k,]=c(start[1:cutnum],exp(coef[2:(cutnum+1)]),coef[(3*cutnum+5):(4*cutnum+4)],LRT$P[2], Waldtest$result$chi2[3],Scoretest$P[2]
                         ,spe[maxroc],sen[maxroc],sen[maxroc]-100*z*sqrt((sen[maxroc]/100)*(1-sen[maxroc]/100)/(sum(na.exclude(outcome)==1)))
                         ,sen[maxroc]+100*z*sqrt((sen[maxroc]/100)*(1-sen[maxroc]/100)/(sum(na.exclude(outcome)==1)))
                         ,sen[maxroc]+spe[maxroc]-1,sqrt((sen[mineuc]-100)^2 + (spe[mineuc]-100)^2),abs(sen[minman]-100)+ abs(100-spe[minman]),as.numeric(roc$auc))

          }
            start[i]=start[i]+cutunit
            k=k+1
            }
          if(cutnum==1){
            return(result)
          }else{
           return(list(start=start,group=group,result=result,k=k))
          }
        }
        else
        {
          while(start[i]<end[i])
          {
            if(i==1){group[i]=sum(userdata$factor<=start[i])} else {group[i]=sum(userdata$factor<=start[i])-sum(group[1:(i-1)])}
            processin=process(userdata,start,end,i+1,cutnum,group,n,cutunit,nmin,result,k)
            group=processin$group
            start=processin$start
            result=processin$result
            k=processin$k;
            start[i] = start[i]+cutunit
            start[i+1] = userdata$factor[sum(userdata$factor<=start[i])+nmin]
          }
          if(i!=1)
          {
            return(list(start=start,group=group,result=result,k=k))
          }
          else
          {
            return(result)
          }
        }
      }

      allcut=na.omit(data.frame(process(userdata,start,end,i=1,cutnum,group,n,cutunit,nmin,result,k)))
if(cutnum==1){
  par(mfrow=c(1,2))
  plotdata=data.frame(cut=allcut[,"Cut1"],OR=allcut[,"OR1"],ORup=allcut[,"OR_up"],ORlow=allcut[,"OR_low"],p=allcut[,"Likelihood.ratio.test"])
  plot(plotdata$cut,plotdata$OR,type="l",lwd=2,xlab="Cutoff point",ylab="Odds ratio",main="Cutoff point and Odds Ratio", ylim=c(min(allcut[,"OR_low"]),max(allcut[,"OR1"])))
  lines(plotdata$cut,plotdata$ORlow,lty=2)
  lines(plotdata$cut,plotdata$ORup,lty=2)
  plot(plotdata$cut,plotdata$p,type="l",lwd=2,xlab="Cutoff point",ylab="p-value",main="Cutoff point and p-value")
  allcut$euclidean=sqrt((100-allcut$Sensitivity)^2 + (100-allcut$Specificity)^2)
  allcut$manhattan=100-allcut$Sensitivity+ 100-allcut$Specificity
     }
    #-----------------------------------------------------------------------------------------------#
      bestcut <- matrix(NA,ncol=cutnum,nrow=1)
      colnames(bestcut)=c(paste("Cut",1:cutnum,sep=""))
      if(cutnum==1){
       cutpoint <- c(allcut[which.min(allcut[,"Likelihood.ratio.test"]),][1:cutnum],allcut[which.min(allcut[,"Waldtest"]),][1:cutnum],
                      allcut[which.min(allcut[,"Scoretest"]),][1:cutnum],allcut[which.min(allcut$Fisher.s.Exact.Test),][1:cutnum],allcut[which.max(allcut[,"Youden"]),][1:cutnum],
                      allcut[which.min(allcut[,"euclidean"]),][1:cutnum],allcut[which.max(allcut[,"manhattan"]),][1:cutnum],allcut[which.max(allcut[,"AUC"]),][1:cutnum])
	}
      else{
       cutpoint <- c(allcut[which.min(allcut[,"Likelihood.ratio.test"]),][1:cutnum],allcut[which.min(allcut[,"Waldtest"]),][1:cutnum],
                      allcut[which.min(allcut[,"Scoretest"]),][1:cutnum],allcut[which.max(allcut[,"Youden"]),][1:cutnum],
                      allcut[which.min(allcut[,"euclidean"]),][1:cutnum],allcut[which.max(allcut[,"manhattan"]),][1:cutnum],allcut[which.max(allcut[,"AUC"]),][1:cutnum])
      }
	cutmatrix <- matrix(as.numeric(cutpoint),ncol=cutnum,byrow=T)
	count=c()
      if(cutnum==1){
      	for(i in 1:8){
		count[i]=paste(cutmatrix[i,],collapse = "")
		}
	}
	else{
	for(i in 1:7){
	count[i]=paste(cutmatrix[i,],collapse = "")
	}
      }
	tablenumber <- data.frame(table(count))
	bestcut[1,] <- cutmatrix[which(count==tablenumber[which.max(tablenumber$Freq),1])[1],]
      if(cutnum==1){
        output=list(allcut=allcut,
                    Likelihood.ratio.test=allcut[which.min(allcut[,"Likelihood.ratio.test"]),],
                    waldtest=allcut[which.min(allcut[,"Waldtest"]),],
                    scoretest=allcut[which.min(allcut[,"Scoretest"]),],
                    ORabsmax=allcut[which.max(abs(allcut$OR1)),],
                    youden=allcut[which.max(allcut$Youden),],
                    Fisher.test=allcut[which.min(allcut$Fisher.s.Exact.Test),],
                    euclidean=allcut[which.min(allcut$euclidean),],
                    manhattan=allcut[which.min(allcut$manhattan),],
                    AUC=allcut[which.max(allcut[,"AUC"]),])
                    #,bestcut=bestcut)
      }
      else
      {
        output=list(allcut=allcut,
                    Likelihood.ratio.test=allcut[which.min(allcut[,"Likelihood.ratio.test"]),],
                    waldtest=allcut[which.min(allcut[,"Waldtest"]),],
                    scoretest=allcut[which.min(allcut[,"Scoretest"]),],
                    youden=allcut[which.max(allcut[,"Youden"]),],
                    euclidean=allcut[which.min(allcut[,"euclidean"]),],
                    manhattan=allcut[which.min(allcut[,"manhattan"]),],
                    AUC=allcut[which.max(allcut[,"AUC"]),])
                    #,bestcut=bestcut)
      }
  }
  return(output)
}