library("survival")
findcutnum = function (factor,outcome,datatype,nmin=20,segment=100)
{
  if (missing(factor))
    stop("The argument factor is missing")
  if (missing(outcome))
    stop("The argument outcome is missing")
  if (missing(datatype))
    stop("The argument datatype is missing")
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
      stop("The argument outcome must be  vector")
  }
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
p=1
coxfit = coxph(Surv(userdata$time,userdata$event) ~ userdata$factor)
originAIC=(-2*coxfit$loglik[2])+2*p
#originAICc=AIC+(2*p*(p+1)/(n-p-1))


cut1AIC=c()
start1=userdata$factor[nmin] #confirmed group1 has 20 person and where cut1 to start cutting
end1=userdata$factor[n-nmin]#cut1 end
while(start1<end1)
         {
	        group1=sum(userdata$factor<=start1) #nummber of group1
		  factor_status=c(rep(0,group1),rep(1,n-group1))
              coxfit = coxph(Surv(userdata$time,userdata$event) ~ factor(factor_status))
              AIC=(-2*coxfit$loglik[2])+2*p
              AICc=AIC+(2*p*(p+1)/(n-p-1))
              cut1AIC=c(cut1AIC,AIC)
          	  start1=start1+cutunit
          }
#-------------------------------------cut2------------------------------------------------------------#
p=2
cut2AIC=c()
start1=userdata$factor[nmin] #confirmed group1 has 20 person and where cut1 to start cutting
end1=userdata$factor[n-nmin*2]#cut1 end

start2=userdata$factor[which(userdata$factor>start1)[nmin]] #confirmed group2 has 20 person and where cut2 to start cutting
end2=userdata$factor[n-nmin] #cut2 end

while(start1<end1)
    {
	group1=sum(userdata$factor<=start1) #nummber of group1
          while(start2<end2)
     		 {
       	     group2=sum(userdata$factor<=start2)-group1
                 factor_status=c(rep(0,group1),rep(1,group2),rep(2,n-group1-group2))
                 coxfit = coxph(Surv(userdata$time,userdata$event) ~ factor(factor_status))
                 AIC=(-2*coxfit$loglik[2])+2*p
                 AICc=AIC+(2*p*(p+1)/(n-p-1))
                 cut2AIC=c(cut2AIC,AIC)
                 start2=start2+cutunit
             }
          start1 = start1+cutunit
          start2 = userdata$factor[sum(userdata$factor<=start1)+nmin]
    }
#-------------------------------------cut3------------------------------------------------------------#
p=3
cut3AIC=c()
start1=userdata$factor[nmin] #confirmed group1 has 20 person and where cut1 to start cutting
end1=userdata$factor[n-nmin*3]#cut1 end

start2=userdata$factor[which(userdata$factor>start1)[nmin]] #confirmed group2 has 20 person and where cut2 to start cutting
end2=userdata$factor[n-nmin*2] #cut2 end

start3=userdata$factor[which(userdata$factor>start2)[nmin]] #confirmed group2 has 20 person and where cut2 to start cutting
end3=userdata$factor[n-nmin]

while(start1<end1)
    {
	group1=sum(userdata$factor<=start1) #nummber of group1
     while(start2<end2)
          {
              group2=sum(userdata$factor<=start2)-group1
           while(start3<end3)
                 {
               	  group3=sum(userdata$factor<=start3)-group1-group2
               	  factor_status=c(rep(0,group1),rep(1,group2),rep(2,group3),rep(3,n-group1-group2-group3))
               	  coxfit=coxph(Surv(userdata$time,userdata$event) ~ factor(factor_status))
                    AIC=(-2*coxfit$loglik[2])+2*p
                    AICc=AIC+(2*p*(p+1)/(n-p-1))
                    start3=start3+cutunit
                    cut3AIC=c(cut3AIC,AIC)
          	   	 }
       		start2 = start2+cutunit
       		start3 = userdata$factor[sum(userdata$factor<=start2)+nmin]
                   }
			start1 = start1+cutunit
			start2 = userdata$factor[sum(userdata$factor<=start1)+nmin]
               }
#-----------------------------------------cut4---------------------------------------------#
p=4
cut4AIC=c()
start1=userdata$factor[nmin] #confirmed group1 has 20 person and where cut1 to start cutting
end1=userdata$factor[n-nmin*4]#cut1 end

start2=userdata$factor[which(userdata$factor>start1)[nmin]] #confirmed group2 has 20 person and where cut2 to start cutting
end2=userdata$factor[n-nmin*3] #cut2 end

start3=userdata$factor[which(userdata$factor>start2)[nmin]] #confirmed group2 has 20 person and where cut2 to start cutting
end3=userdata$factor[n-nmin*2]

start4=userdata$factor[which(userdata$factor>start3)[nmin]] #confirmed group2 has 20 person and where cut2 to start cutting
end4=userdata$factor[n-nmin]

while(start1<end1)
    {
	group1=sum(userdata$factor<=start1) #nummber of group1
      while(start2<end2)
          {
              group2=sum(userdata$factor<=start2)-group1
              while(start3<end3)
                 {
                   group2=sum(userdata$factor<=start2)-group1
                       while(start4<end4)
                      {
                          group4=sum(userdata$factor<=start4)-group1-group2-group3
               	         factor_status=c(rep(0,group1),rep(1,group2),rep(2,group3),rep(3,group4),rep(4,n-group1-group2-group3-group4))
               	         coxfit=coxph(Surv(userdata$time,userdata$event) ~ factor(factor_status))
                           AIC=(-2*coxfit$loglik[2])+2*p
                           AICc=AIC+(2*p*(p+1)/(n-p-1))
                           start4=start4+cutunit
                           cut4AIC=c(cut4AIC,AIC)
          	   	     }
                      start3 = start3+cutunit
       		    start4 = userdata$factor[sum(userdata$factor<=start3)+nmin]
                    }
       		start2 = start2+cutunit
       		start3 = userdata$factor[sum(userdata$factor<=start2)+nmin]
                   }
			start1 = start1+cutunit
			start2 = userdata$factor[sum(userdata$factor<=start1)+nmin]
               }

}#type=sur


#-------------------------------------------------------------------------------------------------------------------#
if (datatype == "logistic"){
	userdata=na.omit(data.frame(outcome,factor))
	index=order(userdata$factor)
	userdata=userdata[index,]

      n=dim(userdata)[1] #number of subject

	range=range(userdata$factor) #range of continous predictor
	cutunit=diff(range)/segment
      logist_red <- glm(userdata$outcome ~userdata$factor, family=binomial())
      originAIC=logist_red$aic

#----------------------------------------cut 1 ---------------------------------------------------------#
cut1AIC=c()
                  start1=userdata$factor[nmin] #confirmed group1 has 20 person and where cut1 to start cutting
			end1=userdata$factor[n-nmin]#cut1 end

                  while(start1<end1)
                    {
			     group1=sum(userdata$factor<=start1) #nummber of group1
                       factor_status=c(rep(0,group1),rep(1,n-group1))
                       logist_fit =glm(userdata$outcome ~factor(factor_status), family=binomial())
                       AIC=logist_fit$aic
          	           start1=start1+cutunit
                       cut1AIC=c(cut1AIC,AIC)
                    }
#----------------------------------------cut 2 ---------------------------------------------------------#
cut2AIC=c()
                  start1=userdata$factor[nmin] #confirmed group1 has 20 person and where cut1 to start cutting
			end1=userdata$factor[n-nmin*2]#cut1 end

			start2=userdata$factor[which(userdata$factor>start1)[nmin]] #confirmed group2 has 20 person and where cut2 to start cutting
			end2=userdata$factor[n-nmin] #cut2 end
  while(start1<end1)
	 {
	      group1=sum(userdata$factor<=start1) #nummber of group1
                 while(start2<end2)
     			    {
       			    group2=sum(userdata$factor<=start2)-group1
                            factor_status=c(rep(0,group1),rep(1,group2),rep(2,n-group1-group2))
               		    logist_fit =glm(userdata$outcome ~factor(factor_status),family=binomial())
                            start2=start2+cutunit
 				    AIC=logist_fit$aic
                            cut2AIC=c(cut2AIC,AIC)
         		     }

       			start1 = start1+cutunit
       			start2 = userdata$factor[sum(userdata$factor<=start1)+nmin]
     				}
#----------------------------------------cut 3 ---------------------------------------------------------#
cut3AIC=c()
                  start1=userdata$factor[nmin] #confirmed group1 has 20 person and where cut1 to start cutting
			end1=userdata$factor[n-nmin*3]#cut1 end

			start2=userdata$factor[which(userdata$factor>start1)[nmin]] #confirmed group2 has 20 person and where cut2 to start cutting
			end2=userdata$factor[n-nmin*2] #cut2 end

			start3=userdata$factor[which(userdata$factor>start2)[nmin]] #confirmed group2 has 20 person and where cut2 to start cutting
			end3=userdata$factor[n-nmin]
	while(start1<end1)
			  {
				group1=sum(userdata$factor<=start1) #nummber of group1

   			while(start2<end2)
     			  {
       			group2=sum(userdata$factor<=start2)-group1

           		while(start3<end3)
            		{
               			group3=sum(userdata$factor<=start3)-group1-group2
               			factor_status=c(rep(0,group1),rep(1,group2),rep(2,group3),rep(3,n-group1-group2-group3))
               			logist_fit =glm(userdata$outcome ~factor(factor_status), family=binomial())
                              AIC=logist_fit$aic
                              cut3AIC=c(cut3AIC,AIC)
          	  			start3=start3+cutunit

         			}

       		start2 = start2+cutunit
       		start3 = userdata$factor[sum(userdata$factor<=start2)+nmin]

     		        }

			start1 = start1+cutunit
			start2 = userdata$factor[sum(userdata$factor<=start1)+nmin]
		        }
#----------------------------------------cut 4 ---------------------------------------------------------#
cut4AIC=c()
                  start1=userdata$factor[nmin] #confirmed group1 has 20 person and where cut1 to start cutting
			end1=userdata$factor[n-nmin*4]#cut1 end

			start2=userdata$factor[which(userdata$factor>start1)[nmin]] #confirmed group2 has 20 person and where cut2 to start cutting
			end2=userdata$factor[n-nmin*3] #cut2 end

			start3=userdata$factor[which(userdata$factor>start2)[nmin]] #confirmed group2 has 20 person and where cut2 to start cutting
			end3=userdata$factor[n-nmin*2]

                  start4=userdata$factor[which(userdata$factor>start3)[nmin]]
			end4=userdata$factor[n-nmin]
	while(start1<end1)
			  {
				group1=sum(userdata$factor<=start1) #nummber of group1

   			while(start2<end2)
     			  {
       			group2=sum(userdata$factor<=start2)-group1

           		while(start3<end3)
            		{
               		  group3=sum(userdata$factor<=start3)-group1-group2
                            while(start4<end4)
            		  {
                              group4=sum(userdata$factor<=start4)-group1-group2-group3
               			factor_status=c(rep(0,group1),rep(1,group2),rep(2,group3),rep(3,group4),rep(4,n-group1-group2-group3-group4))
               			logist_fit =glm(userdata$outcome ~factor(factor_status), family=binomial())
                              AIC=logist_fit$aic
                              cut4AIC=c(cut4AIC,AIC)
          	  			start4=start4+cutunit

         			   }
                          start3 = start3+cutunit
       		       start4 = userdata$factor[sum(userdata$factor<=start3)+nmin]

                        }

       		start2 = start2+cutunit
       		start3 = userdata$factor[sum(userdata$factor<=start2)+nmin]

     		        }

			start1 = start1+cutunit
			start2 = userdata$factor[sum(userdata$factor<=start1)+nmin]
		        }

} #log end
AICmin=which.min(c(min(originAIC),min(cut1AIC),min(cut2AIC),min(cut3AIC),min(cut4AIC)))
cat(paste("AIC have a minimum value ","When cutnum = ",AICmin-1,"\n",sep=""))
return(list(min=c(min(originAIC),min(cut1AIC),min(cut2AIC),min(cut3AIC),min(cut4AIC))))
}