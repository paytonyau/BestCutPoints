# BestSurvCutPoints
### Finding the optimal number and locations of cutpoints for survival analysis

Finding the optimal numbers and locations of cutpoints for survival analysis is a statistical technique used in medical research and other related fields to categorize a continuous variable into groups or intervals. The aim is to find the optimal number and location of cutpoints that provide the best discrimination between the risk of an event (such as death or disease) and non-occurrence of an event. The optimal number and location of cutpoints can be determined using various statistical methods, such as likelihood ratio tests, log-rank tests, or areas under the receiver operating characteristic (ROC) curve. To determine the optimal number of cutpoints, researchers can use methods such as the Akaike information criterion (AIC) or cross-validation. To determine the optimal location of cutpoints, researchers can use genetic algorithms or other optimization techniques. The optimal number and location of cutpoints can help researchers better understand the relationship between the continuous variable and the event of interest, leading to improved decision-making and clinical outcomes.

----------------------------
(A1) `findcutnum` - find **the optimal number of cut-off points** for a continuous risk factor. It minimizes the akaike information criterion (AIC) and handle both survival and binary outcome

**findcutnum (factor,outcome,datatype,nmin=20,segment=100)**

**factor**: continuous risk factor for which to find the optimal cut-offs (Nx1 vector)
**outcome**: a matrix of event and time
**datatype**: specify "survival" for survival outcome or "binary" for binary outcome (string)
**nmin**: minimum number of individuals in each group(positive integer, default=20)
**segment**: total number of pieces(integer, default=100)

(A2) `findcut`- determines **the best location of cut-off points** for a continuous risk factor use contingency tables (X^2^) approach, and it can handle both survival and binary data. For **survival**, the recommended criteria to use are "likelihood ratio test" and "logrank test", while for **binary**, "AUC" and "Likelihood ratio test" are suggested

**findcut(factor,outcome,cutnum,datatype,nmin,segment)**

**factor**: continuous risk factor for which to find the optimal cut-offs (Nx1 vector)
**cutnum**: number of cut-offs (positive integer)
**datatype**: specify "survival" for survival outcome or "binary" for binary outcome (string)
**nmin**: minimum number of individuals in each group(positive integer, default=20)
**segment**: total number of pieces(integer, default=100)

---------------------------------------

#### A demonstration of using `findcutnum` and  `findcut` in the R
```
library("survival")
library("KMsurv")
library("xtable")
library("splines")
library("pROC")
library("aod")

# Survival data
BMI=c(30,16,29,29,21,29,27,24,17,27,22,27,26,16,21,21,23,20,25,23,28,20,22,22,37,23,25,34,31,26)
event=c(1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,1,1,0)
OS=c(138,92,64,15,62,235,214,197,41,33,257,115,123,44,154,71,61,182,75,214,25,217,113,200,175,117,166,0,57,186)
invasion=c(1.1,0.1,1.0,0.8,1.2,1.0,0.3,1.0,0.4,0.6,0.4,0.8,1.0,1.0,1.1,0.5,0.9,1.0,0.6,0.1,1.1,0.4,0.4,1.1,1.1,0.4,1.1,1.2,0.9,0.8)
LVSI=c(0,0,0,1,0,1,0,0,1,1,0,1,1,1,0,0,0,1,1,1,1,0,1,0,0,0,0,1,0,1)
```

###### (1) Analyzing data from survival information 
```
findcutnum(factor=BMI, outcome=cbind(event,OS), datatype="survival", nmin=5, segment=100)

findcut(factor=BMI, outcome=cbind (event,OS), cutnum=2, datatype = "survival", nmin=5, segment=100)
```
###### (2) Analyzing data from dichotomised data
```
findcutnum(factor=invasion, outcome=LVSI, datatype="logistic", nmin=5, segment=100)

findcut(factor=invasion, outcome=LVSI, cutnum=2, datatype = "logistic", nmin=5, segment=100)
```

-------------------------------------
(B1) `findcutnumCox` - find the optimal **number** of cutpoints

**findnumCox(target,event,time,confound,totalcut=3,initial_rr=NULL,initial_cut=NULL,initial_domain=NULL,numgen=10,numcross=20,gap=NULL)**

**target**: a continuous target variable to be categorized (an nx1 vector)
**event**: failure indicator (1: event occurs; 0: right censored)
**time**: observed time
**confound**: an nxq data.frame including all the confounding covariates
**totalcut**: maximum number of cutpoints (default is 3)
**initial_rr**: initial values for relative risk (default is NULL)
**initial_cut**: initial values for the locations of cutpoints (default is NULL)
**initial_domain**: upper and lower bounds for cut points (default is NULL)
**numgen**: maximum number of iterations for genetic algorithms
**numcross**: number of cross-validations (default is NULL)
**gap**: minimum gap between two consecutive cutpoints (default is 0.03)

(B2) `findcutCox` -  find the optimal **location** of cutpoints for a continuous variable

**findcutCox(target,event,time,numbercross,numcut,initial_rr=NULL,initial_cut=NULL,initial_domain=NULL,numgen,gap=0.03)**

**target**: A continuous variable to be categorized (an nx1 vector)
**event**: Failure indicator (1: event occurs; 0: right censored)
**time**: Observed time
**confound**: An nxq data.frame including all confounding covariates
**numcross**: Number of cross validation (ie B in the paper)
**numcut**: Number of cutpoints
**initial_rr**: Initial values for relative risk; Type: list; Default is NULL
**initial_cut**: Initial values for the locations of cutpoints; Type: list; Default is NULL
**initial_domain**: Upper and lower bounds for cut points; Type: a kx2 matrix; Default is NULL; each row of the matrix (a 1x2 vector) represents the lower and upper bound for one cut point.
**numgen**: Maximum number of iterations for genetic algorithms
**gap**: Minimum gap between two consecutive cutpoints, default is 0.03

---------------------------------------
#### An example run using  `findnumCox` and  `findcutCox`
```
library("rgenoud")
library("survival")
library("foreach")
library("doParallel")
library("doRNG")
library("xtable")

data=read.csv("toydata.csv") # in the Tuturial folder
attach(data)
set.seed(30)
ptm <- proc.time()

result<-findnumCox(BMI,Death,Death_surtime,
                   confound = stage3,numcross=20,
                   totalcut=3,initial_rr=NULL,
                   initial_cut=NULL,initial_domain=NULL,
                   numgen=10,gap=NULL)

proc.time()-ptm
result$aic #corrected AIC value
which.min(result$aic) ## number of optimal cut points
result$HR ## corrected hazard ratio (or corrected relative risk)
#Cutpvalue # corrected p-value for each coefficient estimator (Note that for each cross validation, the p-value would be different)
```

###### confound example
```
set.seed(2019)
target=data$BMI;event=data$Death;time=data$Death_surtime;confound=data$stage3
userdata<-na.omit(data.frame(target,event,time,confound))
N<-dim(userdata)[1]
numboot<-sample(N)
initial_rr <- list(c(3,2),c(3,6,2),c(3,4,5,2))
initial_cut <- list(c(19),c(19,30),c(15,25,30))
initial_domain <- list(matrix(c(12,35,0,5,0,5),
                       ncol = 2,byrow = TRUE),
                       matrix(c(15,35,15,35,0,5,0,5,0,5),
                       ncol = 2,byrow = TRUE),
                       matrix(c(15,35,15,35,15,35,0,5,0,5,0,5,0,5),
                       ncol = 2,byrow = TRUE))

## example for no initial values
E <- findcutCox(BMI,Death,Death_surtime,stage3,numcut=3,initial_rr=NULL,
                initial_cut=NULL,initial_domain=NULL,numgen=15,gap=NULL)
E

## example for giving initial values
G <- findcutCox(BMI, Death, Death_surtime, stage3, numcut=3, initial_rr=initial_rr,
                initial_cut=initial_cut, initial_domain=initial_domain, numgen=15, gap=NULL)
G
```

###### non confound example
```
set.seed(2019)
target=data$BMI;event=data$Death;time=data$Death_surtime
userdata<-na.omit(data.frame(target,event,time))
N<-dim(userdata)[1]
numboot<-sample(N)
initial_rr <- list(c(3),c(6,3),c(3,4,5))
initial_cut <- list(c(19),c(19,30),c(15,25,30))
initial_domain <- list(matrix(c(13,48,0,5),ncol = 2,byrow = TRUE),
                       matrix(c(13,35,13,35,0,5,0,5),ncol = 2,byrow = TRUE),
                       matrix(c(13,35,13,35,13,35,0,5,0,5,0,5),ncol = 2,byrow = TRUE))

## example for no initial values (numcut=2，repercenting finding the best 2 cutoffs)
H <- findcutCox(BMI,Death,Death_surtime,confound=NULL,numcut=2,
                initial_rr=NULL,initial_cut=NULL,
                initial_domain=NULL,numgen=15,gap=NULL)
H

## example for giving initial values
I <- findcutCox(BMI,Death, Death_surtime, confound=NULL, numcut=2, initial_rr=initial_rr,
                initial_cut=initial_cut,initial_domain=initial_domain,numgen=15,gap=NULL)
I
```

The [original publication](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0176231), along with its associated scripts (https://osf.io/ef7na/ and http://www.math.nsysu.edu.tw/~cchang/) , were obtained through the Attribution 4.0 International (CC BY 4.0) license. The authors  have provided a step-by-step [user manual (in pdf)](https://github.com/paytonyau/BestSurvCutPoints/blob/main/Tuturial/User_Manual_Fundcut.pdf) to assist users in using the script for `findcutnum` and  `findcut`.  The example to run for `findnumCox` and  `findcutCox` can be found in [the tutorial](https://github.com/paytonyau/BestSurvCutPoints/blob/main/Tuturial/).

##### REFERENCE
Chang, C., Hsieh, M. K., Chang, W. Y., Chiang, A. J., & Chen, J. (2017). Determining the optimal number and location of cutoff points with application to data of cervical cancer. PloS one, 12(4), e0176231.
