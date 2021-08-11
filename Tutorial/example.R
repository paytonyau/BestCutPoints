##Before running this example, first save the files findnumCox.R, findcutCox_gen.R, example.R and toydata.csv in the same directory.
## source files
gc();rm()
setwd('D:/cutpoints') ## please change the directory in which you save the files findnumCox.R, findcutCox_gen.R, example.R and toydata.csv. 
source("findnumCox.R")
source("findcutCox.R")


################################################################################
#### tutorial for using the function "findnumCox"########################################################
################################################################################
#### Aim of function findnumCox: find optimal number of cutpoints
#### Before using findnumCox, first source functions findnumCox.R and "findcutCox.R
#### set up the input arguments
#### target:a target continous variable to be categorized (an nx1 vector)
#### event: failure indicator (1: event occurs; 0: right censored)
#### time: observed time
#### confound: an nxq data.frame including all the confounding covariates
#### totalcut: maximum number of cutpoints(i.e. K in the paper); default=3.
#### initial_rr: initial values for relative risk; type: list; defalut is NULL
#### initial_cut: initial values for the locations of cutpoints; type: list; default is NULL
#### initial_domain: upper and lower bounds for cut points; type: a kx2 matrix; default is NULL; each row of the matrix (a 1x2 vector) 
#### is the lower and upper bound for one cut point
#### numgen: maximum number of iterations for genetic algorithms
#### numcross: number of cross validation (ie B in the paper)
#### gap: minimum gap between two consecutive cutpoints, default is 0.03
################################################################################

################################################################################
#### confound example (without setting the intial values) ###############################
################################################################################
data=read.csv("toydata.csv")
attach(data)
set.seed(30)
ptm <- proc.time()
result<-findnumCox(BMI,                 
           Death,               
           Death_surtime,       
           confound = stage3,   
           numcross=20,          
           totalcut=3,           
           initial_rr=NULL,   
           initial_cut=NULL,    
           initial_domain=NULL, 
           numgen=10,            
           gap=NULL             
)
proc.time()-ptm
result$aic #corrected AIC value
which.min(result$aic) ## number of optimal cut points
result$HR ## corrected hazard ratio (or corrected relative risk)
#Cutpvalue # corrected p-value for each coefficient estimator (Note that for each cross validation, the p-value would be different)
################################################################################
#### confound example (with intial values)####################################
################################################################################
set.seed(3)
result<-findnumCox(BMI,
           Death,
           Death_surtime,
           confound = stage3,
           numcross=20,
           totalcut=3,
           initial_rr=list(c(3,2),c(3,6,2),c(3,4,5,2)),
           initial_cut=list(c(19),c(19,30),c(15,25,30)),
           initial_domain=list(matrix(c(15,35,0,5,0,5),ncol = 2,byrow = TRUE),
                               matrix(c(15,35,15,35,0,5,0,5,0,5),ncol = 2,byrow = TRUE),
                               matrix(c(15,35,15,35,15,35,0,5,0,5,0,5,0,5),ncol = 2,byrow = TRUE)),
           numgen=10,
           gap=0.03)
result$aic # corrected AIC value
which.min(result$aic) # choose optimal number of cut points
result$HR # corrected hazard ratio
#Cutpvalue



################################################################################
#### tutorial for using function findcutCox########################################################
################################################################################
#### Aim of function findcutCox: to find optimal locations of cutpoints
#### Before using findcutCox, first source functions findnumCox.R and "findcutCox.R
#### set up the input arguments: 
#### target:a target continous variable to be categorized (an nx1 vector)
#### event: failure indicator (1: event occurs; 0: right censored)
#### time: observed time
#### confound: an nxq data.frame including all the confounding covariates
#### numcross: number of cross validation (ie B in the paper)
#### numcut: number of cutpoints
#### initial_rr: initial values for relative risk; type: list; defalut is NULL
#### initial_cut: initial values for the locations of cutpoints; type: list; default is NULL
#### initial_domain: upper and lower bounds for cut points; type: a kx2 matrix; default is NULL; each row of the matrix (a 1x2 vector) 
#### is the lower and upper bound for one cut point
#### numgen: maximum number of iterations for genetic algorithms
#### gap: minimum gap between two consecutive cutpoints, default is 0.03
################################################################################

############################################################################
######## confound example ##################################################
############################################################################
set.seed(2019)
target=data$BMI;event=data$Death;time=data$Death_surtime;confound=data$stage3
userdata<-na.omit(data.frame(target,event,time,confound))
N<-dim(userdata)[1]
numboot<-sample(N)
initial_rr <- list(c(3,2),c(3,6,2),c(3,4,5,2))
initial_cut <- list(c(19),c(19,30),c(15,25,30))
initial_domain <- list(matrix(c(12,35,0,5,0,5),ncol = 2,byrow = TRUE),
                       matrix(c(15,35,15,35,0,5,0,5,0,5),ncol = 2,byrow = TRUE),
                       matrix(c(15,35,15,35,15,35,0,5,0,5,0,5,0,5),ncol = 2,byrow = TRUE))
## example for no initial values
E <- findcutCox(BMI,Death,Death_surtime,stage3,numcut=3,initial_rr=NULL,initial_cut=NULL,initial_domain=NULL,numgen=15,gap=NULL)
E
## example for giving initial values
G <- findcutCox(BMI,Death,Death_surtime,stage3,numcut=3,initial_rr=initial_rr,initial_cut=initial_cut,initial_domain=initial_domain,numgen=15,gap=NULL)
G
############################################################################
######## no_confound example ###############################################
############################################################################
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
## example for no initial values (numcut=2ï¼Œè¡¨ç¤ºåªå°‹æ‰¾?…©?€‹å?‡å?†é??)
H <- findcutCox(BMI,Death,Death_surtime,confound=NULL,numcut=2,initial_rr=NULL,initial_cut=NULL,initial_domain=NULL,numgen=15,gap=NULL)
H
## example for giving initial values
I <- findcutCox(BMI,Death,Death_surtime,confound=NULL,numcut=2,initial_rr=initial_rr,initial_cut=initial_cut,initial_domain=initial_domain,numgen=15,gap=NULL)
I
#### excute and save the results##
# K=NULL
# for(i in 1:10){
#   numboot<-sample(N)
#   H <- findcutCox(BMI,Death,Death_surtime,confound=NULL,numcut=3,initial_rr=NULL,initial_cut=NULL,initial_domain=NULL,numgen=5,gap=NULL)
#   K <- rbind(H$cut[1:3],K)
# }
# K
# colMeans(K)


