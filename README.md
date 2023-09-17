
# BestSurvCutPoints

### Finding the optimal number and locations of cutpoints for survival analysis

BestSurvCutPoints is an R package that provides functions for finding optimal cutpoints in survival analysis. It includes four core functions: `findcutnum`, `findcut`, `findnumCox`, and `findcutCox`. It involves categorizing a continuous variable into intervals to best distinguish event occurrence (e.g., death or disease) from non-occurrence. This package utilises AIC or cross-validation to find the ideal number of cutpoints, and genetic algorithms for their precise locations and help researchers and data analysts determine the most suitable cutpoints for survival data analysis.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Function Documentation](#function-documentation)
- [Examples](#examples)
- [References](#References)
- [License](#license)

## Installation

you can install it from GitHub using the `devtools` package:

`devtools::install_github("paytonyau/BestSurvCutPoints")` 

## Usage

To use the BestSurvCutPoints package, load it using the `library` function:

`library(BestSurvCutPoints)` 

Once loaded, you can use the following core functions:

## Function Documentation

### `findcutnum`

The `findcutnum` function helps find **optimal number of cutpoints** for a continuous risk factor. It minimizes the akaike information criterion (AIC) and handles both survival and binary outcome

***findcutnum (factor, outcome, datatype, nmin = 20, segment = 100)***

**Arguments:**

- `factor`: continuous risk factor for which to find the optimal cut-offs (Nx1 vector)
- `outcome`: a matrix of event and time
- `datatype`: specify "survival" for survival outcome or "binary" for binary outcome (string)
- `nmin`: minimum number of individuals in each group(positive integer, default=20)
- `segment`: total number of pieces (integer, default=100)

**Returns:**

A list containing information about the optimal cutpoints and the corresponding test statistics.

### `findcut`

The `findcut` function helps find **the optimal cutpoints** for a continuous risk factor using contingency tables (X^2^). it can handle both survival and binary data. 
For **survival**, the recommended criteria to use are "likelihood ratio test" and "logrank test", while for **binary**, "AUC" and "Likelihood ratio test" are suggested.

***findcut(factor, outcome, cutnum, datatype, nmin, segment)***

**Arguments:**

- `factor`: continuous risk factor for which to find the optimal cut-offs (Nx1 vector)
- `cutnum`: number of cut-offs (positive integer)
- `datatype`: specify "survival" for survival outcome or "binary" for binary outcome (string)
- `nmin`: minimum number of individuals in each group(positive integer, default=20)
- `segment`: total number of pieces(integer, default=100)

**Returns:**

A numeric value representing the optimal cutpoint.

### `findnumCox`

The `findnumCox` function finds the optimal number of cutpoints for survival analysis.

***findnumCox(target, event, time, confound, totalcut = 3, initial_rr = NULL, initial_cut = NULL, initial_domain = NULL, numgen = 10, numcross = 20,gap = NULL)***

**Arguments:**

- `target`: a continuous target variable to be categorized (an nx1 vector)
- `event`: failure indicator (1: event occurs; 0: right censored)
- `time`: observed time
- `confound`: an nxq data.frame including all the confounding covariates
- `totalcut`: maximum number of cutpoints (default is 3)
- `initial_rr`: initial values for relative risk (default is NULL)
- `initial_cut`: initial values for the locations of cutpoints (default is NULL)
- `initial_domain`: upper and lower bounds for cut points (default is NULL)
- `numgen`: maximum number of iterations for genetic algorithms
- `numcross`: number of cross-validations (default is NULL)
- `gap`: minimum gap between two consecutive cutpoints (default is 0.03)


**Returns:**

A list containing AIC values, hazard ratios, and cutpoint p-values for different models with varying numbers of cutpoints.

### `findcutCox`


The `findcutCox` function finds the optimal **location** cutpoints for survival analysis (a continuous variable).

***findcutCox(target,event,time,numbercross,numcut,initial_rr=NULL,initial_cut=NULL,initial_domain=NULL,numgen,gap=0.03)***

**Arguments:**

- `target`: A continuous variable to be categorized (an nx1 vector)
- `event`: Failure indicator (1: event occurs; 0: right censored)
- `time`: Observed time
- `confound`: An nxq data.frame including all confounding covariates
- `numcross`: Number of cross validation (ie B in the paper)
- `numcut`: Number of cutpoints
- `initial_rr`: Initial values for relative risk; Type: list; Default is NULL
- `initial_cut`: Initial values for the locations of cutpoints; Type: list; Default is NULL
- `initial_domain`: Upper and lower bounds for cut points; Type: a kx2 matrix; Default is NULL; each row of the matrix (a 1x2 vector) represents the lower and upper bound for one cut point.
- `numgen`: Maximum number of iterations for genetic algorithms
- `gap`: Minimum gap between two consecutive cutpoints, default is 0.03

**Returns:**

A numeric value representing the optimal cutpoint.

## Examples

Here are some examples of how to use the package:
#### A demonstration of using `findcutnum` and  `findcut` in the R
```
# Load necessary libraries
library("survival")   # Survival analysis library
library("KMsurv")     # Kaplan-Meier survival curves
library("xtable")     # Table generation for reports
library("splines")    # Basis splines for modeling
library("pROC")       # Receiver Operating Characteristic (ROC) analysis
library("aod")        # Analysis of Overdispersed Data

# Survival data
BMI = c(30,16,29,29,21,29,27,24,17,27,22,27,26,16,21,21,23,20,25,23,28,20,22,22,37,23,25,34,31,26)
event = c(1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,1,1,0)
OS = c(138,92,64,15,62,235,214,197,41,33,257,115,123,44,154,71,61,182,75,214,25,217,113,200,175,117,166,0,57,186)
invasion = c(1.1,0.1,1.0,0.8,1.2,1.0,0.3,1.0,0.4,0.6,0.4,0.8,1.0,1.0,1.1,0.5,0.9,1.0,0.6,0.1,1.1,0.4,0.4,1.1,1.1,0.4,1.1,1.2,0.9,0.8)
LVSI = c(0,0,0,1,0,1,0,0,1,1,0,1,1,1,0,0,0,1,1,1,1,0,1,0,0,0,0,1,0,1)

```

##### (1) Analysing data from survival information 
```
# Function to find optimal cutpoints based on BMI for survival analysis
findcutnum(
  factor = BMI,               # Predictor variable
  outcome = cbind(event, OS), # Survival outcome data
  datatype = "survival",      # Data type for analysis
  nmin = 5,                   # Minimum group size
  segment = 100               # Number of segments
)

# Function to find optimal cutpoints based on BMI for survival analysis with 2 cutpoints
findcut(
  factor = BMI,               # Predictor variable
  outcome = cbind(event, OS), # Survival outcome data
  cutnum = 2,                 # Number of cutpoints
  datatype = "survival",      # Data type for analysis
  nmin = 5,                   # Minimum group size
  segment = 100               # Number of segments
)

```
##### (2) Analysing data from dichotomised data
```
# Function to find optimal cutpoints based on invasion for logistic regression analysis
findcutnum(
  factor = invasion,   # Predictor variable (invasion)
  outcome = LVSI,      # Outcome variable (LVSI)
  datatype = "logistic",  # Data type for logistic regression
  nmin = 5,            # Minimum group size
  segment = 100        # Number of segments
)

# Function to find optimal cutpoints based on invasion for logistic regression analysis with 2 cutpoints
findcut(
  factor = invasion,   # Predictor variable (invasion)
  outcome = LVSI,      # Outcome variable (LVSI)
  cutnum = 2,          # Number of cutpoints
  datatype = "logistic",  # Data type for logistic regression
  nmin = 5,            # Minimum group size
  segment = 100        # Number of segments
)

```
-------------------------------------
#### An example run using  `findnumCox` and  `findcutCox`
```
# Load necessary libraries
library("rgenoud")
library("survival")
library("foreach")
library("doParallel")
library("doRNG")
library("xtable")

# Load your dataset from a CSV file (toydata.csv should be in the current working directory)
data <- read.csv("toydata.csv")

# Attach the data for easy access to columns
attach(data)

# Set a random seed for reproducibility
set.seed(30)

# Record the starting time
ptm <- proc.time()

# Call the findnumCox function to find optimal cutpoints
result <- findnumCox(
  BMI,                # Predictor variable (BMI)
  Death,              # Event indicator variable (Death)
  Death_surtime,      # Survival time variable (Death_surtime)
  confound = stage3,  # Confounding variable (stage3)
  numcross = 20,      # Number of cross-validations
  totalcut = 3,       # Total number of cutpoints to consider
  initial_rr = NULL,  # Initial relative risk values (set to NULL)
  initial_cut = NULL, # Initial cutpoint values (set to NULL)
  initial_domain = NULL, # Initial domain values (set to NULL)
  numgen = 10,        # Number of generations for genetic algorithm
  gap = NULL          # Gap value (set to NULL)
)

# Calculate the time taken for the analysis
proc.time() - ptm

# Print the corrected AIC values
result$aic

# Find the number of optimal cutpoints (the index with the minimum AIC)
which.min(result$aic)

# Print the corrected hazard ratios (or relative risks)
result$HR

# Note: Corrected p-values for each coefficient estimator (Cutpvalue) are available in the result but are not printed here.
```

##### confound example
```
# Set a random seed for reproducibility
set.seed(2019)

# Extract data from your dataset
target <- data$BMI
event <- data$Death
time <- data$Death_surtime
confound <- data$stage3

# Create user data based on input and handle missing values
userdata <- na.omit(data.frame(target, event, time, confound))

# Determine the number of observations
N <- dim(userdata)[1]

# Generate random bootstrap samples
numboot <- sample(N)

# Define initial values for relative risks, cutpoints, and domains
initial_rr <- list(c(3, 2), c(3, 6, 2), c(3, 4, 5, 2))
initial_cut <- list(c(19), c(19, 30), c(15, 25, 30))
initial_domain <- list(
  matrix(c(12, 35, 0, 5, 0, 5), ncol = 2, byrow = TRUE),
  matrix(c(15, 35, 15, 35, 0, 5, 0, 5, 0, 5), ncol = 2, byrow = TRUE),
  matrix(c(15, 35, 15, 35, 15, 35, 0, 5, 0, 5, 0, 5, 0, 5), ncol = 2, byrow = TRUE)
)

# Example for no initial values provided
E <- findcutCox(BMI, Death, Death_surtime, stage3, numcut = 3, initial_rr = NULL,
                initial_cut = NULL, initial_domain = NULL, numgen = 15, gap = NULL)

# Print the result
E

# Example for providing initial values
G <- findcutCox(BMI, Death, Death_surtime, stage3, numcut = 3, initial_rr = initial_rr,
                initial_cut = initial_cut, initial_domain = initial_domain, numgen = 15, gap = NULL)

# Print the result
G
```

##### non confound example
```
# Set a random seed for reproducibility
set.seed(2019)

# Extract data from your dataset
target <- data$BMI
event <- data$Death
time <- data$Death_surtime

# Create user data based on input and handle missing values
userdata <- na.omit(data.frame(target, event, time))

# Determine the number of observations
N <- dim(userdata)[1]

# Generate random bootstrap samples
numboot <- sample(N)

# Define initial values for relative risks, cutpoints, and domains
initial_rr <- list(c(3), c(6, 3), c(3, 4, 5))
initial_cut <- list(c(19), c(19, 30), c(15, 25, 30))
initial_domain <- list(
  matrix(c(13, 48, 0, 5), ncol = 2, byrow = TRUE),
  matrix(c(13, 35, 13, 35, 0, 5, 0, 5), ncol = 2, byrow = TRUE),
  matrix(c(13, 35, 13, 35, 13, 35, 0, 5, 0, 5, 0, 5), ncol = 2, byrow = TRUE)
)

# Example for no initial values (finding the best 2 cutoffs)
H <- findcutCox(BMI, Death, Death_surtime, confound = NULL, numcut = 2,
                initial_rr = NULL, initial_cut = NULL,
                initial_domain = NULL, numgen = 15, gap = NULL)

# Print the result
H

# Example for providing initial values
I <- findcutCox(BMI, Death, Death_surtime, confound = NULL, numcut = 2, initial_rr = initial_rr,
                initial_cut = initial_cut, initial_domain = initial_domain, numgen = 15, gap = NULL)

# Print the result
I
```


The [original publication](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0176231), along with its associated scripts (https://osf.io/ef7na/ and http://www.math.nsysu.edu.tw/~cchang/) , were obtained through the Attribution 4.0 International (CC BY 4.0) license. The authors  have provided a step-by-step [user manual (in pdf)](https://github.com/paytonyau/BestSurvCutPoints/blob/main/Tuturial/User_Manual_Fundcut.pdf) to assist users in using the script for `findcutnum` and  `findcut`.  
The example to run for `findnumCox` and  `findcutCox` can be found in [the tutorial](https://github.com/paytonyau/BestSurvCutPoints/blob/main/Tuturial/).

## References
Chang, C., Hsieh, M. K., Chang, W. Y., Chiang, A. J., & Chen, J. (2017). Determining the optimal number and location of cutoff points with application to data of cervical cancer. PloS one, 12(4), e0176231.