#' findnumCox
#'
#' find the optimal number of cutpoints
#' @param findnumCox
#' @return
#' @examples
#' findnumCox(BMI,Death,Death_surtime,confound = stage3,numcross=20,totalcut=3,initial_rr=NULL,initial_cut=NULL,initial_domain=NULL,numgen=10,gap=NULL)
#' @export

# Load necessary libraries
findnumCox <- function(target, event, time, confound, numcross, totalcut = 3, initial_rr = NULL,
                       initial_cut = NULL, initial_domain = NULL, numgen, gap = 0.03) {
  libraries1 <- c("rgenoud", "survival", "foreach", "doParallel", "doRNG", "xtable")
  lapply(libraries1, library, quietly = TRUE, character.only = TRUE)

  # Set global variables
  confound <<- confound

  # Create user data based on input and handle confound variable
  if (is.null(confound)) {
    userdata <- na.omit(data.frame(target, event, time))
  } else {
    userdata <- na.omit(data.frame(target, event, time, confound))
    confound <- data.frame(userdata[, 4:dim(userdata)[2]])
  }

  # Attach user data for easy access
  attach(userdata)
  N <<- dim(userdata)[1]
  numboot <<- sample(N)

  # Perform cross-validation using foreach
  aictest.apply <- foreach(i = 1:1, .combine = cbind) %do% {
    attach(userdata)
    target <- userdata$target
    event <- userdata$event
    time <- userdata$time

    if (is.null(confound)) {
      counfound <- NULL
    } else {
      confound <- userdata$confound
    }

    mt <- replicate(numcross, sample(N))
    result <- cbind(apply(mt, 2, aictest, target = target, event = event, time = time,
                          confound = confound, totalcut = totalcut, initial_rr = initial_rr,
                          initial_cut = initial_cut, initial_domain = initial_domain,
                          numgen = numgen, gap = gap))
  }
  aictest.apply <<- aictest.apply

  # Initialize result variables
  lik_cut1 <- NULL
  hr1 <- NULL
  lik_cut2 <- NULL
  hr2 <- NULL
  lik_cut3 <- NULL
  hr3 <- NULL
  pvalue1 <- NULL
  pvalue2 <- NULL
  pvalue3 <- NULL
  cutpvalue1 <- NULL
  cutpvalue2 <- NULL
  cutpvalue3 <- NULL

  # Calculate AIC, HR, p-values, and cut p-values
  RESULT <- for (i in 1:numcross) {
    A <- aictest.apply[[i]]
    lik_cut1[i] <- aictest.apply[[i]][1, 1]
    lik_cut2[i] <- aictest.apply[[i]][1, 2]
    lik_cut3[i] <- aictest.apply[[i]][1, 3]
    hr1 <- rbind(hr1, A[, 1]$beta)
    hr2 <- rbind(hr2, A[, 2]$beta)
    hr3 <- rbind(hr3, A[, 3]$beta)
    pvalue1 <- rbind(pvalue1, A[, 1]$overallstatistic.cross)
    pvalue2 <- rbind(pvalue2, A[, 2]$overallstatistic.cross)
    pvalue3 <- rbind(pvalue3, A[, 3]$overallstatistic.cross)
    cutpvalue1 <- rbind(cutpvalue1, A[, 1]$cutpvalue)
    cutpvalue2 <- rbind(cutpvalue2, A[, 2]$cutpvalue)
    cutpvalue3 <- rbind(cutpvalue3, A[, 3]$cutpvalue)
  }

  # Calculate AIC
  aic <- list(aic1 = mean((-2) * unlist(lik_cut1) + dim(userdata)[2] - 2),
              aic2 = mean((-2) * unlist(lik_cut2) + dim(userdata)[2] - 1),
              aic3 = mean((-2) * unlist(lik_cut3) + dim(userdata)[2]))

  # Calculate overall p-values
  overall_pvalue <- list(overallpvalue1 = sort(1 - pchisq(2 * pvalue1, 1)),
                         overallpvalue2 = sort(1 - pchisq(2 * pvalue2, 2)),
                         overallpvalue3 = sort(1 - pchisq(2 * pvalue3, 3)))

  # Calculate HR
  hr <- list(hr1 = exp(median(hr1[, 1])),
             hr2 = c(exp(median(hr2[, 1])), exp(median(hr2[, 2]))),
             hr3 = c(exp(median(hr3[, 1])), exp(median(hr3[, 2])), exp(median(hr3[, 3]))))

  # Store cut p-values
  cutpvalue <- list(cutpvalue1 = cutpvalue1, cutpvalue2 = cutpvalue2, cutpvalue3 = cutpvalue3)

  # Assign results to global variables
  aic <<- aic
  HR <<- hr
  Overall_pvalue <<- overall_pvalue
  Cutpvalue <<- cutpvalue

  return(list(aic = aic, HR = HR, Cutpvalue = cutpvalue))
}

# Function to perform AIC test
aictest <- function(target, event, time, confound = NULL, totalcut = totalcut, initial_rr = NULL,
                    initial_cut = NULL, initial_domain = NULL, numgen = numgen, gap = NULL, numboot) {
  confound <<- confound

  # Create user data based on input and handle confound variable
  if (is.null(confound)) {
    userdata <- na.omit(data.frame(target, event, time))
  } else {
    userdata <- na.omit(data.frame(target, event, time, confound))
    confound <- userdata[, 4:dim(userdata)[2]]
  }
  target <- userdata$target
  censor <- userdata$event
  Ti <- userdata$time
  N <- dim(userdata)[1]

  # Define aicfunction for cross-validation
  aicfunction <- function(numcut, numgen, numboot) {
    set1 <- numboot[1:round(N / 2)]
    set2 <- numboot[(round(N / 2) + 1):N] # Cross-validation

    numcut <<- numcut

    generate1 <- maxloglik(target[set1], numcut, Ti[set1], censor[set1], confound[set1], numgen,
                           domain_range(target[set1], initial_domain, numcut, confound[set1]),
                           initial(target[set1], initial_rr, initial_cut, numcut), gap)

    generate2 <- maxloglik(target[set2], numcut, Ti[set2], censor[set2], confound[set2], numgen,
                           domain_range(target[set2], initial_domain, numcut, confound[set2]),
                           initial(target[set2], initial_rr, initial_cut, numcut), gap)

    cut1 <- generate1$par_corr
    cut2 <- generate2$par_corr
    group1 <- cut(target[set1], breaks = c(0, cut2[1:numcut], max(target[set1])),
                  labels = c(0:numcut), quantile = FALSE)
    group2 <- cut(target[set2], breaks = c(0, cut1[1:numcut], max(target[set2])),
                  labels = c(0:numcut), quantile = FALSE)

    if (is.null(confound)) {
      fit <- coxph(Surv(time[c(set1, set2)], censor[c(set1, set2)]) ~
                     factor(c(group1, group2)) + strata(rep(c(1, 2), c(round(N / 2), N - round(N / 2)))))
      overallstatistic.cross <- as.numeric(summary(fit)[["logtest"]]["test"])
    } else {
      fit <- coxph(Surv(time[c(set1, set2)], censor[c(set1, set2)]) ~
                     factor(c(group1, group2)) + strata(rep(c(1, 2), c(round(N / 2), N - round(N / 2)))) +
                     as.matrix(confound[c(set1, set2)]))
      fit.re <- coxph(Surv(time[c(set1, set2)], censor[c(set1, set2)]) ~
                        strata(rep(c(1, 2), c(round(N / 2), N - round(N / 2)))) +
                        as.matrix(confound[c(set1, set2)]))
      overallstatistic.cross <- fit$loglik[2] - fit.re$loglik[2]
    }

    # Wald Test
    varcov = fit$var
    cutpvalue <- as.numeric(summary(fit)$coefficients[, 5])
    beta <- as.numeric(fit$coefficients)

    return(list(aic = (-2 * fit$loglik[2] + 2 * numcut), loglik = fit$loglik[2], cutpvalue = cutpvalue,
                beta = beta, cut1 = cut1, cut2 = cut2, overallstatistic.cross = overallstatistic.cross,
                chiWald = varcov))
  }

  # Function to handle exceptions during AIC calculation
  catchwrongaic <- function(numcut, baby, rand) {
    out <- tryCatch(aicfunction(numcut, baby, rand),
                    error = function(e) mget(c("aic", "loglik", "pvalue", "beta", "cut1", "cut2", "overall", "wald"),
                                             new.env(),
                                             ifnotfound = as.list(list(c("wrong"), c("wrong"), rep(c("wrong"), numcut),
                                                                       rep(c("wrong"), numcut),
                                                                       rep(c("wrong"), 2 * numcut),
                                                                       c("wrong"), c("wrong"))))
    )
    return(out)
  }

  # Apply aicfunction to totalcut number of times
  return(sapply(c(1:totalcut), aicfunction, numgen = numgen, numboot = numboot))
}

# Function to perform maximum log-likelihood estimation
maxloglik <- function(target, numcut, time, censor, confound = NULL, baby, domain, initial, gap = NULL) {
  # Set global variables
  time_global <<- time
  censor_global <<- censor
  target_global <<- target
  nvars <<- 2 * numcut + dim(data.frame(confound))[2]

  if (is.null(confound)) {
    confound_global <<- NULL
  } else {
    confound_global <<- as.matrix(confound)
  }

  numcut_global <<- numcut
  target_max <<- max(target)

  if (is.null(gap)) {
    gap_global <<- quantile(sort(diff(sort(na.omit(target)))), probs = 0.5)
  } else {
    gap_global <<- gap
  }

  # Perform genetic algorithm optimization
  ccc <- genoud(obj, nvars, max = TRUE, pop.size = 100, max.generations = baby, wait.generations = 10,
                hard.generation.limit = TRUE, starting.values = initial, MemoryMatrix = TRUE,
                Domains = domain, solution.tolerance = 0.001, print.level = 0,
                gr = NULL, boundary.enforcement = 2, lexical = FALSE, gradient.check = TRUE)

  # Update par_corr with sorted cut points
  ccc$par_corr <- ccc$par
  ccc$par_corr[1:numcut] <- sort(ccc$par[1:numcut]) + seq(0, gap_global * (numcut - 1), by = gap_global)

  return(ccc)
}

# Function to initialize cut points
initial <- function(target, initial_rr = NULL, initial_cut = NULL, numcut) {
  if (is.null(initial_cut) || is.null(initial_rr)) {
    incut <- quantile(target, probs = seq(0, 1, 1 / (numcut + 1)))
    initial <- c(incut[2:(numcut + 1)], initial_rr)
  } else {
    incut <- array(matrix(c(initial_cut), ncol = 3, byrow = TRUE), dim = c(1, 3))
    inrr <- array(matrix(c(initial_rr), ncol = 3, byrow = TRUE), dim = c(1, 3))
    initial <- c(as.numeric(unlist(incut[1, numcut])), as.numeric(unlist(inrr[1, numcut])))
  }
  return(initial)
}

# Function to define parameter domains
domain_range <- function(target, initial_domain = NULL, numcut, confound = NULL) {
  if (is.null(confound)) {
    vars <- numcut
  } else {
    vars <- numcut + dim(data.frame(confound))[2]
  }

  if (is.null(initial_domain)) {
    cutdomain <- c(quantile(target, probs = 0.1), quantile(target, probs = 0.99))
    rrdomain <- rep(c(-10, 10), vars)
    ran <- if (numcut == 1) {
      c(1, 2)
    } else if (numcut == 2) {
      c(1, 2, 1, 2)
    } else {
      rep(1:2, times = numcut)
    }
    indomain <- matrix(c(cutdomain[ran], rrdomain), ncol = 2, byrow = TRUE)
  } else {
    indomain <- initial_domain[[numcut]]
  }
  return(indomain)
}

# Objective function for genetic algorithm
obj <- function(xx) {
  cutoff <- xx[1:numcut_global] # Cut points
  cut_design <- cut(target_global, breaks = c(0, sort(cutoff) + seq(0, gap_global * (length(cutoff) - 1), by = gap_global), target_max), quantile = FALSE, labels = c(0:numcut_global))
  beta <- xx[(numcut + 1):nvars] # Coefficients of parameters

  if (is.null(confound)) {
    beta1 <- coxph(Surv(time_global, censor_global) ~ cut_design, init = beta, iter.max = 0) # Iteration is zero
  } else {
    beta1 <- coxph(Surv(time_global, censor_global) ~ cut_design + confound_global, init = beta, iter.max = 0) # Iteration is zero
  }

  return(beta1$loglik[2])
}
