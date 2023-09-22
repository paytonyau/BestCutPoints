#' findcut
#'
#' Determines the best location of cut-off points for a continuous risk factor using contingency tables (X^2) approach.
#'
#' @param factor The continuous risk factor.
#' @param outcome The outcome variable. For survival analysis, it should be a matrix with two columns: event and time.
#'                For logistic regression, it should be a binary outcome vector.
#' @param cutnum The number of cut-off points to find.
#' @param datatype The type of analysis to perform, either "survival" or "logistic".
#' @param nmin The minimum number of observations in each group after cutting.
#' @param segment The number of segments to divide the continuous risk factor.
#'
#' @return A list of results depending on the analysis type:
#'   - For "survival" analysis: A list containing various statistics related to cut-off point selection.
#'   - For "logistic" analysis: A list containing various statistics related to cut-off point selection.
#'
#' @examples
#' findcut(factor = BMI, outcome = cbind(event, OS), cutnum = 2, datatype = "survival", nmin = 5, segment = 100)
#' findcut(factor = invasion, outcome = LVSI, cutnum = 2, datatype = "logistic", nmin = 5, segment = 100)
#'
#' @export
#'
#' @seealso [Other functions or packages that are related to this one.]
#'
#' @references
#' Chang, C., Hsieh, M. K., Chang, W. Y., Chiang, A. J., & Chen, J. (2017). Determining the optimal number and location of cutoff points with application to data of cervical cancer. PloS one, 12(4), e0176231.
#'
#' @keywords survival analysis," "logistic regression
#'
#' @family [The package family, if applicable.]
#'
#' @note [Any additional notes or comments.]
#'
#' @author Payton Yau (Package Development)
#'
#' @copyright [Copyright information, if applicable.]
#'
#' @license [License information, if applicable.]
#'

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
    # Prepare the survival data
    delta <- data.frame(outcome)
    colnames(delta) <- c("event", "time")
    userdata <- na.omit(data.frame(delta$event, delta$time, factor))  # Collating survival data
    colnames(userdata) <- c("event", "time", "factor")
    index <- order(userdata$factor)
    userdata <- userdata[index, ]
    n <- dim(userdata)[1]  # Number of subjects
    range <- range(userdata$factor)  # Range of continuous predictor
    cutunit <- diff(range) / segment
    k <- 1  # Count for results
    start <- c()
    end <- c()
    group <- rep(0, cutnum)  # Dummy variable
    ###############################
    # Create result matrix based on the number of cutpoints
    if (cutnum == 1) {
      start[1] <- userdata$factor[nmin]
      end[1] <- userdata$factor[n - nmin * cutnum]
      result <- matrix(NA, ncol = 11, nrow = 100000)
      colnames(result) <- c("Cut1", "Log-rank test", "Gehan-Wilcoxon test",
                            "PH_assumption1", "HR1", "HRlow", "HRup", "P1", "Likelihood ratio test", "Waldtest", "Scoretest")
    } else {
      for (i in 2:cutnum) {
        start[1] <- userdata$factor[nmin]
        end[1] <- userdata$factor[n - nmin * cutnum]
        start[i] <- userdata$factor[which(userdata$factor > start[i - 1])[nmin]]  # Group confirmation
        end[i] <- userdata$factor[n - nmin * (cutnum - i + 1)]
        result <- matrix(NA, ncol = 4 * cutnum + 5, nrow = 100000)
        colnames(result) <- c(paste("Cut", 1:cutnum, sep = ""),
                              "Log-rank test", "Gehan-Wilcoxon test",
                              paste("PH_assumption", 1:cutnum, sep = ""),
                              paste("HR", 1:cutnum, sep = ""), paste("P", 1:cutnum, sep = ""),
                              "Likelihood ratio test", "Waldtest", "Scoretest")
      }
    }

    # Define the process function for recursive calculations
    process <- function(userdata, start, end, i, cutnum, group, n, cutunit, nmin, result, k) {
      if (i == cutnum) {
        while (start[i] < end[i]) {
          if ((i - 1) == 0) {
            group[i] <- sum(userdata$factor <= start[i])
          } else {
            group[i] <- sum(userdata$factor <= start[i]) - sum(group[1:(i - 1)])
          }
          j <- 0:cutnum
          factor_status <- c(rep(j[1:(tail(j + 1, 1) - 1)], group[1:(tail(j + 1, 1) - 1)]), rep(tail(j, 1), n - sum(group[1:tail(j, 1)])))
          coxfit <- coxph(Surv(userdata$time, userdata$event) ~ factor(factor_status))  # Cox proportional hazards model
          model <- summary(coxfit)
          coef <- model$coefficients
          z <- qnorm(1 - (1 - 0.95) / 2)
          lrtest_fit <- survdiff(Surv(userdata$time, userdata$event) ~ factor(factor_status))  # Log-rank test
          wilcoxon_fit <- survdiff(Surv(userdata$time, userdata$event) ~ factor(factor_status), rho = 1)  # Wilcoxon test
          if (i == 1) {
            result[k, ] <- c(start[1], 1 - pchisq(lrtest_fit$chisq, 1), 1 - pchisq(wilcoxon_fit$chisq, 1),
                             cox.zph(coxfit)$table[, 3], coef[, 2], exp(coef[1] - z * coef[3]),
                             exp(coef[1] + z * coef[3]), coef[5], model[["logtest"]]["pvalue"],
                             model[["waldtest"]]["pvalue"], model[["sctest"]]["pvalue"])
          } else {
            result[k, ] <- c(start[1:cutnum], 1 - pchisq(lrtest_fit$chisq, cutnum), 1 - pchisq(wilcoxon_fit$chisq, cutnum),
                             cox.zph(coxfit)$table[, 3][1:cutnum], coef[, 2], coef[, 5],
                             model[["logtest"]]["pvalue"], model[["waldtest"]]["pvalue"], model[["sctest"]]["pvalue"])
          }
          start[i] <- start[i] + cutunit
          k <- k + 1
        }
        if (cutnum == 1) {
          return(result)
        } else {
          return(list(start = start, group = group, result = result, k = k))
        }
      } else {
        while (start[i] < end[i]) {
          if (i == 1) {
            group[i] <- sum(userdata$factor <= start[i])
          } else {
            group[i] <- sum(userdata$factor <= start[i]) - sum(group[1:(i - 1)])
          }
          processin <- process(userdata, start, end, i + 1, cutnum, group, n, cutunit, nmin, result, k)
          group <- processin$group
          start <- processin$start
          result <- processin$result
          k <- processin$k
          start[i] <- start[i] + cutunit
          start[i + 1] <- userdata$factor[sum(userdata$factor <= start[i]) + nmin]
        }
        if (i != 1) {
          return(list(start = start, group = group, result = result, k = k))
        } else {
          return(result)
        }
      }
    }

    # Perform the survival analysis
    allcut <- na.omit(data.frame(process(userdata, start, end, i = 1, cutnum, group, n, cutunit, nmin, result, k)))

    # Create plots and calculate best cutpoints
    if (cutnum == 1) {
      par(mfrow = c(1, 2))
      plotdata <- data.frame(cut = allcut[, "Cut1"], HR = allcut[, "HR1"], HRlow = allcut[, "HRlow"], HRup = allcut[, "HRup"], p = allcut[, "Log.rank.test"])
      plot(plotdata$cut, plotdata$HR, type = "l", lwd = 2, xlab = "Cutoff point", ylab = "Hazard ratio", main = "Cutoff point and Hazard Ratio", ylim = c(min(allcut[, "HRlow"]), max(allcut[, "HR1"])))
      lines(plotdata$cut, plotdata$HRlow, lty = 2)
      lines(plotdata$cut, plotdata$HRup, lty = 2)
      plot(plotdata$cut, plotdata$p, type = "l", lwd = 2, xlab = "Cutoff point", ylab = "p-value", main = "Cutoff point and p-value")
      phokin <- which(allcut[, "PH_assumption1"] > 0.05)
    } else {
      ph <- paste("PH_assumption", 1:cutnum, sep = "")
      phok <- matrix(data = NA, ncol = 10000, nrow = 4)
      for (i in 1:cutnum) {
        maxph <- length(which(allcut[, ph[i]] > 0.05))
        for (x in 1:maxph) {
          phok[i, x] <- which(allcut[, ph[i]] > 0.05)[x]
        }
      }
      phokin <- as.vector(phok[1:cutnum, ])
      for (i in 1:(cutnum - 1)) {
        phokin <- na.exclude(intersect(intersect(phok[i, ], phokin), phok[i + 1, ]))
      }
    }
    ################################################################################
    # Calculate the best cutpoints based on different criteria
    bestcut <- matrix(NA, ncol = cutnum, nrow = 1)
    colnames(bestcut) <- c(paste("Cut", 1:cutnum, sep = ""))
    cutpoint <- c(allcut[which.min(allcut[, "Log.rank.test"]), ][1:cutnum], allcut[which.min(allcut[, "Gehan.Wilcoxon.test"]), ][1:cutnum],
                  allcut[phokin[which.min(allcut[phokin, "Likelihood.ratio.test"])], ][1:cutnum], allcut[phokin[which.min(allcut[phokin, "Waldtest"])], ][1:cutnum],
                  allcut[phokin[which.min(allcut[phokin, "Scoretest"])], ][1:cutnum])
    cutmatrix <- matrix(as.numeric(cutpoint), ncol = cutnum, byrow = TRUE)
    count <- c()
    for (i in 1:5) {
      count[i] <- paste(cutmatrix[i, ], collapse = "")
    }
    tablenumber <- data.frame(table(count))
    bestcut[1, ] <- cutmatrix[which(count == tablenumber[which.max(tablenumber$Freq), 1])[1], ]

    # Prepare the output based on the number of cutpoints
    if (cutnum == 1) {
      output <- list(allcut = allcut,
                     logranktest = allcut[which.min(allcut[, "Log.rank.test"]), ],
                     wilcoxon = allcut[which.min(allcut[, "Gehan.Wilcoxon.test"]), ],
                     logtest = allcut[phokin[which.min(allcut[phokin, "Likelihood.ratio.test"])], ],
                     waldtest = allcut[phokin[which.min(allcut[phokin, "Waldtest"])], ],
                     scoretest = allcut[phokin[which.min(allcut[phokin, "Scoretest"])], ],
                     HRabsmax = allcut[phokin[which.max(abs(allcut[phokin, "HR1"] - 1))], ])
    } else {
      output <- list(allcut = allcut,
                     logranktest = allcut[which.min(allcut[, "Log.rank.test"]), ],
                     wilcoxon = allcut[which.min(allcut[, "Gehan.Wilcoxon.test"]), ],
                     logtest_likelihood.ratio.test = allcut[phokin[which.min(allcut[phokin, "Likelihood.ratio.test"])], ],
                     waldtest = allcut[phokin[which.min(allcut[phokin, "Waldtest"])], ],
                     scoretest = allcut[phokin[which.min(allcut[phokin, "Scoretest"])], ])
    }

    return(output)
    }
  #datatype=survival end
  #----------------------------------------logit regression---------------------------------------------------------#
  # Check if the analysis is logistic regression
  if (datatype == "logistic") {
    # Prepare the logistic regression data
    userdata <- na.omit(data.frame(outcome, factor))
    index <- order(userdata$factor)
    userdata <- userdata[index, ]
    n <- dim(userdata)[1]  # Number of subjects
    range <- range(userdata$factor)  # Range of continuous predictor
    cutunit <- diff(range) / segment

    # Fit a basic logistic regression model
    logist_red <- glm(userdata$outcome ~ 1, family = binomial())  # No other variable, so reduce the model
    k <- 1  # Count for results
    start <- c()
    end <- c()
    group <- rep(0, cutnum)  # Dummy variable
    #----------------------------------------cut---------------------------------------------------------#
    # Create result matrix based on the number of cutpoints
    if (cutnum == 1) {
      start[1] <- userdata$factor[nmin]
      end[1] <- userdata$factor[n - nmin * cutnum]
      result <- matrix(NA, ncol = 13, nrow = 100000)
      colnames(result) <- c("Cut1", "OR1", "OR_up", "OR_low", "P1", "Likelihood ratio test", "Waldtest", "Scoretest", "Specificity", "Sensitivity", "Youden",
                            "Fisher's Exact Test", "AUC")
    } else {
      for (i in 2:cutnum) {
        start[1] <- userdata$factor[nmin]
        end[1] <- userdata$factor[n - nmin * cutnum]
        start[i] <- userdata$factor[which(userdata$factor > start[i - 1])[nmin]]  # Confirm group2 has 20 people and where cut2 to start cutting
        end[i] <- userdata$factor[n - nmin * (cutnum - i + 1)]
        result <- matrix(NA, ncol = 3 * cutnum + 11, nrow = 100000)
        colnames(result) <- c(paste("Cut", 1:cutnum, sep = ""), paste("OR", 1:cutnum, sep = ""), paste("P", 1:cutnum, sep = ""),
                              "Likelihood ratio test", "Waldtest", "Scoretest", "Specificity", "Sensitivity", "Lower", "Higher", "Youden",
                              "euclidean", "manhattan", "AUC")
      }
    }


    # Define the process function for recursive calculations
    process <- function(userdata, start, end, i, cutnum, group, n, cutunit, nmin, result, k) {
      # Check if we have reached the last cutpoint
      if (i == cutnum) {
        while (start[i] < end[i]) {
          # Calculate group based on cutoff values
          if (i == 1) {
            group[i] <- sum(userdata$factor <= start[i])
          } else {
            group[i] <- sum(userdata$factor <= start[i]) - sum(group[1:(i - 1)])
          }
          # Create a factor_status variable
          j <- 0:cutnum
          factor_status <- c(rep(j[1:(tail(j + 1, 1) - 1)], group[1:cutnum]), rep(tail(j, 1), n - sum(group[1:tail(j, 1)])))
          # Fit a logistic regression model
          logist_fit <- glm(userdata$outcome ~ factor(factor_status), family = binomial())
          model <- summary(logist_fit)
          coef <- model$coefficients
          table <- table(userdata$outcome, factor_status)
          Fisher.test <- fisher.test(table)$p.value
          roc <- roc(userdata$outcome, factor_status)
          spe <- roc$specificities * 100
          sen <- roc$sensitivities * 100
          euc <- sqrt((rep(100, cutnum + 2) - spe)^2 + (rep(100, cutnum + 2) - sen)^2)
          manhan <- abs(rep(100, cutnum + 2) - spe) + abs(sen - rep(100, cutnum + 2))
          LRT <- anova(logist_fit, logist_red, test = "Chisq")
          Scoretest <- anova(logist_fit, logist_red, test = "Rao")
          z <- qnorm(1 - (1 - 0.95) / 2)

          if (i == 1) {
            # Handle the first cutpoint differently
            Waldtest <- wald.test(b = coef(logist_fit), Sigma = vcov(logist_fit), Terms = 2)
            result[k, ] <- c(start[1], exp(coef[2]), exp((coef[2] + z * coef[4])), exp((coef[2] - z * coef[4])), coef[8], LRT$P[2], Waldtest$result$chi2[3], Scoretest$P[2],
                             spe[2], sen[2], sen[2] + spe[2] - 1, Fisher.test, as.numeric(roc$auc))
          } else {
            # Handle subsequent cutpoints differently
            Waldtest <- wald.test(b = coef(logist_fit), Sigma = vcov(logist_fit), Terms = 2:(cutnum + 1))
            maxroc <- which.max(sen + spe)
            mineuc <- which.min(euc)
            minman <- which.min(manhan)
            result[k, ] <- c(start[1:cutnum], exp(coef[2:(cutnum + 1)]), coef[(3 * cutnum + 5):(4 * cutnum + 4)], LRT$P[2], Waldtest$result$chi2[3], Scoretest$P[2],
                             spe[maxroc], sen[maxroc], sen[maxroc] - 100 * z * sqrt((sen[maxroc] / 100) * (1 - sen[maxroc] / 100) / (sum(na.exclude(outcome) == 1))),
                             sen[maxroc] + 100 * z * sqrt((sen[maxroc] / 100) * (1 - sen[maxroc] / 100) / (sum(na.exclude(outcome) == 1))),
                             sen[maxroc] + spe[maxroc] - 1, sqrt((sen[mineuc] - 100)^2 + (spe[mineuc] - 100)^2),
                             abs(sen[minman] - 100) + abs(100 - spe[minman]), as.numeric(roc$auc))
          }

          start[i] <- start[i] + cutunit
          k <- k + 1
        }

        # Return the result matrix or a list with updated variables
        if (cutnum == 1) {
          return(result)
        } else {
          return(list(start = start, group = group, result = result, k = k))
        }
      } else {
        # Recursive case for handling multiple cutpoints
        while (start[i] < end[i]) {
          if (i == 1) {
            group[i] <- sum(userdata$factor <= start[i])
          } else {
            group[i] <- sum(userdata$factor <= start[i]) - sum(group[1:(i - 1)])
          }
          # Recursive call to process for the next cutpoint
          processin <- process(userdata, start, end, i + 1, cutnum, group, n, cutunit, nmin, result, k)
          group <- processin$group
          start <- processin$start
          result <- processin$result
          k <- processin$k
          start[i] <- start[i] + cutunit
          start[i + 1] <- userdata$factor[sum(userdata$factor <= start[i]) + nmin]
        }

        # Return a list with updated variables or the final result matrix
        if (i != 1) {
          return(list(start = start, group = group, result = result, k = k))
        } else {
          return(result)
        }
      }
    }


    # Perform the recursive analysis and store the results
    allcut <- na.omit(data.frame(process(userdata, start, end, i = 1, cutnum, group, n, cutunit, nmin, result, k)))

    # Generate plots and calculate best cutpoints
    if (cutnum == 1) {
      par(mfrow = c(1, 2))
      plotdata <- data.frame(cut = allcut[, "Cut1"], OR = allcut[, "OR1"], ORup = allcut[, "OR_up"], ORlow = allcut[, "OR_low"], p = allcut[, "Likelihood.ratio.test"])
      plot(plotdata$cut, plotdata$OR, type = "l", lwd = 2, xlab = "Cutoff point", ylab = "Odds ratio", main = "Cutoff point and Odds Ratio", ylim = c(min(allcut[, "OR_low"]), max(allcut[, "OR1"])))
      lines(plotdata$cut, plotdata$ORlow, lty = 2)
      lines(plotdata$cut, plotdata$ORup, lty = 2)
      plot(plotdata$cut, plotdata$p, type = "l", lwd = 2, xlab = "Cutoff point", ylab = "p-value", main = "Cutoff point and p-value")
      allcut$euclidean <- sqrt((100 - allcut$Sensitivity)^2 + (100 - allcut$Specificity)^2)
      allcut$manhattan <- 100 - allcut$Sensitivity + 100 - allcut$Specificity
    }
    #-----------------------------------------------------------------------------------------------#
    # Calculate the best cutpoints based on different criteria
    bestcut <- matrix(NA, ncol = cutnum, nrow = 1)
    colnames(bestcut) <- c(paste("Cut", 1:cutnum, sep = ""))
    if (cutnum == 1) {
      cutpoint <- c(allcut[which.min(allcut[, "Likelihood.ratio.test"]), ][1:cutnum], allcut[which.min(allcut[, "Waldtest"]), ][1:cutnum],
                    allcut[which.min(allcut[, "Scoretest"]), ][1:cutnum], allcut[which.min(allcut$Fisher.s.Exact.Test), ][1:cutnum],
                    allcut[which.max(allcut[, "Youden"]), ][1:cutnum], allcut[which.min(allcut[, "euclidean"]), ][1:cutnum],
                    allcut[which.max(allcut[, "manhattan"]), ][1:cutnum], allcut[which.max(allcut[, "AUC"]), ][1:cutnum])
    } else {
      cutpoint <- c(allcut[which.min(allcut[, "Likelihood.ratio.test"]), ][1:cutnum], allcut[which.min(allcut[, "Waldtest"]), ][1:cutnum],
                    allcut[which.min(allcut[, "Scoretest"]), ][1:cutnum], allcut[which.max(allcut[, "Youden"]), ][1:cutnum],
                    allcut[which.min(allcut[, "euclidean"]), ][1:cutnum], allcut[which.max(allcut[, "manhattan"]), ][1:cutnum],
                    allcut[which.max(allcut[, "AUC"]), ][1:cutnum])
    }
    cutmatrix <- matrix(as.numeric(cutpoint), ncol = cutnum, byrow = TRUE)
    count <- c()
    if (cutnum == 1) {
      for (i in 1:8) {
        count[i] <- paste(cutmatrix[i, ], collapse = "")
      }
    } else {
      for (i in 1:7) {
        count[i] <- paste(cutmatrix[i, ], collapse = "")
      }
    }
    tablenumber <- data.frame(table(count))
    bestcut[1, ] <- cutmatrix[which(count == tablenumber[which.max(tablenumber$Freq), 1])[1], ]

    # Prepare the output based on the number of cutpoints
    if (cutnum == 1) {
      output <- list(allcut = allcut,
                     Likelihood.ratio.test = allcut[which.min(allcut[, "Likelihood.ratio.test"]), ],
                     waldtest = allcut[which.min(allcut[, "Waldtest"]), ],
                     scoretest = allcut[which.min(allcut[, "Scoretest"]), ],
                     ORabsmax = allcut[which.max(abs(allcut$OR1)), ],
                     youden = allcut[which.max(allcut$Youden), ],
                     Fisher.test = allcut[which.min(allcut$Fisher.s.Exact.Test), ],
                     euclidean = allcut[which.min(allcut$euclidean), ],
                     manhattan = allcut[which.min(allcut$manhattan), ],
                     AUC = allcut[which.max(allcut[, "AUC"]), ])
    } else {
      output <- list(allcut = allcut,
                     Likelihood.ratio.test = allcut[which.min(allcut[, "Likelihood.ratio.test"]), ],
                     waldtest = allcut[which.min(allcut[, "Waldtest"]), ],
                     scoretest = allcut[which.min(allcut[, "Scoretest"]), ],
                     youden = allcut[which.max(allcut$Youden), ],
                     euclidean = allcut[which.min(allcut[, "euclidean"]), ],
                     manhattan = allcut[which.min(allcut[, "manhattan"]), ],
                     AUC = allcut[which.max(allcut[, "AUC"]), ])
    }
    # Return the output
    return(output)
  }
}
