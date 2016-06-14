#' Feature Probability
#'
#' Classifies data based on feature probabilities
#'
#' @param x Complete cleaned dataset in matrix form.
#' @param y Column of data to be classified.
#' @param threshold Sets the correlation threshold for independent variables.  Defaults to 0.
#' @return Returns two variables, mean squared error "\code{MSE}" and predicted values "\code{Predictions}" as well as Prediction Accuracy measured by percentage of exact classifications.
#' @keywords classifier
#' @author Fred Viole, OVVO Financial Systems
#' @references Annotated code is available at
#' \url{https://github.com/OVVO-Financial/NNS/blob/Prelim/Feature_Probability}
#' @examples
#' ## Using 'iris' dataset where predictive attributes are columns 1:4, and the class is column 5.
#' \dontrun{Feature.probability(iris,5)}
#'
#' ## To call mean squared error
#' \dontrun{Feature.probability(iris,5)$MSE}
#'
#' ## To call predicted values
#' \dontrun{Feature.probability(iris,5)$Predictions}
#' @export




Feature.probability = function (x, y,threshold = 0) {

  original.columns  <-  ncol(x)
  original.variable <-  x
  new.variable <-  matrix(nrow=nrow(x))
  ###  Turn each column into numeric values
  for (i in 1:ncol(original.variable)) {
    new.variable = cbind(new.variable, as.numeric(original.variable[, i]))
  }
  x <- new.variable[, c(-1, -(y + 1))]
  y <- new.variable[, (y + 1)]

  z <- cbind(x, y)
  corr <- numeric()
  ###  Find correlations for all features and output variable
  for(i in 1:ncol(x)){
    corr[i] = (VN.dep(x[,i],y,1,print.map = FALSE)[1])
    if(abs(corr[i])<threshold){corr[i]=0}
  }

  kurt.feat = numeric()
  ###  Find kurtosis for all features - uniform means feature present in all output variable instances
  for(i in 1:ncol(x)){
  kurt.feat[i] <- (UPM(4,mean(x[,i]),x[,i])+LPM(4,mean(x[,i]),x[,i]))/
    ((UPM(2,mean(x[,i]),x[,i])+LPM(2,mean(x[,i]),x[,i]))^2)}
  ###  Determine majority output
  majority = as.data.frame(table(y))
  majority = majority[order(majority[,2],decreasing=FALSE),]
  y.fitted <- numeric()
  sortedUnqY <- sort(unique(y))
  ### Probability a feature observation is in the feature distribution
  for (k in 1:length(y)) {
    probs.cor <- matrix(ncol=4)
    probs.kurt <- matrix(ncol=4)
    for(item in sortedUnqY) {
      majority.rank = which(majority[,1] == item)
      for(j in 1:ncol(x)) {
        prob.cor  <- numeric()
        prob.kurt <- numeric()
        sub_x   <- z[y == item, j]
        ###  VN.ANOVA adaptation for a single point probability
        intm    <- abs(LPM(1, x[k, j], sub_x)/(LPM(1, x[k, j], sub_x) + UPM(1, x[k, j], sub_x)) - 0.5)
        prob.cor[j] <- corr[j]*((0.5- intm) / 0.5)
        ###  Determine sub.feature kurtosis - higher is more predictive
        kurt.sub<- (UPM(4,mean(sub_x),sub_x)+LPM(4,mean(sub_x),sub_x))/
          ((UPM(2,mean(sub_x),sub_x)+LPM(2,mean(sub_x),sub_x))^2)
        prob.kurt[j] <- (kurt.sub/kurt.feat[j])*((0.5- intm) / 0.5)
        probs.cor   <- rbind(probs.cor, cbind(x[k, j], prob.cor[j], item, majority.rank))
        probs.kurt  <- rbind(probs.kurt, cbind(x[k, j], prob.kurt[j], item, majority.rank))
      }
    }
    probs.cor[is.na(probs.cor)] <-0
    probs.cor <- probs.cor[-1,]
    probs.kurt[is.na(probs.kurt)] <-0
    probs.kurt <- probs.kurt[-1,]
    vote <- numeric()
    ### Each feature then votes for a classification
    for (v in 1:ncol(x)) {
      reduced.probs <- matrix(ncol=4)
      row.seq <- seq(v, length(probs.cor[, 1]), ncol(x))
      reduced.probs <- rbind(probs.cor[(row.seq), ],probs.kurt[(row.seq), ])

      reduced.probs <-(reduced.probs[order(reduced.probs[,2],reduced.probs[,4],decreasing=TRUE),])

      if (length(unique(reduced.probs[, 2])) > 1) {
        vote[v] <- reduced.probs[which.max(reduced.probs[, 2]), 3]
      }
    }
    outcome = as.data.frame(table(vote))
    #return(length(outcome[,2]))
    if (length(outcome[, 2]) > 1) {
      ties <- (sum((outcome[, 2]) == max(outcome[, 2])))
    } else {
      ties <- 0
    }
    if (ties > 1) {
      tie.breaker <-  reduced.probs[which.max(reduced.probs[, 2]), 3]
      y.fitted[k] <-  as.numeric(as.character(tie.breaker))
    }
    else  {
      y.fitted[k] <- as.numeric(as.character(outcome[which.max(outcome[, 2]), 1]))
    }
  }

  MSE = mean((y.fitted-y)^2)
  return(list("MSE"=MSE,"Predictions"=y.fitted))
}
