#' Feature Probability
#'
#' Classifies data based on feature probabilities
#'
#' @param x Complete cleaned dataset in matrix form.
#' @param y Column of data to be classified.
#' @param threshold Sets the correlation threshold for independent variables.  Defaults to 0.
#' @param point.est IV data point(s) to be classified, in matrix form.
#' @return Returns variables, \code{"MSE"} mean squared error, \code{"Fitted"} for only the fitted values of the DV, and  \code{"Point.est"} for predicted values.
#' @keywords classifier
#' @author Fred Viole, OVVO Financial Systems
#' @examples
#' ## Using 'iris' dataset where predictive attributes are columns 1:4, and the class is column 5.
#' Feature.probability(iris,5)
#'
#' ## To call mean squared error
#' Feature.probability(iris,5)$MSE
#'
#' ## To call fitted values
#' Feature.probability(iris,5)$Fitted
#'
#' ## To generate a single predicted value
#' Feature.probability(iris,5, point.est=cbind(5.1,3.5,1.4,0.2))$Point.est
#'
#' ## To generate multiple predicted values
#' Feature.probability(iris,5, point.est=(iris[1:10,1:4]))$Point.est
#' @export




Feature.probability = function (x, y,threshold = 0,point.est=NULL) {

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
    corr[i] = VN.dep(x[,i],y,degree=0,print.map = FALSE)$Correlation
   # corr[i]=cor(x[,i],y)
    if(abs(corr[i])<threshold){corr[i]=0}
  }

  kurt.feat = numeric()
  ###  Find kurtosis for all features - uniform distribution translates to the feature present in all output variable instances
  for(i in 1:ncol(x)){
    kurt.feat[i] <- (UPM(4,mean(x[,i]),x[,i])+LPM(4,mean(x[,i]),x[,i]))/
      ((UPM(2,mean(x[,i]),x[,i])+LPM(2,mean(x[,i]),x[,i]))^2)}
  ###  Determine majority output
  majority = as.data.frame(table(y))
  majority = majority[order(majority[,2],decreasing=FALSE),]
  y.fitted <- numeric()
  y.estimate <- numeric()
  sortedUnqY <- sort(unique(y))
  ### Probability a feature observation is in the feature distribution
  for (k in 1:length(y)) {
    probs.cor <- matrix(ncol=4)
    probs.kurt <- matrix(ncol=4)
    point.probs.cor <- matrix(ncol=4)
    point.probs.kurt <- matrix(ncol=4)

    for(item in sortedUnqY) {
      majority.rank = which(majority[,1] == item)
      for(j in 1:ncol(x)) {
        prob.cor  <- numeric()
        prob.kurt <- numeric()
        point.prob.cor  <- numeric()
        point.prob.kurt <- numeric()

        sub_x   <- z[z[,(ncol(x)+1)] == item, j]
        ###  VN.ANOVA adaptation for a single point probability
        intm    <- abs(LPM(1, x[k, j], sub_x)/(LPM(1, x[k, j], sub_x) + UPM(1, x[k, j], sub_x)) - 0.5)
        prob.cor[j] <- corr[j]*((0.5- intm) / 0.5)

        ###  Determine sub.feature kurtosis - higher is more predictive
        kurt.sub<- (UPM(4,mean(sub_x),sub_x)+LPM(4,mean(sub_x),sub_x))/
          ((UPM(2,mean(sub_x),sub_x)+LPM(2,mean(sub_x),sub_x))^2)
        prob.kurt[j] <- (kurt.sub/kurt.feat[j])*((0.5- intm) / 0.5)


        probs.cor   <- rbind(probs.cor, cbind(x[k, j], prob.cor[j], item, majority.rank))
        probs.kurt  <- rbind(probs.kurt, cbind(x[k, j], prob.kurt[j], item, majority.rank))


        if(!is.null(point.est)){
          q=min(k,length(point.est[,1]))
          pointm <- abs(LPM(1, point.est[q,j], sub_x)/(LPM(1, point.est[q,j], sub_x) + UPM(1, point.est[q,j], sub_x)) - 0.5)
          point.prob.cor[j] <- corr[j]*((0.5- pointm) / 0.5)
          point.prob.kurt[j] <- (kurt.sub/kurt.feat[j])*((0.5- pointm) / 0.5)
          point.probs.cor <- rbind(point.probs.cor, cbind(point.est[q,j], point.prob.cor[j], item, majority.rank))
          point.probs.kurt <- rbind(point.probs.kurt, cbind(point.est[q,j], point.prob.kurt[j], item, majority.rank))
        }



      }
    }



    probs.cor[is.na(probs.cor)] <-0
    probs.cor <- probs.cor[-1,]
    probs.kurt[is.na(probs.kurt)] <-0
    probs.kurt <- probs.kurt[-1,]
    if(!is.null(point.est)){
      point.probs.cor[is.na(point.probs.cor)] <-0
      point.probs.cor <- point.probs.cor[-1,]
      point.probs.kurt[is.na(point.probs.kurt)] <-0
      point.probs.kurt <- point.probs.kurt[-1,]
    }

    vote <- numeric()
    point.vote<- numeric()
    ### Each feature then votes for a classification
    for (v in 1:ncol(x)) {
      reduced.probs <- matrix(ncol=4)
      row.seq <- seq(v, length(probs.cor[, 1]), ncol(x))
      reduced.probs <- rbind(probs.cor[(row.seq), ],probs.kurt[(row.seq), ])
      if(!is.null(point.est)){
        reduced.point.probs <- matrix(ncol=4)
        point.row.seq <- seq(v, length(point.probs.cor[, 1]), ncol(x))
        reduced.point.probs <- rbind(point.probs.cor[(point.row.seq), ],point.probs.kurt[(point.row.seq), ])
        reduced.point.probs <-(reduced.point.probs[order(reduced.point.probs[,2],reduced.point.probs[,4],decreasing=TRUE),])
      }
      reduced.probs <-(reduced.probs[order(reduced.probs[,2],reduced.probs[,4],decreasing=TRUE),])


      if (length(unique(reduced.probs[, 2])) >= 1) {
        vote[v] <- reduced.probs[which.max(reduced.probs[, 2]), 3]
        if(!is.null(point.est)){
          point.vote[v] <- reduced.point.probs[which.max(reduced.point.probs[, 2]), 3]
        }
      }
    }

    outcome = as.data.frame(table(vote))
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

    if(!is.null(point.est)){
      qq=min(k,length(point.est[,1]))
      point.outcome = as.data.frame(table(point.vote))
      if (length(point.outcome[, 2]) > 1) {
        point.ties <- (sum((point.outcome[, 2]) == max(point.outcome[, 2])))
      } else {
        point.ties <- 0
      }

      if (point.ties > 1) {
        point.tie.breaker <-  reduced.point.probs[which.max(reduced.point.probs[, 2]), 3]
        y.estimate[qq] <-  as.numeric(as.character(point.tie.breaker))
      }
      else  {
        y.estimate[qq] <- as.numeric(as.character(point.outcome[which.max(point.outcome[, 2]), 1]))
      }
    }




  }


  MSE = mean((y.fitted-y)^2)
  if(!is.null(point.est)){
    return(list("MSE"=MSE,"Fitted"=y.fitted,'Point.est'=y.estimate))} else {
      return(list("MSE"=MSE,"Fitted"=y.fitted))
    }
}
