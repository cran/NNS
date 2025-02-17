#' NNS ARMA
#'
#' Autoregressive model incorporating nonlinear regressions of component series.
#'
#' @param variable a numeric vector.
#' @param h integer; 1 (default) Number of periods to forecast.
#' @param training.set numeric; \code{NULL} (default) Sets the number of variable observations
#'
#'  \code{(variable[1 : training.set])} to monitor performance of forecast over in-sample range.
#' @param seasonal.factor logical or integer(s); \code{TRUE} (default) Automatically selects the best seasonal lag from the seasonality test.  To use weighted average of all seasonal lags set to \code{(seasonal.factor = FALSE)}.  Otherwise, directly input known frequency integer lag to use, i.e. \code{(seasonal.factor = 12)} for monthly data.  Multiple frequency integers can also be used, i.e. \code{(seasonal.factor = c(12, 24, 36))}
#' @param modulo integer(s); NULL (default) Used to find the nearest multiple(s) in the reported seasonal period.
#' @param mod.only logical; \code{TRUE} (default) Limits the number of seasonal periods returned to the specified \code{modulo}.
#' @param weights numeric or \code{"equal"}; \code{NULL} (default) sets the weights of the \code{seasonal.factor} vector when specified as integers.  If \code{(weights = NULL)} each \code{seasonal.factor} is weighted on its \link{NNS.seas} result and number of observations it contains, else an \code{"equal"} weight is used.
#' @param best.periods integer; [2] (default) used in conjunction with \code{(seasonal.factor = FALSE)}, uses the \code{best.periods} number of detected seasonal lags instead of \code{ALL} lags when
#' \code{(seasonal.factor = FALSE, best.periods = NULL)}.
#' @param negative.values logical; \code{FALSE} (default) If the variable can be negative, set to
#' \code{(negative.values = TRUE)}.  If there are negative values within the variable, \code{negative.values} will automatically be detected.
#' @param method options: ("lin", "nonlin", "both", "means"); \code{"nonlin"} (default)  To select the regression type of the component series, select \code{(method = "both")} where both linear and nonlinear estimates are generated.  To use a nonlinear regression, set to
#' \code{(method = "nonlin")}; to use a linear regression set to \code{(method = "lin")}.  Means for each subset are returned with \code{(method = "means")}.
#' @param dynamic logical; \code{FALSE} (default) To update the seasonal factor with each forecast point, set to \code{(dynamic = TRUE)}.  The default is \code{(dynamic = FALSE)} to retain the original seasonal factor from the inputted variable for all ensuing \code{h}.
#' @param shrink logical; \code{FALSE} (default) Ensembles forecasts with \code{method = "means"}.
#' @param plot logical; \code{TRUE} (default) Returns the plot of all periods exhibiting seasonality and the \code{variable} level reference in upper panel.  Lower panel returns original data and forecast.
#' @param seasonal.plot logical; \code{TRUE} (default) Adds the seasonality plot above the forecast.  Will be set to \code{FALSE} if no seasonality is detected or \code{seasonal.factor} is set to an integer value.
#' @param pred.int numeric [0, 1]; \code{NULL} (default) Plots and returns the associated prediction intervals for the final estimate.  Constructed using the maximum entropy bootstrap \link{NNS.meboot} on the final estimates.
#' @return Returns a vector of forecasts of length \code{(h)} if no \code{pred.int} specified.  Else, returns a \code{data.table} with the forecasts as well as lower and upper prediction intervals per forecast point.
#' @note
#' For monthly data series, increased accuracy may be realized from forcing seasonal factors to multiples of 12.  For example, if the best periods reported are: \{37, 47, 71, 73\}  use
#' \code{(seasonal.factor = c(36, 48, 72))}.
#'
#' \code{(seasonal.factor = FALSE)} can be a very computationally expensive exercise due to the number of seasonal periods detected.
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#'
#' Viole, F. (2019) "Forecasting Using NNS"  \doi{10.2139/ssrn.3382300}
#'
#' @examples
#'
#' ## Nonlinear NNS.ARMA using AirPassengers monthly data and 12 period lag
#' \dontrun{
#' NNS.ARMA(AirPassengers, h = 45, training.set = 100, seasonal.factor = 12, method = "nonlin")
#'
#' ## Linear NNS.ARMA using AirPassengers monthly data and 12, 24, and 36 period lags
#' NNS.ARMA(AirPassengers, h = 45, training.set = 120, seasonal.factor = c(12, 24, 36), method = "lin")
#'
#' ## Nonlinear NNS.ARMA using AirPassengers monthly data and 2 best periods lag
#' NNS.ARMA(AirPassengers, h = 45, training.set = 120, seasonal.factor = FALSE, best.periods = 2)
#' }
#' @export



# Autoregressive Model
NNS.ARMA <- function(variable,
                     h = 1,
                     training.set = NULL,
                     seasonal.factor = TRUE,
                     weights = NULL,
                     best.periods = 1,
                     modulo = NULL,
                     mod.only = TRUE,
                     negative.values = FALSE,
                     method = "nonlin",
                     dynamic = FALSE,
                     shrink = FALSE,
                     plot = TRUE,
                     seasonal.plot = TRUE,
                     pred.int = NULL){
  
  
  if(is.numeric(seasonal.factor) && dynamic) stop('Hmmm...Seems you have "seasonal.factor" specified and "dynamic = TRUE".  Nothing dynamic about static seasonal factors!  Please set "dynamic = FALSE" or "seasonal.factor = FALSE"')
  
  if(any(class(variable)%in%c("tbl","data.table"))) variable <- as.vector(unlist(variable))
  
  if(sum(is.na(variable)) > 0) stop("You have some missing values, please address.")
  
  method <- tolower(method)
  if(method == "means") shrink <- FALSE
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  if(!is.null(best.periods) && !is.numeric(seasonal.factor)) seasonal.factor <- FALSE
  
  label <- deparse(substitute(variable))
  variable <- as.numeric(variable)
  OV <- variable
  
  if(min(variable) < 0) negative.values <- TRUE
  
  if(!is.null(training.set)){
    variable <- variable[1 : training.set]
    FV <- variable[1 : training.set]
  } else {
    training.set <- length(variable)
    variable <- variable
    FV <- variable
  }
  
  Estimates <- numeric(length = h)
  
  
  if(is.numeric(seasonal.factor)){
    seasonal.plot = FALSE
    M <- matrix(seasonal.factor, ncol=1)
    colnames(M) <- "Period"
    lag <- seasonal.factor
    output <- numeric(length(seasonal.factor))
    for(i in 1 : length(seasonal.factor)){
      rev.var <- variable[seq(length(variable), 1, -i)]
      output[i] <- abs(sd(rev.var) / mean(rev.var))
    }
    
    if(is.null(weights)){
      Relative.seasonal <- output / abs(sd(variable)/mean(variable))
      Seasonal.weighting <- 1 / Relative.seasonal
      Observation.weighting <- 1 / sqrt(seasonal.factor)
      Weights <- (Seasonal.weighting * Observation.weighting) / sum(Observation.weighting * Seasonal.weighting)
      seasonal.plot <- FALSE
    } else {
      Weights <- weights
    }
    
  } else {
    M <- NNS.seas(variable, plot=FALSE, modulo = modulo, mod.only = mod.only)
    if(!is.list(M)){
      M <- t(1)
    } else {
      if(is.null(best.periods)){
        M <- M$all.periods
      } else {
        if(!seasonal.factor && is.numeric(best.periods) && (length(M$all.periods$Period) < best.periods)){
          best.periods <- length(M$all.periods$Period)
        }
        if(!seasonal.factor && is.null(best.periods)){
          best.periods <- length(M$all.periods$Period)
        }
        M <- M$all.periods[1 : best.periods, ]
      }
    }
    
    ASW <- ARMA.seas.weighting(seasonal.factor, M)
    lag <- ASW$lag
    
    if(is.null(weights)) Weights <- ASW$Weights else Weights <- weights
    
    if(is.character(weights)) Weights <- rep(1/length(lag), length(lag))
    
  }
  
  # Vectorized linear
  Lin.Reg.Estimates <- list()
  Regression.Estimates_means <- list()
  
  if (method == "lin" && is.numeric(seasonal.factor) && length(seasonal.factor) == 1) {
    for(k in lag) {
      lag.idx <- which(k == lag)
      GV.lin <- generate.lin.vectors(variable, lag[lag.idx], h)
      
      # Generate linear regression estimates for each lag
      Lin.Regression.Estimates <- lapply(1:min(h, lag[lag.idx]), function(i) {
        last.xs <- tail(GV.lin$Component.index[[i]], 1)
        lin.reg <- fast_lm(GV.lin$Component.index[[i]], GV.lin$Component.series[[i]])
        coefs <- lin.reg$coef
        
        return(as.numeric(coefs[1] + coefs[2] * unlist(GV.lin$forecast.values[[i]])))
      })
      
      Lin.Reg.Estimates[[lag.idx]] <- unlist(Lin.Regression.Estimates)[order(unlist(GV.lin$forecast.index))] * Weights[lag.idx]
      if((method=="means") || shrink){
        Regression.Estimates_means <- unlist(lapply(GV.lin$Component.series, function(x) mean(x) * Weights[lag.idx]) )
        if(shrink) Lin.Reg.Estimates[[lag.idx]] <- (Lin.Reg.Estimates[[lag.idx]] + Regression.Estimates_means) / 2 else Lin.Reg.Estimates <- Regression.Estimates_means
      }
    }
    
    # Calculate weighted sum of regression estimates for each lag
    Lin.estimates <- Reduce(`+`, Lin.Reg.Estimates)
    
    if(!negative.values) Lin.estimates <- pmax(0, Lin.estimates)
    
    Estimates <- Lin.estimates
    variable <- c(variable, Estimates)
    FV <- variable
  } else {
    
    # Regression for each estimate in h
    for (j in 1:h) {
      # Regenerate seasonal.factor if dynamic
      if (dynamic) {
        seas.matrix <- NNS.seas(variable, plot = FALSE)
        if (!is.list(seas.matrix)) {
          M <- t(1)
        } else {
          if (is.null(best.periods)) {
            M <- seas.matrix$all.periods
            best.periods <- length(M$all.periods$Period)
          } else {
            if (length(M$all.periods$Period) < best.periods) {
              best.periods <- length(M$all.periods$Period)
            }
            M <- seas.matrix$all.periods[1:best.periods, ]
          }
        }
        
        ASW <- ARMA.seas.weighting(seasonal.factor, M)
        lag <- ASW$lag
        Weights <- ASW$Weights
      }
      
      # Re-Generate vectors for 1:lag if dynamic
      GV <- generate.vectors(variable, lag)
      Component.index <- GV$Component.index
      Component.series <- GV$Component.series
      
      # Regression on Component Series
      ## Regression on Component Series
      if (method %in% c("nonlin", "both")) {
        Regression.Estimates <- sapply(seq_along(lag), function(i) {
          x <- Component.index[[i]]
          y <- Component.series[[i]]
          
          last.y <- tail(y, 1)
          
          reg.points <- NNS.reg(x, y, return.values = FALSE, plot = FALSE, multivariate.call = TRUE)
          reg.points <- reg.points[complete.cases(reg.points), ]
          
          xs <- tail(reg.points$x, 1) - reg.points$x
          ys <- tail(reg.points$y, 1) - reg.points$y
          
          xs <- head(xs, -1)
          ys <- head(ys, -1)
          
          run <- mean(rep(xs, (1:length(xs))^2))
          rise <- mean(rep(ys, (1:length(ys))^2))
          
          last.y + (rise / run)
        })
        
        Regression.Estimates <- pmax(0, Regression.Estimates)
        Nonlin.estimates <- sum(Regression.Estimates * Weights)
      }
      
      if (method %in% c("lin", "both", "means")) {
        Lin.Regression.Estimates <- sapply(seq_along(lag), function(i) {
          last.x <- tail(Component.index[[i]], 1)
          lin.reg <- fast_lm(Component.index[[i]], Component.series[[i]])
          coefs <- lin.reg$coef
          return(as.numeric(coefs[1] + coefs[2] * (last.x + 1)))
        })
        
        Lin.Regression.Estimates <- unlist(Lin.Regression.Estimates)
        
        if (method %in% c("means", "shrink")) {
          Regression.Estimates_means <- sapply(Component.series, mean)
          if (shrink) Lin.Regression.Estimates <- (Lin.Regression.Estimates + Regression.Estimates_means) / 2 else Lin.Regression.Estimates <- Regression.Estimates_means
        }
        
        Lin.estimates <- sum(Lin.Regression.Estimates * Weights)
        if(!negative.values)  Lin.estimates <- pmax(0, Lin.estimates)
      }
      
      if (method == "lin") Estimates[j] <- sum(Lin.estimates * Weights)
      if (method == 'both') Estimates[j] <- mean(c(Lin.estimates, Nonlin.estimates))
      if (method == "nonlin")  Estimates[j] <- sum(Nonlin.estimates * Weights)
      
      variable <- c(variable, Estimates[j])
      FV <- variable
    } # j loop
  }
  
  if(!is.null(pred.int)){
    if (method != "means") lin.resid <- mean(abs(Lin.Regression.Estimates - mean(Lin.Regression.Estimates)))
    PIs <- do.call(cbind, NNS.MC(Estimates, lower_rho = 0, upper_rho = 1, by = .2, exp = 2)$replicates)
    lin.resid <- mean(unlist(lin.resid))
    lin.resid[is.na(lin.resid)] <- 0
    
    upper_lower <- apply(PIs, 1, function(z) list(UPM.VaR((1-pred.int)/2, 0, z), abs(LPM.VaR((1-pred.int)/2, 0, z)))) 
    upper_PIs <- as.numeric(lapply(upper_lower, `[[`, 1)) + lin.resid
    lower_PIs <- as.numeric(lapply(upper_lower, `[[`, 2)) - lin.resid
  } else lin.resid <- 0
  
  #### PLOTTING
  if(plot){
    original.par = par(no.readonly = TRUE)
    if(seasonal.plot){
      par(mfrow = c(2, 1))
      if(ncol(M) > 1){
        plot(unlist(M[, 1]), unlist(M[, 2]),
             xlab = "Period", ylab = "Coefficient of Variation", main = "Seasonality Test", ylim = c(0, 1.5 * unlist(M[, 3])[1]))
        points(unlist(M[ , 1]), unlist(M[ , 2]), pch = 19, col = 'red')
        abline(h = unlist(M[, 3])[1], col = "red", lty = 5)
        text((min(unlist(M[ , 1])) + max(unlist(M[ , 1]))) / 2, unlist(M[, 3])[1], pos = 3, "Variable Coefficient of Variation", col = 'red')
      } else {
        plot(1,1, pch = 19, col = 'blue', xlab = "Period", ylab = "Coefficient of Variation", main = "Seasonality Test",
             ylim = c(0, 2 * abs(sd(FV) / mean(FV))))
        text(1, abs(sd(FV) / mean(FV)), pos = 3, "NO SEASONALITY DETECTED", col = 'red')
      }
    }
    
    
    if(is.null(label)) label <- "Variable"
    
    
    if(!is.null(pred.int)){
      plot(OV, type = 'l', lwd = 2, main = "NNS.ARMA Forecast", col = 'steelblue',
           xlim = c(1, max((training.set + h), length(OV))),
           ylab = label, ylim = c(min(Estimates, OV,  unlist(PIs) ), max(OV, Estimates, unlist(PIs) )) )
      
      
      polygon(c((training.set+1) : (training.set+h), rev((training.set+1) : (training.set+h))),
              c(lower_PIs, rev(upper_PIs)),
              col = rgb(1, 192/255, 203/255, alpha = 0.5),
              border = NA)
      
      
      lines(OV, type = 'l', lwd = 2, col = 'steelblue')
      
      lines((training.set + 1) : (training.set + h), Estimates, type = 'l', lwd = 2, lty = 1, col = 'red')
      segments(training.set, FV[training.set], training.set + 1, Estimates[1],lwd = 2,lty = 1,col = 'red')
      legend('topleft', bty = 'n', legend = c("Original", paste0("Forecast ", h, " period(s)")), lty = c(1, 1), col = c('steelblue', 'red'), lwd = 2)
    } else {
      plot(OV, type = 'l', lwd = 2, main = "NNS.ARMA Forecast", col = 'steelblue',
           xlim = c(1, max((training.set + h), length(OV))),
           ylab = label, ylim = c(min(Estimates, OV), max(OV, Estimates)))
      
      if(training.set[1] < length(OV)){
        lines((training.set + 1) : (training.set + h), Estimates, type = 'l',lwd = 2, lty = 3, col = 'red')
        segments(training.set, FV[training.set], training.set + 1, Estimates[1], lwd = 2, lty = 3, col = 'red')
        legend('topleft', bty = 'n', legend = c("Original", paste0("Forecast ", h, " period(s)")), lty = c(1, 2), col = c('steelblue', 'red'), lwd = 2)
      } else {
        lines((training.set + 1) : (training.set + h), Estimates, type = 'l', lwd = 2, lty = 1, col = 'red')
        segments(training.set, FV[training.set], training.set + 1, Estimates[1], lwd = 2, lty = 1, col = 'red')
        legend('topleft', bty = 'n', legend = c("Original", paste0("Forecast ", h, " period(s)")),lty = c(1, 1), col = c('steelblue', 'red'), lwd = 2)
      }
      
      
    }
    points(training.set, OV[training.set], col = "green", pch = 18)
    points(training.set + h, tail(FV, 1), col = "green", pch = 18)
    
    par(original.par)
  }
  
  
  options(warn = oldw)
  if(!is.null(pred.int)){
    
    results <- cbind.data.frame(Estimates,  pmin(Estimates, lower_PIs),  pmax(Estimates, upper_PIs))
    colnames(results) = c("Estimates",
                          paste0("Lower ", round(pred.int*100,2), "% pred.int"),
                          paste0("Upper ", round(pred.int*100,2), "% pred.int"))
    return(data.table::data.table(results))
  } else {
    return(Estimates)
  }
}