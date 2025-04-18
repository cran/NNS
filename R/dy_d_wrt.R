#' Partial Derivative dy/d_[wrt]
#'
#' Returns the numerical partial derivative of \code{y} with respect to [wrt] any regressor for a point of interest.  Finite difference method is used with \link{NNS.reg} estimates as \code{f(x + h)} and \code{f(x - h)} values.
#'
#' @param x a numeric matrix or data frame.
#' @param y a numeric vector with compatible dimensions to \code{x}.
#' @param wrt integer; Selects the regressor to differentiate with respect to (vectorized).
#' @param eval.points numeric or options: ("obs", "apd", "mean", "median", "last"); Regressor points to be evaluated.
#' \itemize{
#' \item Numeric values must be in matrix or data.frame form to be evaluated for each regressor, otherwise, a vector of points will evaluate only at the \code{wrt} regressor.  See examples for use cases.
#' \item Set to \code{(eval.points = "obs")} (default) to find the average partial derivative at every observation of the variable with respect to \emph{for specific tuples of given observations.}
#' \item Set to \code{(eval.points = "apd")} to find the average partial derivative at every observation of the variable with respect to \emph{over the entire distribution of other regressors.}
#' \item Set to \code{(eval.points = "mean")} to find the partial derivative at the mean of value of every variable.
#' \item Set to \code{(eval.points = "median")} to find the partial derivative at the median value of every variable.
#' \item Set to \code{(eval.points = "last")} to find the partial derivative at the last observation of every value (relevant for time-series data).
#' }
#' @param mixed logical; \code{FALSE} (default) If mixed derivative is to be evaluated, set \code{(mixed = TRUE)}.
#' @param messages logical; \code{TRUE} (default) Prints status messages.
#' @return Returns column-wise matrix of wrt regressors:
#' \itemize{
#' \item{\code{dy.d_(...)[, wrt]$First}} the 1st derivative
#' \item{\code{dy.d_(...)[, wrt]$Second}} the 2nd derivative
#' \item{\code{dy.d_(...)[, wrt]$Mixed}} the mixed derivative (for two independent variables only).
#' }
#'
#'
#' @note For binary regressors, it is suggested to use \code{eval.points = seq(0, 1, .05)} for a better resolution around the midpoint.
#'
#' @author Fred Viole, OVVO Financial Systems
#'
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#'
#' Vinod, H. and Viole, F. (2020) "Comparing Old and New Partial Derivative Estimates from Nonlinear Nonparametric Regressions"  \doi{10.2139/ssrn.3681104}
#'
#' @examples
#' \dontrun{
#' set.seed(123) ; x_1 <- runif(1000) ; x_2 <- runif(1000) ; y <- x_1 ^ 2 * x_2 ^ 2
#' B <- cbind(x_1, x_2)
#'
#' ## To find derivatives of y wrt 1st regressor for specific points of both regressors
#' dy.d_(B, y, wrt = 1, eval.points = t(c(.5, 1)))
#'
#' ## To find average partial derivative of y wrt 1st regressor,
#' only supply 1 value in [eval.points], or a vector of [eval.points]:
#' dy.d_(B, y, wrt = 1, eval.points = .5)
#'
#' dy.d_(B, y, wrt = 1, eval.points = fivenum(B[,1]))
#'
#'
#' ## To find average partial derivative of y wrt 1st regressor,
#' for every observation of 1st regressor:
#' apd <- dy.d_(B, y, wrt = 1, eval.points = "apd")
#' plot(B[,1], apd[,1]$First)
#'
#' ## 95% Confidence Interval to test if 0 is within
#' ### Lower CI
#' LPM.VaR(.025, 0, apd[,1]$First)
#'
#' ### Upper CI
#' UPM.VaR(.025, 0, apd[,1]$First)
#' }
#' @export



dy.d_ <- function(x, y, wrt,
                  eval.points = "obs",
                  mixed = FALSE,
                  messages = TRUE){
  
  
  
  n <- nrow(x)
  l <- ncol(x)
  
  if(is.null(l)) stop("Please ensure (x) is a matrix or data.frame type object.")
  if(l < 2) stop("Please use NNS::dy.dx(...) for univariate partial derivatives.")
  
  dummies <- list()
  for(i in 1:l){
    dummies[[i]] <- factor_2_dummy_FR(x[,i])
    if(!is.null(ncol(dummies[i][[1]]))) colnames(dummies[i][[1]]) <- paste0(colnames(x)[i], "_", colnames(dummies[i][[1]]))
  }
  x <- do.call(cbind, dummies)
  
  if(messages) message("Currently generating NNS.reg finite difference estimates...Regressor ", wrt,"\r", appendLF=TRUE)
  
  
  if(is.null(colnames(x))){
    colnames.list <- lapply(1 : l, function(i) paste0("X", i))
    colnames(x) <- as.character(colnames.list)
  }
  
  if(any(class(x)%in%c("tbl","data.table")))  x <- as.data.frame(x)
  if(!is.null(y) && any(class(y)%in%c("tbl","data.table"))) y <- as.vector(unlist(y))
  
  if(l != 2) mixed <- FALSE
  
  if(is.character(eval.points)){
    eval.points <- tolower(eval.points)
    if(eval.points == "median"){
      eval.points <- t(apply(x, 2, median))
    } else {
      if(eval.points == "last"){
        eval.points <- tail(x, 1)
      } else {
        if(eval.points == "mean"){
          eval.points <- t(apply(x, 2, mean))
        } else {
          if(eval.points == "apd"){
            eval.points <- as.vector(x[ , wrt, drop = FALSE])
          } else {
            eval.points <- x
          }
        }
      }
    }
  }
  
  original.eval.points.min <- eval.points
  original.eval.points.max <- eval.points
  original.eval.points <- eval.points
  
  norm.matrix <- apply(x, 2, function(z) NNS.rescale(z, 0, 1))
 
  zz <- max(NNS.dep(x[,wrt], y, asym = TRUE)$Dependence, NNS.copula(cbind(x[,wrt],x[,wrt],y)), NNS.copula(cbind(norm.matrix[,wrt], norm.matrix[,wrt], y)))
 
  h_s <- seq(2, 10, 2)

  results <- vector(mode = "list", length(h_s))
  
  for(h in h_s){
    index <- which(h == h_s)
    if(is.vector(eval.points) || ncol(eval.points) == 1){
      eval.points <- unlist(eval.points)

      h_step <- gravity(abs(diff(x[,wrt]))) * h_s[index]
      
      if(h_step==0) h_step <- ((abs((max(x[,wrt]) - min(x[,wrt])) ))/length(x[,wrt])) * h_s[index]
      
      original.eval.points.min <- original.eval.points.min - h_step
      original.eval.points.max <- h_step + original.eval.points.max
      
      seq_by <- max(.01, (1 - zz)/2)
      
      deriv.points <- apply(x, 2, function(z) LPM.VaR(seq(0,1,seq_by), 1, z))
      
      sampsize <- length(seq(0, 1, seq_by))
      
      if(ncol(deriv.points)!=ncol(x)){
        deriv.points <- matrix(deriv.points, ncol = l, byrow = FALSE)
      }
      
      
      
      deriv.points <- data.table::data.table(do.call(rbind, replicate(3*length(eval.points), deriv.points, simplify = FALSE)))
     
      data.table::set(deriv.points, i = NULL, j = as.integer(wrt), value = rep(unlist(rbind(original.eval.points.min,
                                                                                            eval.points,
                                                                                            original.eval.points.max))
                                                                               , each = sampsize, length.out = nrow(deriv.points) ))
      
      
      colnames(deriv.points) <- colnames(x)
      
      distance_wrt <- h_step
      
      
      position <- rep(rep(c("l", "m", "u"), each = sampsize), length.out = nrow(deriv.points))
      id <- rep(1:length(eval.points), each = 3*sampsize, length.out = nrow(deriv.points))
      
      
      if(messages) message(paste("Currently evaluating the ", nrow(deriv.points), " required points "  ), index, " of ", length(h_s),"\r", appendLF=FALSE)
      
      
      
      estimates <- NNS.reg(x, y, point.est = deriv.points, dim.red.method = "equal", plot = FALSE, threshold = 0, order = NULL, point.only = TRUE, ncores = 1)$Point.est
      
      estimates <- data.table::data.table(cbind(estimates = estimates,
                                                position = position,
                                                id = id))
      
      lower_msd <- estimates[position=="l", sapply(.SD, function(x) list(mean=gravity(as.numeric(x)), sd=sd(as.numeric(x)))), .SDcols = "estimates", by = id]
      lower <- lower_msd$V1
      lower_sd <- lower_msd$V2
      
      fx_msd <- estimates[position=="m", sapply(.SD, function(x) list(mean=gravity(as.numeric(x)), sd=sd(as.numeric(x)))), .SDcols = "estimates", by = id]
      f.x <- fx_msd$V1
      f.x_sd <- fx_msd$V2
      
      upper_msd <- estimates[position=="u", sapply(.SD, function(x) list(mean=gravity(as.numeric(x)), sd=sd(as.numeric(x)))), .SDcols = "estimates", by = id]
      upper <- upper_msd$V1
      upper_msd <- upper_msd$V2
      
      rise_1 <- upper - f.x 
      rise_2 <- f.x - lower
      
    } else {
      
      n <- nrow(eval.points)
      original.eval.points <- eval.points

      h_step <- gravity(abs(diff(x[,wrt]))) * h_s[index]
      
      if(h_step==0) h_step <- ((abs((max(x[,wrt]) - min(x[,wrt])) ))/length(x[,wrt])) * h_s[index]
      
      original.eval.points.min[ , wrt] <- original.eval.points.min[ , wrt] - h_step
      original.eval.points.max[ , wrt] <- h_step + original.eval.points.max[ , wrt]
      
      deriv.points <- rbind(original.eval.points.min,
                            original.eval.points,
                            original.eval.points.max)
      
      if(messages) message("Currently generating NNS.reg finite difference estimates...bandwidth ", index, " of ", length(h_s),"\r" ,appendLF=FALSE)
      
      
      estimates <- NNS.reg(x, y, point.est = deriv.points, dim.red.method = "equal", plot = FALSE, threshold = 0, order = NULL, point.only = TRUE, ncores = 1)$Point.est

      lower <- head(estimates,n)
      f.x <- estimates[(n+1):(2*n)]
      upper <- tail(estimates,n)
      
      rise_1 <- upper - f.x 
      rise_2 <- f.x - lower
      
      distance_wrt <- h_step
    }
    
    
    if(mixed){
      if(is.null(dim(eval.points))){
        if(length(eval.points)!=2) stop("Mixed Derivatives are only for 2 IV")
      } else {
        if(ncol(eval.points) != 2) stop("Mixed Derivatives are only for 2 IV")
      }
      
      if(!is.null(dim(eval.points))){
        h_step_1 <- gravity(abs(diff(x[,1]))) * h_s[index]
        if(h_step_1==0) h_step_1 <- ((abs((max(x[,1]) - min(x[,1])) ))/length(x[,1])) * h_s[index]
        
        
        h_step_2 <- gravity(abs(diff(x[,2]))) * h_s[index]
        if(h_step_2==0) h_step_2 <- ((abs((max(x[,2]) - min(x[,2])) ))/length(x[,2])) * h_s[index]
        
        mixed.deriv.points <- matrix(c(h_step_1 + eval.points[,1], h_step_2 + eval.points[,2],
                                       eval.points[,1] - h_step_1, h_step_2 + eval.points[,2],
                                       h_step_1 + eval.points[,1], eval.points[,2] - h_step_2,
                                       eval.points[,1] - h_step_1, eval.points[,2] - h_step_2), ncol = 2, byrow = TRUE)
        
        mixed.distances <- 4 * (h_step_1  * h_step_2)
        
      } else {
        mixed.deriv.points <- matrix(c(h_step + eval.points,
                                       eval.points[1] - h_step, h_step + eval.points[2],
                                       h_step + eval.points[1], eval.points[2] - h_step,
                                       eval.points - h_step), ncol = 2, byrow = TRUE)
        
        mixed.distances <- 4 * (h_step^2)
      }
      
      
      mixed.estimates <- NNS.reg(x, y, point.est = mixed.deriv.points, dim.red.method = "equal", plot = FALSE, threshold = 0, order = NULL, point.only = TRUE, ncores = 1)$Point.est
      
      
      z <- matrix(mixed.estimates, ncol=4, byrow=TRUE)
      z <- z[,1] + z[,4] - z[,2] - z[,3]
      mixed <- (z / mixed.distances)
      
      results[[index]] <- list("First" = (rise_1 + rise_2)/(2 * distance_wrt),
                               "Second" = (upper - f.x + lower) / ((distance_wrt) ^ 2),
                               "Mixed" = mixed)
      
    } else {
      results[[index]] <- list("First" = (rise_1 + rise_2)/(2 * distance_wrt),
                               "Second" = (upper - f.x + lower) / ((distance_wrt) ^ 2) )
    }
    
    
  }

  if(mixed){
    final_results <- list("First" = apply((do.call(cbind, (lapply(results, `[[`, 1)))), 1, function(x) mean(rep(x, length(x):1))),
                          "Second" = apply((do.call(cbind, (lapply(results, `[[`, 2)))), 1, function(x) mean(rep(x, length(x):1))),
                          "Mixed" = apply((do.call(cbind, (lapply(results, `[[`, 3)))), 1, function(x) mean(rep(x, length(x):1))))
  } else {
    final_results <- list("First" = apply((do.call(cbind, (lapply(results, `[[`, 1)))), 1, function(x) mean(rep(x, length(x):1))),
                          "Second" = apply((do.call(cbind, (lapply(results, `[[`, 2)))), 1, function(x) mean(rep(x, length(x):1))))
    
  }
  if(messages) message("","\r", appendLF=TRUE)
  return(final_results)
  
}

dy.d_ <- Vectorize(dy.d_, vectorize.args = c("wrt"))