#' NNS Causation
#'
#' Returns the causality from observational data between two variables.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y \code{NULL} (default) or a numeric vector with compatible dimensions to \code{x}.
#' @param factor.2.dummy logical; \code{FALSE} (default) Automatically augments variable matrix with numerical dummy variables based on the levels of factors.  Includes dependent variable \code{y}.
#' @param tau options: ("cs", "ts", integer); 0 (default) Number of lagged observations to consider (for time series data).  Otherwise, set \code{(tau = "cs")} for cross-sectional data.  \code{(tau = "ts")} automatically selects the lag of the time series data, while \code{(tau = [integer])} specifies a time series lag.
#' @param plot logical; \code{FALSE} (default) Plots the raw variables, tau normalized, and cross-normalized variables.
#' @return Returns the directional causation (x ---> y) or (y ---> x) and net quantity of association.  For causal matrix, directional causation is returned as ([column variable] ---> [row variable]).  Negative numbers represent causal direction attributed to [row variable].
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#' @examples
#'
#' \dontrun{
#' ## x causes y...
#' set.seed(123)
#' x <- rnorm(1000) ; y <- x ^ 2
#' NNS.caus(x, y, tau = "cs")
#'
#' ## Causal matrix without per factor causation
#' NNS.caus(iris, tau = 0)
#'
#' ## Causal matrix with per factor causation
#' NNS.caus(iris, factor.2.dummy = TRUE, tau = 0)
#' }
#' @export

NNS.caus <- function(x, y = NULL,
                     factor.2.dummy = FALSE,
                     tau = 0,
                     plot = FALSE){

  if(!is.null(y))  if(sum(is.na(cbind(x,y))) > 0) stop("You have some missing values, please address.")
  if(is.null(y))  if(sum(is.na(x)) > 0) stop("You have some missing values, please address.")

  orig.tau <- tau
  orig.plot <- plot

  if(any(class(x)%in%c("tbl","data.table")) && dim(x)[2]==1) x <- as.vector(unlist(x))
  if(any(class(x)%in%c("tbl","data.table"))) x <- as.data.frame(x)
  if(!is.null(y) && any(class(y)%in%c("tbl","data.table"))) y <- as.vector(unlist(y))


  if(factor.2.dummy){
    if(!is.null(dim(x))){
      if(!is.numeric(x)){
        x <- do.call(cbind, lapply(x, factor_2_dummy_FR))
      } else {
        x <- apply(x, 2, as.double)
      }
      if(is.list(x)){
        x <- do.call(cbind, x)
        x <- apply(x, 2, as.double)
      }

    } else {
      x <- factor_2_dummy(x)
      if(is.null(dim(x))){
        x <- as.double(x)
      } else {
        x <- apply(x, 2, as.double)
      }
    }
  }

  if(!is.null(y)){
    if(is.factor(y)) y <- as.numeric(y)

    if(is.numeric(tau)){
      Causation.x.given.y <- Uni.caus(x, y, tau = tau, plot = FALSE)
      Causation.y.given.x <- Uni.caus(y, x, tau = tau, plot = FALSE)

      Causation.x.given.y[is.na(Causation.x.given.y)] <- 0
      Causation.y.given.x[is.na(Causation.y.given.x)] <- 0

      if(Causation.x.given.y == Causation.y.given.x ||
         Causation.x.given.y == 0 || Causation.y.given.x == 0){
        Causation.x.given.y <- Uni.caus(x, y, tau = tau, plot = FALSE)
        Causation.y.given.x <- Uni.caus(y, x, tau = tau, plot = FALSE)
        Causation.x.given.y[is.na(Causation.x.given.y)] <- 0
        Causation.y.given.x[is.na(Causation.y.given.x)] <- 0
      }
    }

    if(tau == "cs"){
      Causation.x.given.y <- Uni.caus(x, y, tau = 0, plot = FALSE)
      Causation.y.given.x <- Uni.caus(y, x, tau = 0, plot = FALSE)

      Causation.x.given.y[is.na(Causation.x.given.y)] <- 0
      Causation.y.given.x[is.na(Causation.y.given.x)] <- 0

      if(Causation.x.given.y == Causation.y.given.x ||
         Causation.x.given.y == 0 || Causation.y.given.x == 0){
        Causation.x.given.y <- Uni.caus(x, y, tau = 0, plot = FALSE)
        Causation.y.given.x <- Uni.caus(y, x, tau = 0, plot = FALSE)

        Causation.x.given.y[is.na(Causation.x.given.y)] <- 0
        Causation.y.given.x[is.na(Causation.y.given.x)] <- 0
      }
    }

    if(tau == "ts"){
      Causation.y.given.x <- Uni.caus(y, x, tau = 3, plot = FALSE)
      Causation.x.given.y <- Uni.caus(x, y, tau = 3, plot = FALSE)

      Causation.x.given.y[is.na(Causation.x.given.y)] <- 0
      Causation.y.given.x[is.na(Causation.y.given.x)] <- 0

    }


      if(abs(Causation.x.given.y) <= abs(Causation.y.given.x)){
        if(plot){
          # For plotting only
          if(tau == "cs") tau <- 0

          if(tau == "ts") tau <- 3

          Uni.caus(y, x, tau = tau, plot = plot)
        }
        return(c(Causation.x.given.y = Causation.x.given.y,
                 Causation.y.given.x = Causation.y.given.x,
                 "C(x--->y)" = sign(Causation.y.given.x) * (abs(Causation.y.given.x - Causation.x.given.y)/2)))
      } else {
        if(plot){
          # For plotting only
          if(tau == "cs") tau <- 0
          if(tau == "ts") tau <- 3

          Uni.caus(x, y, tau = tau, plot = plot)
        }
        return(c(Causation.x.given.y = Causation.x.given.y,
                 Causation.y.given.x = Causation.y.given.x,
                 "C(y--->x)" = sign(Causation.x.given.y) * (abs(Causation.x.given.y - Causation.y.given.x)/2)))
      }


  } else {

    NNS.caus.matrix(x, tau = orig.tau)
  }


}
