#' LPM VaR
#'
#' Generates a value at risk (VaR) quantile based on the Lower Partial Moment ratio.
#'
#' @param percentile numeric [0, 1]; The percentile for left-tail VaR (vectorized).
#' @param degree integer; \code{(degree = 0)} for discrete distributions, \code{(degree = 1)} for continuous distributions.
#' @param x a numeric vector.
#' @return Returns a numeric value representing the point at which \code{"percentile"} of the area of \code{x} is below.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100)
#'
#' ## For 5th percentile, left-tail
#' LPM.VaR(0.05, 0, x)
#' }
#' @export

LPM.VaR <- function(percentile, degree, x){

    if(any(class(x)%in%c("tbl","data.table"))) x <- as.vector(unlist(x))

    percentile <- pmax(pmin(percentile, 1), 0)

    if(degree == 0){
        return(quantile(x, percentile, na.rm = TRUE))
    } else {
        func <- function(b){
            abs(LPM.ratio(degree, b, x) - percentile)
        }
        if(min(x)!=max(x)) return(optimize(func, c(min(x),max(x)))$minimum) else return(min(x))
    }
}

LPM.VaR <- Vectorize(LPM.VaR, vectorize.args = "percentile")


#' UPM VaR
#'
#' Generates an upside value at risk (VaR) quantile based on the Upper Partial Moment ratio
#' @param percentile numeric [0, 1]; The percentile for right-tail VaR (vectorized).
#' @param degree integer; \code{(degree = 0)} for discrete distributions, \code{(degree = 1)} for continuous distributions.
#' @param x a numeric vector.
#' @return Returns a numeric value representing the point at which \code{"percentile"} of the area of \code{x} is above.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#'
#' ## For 5th percentile, right-tail
#' UPM.VaR(0.05, 0, x)
#' @export

UPM.VaR <- function(percentile, degree, x){

    if(any(class(x)%in%c("tbl","data.table"))) x <- as.vector(unlist(x))

    percentile <- pmax(pmin(percentile, 1), 0)

    if(degree==0){
        return(quantile(x, 1 - percentile, na.rm = TRUE))
    } else {
        func <- function(b){
            abs(LPM.ratio(degree, b, x) - (1 - percentile))
        }
        if(min(x)!=max(x)) return(optimize(func, c(min(x),max(x)))$minimum) else return(min(x))
    }

}

UPM.VaR <- Vectorize(UPM.VaR, vectorize.args = "percentile")

