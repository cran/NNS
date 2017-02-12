#' NNS Causation
#'
#' Returns the causality from observational data between two variables
#'
#' @param x Variable
#' @param y Variable
#' @param tau Number of lagged observations to consider
#' @return Returns the directional causation and quantity of association.
#' @keywords causation
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{NNS.caus(x,y,3)}
#' @export

NNS.caus <- function(x,y,tau){

  if(abs(Uni.caus(x,y,tau))<abs(Uni.caus(y,x,tau))){
      return(c(Causation.x.given.y = Uni.caus(x,y,tau),
                Causation.y.given.x = Uni.caus(y,x,tau),
                "C(y--->x)" = Uni.caus(y,x,tau)-Uni.caus(x,y,tau)))
  } else {
      return(c(Causation.x.given.y = Uni.caus(x,y,tau),
                Causation.y.given.x = Uni.caus(y,x,tau),
                "C(x--->y)" = Uni.caus(x,y,tau)-Uni.caus(y,x,tau)))
    }
}
