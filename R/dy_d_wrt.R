#' Partial Derivative dy/d[wrt]
#'
#' Returns the numerical partial derivate of y with respect to [wrt] any independent variable for a point of interest.
#'
#' @param B Complete dataset of independent variables (IV) in matrix form.
#' @param y Dependent Variable
#' @param wrt Selects the IV to differentiate with respect to.
#' @param order VN.reg order, defaults to 'max'.
#' @param local.points IV points to be evaluated.
#' @param h Percentage step used for finite step method.  Defaults to \code{h=.1} representing a 10 percent step from the value of the IV.
#' @param n.best Sets the number of closest regression points to use in kernel weighting.  Defaults to number of independent variables.
#' @param mixed If mixed derivative is to be evaluated, set \code{mixed=TRUE}.  Defaults to FALSE.
#' @keywords partial derivative, nonlinear regression
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123);x_1<-runif(100);x_2<-runif(100); y<-x_1^2*x_2^2
#' B=cbind(x_1,x_2)
#' ## To find derivative of y wrt 1st independent variable
#' dy.d_(B,y,wrt=1,local.points=c(.5,.5))
#' @export


dy.d_<- function(B,y,wrt,local.points,order='max',h=.1,n.best=NULL,mixed=FALSE){
  original.local.points.min=local.points
  original.local.points.max=local.points

  original.local.points.min[wrt] = (1-h)*original.local.points.min[wrt]
  original.local.points.max[wrt] = (1+h)*original.local.points.max[wrt]

  deriv.points = matrix(c(original.local.points.min,local.points,original.local.points.max),ncol=length(local.points),byrow = TRUE)

  estimates = VN.reg(B,y,order=order,point.est = deriv.points,n.best=n.best)$prediction
  lower=estimates[1]
  two.f.x = 2*estimates[2]
  upper=estimates[3]

  rise = upper-lower
  run = original.local.points.max[wrt]-local.points[wrt]

  if(mixed==TRUE){
  if(length(local.points)!=2){return("Mixed Derivatives are only for 2 IV")}
  mixed.deriv.points = matrix(c((1+h)*local.points,
                                (1-h)*local.points[1],(1+h)*local.points[2],
                                (1+h)*local.points[1],(1-h)*local.points[2],
                                (1-h)*local.points),ncol=2,byrow = TRUE)


  mixed.estimates = VN.reg(B,y,order=order,point.est=mixed.deriv.points,n.best = n.best)$prediction
  mixed.first = mixed.estimates[1]

  mixed.second = mixed.estimates[2]

  mixed.third = mixed.estimates[3]

  mixed.fourth = mixed.estimates[4]


  return(list("First Derivative"=rise/(2*run),"Second Derivative"=(upper - two.f.x + lower)/(run^2),"Mixed Derivative"=(mixed.first-mixed.second-mixed.third+mixed.fourth)/(4*(run^2))))
  } else
  { return(list("First Derivative"=rise/(2*run),"Second Derivative"=(upper - two.f.x + lower)/(run^2)))}

}
