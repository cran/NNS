#' Partial Derivative dy/d[wrt]
#'
#' Returns the numerical partial derivate of y with respect to [wrt] any regressor for a point of interest.  Finite difference method is used with \link{VN.reg} estimates as f(x+h) and f(x-h) values.
#'
#' @param B Complete dataset of regressors in matrix form.
#' @param y Dependent Variable
#' @param wrt Selects the regressor to differentiate with respect to.
#' @param order VN.reg order, defaults to 1 for multivariate regressions.  If error, make sure \code{order=1}.
#' @param s.t.n Signal to noise parameter, sets the threshold of \code{VN.dep} which reduces \code{"order"} when \code{order=NULL}.  Defaults to 0.9 to ensure high dependence for higher \code{"order"} and endpoint determination.
#' @param eval.points Regressor points to be evaluated.  Set to \code{eval.points="median"} to find partial derivatives at the median of every variable.  Set to \code{eval.points="last"} to find partial derivatives at the last value of every variable.
#' @param h Percentage step used for finite step method.  Defaults to \code{h=.1} representing a 10 percent step from the value of the regressor.
#' @param n.best Sets the number of closest regression points to use in kernel weighting.  Defaults to 2.
#' @param mixed If mixed derivative is to be evaluated, set \code{mixed=TRUE}.  Defaults to FALSE.
#' @param plot Set to \code{plot=TRUE} to view plot, defaults to FALSE.
#' @param precision Sets the number of regression points for estimates.  Set to \code{"HIGH"} where the limit condition of every observation as a regression point. Defaults to \code{"LOW"}.
#' @param norm Normalizes regressors between 0 and 1 for multivariate regression when set to \code{norm="std"}, or normalizes regressors according to \link{VN.norm} when set to \code{norm="VN"}. Defaults to NULL.
#' @param noise.reduction In low signal:noise situations, \code{noise.reduction="median"} uses medians instead of means for partitions, while \code{noise.reduction="mode"} (Default setting) uses modes instead of means for partitions.  \code{noise.reduction=NULL} allows for maximum possible fit.
#' @return Returns the 1st derivative \code{"First Derivative"}, 2nd derivative \code{"Second Derivative"}, and mixed derivative \code{"Mixed Derivative"} (for two independent variables only).
#' @keywords partial derivative, nonlinear regression
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123);x_1<-runif(100);x_2<-runif(100); y<-x_1^2*x_2^2
#' B=cbind(x_1,x_2)
#' ## To find derivatives of y wrt 1st regressor
#' dy.d_(B,y,wrt=1,eval.points=c(.5,.5))
#' @export


dy.d_<- function(B,y,wrt,eval.points="median",order=NULL,s.t.n=0.9,h=.1,n.best=2,mixed=FALSE,plot=FALSE,precision="LOW",norm=NULL,noise.reduction=NULL){
  if(eval.points[1]=="median"){
    eval.points=numeric()
    eval.points=apply(B,2,median)}
  if(eval.points[1]=="last"){
    eval.points=numeric()
    eval.points=as.numeric(B[length(B[,1]),])}

  original.eval.points.min=eval.points
  original.eval.points.max=eval.points

  original.eval.points.min[wrt] = (1-h)*original.eval.points.min[wrt]
  original.eval.points.max[wrt] = (1+h)*original.eval.points.max[wrt]

  deriv.points = matrix(c(original.eval.points.min,eval.points,original.eval.points.max),ncol=length(eval.points),byrow = TRUE)


  estimates = VN.reg(B,y,order=order,point.est = deriv.points,n.best=n.best,s.t.n = s.t.n,plot=plot,precision = precision,norm=norm,noise.reduction=noise.reduction)$Point.est

  lower=estimates[1]
  two.f.x = 2*estimates[2]
  upper=estimates[3]

  rise = upper-lower

  distance.1 = sqrt(sum(sweep(t(c(original.eval.points.max)),2,t(c(eval.points)))^2))
  distance.2 = sqrt(sum(sweep(t(c(original.eval.points.min)),2,t(c(eval.points)))^2))
  run=distance.1+distance.2

  if(mixed==TRUE){
  if(length(eval.points)!=2){return("Mixed Derivatives are only for 2 IV")}
  mixed.deriv.points = matrix(c((1+h)*eval.points,
                                (1-h)*eval.points[1],(1+h)*eval.points[2],
                                (1+h)*eval.points[1],(1-h)*eval.points[2],
                                (1-h)*eval.points),ncol=2,byrow = TRUE)


  mixed.estimates = VN.reg(B,y,order=order,point.est=mixed.deriv.points,n.best = n.best,s.t.n = s.t.n,plot=plot,precision = precision,noise.reduction=noise.reduction)$Point.est
  mixed.first = mixed.estimates[1]

  mixed.second = mixed.estimates[2]

  mixed.third = mixed.estimates[3]

  mixed.fourth = mixed.estimates[4]



  return(list("First Derivative"=rise/(run),"Second Derivative"=(upper - two.f.x + lower)/((.5*run)^2),"Mixed Derivative"=(mixed.first-mixed.second-mixed.third+mixed.fourth)/(4*((.5*run)^2))))
  } else
  { return(list("First Derivative"=rise/(run),"Second Derivative"=(upper - two.f.x + lower)/((.5*run)^2)))}

}
