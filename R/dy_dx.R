#' Partial Derivative dy/dx
#'
#' Returns the numerical partial derivate of y wrt x for a point of interest.
#'
#' @param x Independent Variable
#' @param y Dependent Variable
#' @param order Controls the number of partial moment quadrant means.  Defaults to NULL to allow \link{NNS.reg} to determine optimal order based on R2 of regression.  \code{order='max'} generates a more accurate derivative for well specified cases.
#' @param s.t.n Signal to noise parameter, sets the threshold of \code{NNS.dep} which reduces \code{"order"} when \code{order=NULL}.  Defaults to 0.9 to ensure high dependence for higher \code{"order"} and endpoint determination.
#' @param eval.point Independent variable point to be evaluated.  Defaults to \code{median(x)}.
#' @param deriv.order For second derivative estimate of \code{f(x)}, set \code{deriv.order=2}.  Defaults to first derivative.
#' @param h Percentage step used for finite step method.  Defaults to \code{h=.05} representing a 5 percent step from the value of the independent variable.
#' @param noise.reduction In low signal:noise situations, \code{noise.reduction="median"} uses medians instead of means for partitions, while \code{noise.reduction="mode"} uses modes instead of means for partitions.  \code{noise.reduction=NULL}  (Default setting) allows for maximum possible fit in \link{NNS.reg}.
#' @return Returns the value of the partial derivative estimate for the given order.
#' @keywords partial derivative, nonlinear regression
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' x<-seq(0,2*pi,pi/100); y<-sin(x)
#' dy.dx(x,y,eval.point=1.75)
#' @export


dy.dx <- function(x,y,order=NULL,s.t.n=0.9,eval.point=median(x),deriv.order=1,h=.05,noise.reduction=NULL){
  IV=x
  DV=y

  if(deriv.order==1){
    reg.output <- NNS.reg(IV,DV,plot = FALSE,return.values = TRUE,order=order,s.t.n = s.t.n,noise.reduction = noise.reduction)
    xonly.output <- NNS.reg(IV,DV,type="XONLY",plot=FALSE,return.values = TRUE,order = order,s.t.n = s.t.n,noise.reduction = noise.reduction)


      if(reg.output$R2 >= xonly.output$R2){

      output<- reg.output$derivative
      if(length(output[,1])==1){return(output[,1])}
      if((output[,3][which(eval.point<output[,3])-1][1])<eval.point){
      return(output[,1][which(eval.point<output[,3])][1])}
      else{
        return(mean(c(output[,1][which(eval.point<output[,3])][1],output[,1][which(eval.point<output[,3])-1][1])))
      }
    }
    else {

      output<- xonly.output$derivative
      if(length(output[,1])==1){return(output[,1])}
      if(output[,3][which(eval.point<output[,3])-1][1]<eval.point){
        return(output[,1][which(eval.point<output[,3])][1])}
      else{

        return(mean(c(output[,1][which(eval.point<output[,3])][1],output[,1][which(eval.point<output[,3])-1][1])))
      }
    }
  }

  else{
    ## Second derivative form:
    # f(x+h) - 2(f(x)) +f(x-h) / h^2

    deriv.points=matrix(c((1+h)*eval.point,eval.point,(1-h)*eval.point),ncol = length(eval.point),byrow = TRUE)
    second.deriv.estimates= NNS.reg(IV,DV,plot = FALSE,return.values = TRUE,order=order,point.est = deriv.points,s.t.n = s.t.n,noise.reduction = noise.reduction)$Point.est
    f.x_h = second.deriv.estimates[1]

    two_f.x = 2*second.deriv.estimates[2]

    f.x__h = second.deriv.estimates[3]

    run = ((1+h)*eval.point) - eval.point
  return((f.x_h - two_f.x + f.x__h)/(run^2))

  }

    }


