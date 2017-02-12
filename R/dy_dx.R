#' Partial Derivative dy/dx
#'
#' Returns the numerical partial derivate of y wrt x for a point of interest.
#'
#' @param x Independent Variable
#' @param y Dependent Variable
#' @param order Controls the number of partial moment quadrant means.  Defaults to \code{order='max'} which generates a more accurate derivative for well specified cases.
#' @param s.t.n Signal to noise parameter, sets the threshold of \code{NNS.dep} which reduces \code{"order"} when \code{order=NULL}.  Defaults to 0.99 to ensure high dependence for higher \code{"order"} and endpoint determination.
#' @param eval.point Independent variable point to be evaluated.  Defaults to \code{eval.point=median(x)}.  Set to \code{eval.points="overall"} to find an overall partial derivative estimate.
#' @param deriv.order For second derivative estimate of \code{f(x)}, set \code{deriv.order=2}.  Defaults to first derivative.
#' @param h Percentage step used for finite step method.  Defaults to \code{h=.01} representing a 5 percent step from the value of the independent variable.
#' @param noise.reduction In low signal to noise situations, \code{noise.reduction="median"} uses medians instead of means for partitions, while \code{noise.reduction="mode"} uses modes instead of means for partitions.  \code{noise.reduction="off"}  allows for maximum possible fit in \link{NNS.reg}. Default setting is \code{noise.reduction="mean"}.
#' @param deriv.method Determines the partial derivative from the coefficient of the \link{NNS.reg} output when \code{deriv.method="NNS"} or generates a partial derivative using the finite step method \code{deriv.method="FS"} (Defualt).
#' @return Returns the value of the partial derivative estimate for the given order.
#' @keywords partial derivative
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' x<-seq(0,2*pi,pi/100); y<-sin(x)
#' dy.dx(x,y,eval.point=1.75)
#' @export

dy.dx <- function(x,y,order=NULL,s.t.n=0.99,eval.point=median(x),deriv.order=1,h=.05,noise.reduction='mean',deriv.method="FS"){

  if(eval.point=='overall'){

  ranges=NNS.reg(x,y,order=order,noise.reduction=noise.reduction,plot=F)$derivative
  range.weights=numeric()
  for(i in 1:length(ranges[,1])){
  range.weights[i]=sum(x>=ranges[i,2]&x<ranges[i,3])/length(x)
  }

  range.weights[length(range.weights)]=range.weights[length(range.weights)]+1/length(x)
  return(sum(ranges[,1]*range.weights))

  } else {

  original.eval.point.min=eval.point
  original.eval.point.max=eval.point

  eval.point.min = (1-h)*original.eval.point.min
  eval.point.max = (1+h)*original.eval.point.max


  run=eval.point.max-eval.point.min

  if(deriv.order==1){

  if(deriv.method=="FS"){
  estimates=NNS.reg(x,y,plot = FALSE,order=order,s.t.n = s.t.n,noise.reduction = noise.reduction,point.est = c(eval.point.min,eval.point.max))$Point.est

    rise=estimates[2]-estimates[1]

    return(rise/run) } else {

     reg.output <- NNS.reg(x,y,plot = FALSE,return.values = TRUE,order=order,s.t.n = s.t.n,noise.reduction = noise.reduction)

     output<- reg.output$derivative
      if(length(output[,1])==1){return(output[,1])}
      if((output[,3][which(eval.point<output[,3])-1][1])<eval.point){
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
    second.deriv.estimates= NNS.reg(x,y,plot = FALSE,return.values = TRUE,order=order,point.est = deriv.points,s.t.n = s.t.n,noise.reduction = noise.reduction)$Point.est
    f.x_h = second.deriv.estimates[1]

    two_f.x = 2*second.deriv.estimates[2]

    f.x__h = second.deriv.estimates[3]

    run = ((1+h)*eval.point) - eval.point
  return((f.x_h - two_f.x + f.x__h)/(run^2))

  }

    }

}
