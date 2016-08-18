#' Partial Derivative dy/dx
#'
#' Returns the numerical partial derivate of y wrt x for a point of interest.
#'
#' @param x Independent Variable
#' @param y Dependent Variable
#' @param order Controls the number of partial moment quadrant means.  Defaults to NULL to allow \link{VN.reg} to determine optimal order based on R2 of regression.  \code{order='max'} generates a more accurate derivative for well specified cases.
#' @param local.point Independent variable point to be evaluated.  Defaults to \code{median(x)}.
#' @keywords partial derivative, nonlinear regression
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' x<-seq(0,2*pi,pi/100); y<-sin(x)
#' dy.dx(x,y,local.point=1.75)
#' @export


dy.dx <- function(x,y,order=NULL,local.point=median(x)){
  IV=x
  DV=y


    reg.output <- VN.reg(IV,DV,plot = FALSE,return.values = TRUE,order=order)
    xonly.output <- VN.reg(IV,DV,type="XONLY",plot=FALSE,return.values = TRUE,order = order)


      if(reg.output$R2 >= xonly.output$R2){

      output<- reg.output$derivative
      if(output[,3][which(local.point<output[,3])-1][1]<local.point){
        return(output[,1][which(local.point<output[,3])][1])}
      else{
        return(mean(c(output[,1][which(local.point<output[,3])][1],output[,1][which(local.point<output[,3])-1][1])))
      }
    }
    else {
      output<- xonly.output$derivative
      if(output[,3][which(local.point<output[,3])-1][1]<local.point){
        return(output[,1][which(local.point<output[,3])][1])}
      else{
        return(mean(c(output[,1][which(local.point<output[,3])][1],output[,1][which(local.point<output[,3])-1][1])))
      }
    }
  }


